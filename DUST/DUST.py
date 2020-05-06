import os

import pandas as pd
import xarray as xr
from xarray import Dataset
import cartopy as cr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from utils.read_output import *
import matplotlib as mpl
import numpy as np
import sys
from IPython import embed
from utils.maps import base_map_func
from utils.utils import _gen_log_clevs, _gen_flexpart_colormap
import matplotlib as mpl

mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
#Note Does not work with nested output yet!


def read_out_directory(output_dir, nclusters=5):
    traj_df = None
    outs = {}
    for file in os.listdir(output_dir):
        if file.endswith('.nc'):
            ncfile = xr.open_dataset(output_dir + file)
            outs['ems_sens'] = ncfile 
        elif 'COMMAND.namelist' in file:
            com_dict = read_command_namelist(output_dir)
            outs['command'] = com_dict
        elif 'RELEASES.namelist' in file:
            rel_df = read_release_namelist(output_dir)
            outs['releases'] = rel_df
        elif 'OUTGRID.namelist' in file:
            out_dict = read_outGrid_namelist(output_dir)
        elif 'trajectories.txt' in file:
            traj_df = read_trajectories(output_dir, nclusters)
            outs['trajectories'] = traj_df
        else:
            continue
    
    return outs

@xr.register_dataset_accessor('fp')
class FLEXPART:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset

    
    def plot_total_column(self,pointspec, **kwargs):
        fig, ax = plot_emission_sensitivity(pointspec,height=None, **kwargs)


    def plot_emission_sensitivity(self,pointspec, height,
                                    plotting_method = 'pcolormesh',
                                    info_loc = 'lower right',
                                    log = True,
                                    vmin = None,
                                    vmax = None,
                                    figure=None,
                                    ax =None,
                                    title=None,
                                    extent = None, 
                                    cmap =None):
        ax = ax
        units = self._obj['spec001_mr'].units
        data = self._obj['spec001_mr'][:,pointspec,:,:,:,:]
        if height == None:
            data = data.sum(dim = 'height')
            data = data.sum(dim = 'time')
            data = data[0,:,:]
 
        else:
            height_index = np.argwhere(self._obj.height.values == height)
            if height_index.shape[1] == 0:
                raise ValueError("height param {} is not a valid one.".format(height))

            height_index = np.reshape(height_index, height_index.shape[1])
            data = data[:,:,height_index,:,:]
            data = data.sum(dim = 'time')
            data = data[0,0,:,:]

        lons = self._obj.longitude
        lats = self._obj.latitude
        lon0 = self._obj.RELLNG1[pointspec]
        lat0 = self._obj.RELLAT1[pointspec]



        if figure == None:
            fig = plt.figure()
        else:
            fig = figure
        
        if ax == None:
            ax = base_map_func()
            
        if extent == None and ax.get_extent() == None:
            ax.set_extent([70,120, 25, 50], crs=ccrs.PlateCarree())
        elif extent ==None:
            pass
        else:
            ax.set_extent(extent)
        if cmap == None:
            cmap = _gen_flexpart_colormap()
        else:
            cmap = cmap

        if vmin  ==None and vmax == None:
            dat_min = data.min()
            dat_max = data.max()
        elif vmin != None and vmax == None:
            dat_min = vmin
            dat_max = data.max()
        elif vmin == None and vmax != None:
            dat_min = data.min()
            dat_max = vmax
        else:
            dat_max = vmax
            dat_min = vmin

        if log:
            levels = _gen_log_clevs(dat_min, dat_max)
            norm = mpl.colors.LogNorm(vmin=levels[0], vmax=levels[-1])
        else:
            levels = list(np.arange(dat_min, dat_max, (dat_max - dat_min) / 100))
            norm = None

        if plotting_method == 'pcolormesh':
            im = ax.pcolormesh(lons, lats, data, transform  = ccrs.PlateCarree(),
                    norm=norm, 
                    cmap = cmap)
        elif plotting_method =='contourf':
            im = ax.contourf(lons,lats, data, transform  = ccrs.PlateCarree(),
                    norm=norm, 
                    cmap = cmap, levels=levels)
        else:
            raise ValueError("`method` param '%s' is not a valid one." % plotting_method)

        ax.scatter(lon0, lat0, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')
            
        im.cmap.set_over(color='k', alpha=0.8)

        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        clabels = list(levels[::10])  # #clevs, by 10 steps
        clabels.append(levels[-1])  # add the last label

        cb = plt.colorbar(im,cax=cax,label = units, extend = 'max')
        cb.set_ticks(clabels)

        cb.set_ticklabels(['%3.2g' % cl for cl in clabels])
        cb.ax.minorticks_on()

        plt.axes(ax)

        # Create simulation info string 
        start_date = pd.to_datetime(self._obj.iedate + self._obj.ietime, yearfirst=True)
        rel_start = start_date + self._obj.RELSTART[pointspec].values
        rel_end = start_date + self._obj.RELEND[pointspec].values
        rel_part = self._obj.RELPART[pointspec].values
        info_str = ('FLEXPART {}\n'.format(self._obj.source[:25].strip()) +
                    'Release start : {}\n'.format(rel_start.strftime(format='%d/%m/%y %H:%M')) +
                    'Release end : {}\n'.format(rel_end.strftime(format = '%d/%m/%y %H:%M')) +
                    'Particles released : {:.2E}'.format(rel_part) 
                    )

        anc_text = AnchoredText(info_str, loc=info_loc ,bbox_transform=ax.transAxes,prop=dict(size=8))
        ax.add_artist(anc_text)

        if title != None:
            ax.set_title(title)
        else:
            if self._obj.ldirect == -1: 
                ax.set_title('FLEXPART backward trajectory simulation starting at {} UTC'.format(start_date.strftime(format= '%b %d %Y %H%M'))
                ,fontsize=16)
            else:
                ax.set_title('FLEXPART forward trajectory simulation starting at {} UTC'.format(start_date.strftime(format= '%b %d %Y %H%M'))
                ,fontsize=16)
        return fig, ax


@xr.register_dataset_accessor('dust')
class FLEXDUST:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset
    
    def plot_emission_field(self):
        pass
        


if __name__ == "__main__":
    dset = xr.open_dataset('/opt/uio/flexpart/Compleated_runs/20190306_15/output/grid_time_20190306150000.nc')
    dset.fp.plot_emission_sensitivity(1, height=100)