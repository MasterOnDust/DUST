

import pandas as pd

import cartopy as cr
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.dates as mdates
import matplotlib as mpl

import numpy as np
import sys
import os

from utils.read_output import *
from utils.maps import base_map_func
from utils.utils import _gen_log_clevs, _gen_flexpart_colormap

import xarray as xr
from xarray import Dataset

from IPython import embed


mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
#Note Does not work with nested output yet!


def read_flexpart_output(output_dir, nclusters=5):
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


def read_flexdust_output(path_output_folder):
    outdir = {}
    for output_file in os.listdir(path_output_folder):
        if output_file.endswith('.nc'):
            ncfile = path_output_folder + output_file
            dset = xr.open_dataset(ncfile, decode_times=False)
            s_date = dset.startdate.values 
            s_hour = dset.starthour.values
            s_dT =  pd.to_timedelta(s_hour,unit='h') 
            sTime = pd.to_datetime(s_date, format='%Y%m%d') + s_dT 
            time_index = np.unique(np.reshape(dset.Date.values,dset.Date.shape[0]*2))
            time_freq = int((time_index[1]- time_index[0])/60/60)
            nTimeSteps = len(time_index)-1
            time_index = pd.date_range(start='{}'.format(sTime.strftime('%Y%m%d %H:%M:%S').values[
                    0]), periods=nTimeSteps, freq='{}h'.format(time_freq))
            dset['time'] = time_index
            outdir['dset'] =dset
        elif output_file =='Summary.txt':
            summary_dir = read_flex_dust_summary(path_output_folder + 'Summary.txt')
            outdir['Summary'] = summary_dir
        else:
            continue

    return outdir

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


    def plot_emission_time_seires(self, start_date=None, 
                                        end_date=None,
                                        x_date_format = None,
                                        figsize=(10,6),
                                        title=None,
                                        unit='kg',
                                             **kwargs):
        time = self._obj['time']
        if unit =='kg':
            emssions = (self._obj['Emission']*self._obj.area.values).sum(dim=('lon','lat'))
            units = 'kg'
        elif unit =='kg/m2':
            emssions = self._obj['Emission'].sum(dim=('lon','lat'))
            units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

        else:
            raise(ValueError("method` param {} is not a valid one. Try 'kg' or kg/m2".format(unit)))
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        ax.plot(time,emssions, **kwargs)
        date0 = np.datetime_as_string(time[0].values, unit='D')
        date_end = np.datetime_as_string(time[-1].values, unit='D')
        lon0 = self._obj.lon.min().values; lon1 = self._obj.lon.max().values
        lat0 = self._obj.lat.min().values; lat1 = self._obj.lat.max().values
        if title == None:
            plt.suptitle('Total dust emissions {} - {}'.format(date0,date_end), fontsize = 18)
        else:
            plt.suptitle(title + ' {} - {}'.format(date0,date_end), fontsize = 18)
        plt.title('lon0 = {:.2f} , lat0 = {:.2f}, lon1 = {:.2f}, lat1 = {:.2f}'.format(lon0,lat0,lon1,lat1),fontsize = 12)
        
        ax.set_ylabel(units)
        
        
        fig.autofmt_xdate()
        ax.fmt_data = mdates.DateFormatter(x_date_format)
        ax.grid(linestyle='-')  
        fig.add_subplot(ax)

    def plot_emission_map(self, ax=None, 
                            plotting_method='pcolormesh',
                            fig=None,
                            cmap =None,
                            vmin=None,
                            vmax=None,
                            title=None,
                            log=False):

        
        if fig == None:
            fig = plt.figure(figsize=(10,8))
        else:
            fig = fig
        if title == None:
            plt.suptitle('FLEXDUST estimated accumulated emissions', fontsize=18,y=0.8)
        else:
            plt.suptitle(title)

        date0 = np.datetime_as_string(self._obj.time[0].values, unit='D')
        date_end = np.datetime_as_string(self._obj.time[-1].values, unit='D')
        
        plt.title('start date: {} - end date {}'.format(date0, date_end), fontsize=12)

        if ax == None:
            ax = base_map_func()
        else:
            ax = ax
        lons = self._obj.lon.values
        lats = self._obj.lat.values
        data = self._obj.Emission.sum(dim='time')

        
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
        
        try:
            if self._obj.Emission.units == 'kg/m2':
                units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

            else:
                units = self._obj.Emission.units
        except AttributeError:
            units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

        im.cmap.set_over(color='k', alpha=0.8)

        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        clabels = list(levels[::10])  # #clevs, by 10 steps
        clabels.append(levels[-1])  # add the last label

        cb = plt.colorbar(im,cax=cax,label = units, extend = 'max')
        cb.set_ticks(clabels)

        cb.set_ticklabels(['%3.2g' % cl for cl in clabels])
        cb.ax.minorticks_on()




        plt.axes(ax)




if __name__ == "__main__":
    dset = xr.open_dataset('/opt/uio/flexpart/Compleated_runs/20190306_15/output/grid_time_20190306150000.nc')
    dset.fp.plot_emission_sensitivity(1, height=100)