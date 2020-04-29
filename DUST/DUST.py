import os

import pandas as pd
import xarray as xr
from xarray import Dataset
import cartopy as cr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from utils.read_output import *
import matplotlib as mpl
import numpy as np
import sys
from IPython import embed
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

@xr.register_dataset_accessor('dust')
class DUST:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset

    def plot_emission_sensitivity(self,pointspec ,
                                    figure=None,
                                    height = None ,
                                    title=None,
                                    extent = None, 
                                    cmap ='YlOrBr'):
        units = self._obj['spec001_mr'].units
        data = self._obj['spec001_mr'][:,pointspec,:,:,:,:]
        if height == None:
            data = data.sum(dim = 'height') 
        else:
            height_index = np.argwhere(self._obj.height.values == height)
            if height_index.shape[1] == 0:
                raise IndexError
            height_index = np.reshape(height_index, height_index.shape[1])
            data = data[:,:,height_index,:,:]
    
        data = data.sum(dim = 'time')

        lons = self._obj.longitude
        lats = self._obj.latitude
        
        if figure == None:
            fig = plt.figure()
        else:
            fig = figure
        ax = plt.axes(projection=ccrs.PlateCarree())

        if extent == None:
            ax.set_extent()
        data = data[0,0,:,:]

        # embed()
        im = ax.contourf(lons, lats, data, transform  = ccrs.PlateCarree(),levels=10, norm=mpl.colors.LogNorm(), cmap = cmap,
                        vmin = 0.001)
        ax.coastlines()
        # ax.add_feature(land_50m)
        ax.add_feature(cr.feature.BORDERS)
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = True, color = 'grey', alpha = 0.6, linestyle = '--')
        gl.xlabels_top = False; gl.ylabels_right = False
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        plt.colorbar(im,cax=cax,label = units, extend='max')
        return fig, ax
        

