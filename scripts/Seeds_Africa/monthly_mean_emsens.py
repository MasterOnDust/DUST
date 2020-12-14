#!/usr/bin/env python

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.gridliner import Gridliner
import cartopy as cr
import cartopy.io.img_tiles as cimgt
import numpy as np
from DUST.plot.utils import _gen_flexpart_colormap
import argparse as ap
import pandas as pd

def plot_monthly_mean(paths, fname):
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul','Aug','Sep', 'Oct', 'Nov', 'Dec']
    fig, axes = plt.subplots(nrows = 3, ncols=4,subplot_kw={'projection':ccrs.PlateCarree()}, figsize=(16,14))
    fig.suptitle('Montly averaged Emission Sensitivity {}'.format(fname))
    df = pd.read_csv('african_mountains.csv')
    for ax, path, month in zip(axes.flatten(), paths, months):
        ds = xr.open_dataset(path)
        with xr.set_options(keep_attrs=True):
            ds_monthly_mean = ds.sum(dim='btime').mean(dim='time')
        ds.close()
        ds = ds_monthly_mean
        ds_plot = xr.where(ds.spec001_mr > 0.1, ds.spec001_mr, np.nan)
        stamen_terrain = cimgt.Stamen('terrain-background')
        cmap = _gen_flexpart_colormap()

        ax.add_image(stamen_terrain, 7)




        
        im = xr.plot.pcolormesh(ds_plot, norm=colors.LogNorm(0.1, 1e3),cmap=cmap,
                           extend='max', add_colorbar=False,ax=ax)
        ax.coastlines()
        ax.add_feature(cr.feature.BORDERS,color='gray')
        ax.add_feature(cr.feature.RIVERS)
        ax.add_feature(cr.feature.LAKES)
        ax.scatter(ds_monthly_mean.RELLNG1, ds_monthly_mean.RELLAT1, marker = '*', 
                   s=40, transform = ccrs.PlateCarree(), color ='black')
        gl = ax.gridlines(transform = ccrs.PlateCarree(), draw_labels = True, linestyle ='--')
        ax.scatter(df['lon'],df['lat'], color='red', s=2.8)
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        ax.set_extent((22.950000000000003, 45.0, -8.70, 15.75))
        ax.set_title(month)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
    fig.colorbar(im ,cax=cbar_ax, label='Sensitivity to emissions [S]')
    plt.savefig('figs/{}.png'.format(fname), dpi=300, bbox_inches='tight')
if __name__=="__main__":
    parser = ap.ArgumentParser(description='Plot monthly mean maps')
    parser.add_argument('paths', nargs='+', help='Paths to merged output folder')
    
    args = parser.parse_args()
    paths = args.paths
    
    #paths = glob.glob(path, recursive=True)
    fname = paths[0].split('/')[-1].split('_')[0]
    paths.sort()
    print(fname)
    plot_monthly_mean(paths, fname)
