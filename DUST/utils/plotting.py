
import cartopy as cr
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation

from DUST.utils.utils import _gen_log_clevs, _gen_flexpart_colormap
from DUST.utils.maps import base_map_func

import plotly.express as px

import pandas as pd

import numpy as np
def mpl_base_map_plot(data,
                    plotting_method = 'pcolormesh',
                    log = True,
                    vmin = None,
                    vmax = None,
                    ax = None,
                    **kwargs):
    print()
    """
    DESCRIPTION
    ===========
        Backbone of DUST ploting functionally, should be as general as possible
        Data should be a 2D xarray.dataarray
    USAGE:
    ======
    """
    default_options = {'cmap': None, 'fig' : None, 'mark_receptor': False,
                         'colorbar': True}

    default_options.update(kwargs)

    
    if ax == None:
        ax = plt.axes(projection=ccrs.PlateCarree())
        base_map_func(ax)
    else: 
        ax = ax
    if default_options['cmap'] == None:
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
        im = ax.pcolormesh(data.lon, data.lat, data.values, transform  = ccrs.PlateCarree(),
                norm=norm, 
                cmap = cmap)
    elif plotting_method =='contourf':
        im = ax.contourf(data.lon,data.lat, data.values, transform  = ccrs.PlateCarree(),
                norm=norm, 
                cmap = cmap, levels=levels)
    else:
        raise ValueError("`method` param '%s' is not a valid one." % plotting_method)
    

    if default_options['mark_receptor']:
        ax.scatter(data.lon0, data.lat0, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')
    if default_options['colorbar']:    
        im.cmap.set_over(color='k', alpha=0.8)

        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        clabels = list(levels[::10])  # #clevs, by 10 steps
        clabels.append(levels[-1])  # add the last label

        cb = plt.colorbar(im,cax=cax,label = data.units, extend = 'max')
        cb.set_ticks(clabels)

        cb.set_ticklabels(['%3.2g' % cl for cl in clabels])
        cb.ax.minorticks_on()

        plt.axes(ax)

    return ax


def make_animation(data ,map_func, 
                        figsize = (12,10), 
                        fps = 20, 
                        extent =[70,120, 25, 50], 
                        intervall = 150,**kwargs):
    """
    DESCRIPTION
    ===========

        Create animation of image sequence made from the a 3D data array
        data is xarray.dataarray, with one temporal dimension and two spatial dimmension
        (lon and lat). Saves animation as an mp4 video file. 
     
    """
    default_options = dict(title = 'DUST animation', comment='Movie for sequence of images'
                )
    default_options.update(kwargs)
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title=default_options['title'], artist='DUST',
                comment=default_options['comment'])
    writer = FFMpegWriter(fps=fps, metadata=metadata)
    fig = plt.figure(figsize=figsize, frameon=False)
    fig.suptitle('{}'.format(data.name_location), fontsize=18,y=0.8)
    ax = fig.add_subplot(1, 1, 1)
    ax =plt.axes(projection=ccrs.PlateCarree())
    frames = [d_i for d_i in data]
    ani = animation.FuncAnimation(fig, _animate, frames=frames, 
            fargs=(ax, map_func, extent), interval=intervall)
    return ani
def _animate(d_i, ax, map_func, extent):
    ax.clear()
    map_func(ax, extent= extent)
   
    date = pd.to_datetime(str(d_i.time.values))
    ax.set_title(date.strftime('%Y%m%d %H%M'))
    next_frame = mpl_base_map_plot(d_i, ax=ax,colorbar=False, mark_receptor=True)