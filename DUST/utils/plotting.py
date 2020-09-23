
import cartopy as cr
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from matplotlib.offsetbox import AnchoredText

from DUST.utils.utils import _gen_log_clevs, _gen_flexpart_colormap
from DUST.utils.maps import base_map_func

import pandas as pd

import numpy as np

from functools import partial
from collections import namedtuple

def create_info_str(ax,info_dict, loc):
    info_str = ''
    for key, item in info_dict.items():
        info_str = info_str + '{} : {}\n'.format(key,item) 

    anc_text = AnchoredText(info_str, loc=loc ,bbox_transform=ax.transAxes,prop=dict(size=8))
    ax.add_artist(anc_text)


def mpl_base_map_plot_xr(dataset, ax,
                    plotting_method = 'pcolormesh',
                    log = True,
                    vmin = None,
                    vmax = None,
                    mark_receptor = False,
                    colorbar =True,
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
    default_options = {'cmap': None}

    default_options.update(kwargs)

    if default_options['cmap'] == None:
        cmap = _gen_flexpart_colormap()
        default_options.pop('cmap')
    else:
        cmap = default_options.pop('cmap')

    if vmin  ==None and vmax == None:
        dat_min = dataset[dataset.varName].min()
        dat_max = dataset[dataset.varName].max()
    elif vmin != None and vmax == None:
        dat_min = vmin
        dat_max = dataset[dataset.varName].max()
    elif vmin == None and vmax != None:
        dat_min = dataset[dataset.varName].min()
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
         dataset[dataset.varName].plot.pcolormesh(ax=ax,
                norm=norm,
                cmap = cmap, add_colorbar=colorbar, levels=levels,extend='max',**default_options)
    elif plotting_method =='contourf':
        ax.add_artist(dataset[dataset.varName].plot.contourf(ax=ax, transform  = ccrs.PlateCarree(),
                norm=norm,
                cmap = cmap, levels=levels, add_colorbar=colorbar,extend='max',**default_options))
    else:
        raise ValueError("`method` param '%s' is not a valid one." % plotting_method)


    if mark_receptor:
        ax.scatter(dataset.RELLNG, dataset.RELLAT, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')

    return ax


def make_animation(data ,map_func, title,
                        figsize=(16,8),
                        fps = 20,
                        extent =[70,120, 25, 50],
                        intervall = 150,
                        vmin = None,
                        vmax = None,
                        **kwargs):
    """
    DESCRIPTION
    ===========

        Create animation of image sequence made from the a 3D data array
        data is xarray.dataarray, with one temporal dimension and two spatial dimmension
        (lon and lat). Saves animation as an mp4 video file.

    """

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
    
    default_options = dict(title = 'DUST animation', comment='Movie for sequence of images'
                )
    default_options.update(kwargs)
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title=default_options['title'], artist='DUST',
                comment=default_options['comment'])
    writer = FFMpegWriter(fps=fps, metadata=metadata)
    fig ,ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
                            
    

    
    Artists = namedtuple("Artists", ("mesh", "time"))
    

    cmap = _gen_flexpart_colormap()
    levels = _gen_log_clevs(dat_min,dat_max)
    norm = mpl.colors.LogNorm(vmin=levels[0], vmax=levels[-1])

    artists = Artists(
    ax.pcolormesh(data.lon.values, data.lat.values, data[0].values, animated=True,
    transform  = ccrs.PlateCarree(),
                norm=norm,
                cmap = cmap),
    ax.text(1, 1, "", fontsize=20, transform=ax.transAxes, 
    horizontalalignment='right', verticalalignment='bottom'),)

    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    clabels = list(levels[::10])  # #clevs, by 10 steps
    clabels.append(levels[-1])  # add the last label
    cb = plt.colorbar(artists.mesh,cax=cax,label = data.units, extend = 'max')
    cb.set_ticks(clabels)

    cb.set_ticklabels(['%3.2g' % cl for cl in clabels])

    frames = [d_i for d_i in data]

    init = partial(_init_fig, ax=ax, fig=fig, artists=artists,extent=extent,map_func=map_func, 
                    title=title, data=data[0])
    update = partial(_update_artist,artists=artists)

    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, init_func=init
             ,interval=intervall, repeat_delay=5000)


    return ani

def _init_fig(fig, ax, artists, extent, map_func , title, data):
    ax.set_title(title, fontsize=22)
    ax = map_func(ax, extent)

    ax.scatter(data.lon0, data.lat0, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')

    artists.mesh.set_array([])
    return artists

def _update_artist(frame, artists):
    # print(frame.values.ravel)
    artists.mesh.set_array(frame.values.ravel())
    date = pd.to_datetime(str(frame.time.values))
    artists.time.set_text(date.strftime('%Y%m%d %H%M'))


def _animate(d_i, ax, extent):

    date = pd.to_datetime(str(d_i.time.values))
    ax.set_title(date.strftime('%Y%m%d %H%M'))
    fig, ax = mpl_base_map_plot(d_i, ax=ax,colorbar=False, mark_receptor=True)
    