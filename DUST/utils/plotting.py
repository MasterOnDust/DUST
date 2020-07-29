import cartopy as cr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl

from DUST.utils.utils import _gen_log_clevs, _gen_flexpart_colormap
from DUST.utils.maps import base_map_func

def mpl_base_map_plot(data,
                    plotting_method = 'pcolormesh',
                    log = True,
                    vmin = None,
                    vmax = None,
                    **kwargs):

    """
    DESCRIPTION
    ===========
        Backbone of DUST ploting functionally, should be as general as possible
        Data should be a 2D xarray.dataarray
    USAGE:
    ======
    """
    default_options = {'cmap': None, 'ax': None, 'fig' : None, 'mark_recptor': False,
                        'figsize': (10,6)}

    default_options.update(kwargs)

    
    if default_options['fig'] == None:
        fig = plt.figure(figsize=default_options['figsize'])
    else:
        fig = figure
    
    if default_options['ax'] == None:
        ax = base_map_func()
    else: 
        ax = default_options['ax']

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
        im = ax.pcolormesh(data.lon, data.lat, data, transform  = ccrs.PlateCarree(),
                norm=norm, 
                cmap = cmap)
    elif plotting_method =='contourf':
        im = ax.contourf(data.lon,data.lat, data, transform  = ccrs.PlateCarree(),
                norm=norm, 
                cmap = cmap, levels=levels)
    else:
        raise ValueError("`method` param '%s' is not a valid one." % plotting_method)
    

    if default_options['mark_receptor'] == True:
        ax.scatter(data.lon0, data.lat0, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')
        
    im.cmap.set_over(color='k', alpha=0.8)

    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    clabels = list(levels[::10])  # #clevs, by 10 steps
    clabels.append(levels[-1])  # add the last label

    cb = plt.colorbar(im,cax=cax,label = data.units, extend = 'max')
    cb.set_ticks(clabels)

    cb.set_ticklabels(['%3.2g' % cl for cl in clabels])
    cb.ax.minorticks_on()

    plt.axes(ax)

    return fig, ax

