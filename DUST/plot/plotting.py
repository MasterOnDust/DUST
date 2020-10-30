
import cartopy as cr
import cartopy.crs as ccrs

from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import matplotlib as mpl

from .utils import _gen_log_clevs, _gen_flexpart_colormap
from .maps import base_map_func

import pandas as pd

import numpy as np

from functools import partial
from collections import namedtuple


mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.labelsize'] = 'medium'
mpl.rcParams['ytick.labelsize'] = 'medium'

def plot_emission_time_series(self, time_slice = None,
                                    x_date_format = None,
                                    title=None,
                                    unit='kg',
                                    subtitle = None,
                                    fig= None,
                                    ax = None,
                                    mark_days = None,
                                    plot_kwargs = {},
                                            **fig_kwargs):
    if time_slice != None:
        _obj = self._obj.sel(time=time_slice)
    else:
        _obj = self._obj
    print(_obj.dims)
    time = _obj.time

    if unit =='kg':
        emissions = self._integrate_area(unit, _obj)
        units = 'kg'
    elif unit =='kg/m2':
        emissions = self._integrate_area(unit, obj)
        units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

    else:
        raise(ValueError("method` param {} is not a valid one. Try 'kg' or kg/m2".format(unit)))
    if fig == None and ax==None:
        fig, ax = plt.subplots(1,1,**fig_kwargs)
    elif ax==None:
        ax = fig.add_axes(ax)
    else:
        raise(ValueError('matplotlib.axes has to have a corresponding figure'))

    ax.plot(time,emissions, **plot_kwargs)
    date0 = np.datetime_as_string(time[0].values, unit='D')
    date_end = np.datetime_as_string(time[-1].values, unit='D')

    if title == None:
        pass
    elif title =='default':
        plt.suptitle('Total dust emissions {} - {}'.format(date0,date_end), fontsize = 18)
    else:
        plt.suptitle(title + ' {} - {}'.format(date0,date_end), fontsize = 18)

    if subtitle == 'default':
        lon0 = _obj.lon.min().values; lon1 = _obj.lon.max().values
        lat0 = _obj.lat.min().values; lat1 = _obj.lat.max().values
        plt.title('lon0 = {:.2f} , lat0 = {:.2f}, lon1 = {:.2f}, lat1 = {:.2f}'.format(lon0,lat0,lon1,lat1),fontsize = 12)
    elif subtitle == None:
        pass
    else:
        plt.title(subtitle, fontsize = 12)
    ax.set_ylabel(units)

    locator = mdates.AutoDateLocator(minticks=5, maxticks=12)

    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    ax.grid(linestyle='-')


def plot_emission_map(self, ax=None,
                        plotting_method='pcolormesh',
                        reduce='sum',
                        freq=None,
                        fig=None,
                        cmap =None,
                        vmin=None,
                        vmax=None,
                        title=None,
                        log=False,
                        time_slice = None,
                        mapfunc = None,
                        unit='kg',**fig_kwargs):


    if time_slice != None:
        _obj = self._obj.sel(time=time_slice)
    else:
        _obj = self._obj

    if fig == None and ax==None:
        fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=ccrs.PlateCarree()),**fig_kwargs)

    elif ax==None:
        ax = fig.add_axes(ax)

    elif fig==None and ax != None:
        raise(ValueError('matplotlib.axes has to have a corresponding figure'))
    else:
        pass

    if freq == None:
        pass
    else:
        _obj = self.resample_data(freq=freq, method='sum', dset = _obj)

    if 'time' not in _obj.dims:
        data = _obj.Emission
        date0 = pd.to_datetime(_obj.time.values).strftime('%y%m%d %H')
        date_end = pd.to_datetime(_obj.time.values).strftime('%y%m%d %H')
    elif reduce == 'mean':
        data = _obj.Emission.mean(dim='time', keep_attrs=True)
        date0 = pd.to_datetime(_obj.time.values[0]).strftime('%y%m%d %H')
        date_end = pd.to_datetime(_obj.time.values[-1]).strftime('%y%m%d %H')
    elif reduce == 'sum':
        data = _obj.Emission.sum(dim='time', keep_attrs=True)
        date0 = np.datetime_as_string(_obj.time.values[0], unit='D')
        date_end = np.datetime_as_string(_obj.time.values[-1], unit='D')
    else:
        raise ValueError("`reduce` param '%s' is not a valid one." % unit)

    if unit == 'kg':
        data = data*_obj.area.values[0]
        data = data.assign_attrs(units = '$\mathrm{kg}$')
    elif unit== 'kg/m2':
        data = data
        data = data.assign_attrs(units ='$\mathrm{kg}\; \mathrm{m}^{-2}$')
    else:
        raise ValueError("`unit` param '%s' is not a valid one." % unit)
    plt.title('start time: {} - end time {}'.format(date0, date_end), fontsize=12)

    if title == None:
        plt.suptitle('FLEXDUST estimated accumulated emissions',y=0.9, fontsize=18)
    else:
        plt.suptitle(title, fontsize=18)
    if mapfunc != None:
        ax = mapfunc(ax)

    fig, ax = mpl_base_map_plot(data,ax, fig,plotting_method,log=log, cmap=cmap, vmin=vmin, vmax=vmax)
    return fig, ax

def plot_emission_sensitivity(dset,
                                data_var,
                                ax=None,
                                plotting_method = 'pcolormesh',
                                map_func = None,
                                info_loc = 'lower right',
                                log = True,
                                vmin = None,
                                vmax = None,
                                title=None,
                                extent = None):
    """
    Description
    ===========

        Main function for plotting flexpart emission sensitivity. Require a 2D data array 
        Subsetting and so one should instead be done by using xarray.

    USAGE
    =====

        fig, ax = dset.fp.plot_emssion_sensitivity(height)

        required argument:
            height : height of flexpart output level
        optional arguments:
            point           : if xarray.dataset contains several receptor points
            plotting_method : either 'pcolormesh' or 'contourf' (default='pcolormesh')
            info_loc        : location of info_box, (default = 'lower right')
            log             : log colorscale True/False, (default = True)
            vmin            : lower bound of colorscale (default = None)
            vmax            : upper bound of colorscale (default = None)
            title           : plot title, if None default title is created.
            extent          : extent of the output map,if not already set, default
                                [70,120, 25, 50]
            btimeRange      : How far back in time do you want to at emission sensitivity,
                                using slice('startdate', 'enddate')
            timeRange       : Forward time slice, only applicable when looking at many flexpart
                                output files
            **fig_kwargs    : figure kwargs when creating matplotlib.figure object

    """
    if ax == None:
        ax = plt.axes(projection=ccrs.PlateCarree)

    if map_func == None:


    ax = mpl_base_map_plot_xr(dset, ax=ax,
                                plotting_method=plotting_method,
                                mark_receptor = True, vmin=vmin, vmax=vmax
                                )
    return ax


def create_info_str(ax,info_dict, loc):
    info_str = ''
    for key, item in info_dict.items():
        info_str = info_str + '{} : {}\n'.format(key,item) 

    anc_text = AnchoredText(info_str, loc=loc ,bbox_transform=ax.transAxes,prop=dict(size=8))
    ax.add_artist(anc_text)


def mpl_base_map_plot_xr(dataset, ax,
                    plotting_method = 'pcolormesh',
                    datavar = None,
                    log = True,
                    vmin = None,
                    vmax = None,
                    mark_receptor = False,
                    colorbar =True,
                    **kwargs):
    """
    DESCRIPTION
    ===========
        Backbone of DUST ploting functionally, should be as general as possible
        Data should be a 2D xarray.dataarray
    USAGE:
    ======
    """
    if datavar == None:
        varName = dataset.varName
    else:
        varName = datavar
    default_options = {'cmap': None}

    default_options.update(kwargs)

    if default_options['cmap'] == None:
        cmap = _gen_flexpart_colormap()
        default_options.pop('cmap')
    else:
        cmap = default_options.pop('cmap')

    if vmin  ==None and vmax == None:
        dat_min = dataset[varName].min()
        dat_max = dataset[varName].max()
    elif vmin != None and vmax == None:
        dat_min = vmin
        dat_max = dataset[varName].max()
    elif vmin == None and vmax != None:
        dat_min = dataset[varName].min()
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
         dataset[varName].plot.pcolormesh(ax=ax,
                norm=norm,
                cmap = cmap, add_colorbar=colorbar, levels=levels,extend='max',**default_options)
    elif plotting_method =='contourf':
        ax.add_artist(dataset[varName].plot.contourf(ax=ax, transform  = ccrs.PlateCarree(),
                norm=norm,
                cmap = cmap, levels=levels, add_colorbar=colorbar,extend='max',**default_options))
    else:
        raise ValueError("`method` param '%s' is not a valid one." % plotting_method)


    if mark_receptor:
        ax.scatter(dataset.RELLNG, dataset.RELLAT, marker = '*', s=40, transform = ccrs.PlateCarree(), color ='black')

    return ax



    