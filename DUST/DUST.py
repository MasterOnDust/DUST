import pandas as pd

import cartopy as cr
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.dates as mdates
import matplotlib as mpl
from matplotlib.ticker import LogFormatter , LogFormatterMathtext, LogFormatterSciNotation

import numpy as np
import sys
import os
import glob

from DUST.utils.read_output import *
from DUST.utils.plotting import mpl_base_map_plot, mpl_base_map_plot_xr
from DUST.utils.maps import base_map_func
from DUST.utils.utils import _gen_log_clevs, _gen_flexpart_colormap, _fix_time_flexdust
from DUST.utils.multiply_emsfield import multi_flexpart_flexdust

import xarray as xr
import xarray

from IPython import embed


mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.labelsize'] = 'medium'
mpl.rcParams['ytick.labelsize'] = 'medium'
#Note Does not work with nested output yet!

def multiply_flexpart_flexdust(flexdust, outpath='./out', locations=None, ncFiles = None, path = None, zlib=True):
    """
    DESCRIPTION
    ===========

        Goes through all subdirectories in path looking for flexpart netcdf files with .nc
        file extension. Then multiply corresponding emission sensitivities with emission field

    USAGE:
    ======

        flexdust           : path to flexdust output folder/output or xarray.Dataset containing flexdust output
        outpath (optional) : directory for where the combined data are going to be stored
        locations(optional): receptor location in flexpart either interger or string containing the name,
                             correspoding to RELCOM in FLEXPART release file. If not provided all locations is used.
        ncFiles (optional) : List of paths to flexpart output files
        path (optional)    : path to top directory of flexpart output which is search recursively looking
                             for FLEXPART netcdf files, either ncFiles or path has to be provided!

        returns: python list containing paths to multiplied flexpart / flexdust output

    """

    if ncFiles:
        ncFiles = ncFiles
    elif  path:
        ncFiles = glob.glob(path + "/**/grid*.nc", recursive=True)
    else:
        raise(NameError("Both path and ncFiles is None"))

    files = []
    if os.path.isdir(outpath):
        pass
    else:
        os.mkdir(outpath)


    d = xr.open_dataset(ncFiles[0])
    for i , com in enumerate(d.RELCOM):
        loc = str(com.values)[2:].strip().split()[0]
        if locations:
            if loc or i in locations:
                files.append(multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib), )
        else:
            files.append(multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib))


    d.close()
    return files

def read_multiple_flexpart_output(path):
    """
    DESCRIPTION
    ===========
        Reads in FLEXPART output files in a parent directory
        recursively looks for netcdf files with *.nc file extension

    USAGE
    =====
        dsets = read_multiple_flexpart_output(path)

        return : python dictionary containing xarray datasets

    """
    def not_usefull(ds):
        essentials = ['spec001_mr']
        return  [v for v in ds.data_vars if v not in essentials]

    def pre(ds):
        ds = ds.rename(dict(longitude = 'lon', time= 'btime', latitude='lat'))
        ds = ds.assign(btime = ds.btime/(3600))
        ds.btime.attrs['units'] = 'time along backtrajectory in hours'

        ds = ds.assign_coords(time=pd.to_datetime(ds.iedate + ds.ietime))



        return ds
    if isinstance(path,list)==True:
        nc_files = path
    else:
        try:
            df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
            nc_files = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
        except FileNotFoundError:
            nc_files = glob.glob(path + "**/output/grid*.nc", recursive=True) #recursively find FLEXPART output files
    
    d0 = xr.open_dataset(nc_files[0])
    relcom = str(d0.RELCOM.values)[2:].strip()[:-1].split()
    dend = xr.open_dataset(nc_files[-1])
    data_vars = {
        'RELINT' : d0.RELSTART - d0.RELEND, 
        'RELCOM' : d0.RELCOM,
        'RELLNG' : d0.RELLNG1,
        'RELLNG2'  : d0.RELLNG2,
        'RELLAT' : d0.RELLAT1,
        'RELLAT2' : d0.RELLAT2,
        'RELZZ1' : d0.RELZZ1,
        'RELZZ2' : d0.RELZZ2,
        'RELKINDZ': d0.RELKINDZ,
        'ORO' : d0.ORO.rename(longitude='lon', latitude='lat'),
        'RELPART' : d0.RELPART
    }

    d0.close()

    dsets = xr.open_mfdataset(nc_files, preprocess=(pre), decode_times=False,
                                combine='nested', concat_dim='time', parallel=True,data_vars='minimal')



    dsets = dsets.drop(not_usefull(dsets))

    dsets = dsets.assign(data_vars)
    dsets.RELINT.attrs['long_name'] = 'Time intervall of particle release'
    edate = pd.to_datetime(dsets.time[-1].values)

    dsets = dsets.assign_attrs(dict(iedate = edate.strftime('%Y%m%d'),
                                ietime = edate.strftime('%H%M%S'),
                                dataVar='spec001_mr',
                                relpart = 'Trajectories released per time step'.format(d0.RELPART),
                                history = """FLEXPART emission sensitvity, btime hours since release,
btime should be summed up before time averaging of the output. FLEXPART model output at each time step in has been concatenated 
along the time dimension. Each FLEXPART simulation has the exact same  model settings."""))

    return dsets


def read_flexpart_output(flexpart_output,ldirect = -1 ,nclusters=5):
    """
    DESCRIPTION
    ===========

        Takes in the path to a flexpart output directory or a FLEXPART output file.
        If a path to a output directory is provided then it will read all the FLEXPART
        output files present in that directory and return it as a python directory. If a full path
        to a FLEXPART output file is provided then it will only read that specific output file.


    USAGE
    =====

        flexpart_output : str containing path to FLEXPART output directory or a FLEXPART
                          output file.

        fp = read_flexpart(flexpart_output)

        return : python directory / pandas.dataframe / xarray.dataset

    """
    def not_usefull(ds):
        essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                  'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
        return  [v for v in ds.data_vars if v not in essentials]
    outs = {}
    if os.path.isfile(flexpart_output):
        path = flexpart_output.split('/')
        files = [path[-1]]
        out_dir = str()
        for s in path[:-1]:
            out_dir += s + '/'
    elif os.path.isdirectory(flexpart_output):
        files = os.listdir(flexpart_output)
        if flexpart_output.endswith('/'):
            out_dir = flexpart_output
        else:
            out_dir = flexpart_output + '/'

    else:
        raise(ValueError('flexpart_output has to be a string not {}'.format(type(flexpart_output))))


    for f in files:
        if f.endswith('.nc'):
            ncfile = xr.open_dataset(out_dir + f, decode_times=False)
            ncfile = ncfile.rename(dict(latitude = 'lat', longitude='lon', time='btime'))
            s_time = pd.to_datetime(ncfile.iedate + ncfile.ietime)
            # ncfile = ncfile.expand_dims({'time':1})
            if ldirect < 0:
                ncfile = ncfile.drop_vars(not_usefull(ncfile))
            for data_var in ncfile.data_vars:
                if 'btime' and 'pointspec' in ncfile[data_var].dims:

                    da = ncfile[data_var]
                    da = da.assign_coords(time=s_time)
                    # print(da.dims)
                    da = da.expand_dims({'time': 1})
                    # print(da.dims)
                    ncfile[data_var] = da
                else:
                    continue


            outs['data'] = ncfile
        elif 'COMMAND.namelist' in f:
            com_dict = read_command_namelist(out_dir)
            outs['command'] = com_dict
        elif 'RELEASES.namelist' in f:
            rel_df = read_release_namelist(out_dir)
            outs['releases'] = rel_df
        elif 'OUTGRID.namelist' in f:
            out_dict = read_outGrid_namelist(out_dir)
        elif 'trajectories.txt' in f:
            outs['trajectories'] = read_trajectories(out_dir, nclusters)
        else:
            continue
    if len(outs.keys()) == 1:
        key = [key for key in outs.keys()][0]
        return outs[key]
    else:
        return outs


def read_flexdust_output(path_output):
    """
    DESCRIPTION
    ===========
        Takes in a path to a FLEXDUST output folder/output file.
        If path directly to the FLEXDUST output file is provided then xarray.dataset is returned
        otherwise return a python dictionary, contaning the information in the summary file and
        a xarray.dataset

    USAGE
    =====
        fd = read_flexdust_output(path_output)

        returns : python dictionart/xarray.dataset

    """
    outdir = {}
    if path_output.endswith('.nc'):
        outdir = _fix_time_flexdust(path_output)
    else:
        for output_file in os.listdir(path_output):
            if output_file.endswith('.nc'):
                ncfile = path_output + output_file

                outdir['dset'] =_fix_time_flexdust(ncfile)
            elif output_file.endswith('.txt'):
                summary_dir = read_flex_dust_summary(path_output + output_file)
                outdir['Summary'] = summary_dir
            else:
                continue
    if len(outdir.keys()) == 1:
        key = [key for key in outs.keys()][0]
        return outdir[key]
    else:
        return outdir
class DUSTBase(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

        self._x_dim = None
        self._y_dim = None
        self._RELCOM = None

        if "lon" in self._obj.dims and "lat" in self._obj.dims:
            self._x_dim = "lon"
            self._y_dim = "lat"
        elif "longitude" in self._obj.dims and "latitude" in self._obj.dims:
            self._x_dim = "longitude"
            self._y_dim = "latitude"
        else:
            raise(ValueError("Spatial dimmensions not recognized try to rename spatial dims"))


    def sum(self,dim,keep_attrs=True, **kwargs):
        if 'RELCOM' in self._obj.data_vars:
            RELCOM = self._RELCOM
            ds_sum = self._obj.sum(dim, keep_attrs=keep_attrs, **kwargs)
            ds_sum = ds_sum.assign(RELCOM = self._RELCOM)
        else:
            ds_sum = self._obj.sum(dim, keep_attrs=keep_attrs, **kwargs)

        return ds_sum


    def mean(self, dim, keep_attrs=True, **kwargs):
        if 'RELCOM'in self._obj.data_vars:

            ds_mean = self._obj.mean(dim, keep_attrs=keep_attrs, **kwargs)
            ds_mean = ds_mean.assign(RELCOM = self._RELCOM)
        else:

            ds_mean = self._obj.mean(dim, keep_attrs=keep_attrs, **kwargs)
        return ds_mean

    def __repr__(self):
        return self._obj.__repr__()

@xr.register_dataset_accessor('fp')
class FLEXPART(DUSTBase):
    def __init__(self, xarray_obj):
        super(FLEXPART, self).__init__(xarray_obj)

        self._RELCOM = self._obj.RELCOM


    def select_receptor_point(self,pointspec, nAgeclass =False, inplace = False):
        """
        DESCRIPTION
        ===========
            Selects receptor point

        USAGE:
        =====
            dset_point = dset.fp.select_receptor_point(pointspec):

        """

        if inplace == True:
            _obj = self._select_receptor_point(self._obj,pointspec,nAgeclass)
            _obj = _obj.sel(nageclass=0)
            self._obj = _obj
        else:
            _obj = self._select_receptor_point(self._obj,pointspec,nAgeclass)
            _obj = _obj.sel(nageclass=0)
        return _obj

    def _select_receptor_point(self,dataset,pointspec, nAgeclass =False):

        return dataset.sel(pointspec=pointspec,  numspec=pointspec, numpoint=pointspec)


    def plot_total_column(self,**kwargs):
        """
        DESCRIPTION
        ===========

            Integrates the height dimmension over all model output layers
            and plots the emission sensitivity

        USAGE
        =====
            fig,ax = dset.fp.plot_emission_sensitivity()

            returns: matplotlib.figure, matplotlib.axes
        """
        fig, ax = self.plot_emission_sensitivity(height=None, **kwargs)
        return fig, ax

    def _set_da_attrs(self,data, _obj=None):
        if _obj !=None:
            _obj = _obj
        else:
            _obj = self._obj


        data = data.assign_attrs(lon0 =_obj.RELLNG1.values,
                        lat0 = _obj.RELLAT1.values,
                        relcom = _obj.RELCOM.values )
        return data

    def make_data_container(self,height=None,
                            btimeRange=None, timeRange=None, data_var='spec001_mr', ds=None):
        """
        DESCRIPTION
        ===========
            Constructs the xarray.DataArray which is expected by the mpl_base_map_plot
            routine.

        """
        if ds != None:
            data = ds[data_var]
        else:
            data = self._obj[data_var]
            print(data.dims)

        if height == None:
            data = data.sum(dim='height')
        else:
            try:
                data = data.sel(height=height)
            except KeyError:
                print('Height = {} is not a valid height, check height defined in OUTGRID'.format(height))
        if 'btime' in data.dims:
            if btimeRange ==None:
                data= data.sum(dim='btime', keep_attrs=True)
            else:
                try:
                    data = data.sel(time=btimeRange).sum(dim='time', keep_attrs=True)
                except KeyError:
                    print("Invalid time range provided check `timeRange`")


        if 'time' in data.dims:
            if timeRange == None:
                data = data.mean(dim='time', keep_attrs=True)
            else:
                try:
                    data = data.sel(time=timeRange).mean(dim='time', keep_attrs=True)
                except KeyError:
                    print("Invalid time range provided check `btimeRange`")


        data = self._set_da_attrs(data, ds)
        return data

    def plot_emission_sensitivity(self, height,
                                    data_var = 'spec001_mr',
                                    point = None,
                                    plotting_method = 'pcolormesh',
                                    info_loc = 'lower right',
                                    log = True,
                                    vmin = None,
                                    vmax = None,
                                    figure=None,
                                    ax =None,
                                    title=None,
                                    mapfunc = None,
                                    extent = None,
                                    btimeRange = None,
                                    timeRange = None,
                                    **fig_kwargs):
        """
        Description
        ===========

            Main function for plotting flexpart emission sensitivity

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
                figure          : matplotlib.figure object (default = None),
                                  if not provided new figure is created
                mapfunc         : function for making cartopy map, mapfunc should only take,
                                  cartopy.GeoAxes as argument and return cartopy.GeoAxes
                title           : plot title, if None default title is created.
                extent          : extent of the output map,if not already set, default
                                  [70,120, 25, 50]
                btimeRange      : How far back in time do you want to at emission sensitivity,
                                  using slice('startdate', 'enddate')
                timeRange       : Forward time slice, only applicable when looking at many flexpart
                                  output files
                **fig_kwargs    : figure kwargs when creating matplotlib.figure object

        """

        if 'pointspec' not in self._obj.dims:
            _obj = self._obj
        elif point in self._obj.pointspec:
            # print(pointspec)
            _obj = self.select_receptor_point(point, inplace=False)

        else:
            raise(ValueError('Pointspec not selected use select_receptor first or set point != None'))
        data = self.make_data_container(height, btimeRange, timeRange, data_var, ds =_obj)


        if 'nageclass' in data.dims:
            data = data.sel(nageclass=0)
        if figure == None and ax==None:
            fig, ax = plt.subplots(1,1, subplot_kw=dict(projection=ccrs.PlateCarree()),**fig_kwargs)
        elif ax==None:
            ax = fig.add_axes(ax)
        elif figure==None and ax != None:
            raise(ValueError('matplotlib.axes has to have a corresponding figure'))
        else:
            fig=figure

        if  mapfunc:
            ax = mapfunc(ax)
        else:
            ax = base_map_func(ax)

        fig, ax = mpl_base_map_plot(data, ax=ax, fig=fig,
                                    plotting_method=plotting_method,
                                    mark_receptor = True, vmin=vmin, vmax=vmax
                                    )
        start_date = pd.to_datetime(_obj.iedate + _obj.ietime, yearfirst=True)
        rel_start = start_date + pd.to_timedelta(_obj.RELSTART.values)
        rel_end = start_date + pd.to_timedelta(_obj.RELEND.values)
        rel_part = _obj.RELPART.values

        info_str = ('FLEXPART {}\n'.format(_obj.source[:25].strip()) +
                    'Release start : {}\n'.format(rel_start.strftime('%d/%m/%y %H:%M')) +
                    'Release end : {}\n'.format(rel_end.strftime('%d/%m/%y %H:%M')) +
                    'Particles released : {:.2E}'.format(rel_part)
                    )

        anc_text = AnchoredText(info_str, loc=info_loc ,bbox_transform=ax.transAxes,prop=dict(size=8))
        ax.add_artist(anc_text)

        if extent == None and ax.get_extent() == None:
            ax.set_extent([70,120, 25, 50], crs=ccrs.PlateCarree())
        elif extent ==None:
            pass
        else:
            ax.set_extent(extent)

        if title != None:
            ax.set_title(title)
        else:
            if _obj.ldirect == -1:
                ax.set_title('FLEXPART backward trajectory simulation starting at {} UTC'.format(start_date.strftime('%b %d %Y %H%M'))
                ,fontsize=16)
            else:
                ax.set_title('FLEXPART forward trajectory simulation starting at {} UTC'.format(start_date.strftime('%b %d %Y %H%M'))
                ,fontsize=16)
        return fig, ax


@xr.register_dataset_accessor('fd')
class FLEXDUST:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset

    def __repr__(self):
        return self._obj.__repr__()


    def _integrate_area(self, unit, dset=None, inplace=False):
        if dset != None:
            _obj=dset
        else:
            _obj=self.obj

        if unit == 'kg':
            emission_series = (_obj['Emission']*_obj.area.values).sum(dim=('lon','lat'))
        elif unit =='kg/m2':
            emission_series = _obj['Emission'].sum(dim=('lon','lat'))
            units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

        else:
            raise(ValueError("method` param {} is not a valid one. Try 'kg' or kg/m2".format(unit)))
        if inplace:
            self._obj = _obj

        return emission_series

    def get_total_emission(self, unit='kg'):
        area_integrated = self._integrate_area(unit=unit)
        tot_emissions = area_integrated.sum(dim='time')
        return tot_emissions

    def emission_time_series_to_df(self):
        emissions_kg = self._integrate_area('kg')
        emissions_kg_m2 = self._integrate_area('kg/m2')
        df = pd.DataFrame(columns=['emissions(kg)','emissions(kg/m2)'], index=emissions_kg.time.values)
        date0 = np.datetime_as_string(self._obj.time[0].values, unit='D')
        date_end = np.datetime_as_string(self._obj.time[-1].values, unit='D')
        df['emissions(kg)'] = emissions_kg.to_dataframe()
        df['emissions(kg/m2)'] = emissions_kg_m2.to_dataframe()
        return df

    def emission_time_series_to_csv(self, filename=None):
        df = self.emission_time_series_to_df()

        if filename == None:
            filename = 'Emissions' + df.index[0].strftime(format='%b_%d_%Y') + df.index[-1].strftime(format='%b_%d_%Y') + '.csv'
        else:
            filename =filename
        df.to_csv(filename)

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


    def resample_data(self, freq, method='mean', dset=None):
        if dset != None:
            _obj = dset
        else:
            _obj = self._obj
        if method == 'mean':
            _obj =  _obj.resample(time=freq).mean()
        elif method =='sum':
            _obj = _obj.resample(time=freq).sum()
        else:
            raise ValueError("`method` param '%s' is not a valid one." % method)

@xr.register_dataset_accessor('srr')
class PROCESS_SRR:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        varNames = ['WetDep', 'DryDep', 'Conc']
        if 'WetDep' in self._obj.data_vars:
            var = 'WetDep'
        elif 'DryDep' in self._obj.data_vars:
            var = 'DryDep'
        elif 'Conc' in self._obj.data_vars:
            var = 'Conc'
        else:
            raise(KeyError('''Dataset does not contain flexpart SRR values {},
                                does not have the correct dimensions'''.format(self._obj.dims)))
        self._obj = self._obj.assign_attrs(varName=var)
        self.var = var

    def make_time_seires(self,timeRange=None, btimeRange=None):
        _obj = self._obj
        if timeRange !=None:
            _obj = _obj.sel(time=timeRange)
        if btimeRange != None:
            _obj= _obj.sel(btime = btimeRange)
        b0 = _obj.btime[0].values
        b_end = _obj.btime[-1].values

        _obj = _obj.sum(dim='btime', keep_attrs=True)

        data = _obj[self.var]

        data = data.sum(dim=['lat','lon'], keep_attrs=True)
        _obj = _obj.drop_dims(['lat','lon'])

        _obj = _obj.assign({self.var : data})
        s_time = pd.to_datetime(_obj.time[0].values).strftime('%Y%m%d %H:%M')
        e_time = pd.to_datetime(_obj.time[-1].values).strftime('%Y%m%d %H:%M')
        _obj = _obj.assign_attrs({'start-date': s_time, 'end-date' : e_time,
                                 'bstart': '{} s'.format(b0), 'bstop': '{} s'.format(b_end)})


        return _obj

    def plot_time_series(self, timeRange=None, btimeRange=None, ax=None,
                        fig=None, **plot_kwargs):
        if 'btime' in self._obj.dims:

            _obj = self.make_time_seires(timeRange, btimeRange)
        
        else:
            _obj = self._obj

        if ax == None: 
            ax=plt.gca()

        if 'label' in plot_kwargs.keys():
            _obj[self.var].plot(ax=ax, **plot_kwargs)
        else:
            _obj[self.var].plot(label = _obj.receptor_name, ax=ax, **plot_kwargs)

        


        # locator = mdates.AutoDateLocator(minticks=minticks, maxticks=maxticks)
        locator = mdates.DayLocator(interval=3)
        locator_minor = mdates.HourLocator(interval=12) 

        formatter = mdates.DateFormatter('%m-%d')
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_minor_locator(locator_minor)
        ax.xaxis.set_major_formatter(formatter)

        ax.grid(linestyle='-')

        return ax

    def make_data_container(self, timeRange=None, btimeRange=None, reduce = 'mean'):
        _obj = self._obj

        if timeRange !=None:
            _obj = _obj.sel(time=timeRange)
        if btimeRange != None:
            _obj= _obj.sel(btime = btimeRange)
        _obj = _obj.sum(dim='btime', keep_attrs=True)
        s_time = pd.to_datetime(_obj.time[0].values).strftime('%Y%m%d %H:%M')
        e_time = pd.to_datetime(_obj.time[-1].values).strftime('%Y%m%d %H:%M')
        if reduce =='mean':
            _obj = _obj.mean(dim='time', keep_attrs=True)
            _obj = _obj.assign_attrs(info='Mean {}'.format(self.var))
        elif reduce == 'cumulative':
            _obj = _obj.sum(dim='time', keep_attrs=True)
            _obj = _obj.assign_attrs(info='cumulative {}'.format(self.var))
        else:
            raise(ValueError('''reduce = {} is not a valid option,
                                shoud be either mean or cumulative'''.format(reduce)))
        _obj = _obj.assign_attrs({'start-date': s_time, 'end-date' : e_time})
        return _obj


    def plot_source_contribution_map(self, timeRange=None,
                                            btimeRange= None,
                                            reduce='mean',
                                            vmin = None,
                                            vmax = None,
                                            fig=None,
                                            ax = None,
                                            **fig_kwargs):
        _obj = self.make_data_container(timeRange, btimeRange, reduce)
        if fig == None and ax==None:
            fig, ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()} ,**fig_kwargs)
        elif ax !=None and fig == None:
            raise(ValueError('matplotlib.axes has to have a corresponding figure'))

        elif ax==None and fig != None:
            ax = fig.add_axes(ax)
        else:
            pass

        ax = mpl_base_map_plot_xr(_obj,ax,vmin=vmin,vmax=vmax, mark_receptor=True)
