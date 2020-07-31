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
from DUST.utils.plotting import mpl_base_map_plot
from DUST.utils.maps import base_map_func
from DUST.utils.utils import _gen_log_clevs, _gen_flexpart_colormap, _fix_time_flexdust
from DUST.utils.multiply_emsfield import multi_flexpart_flexdust

import xarray as xr
import xarray

from IPython import embed


mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
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

def read_multiple_flexpart_output(path, ldirect=-1):
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
        essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                  'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

    def pre(ds):
        ds = ds.rename(rename(dict(longitude = 'lon', time= 'btime', latitude='lat')))
        ds = ds.assign_coords(time=pd.to_datetime(ds.iedate + ds.ietime))
        return ds
    

    nc_files = glob.glob(path + "/**/grid*.nc", recursive=True)
    dsets = xr.open_mfdataset(nc_files, preprocess=pre, decode_times=False, 
                            combine='nested', concat_dim='time', parallel=True)
    if ldirect < 0:
        dsets = dsets.drop(not_usefull(dsets))

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
            ncfile = xr.open_dataset(out_dir + f)
            ncfile = ncfile.rename(dict(latitude = 'lat', longitude='lon'))
            if ldirect < 0:
                outs['data'] = ncfile.drop(not_usefull(ncfile))
            else:
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
    if len(outs.keys()) == 1:
        key = [key for key in outs.keys()][0]
        return outdir[key]
    else:
        return outdir

@xr.register_dataset_accessor('fp')
class FLEXPART:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset

    def select_receptor_point(self,pointspec, nAgeclass =False):
        """
        DESCRIPTION
        ===========
            Selects receptor point

        USAGE:
        =====
            dset_point = dset.fp.select_receptor_point(pointspec): 
        """
        if nAgeclass == False:
            _obj = self._obj.sel(nageclass=0)

        return _obj.sel(pointspec=pointspec,  numspec=pointspec, numpoint=pointspec)
    def plot_total_column(self,**kwargs):
        """
        DESCRIPTION
        ===========
            
            Integrates the height dimmension over all model output layers 
            and plots the emission sensitivity 
        
        USAGE
        =====
            fig,ax = dset.fp.plot_emission_sensitivity

            returns: matplotlib.figure, matplotlib.axes
        """
        fig, ax = plot_emission_sensitivity(height=None, **kwargs)
        return fig, ax

    def plot_emission_sensitivity(self, height,
                                    plotting_method = 'pcolormesh',
                                    info_loc = 'lower right',
                                    log = True,
                                    vmin = None,
                                    vmax = None,
                                    figure=None,
                                    ax =None,
                                    title=None,
                                    extent = None, 
                                    btimeRange = None,
                                    **kwargs):
        
  
        self._obj = self._obj.sel(nageclass=0) #I'm not using ageclass at the moment, this is just to remove the dimmension
        data = self._obj.spec001_mr
        
        data = data.assign_attrs(lon0 =self._obj.RELLNG1.values,
                        lat0 = self._obj.RELLAT1.values,
                        relcom = self._obj.RELCOM.values )
        

        

        if btimeRange == None:
            data = data.sum(dim='time', keep_attrs=True)
        else:
            try:
                data = data.sel(time=btimeRange).sum(dim='time', keep_attrs=True)
            except KeyError:
                print("Invalid time range provided check `btimeRange`")
        
        if height == None:
            data = data.sum(dim='height', )
        else:
            try:
                data = data.sel(height=height)
            except KeyError:
                print('Height = {} is not a valid height, check height defined in OUTGRID'.format(height))


        fig, ax = mpl_base_map_plot(data, 
                                    plotting_method=plotting_method,
                                    mark_receptor = True
                                    )
        start_date = pd.to_datetime(self._obj.iedate + self._obj.ietime, yearfirst=True)
        rel_start = start_date + self._obj.RELSTART.values
        rel_end = start_date + self._obj.RELEND.values
        rel_part = self._obj.RELPART.values
        info_str = ('FLEXPART {}\n'.format(self._obj.source[:25].strip()) +
                    'Release start : {}\n'.format(rel_start.strftime(format='%d/%m/%y %H:%M')) +
                    'Release end : {}\n'.format(rel_end.strftime(format = '%d/%m/%y %H:%M')) +
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
            if self._obj.ldirect == -1: 
                ax.set_title('FLEXPART backward trajectory simulation starting at {} UTC'.format(start_date.strftime(format= '%b %d %Y %H%M'))
                ,fontsize=16)
            else:
                ax.set_title('FLEXPART forward trajectory simulation starting at {} UTC'.format(start_date.strftime(format= '%b %d %Y %H%M'))
                ,fontsize=16)
        return fig, ax



@xr.register_dataset_accessor('fd')
class FLEXDUST:
    def __init__(self, xarray_dset):
        self._obj = xarray_dset
        

    def _integrate_area(self, unit):
        if unit == 'kg':
            emission_series = (self._obj['Emission']*self._obj.area.values).sum(dim=('lon','lat'))
        elif unit =='kg/m2':
            emission_series = self._obj['Emission'].sum(dim=('lon','lat'))
            units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

        else:
            raise(ValueError("method` param {} is not a valid one. Try 'kg' or kg/m2".format(unit)))
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

    def plot_emission_time_series(self, start_date=None, 
                                        end_date=None,
                                        x_date_format = None,
                                        title=None,
                                        unit='kg',
                                        subtitle = None,
                                        ax = None,
                                        mark_days = None,

                                             **kwargs):
        time = self._obj['time']
        if unit =='kg':
            emissions = self._integrate_area(unit)
            units = 'kg'
        elif unit =='kg/m2':
            emissions = self._integrate_area(unit)
            units = '$\mathrm{kg}\; \mathrm{m}^{-2}$'

        else:
            raise(ValueError("method` param {} is not a valid one. Try 'kg' or kg/m2".format(unit)))
        if ax == None:
            ax = plt.axes()
        else:
            ax = ax
        ax.plot(time,emissions, **kwargs)
        date0 = np.datetime_as_string(time[0].values, unit='D')
        date_end = np.datetime_as_string(time[-1].values, unit='D')

        if title == None:
            pass
        elif title =='default':
            plt.suptitle('Total dust emissions {} - {}'.format(date0,date_end), fontsize = 18)
        else:
            plt.suptitle(title + ' {} - {}'.format(date0,date_end), fontsize = 18)
        
        if subtitle == 'default':
            lon0 = self._obj.lon.min().values; lon1 = self._obj.lon.max().values
            lat0 = self._obj.lat.min().values; lat1 = self._obj.lat.max().values
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
                            unit='kg'):

        


        date0 = np.datetime_as_string(self._obj.time[0].values, unit='D')
        date_end = np.datetime_as_string(self._obj.time[-1].values, unit='D')
        
        plt.title('start date: {} - end date {}'.format(date0, date_end), fontsize=12)

        if ax == None:
            ax = base_map_func()
        else:
            ax = ax

        if freq == None:
            pass
        else:
            self._obj = self.resample_data(freq=freq, method='sum')

        if reduce == 'mean':
            data = self._obj.Emission.mean(dim='time', keep_attrs=True)
        elif reduce == 'sum': 
            data = self._obj.Emission.sum(dim='time', keep_attrs=True)
        else:
            raise ValueError("`reduce` param '%s' is not a valid one." % unit)

        if unit == 'kg':
            data = data*self._obj.area.values[0]
            data = data.assign_attrs(units = '$\mathrm{kg}$')
        elif unit== 'kg/m2':
            data = data
            data = data.assign_attrs(units ='$\mathrm{kg}\; \mathrm{m}^{-2}$')
        else:
            raise ValueError("`unit` param '%s' is not a valid one." % unit)

        if title == None:
            plt.suptitle('FLEXDUST estimated accumulated emissions', fontsize=18,y=0.8)
        else:
            plt.suptitle(title, fontsize=18,y=0.8)
        fig, ax = mpl_base_map_plot(data,log=log, cmap=cmap)
        return fig, ax


    def resample_data(self, freq, method='mean'):
        if method == 'mean':
            return self._obj.resample(time=freq).mean()
        elif method =='sum':
            return self._obj.resample(time=freq).sum()
        else:
            raise ValueError("`method` param '%s' is not a valid one." % method)



if __name__ == "__main__":
    dset = xr.open_dataset('/opt/uio/flexpart/Compleated_runs/20190306_15/output/grid_time_20190306150000.nc')
    dset.fp.plot_emission_sensitivity(1, height=100)