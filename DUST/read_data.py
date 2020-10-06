
import os
import pandas as pd
import numpy as np
import glob
import xarray as xr

from .utils.read_output import read_command_namelist, read_outGrid_namelist, read_release_namelist, read_flex_dust_summary
from .utils.utils import _fix_time_flexdust

"""
This file contain functions for preparing data for analysis

Standard conventions:
    "time" : allways forward time
    "btime" : time along backward trajectory


"""



def read_multiple_flexpart_outputs(path):
    """
    DESCRIPTION
    ===========
        Reads in FLEXPART output files in a parent directory, 
        if there is an AVAILABLE_OUTPUT file it reads the path from there
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
                                varName='spec001_mr',
                                relpart = 'Trajectories released per time step'.format(d0.RELPART),
                                history = """FLEXPART emission sensitvity, btime hours since release,
btime should be summed up before time averaging of the output. FLEXPART model output at each time step in has been concatenated 
along the time dimension. Each FLEXPART simulation has the exact same  model settings."""))

    return dsets

def read_flexpart_metadata(path_output_folder):
    

    """
    DESCRIPTION
    ===========

        This should be rewritten so that it also can concatinated  or maybe I should 
        just remove this since xarray.open_dataset does better? I could have one function 
        called read metadata and one called read trajectories. 

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
    outs = {}
    if os.path.isfile(path_output_folder):
        path = path_output_folder.split('/')
        files = [path[-1]]
        out_dir = str()
        for s in path[:-1]:
            out_dir += s + '/'
    elif os.path.isdirectory(path_output_folder):
        files = os.listdir(path_output_folder)
        if path_output_folder.endswith('/'):
            out_dir = path_output_folder
        else:
            out_dir = path_output_folder + '/'

    else:
        raise(ValueError('flexpart_output has to be a string not {}'.format(type(path_output_folder))))


    for f in files:
        if 'COMMAND.namelist' in f:
            com_dict = read_command_namelist(out_dir)
            outs['command'] = com_dict
        elif 'RELEASES.namelist' in f:
            rel_df = read_release_namelist(out_dir)
            outs['releases'] = rel_df
        elif 'OUTGRID.namelist' in f:
            out_dict = read_outGrid_namelist(out_dir)
            continue
    if len(outs.keys()) == 1:
        key = [key for key in outs.keys()][0]
        return outs[key]
    else:
        return outs

def read_flexpart_trajectories(path_to_available_output, sdate=None, edate=None):
    """
    DESCRIPTION
    ===========
        Assume that AVAILABLE_OUTPUT is present, AVAILABLE_OUTPUT generated by the 
        list_available_output.py in the scripts folder 

    """
   
    path = path_to_available_output[:-16]
    df = pd.read_csv(path_to_available_output, index_col=0)
    trajectories = [path+row['dir_paths'] + '/'+ 'trajectories.txt' for index,row in df.iterrows()]
    return trajectories

def read_flexdust_output(path_output, **xarray_kwargs):
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
        outdir = _fix_time_flexdust(path_output, **xarray_kwargs)
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
        key = [key for key in outdir.keys()][0]
        return outdir[key]
    else:
        return outdir

def read_flexpart_output(path_output, dataVars='spec001_mr',ldirect=-1,**dset_kwargs):
    """
    DESCRIPTION
    ===========

        Read flexpart netcdf output file and prepare the dataset for further analysis.
        
    """



    
    if ldirect == -1:
        dset = xr.open_dataset(path_output,decode_times=False,**dset_kwargs)
        usefull = ['RELCOM', 'RELLNG1', 'RELLNG2', 'RELLAT1','RELLAT2', 'RELZZ1', 'RELZZ2'
          ,'RELKINDZ', 'RELSTART', 'RELEND', 'RELPART','ORO']
        usefull.append(dataVars)
        print(usefull)
        not_usefull = [v for v in dset.data_vars if v not in usefull]
        dset = dset.drop(not_usefull)
        dset = dset.rename(time='btime')
        dset['btime'] = dset.btime.assign_attrs({'long_name':'time along back trajectory', 'units':'hours'})
        dset = dset.assign_attrs({'varName':dataVars})
        dset = dset.assign_coords(btime=(dset.btime/3600).astype(np.short))
    else:
        dset = xr.open_dataset(path_output, **dset_kwargs)

    dset = dset.swap_dims({'numpoint':'pointspec'})
    dset.re
    relcoms = dset.RELCOM.str.strip().str.decode('utf-8').values

    dset = dset.assign(RELCOM = xr.DataArray(relcoms, dims=('pointspec'),
                    coords={'pointspec' : relcoms},
                  attrs={'long_name':'release point name'}))
    dset = dset.rename(pointspec='point')
    dset['point'] = dset.point.assign_attrs({'long_name': 'Name of release location'})
    dset = dset.rename({'longitude':'lon', 'latitude':'lat'})
    

    dset = dset.squeeze()    


    return dset

def read_data(path_output, dataVar='spec001_mr',**dset_kwargs):
    """
    DESCRIPTION
    ===========

        Base function from reading datafiles
    """
    
    pass

def read_grib_data(path, shortName=None,**dset_kwargs):
    """
    DESCRIPTION
    ===========
        Reading grib files for analysing model forcing. 

    """
    if shortName==None:
        ds_grib = xr.open_dataset(path, 
                          engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'heightAboveGround', 'edition':1}})
    else:
        ds_grib = xr.open_dataset(path, 
                        engine='cfgrib', backend_kwargs={'filter_by_keys':{'shortName':shortName}})
    
    return ds_grib