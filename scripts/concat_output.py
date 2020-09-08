import time
import xarray as xr

import pandas as pd
from functools import partial
import os
import shutil
import glob

import argparse as ap

from dask.distributed import Client, LocalCluster

def not_usefull(ds):
    essentials = ['RELPART','RELXMASS', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

def pre_sum(ds, pointspec):
    ds = ds.rename(dict(longitude = 'lon', time= 'btime', latitude='lat'))
    ds = ds.sel(pointspec=pointspec, numpoint=pointspec, numspec=pointspec, nageclass=0)
    ds = ds.assign_coords(time=pd.to_datetime(ds.iedate + ds.ietime))
    ds = ds.drop(not_usefull(ds))
    ds = ds.sum(dim='btime')


    return ds


def concat_output(ncfiles,pointspec,outpath='',client=None,cluster_kwargs={}):
    """
    DESCRIPTION:
    ===========
        Concat the FLEXPART output netcdf files along a single time dimmension and 
        then store the output to a single netcdf file. It will create an new netcdf 
        file for every receptor location. NB this script is make for concatination
        of backward simulations. 
    
    USEAGE:
    ======
        
        conc_output(ncfiles)


        ncfiles : List containing path to flexpart output files.
    """
    # if client == None:
    #     cluster = LocalCluster(**cluster_kwargs)
    #     client = Client(cluster)
    
    
    def not_usefull(ds):
        essentials = ['RELPART','RELXMASS', 'spec001_mr']
        return  [v for v in ds.data_vars if v not in essentials]

    def pre_sum(ds, pointspec):
        ds = ds.rename(dict(longitude = 'lon', time= 'btime', latitude='lat'))
        ds = ds.sel(pointspec=pointspec, numpoint=pointspec, numspec=pointspec, nageclass=0)
        ds = ds.assign_coords(time=pd.to_datetime(ds.iedate + ds.ietime))
        ds = ds.drop(not_usefull(ds))


        return ds

    

    d0 = xr.open_dataset(ncfiles[0]).sel(pointspec=0, numpoint=0, numspec=0)
    relcom = str(d0.RELCOM.values)[2:].strip()[:-1].split()
    dend = xr.open_dataset(ncfiles[-1]).sel(pointspec=0, numpoint=0, numspec=0)
    ind_receptor =d0.ind_receptor
    if ind_receptor == 1:
        f_name = 'grid_time'
    elif ind_receptor == 3:
        f_name = 'grid_wetdep'
    elif ind_receptor ==4:
        f_name = 'grid_drydep'
    else:
        f_name = 'Unknown'
    
    data_vars = {
        'RELINT' : d0.RELEND - d0.RELSTART, 
        'START' : d0.RELSTART,
        'END' :  dend.RELSTART,
        'RELCOM' : d0.RELCOM,
        'RELLNG' : d0.RELLNG1,
        'RELLNG2'  : d0.RELLNG2,
        'RELLAT' : d0.RELLAT1,
        'RELLAT2' : d0.RELLAT2,
        'RELZZ1' : d0.RELZZ1,
        'RELZZ2' : d0.RELZZ2,
        'RELKINDZ': d0.RELKINDZ
    }
    attrsd0 = d0.attrs
    attrs_new = {'history': "Created " + time.ctime(time.time()),
                           'title' : 'FLEXPART backward simulation model output',
                           'iedate' : dend.iedate,
                           'ietime' : dend.ietime,
                           'varName' : 'spec001_mr'}
    d0.close()
    dend.close()
    outFileName =  '_'.join(relcom) +'_{}_{}_{}'.format(f_name,attrsd0['ibdate'], attrs_new['iedate']) + '.nc'
    pre = partial(pre_sum,pointspec=0)
    dsets = xr.open_mfdataset(ncfiles, preprocess=pre, decode_times=False,
                            combine='nested', concat_dim='time', parallel=True,data_vars='minimal')
    dsets = dsets.assign(data_vars)
    dsets = dsets.assign_attrs(attrsd0)
    dsets = dsets.assign_attrs(attrs_new)
    dsets.to_netcdf(outpath+outFileName,encoding={'spec001_mr' : {'zlib': True, 'complevel': 6}}, unlimited_dims='time')

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Concat FLEXPART output from backward simulation along a single time dimmension')
    parser.add_argument('path', help='path to top directory containing flexpart output')
    parser.add_argument('--outpath', '--op', help='where the concatinated output should be stored', default='.')
    parser.add_argument('--locations', '--loc', help='number of location to concatinate output for', default='ALL')
    parser.add_argument('--memory_limit', default='4GB', help='memory limit local dask cluster')
    parser.add_argument('--n_worker', default=8, type=int, help='Number of dask workers')
    parser.add_argument('--use_cluster', '--uc', action='store_true', help='Weather to use cluster or not')
    args = parser.parse_args()
    outpath = args.outpath
    path = args.path
    locations = args.locations
    n_workers=args.n_worker
    memory_limit=args.memory_limit
    uc = args.use_cluster
    #IF AVAILABLE_OUPUT file is created, use that instead of recursive search, slow on mounted system 
    if path.endswith('/') == False:
        path = path +'/'
    #IF AVAILABLE_OUPUT file is created, use that before recursive search, slow on mounted system 
    try:
        df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
        ncFiles = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
    except FileNotFoundError:
        ncFiles = glob.glob(path + "**/output/grid*.nc", recursive=True) #recursively find FLEXPART output files
    if uc:
        cluster = LocalCluster(n_workers=n_workers, memory_limit=args.memory_limit,dashboard_address=':4444')
        client = Client(cluster)
    else:
        client=None
    
    d = xr.open_dataset(ncFiles[0])
    relCOMS = d.RELCOM
    ind_receptor = d.ind_receptor
    d.close()
    

    if ind_receptor == 1:
        f_name = 'Conc'
    elif ind_receptor == 3:
        f_name = 'WetDep'
    elif ind_receptor ==4:
        f_name = 'DryDep'
    else:
        f_name = 'Unknown'
    
    

    e_time = pd.to_datetime(ncFiles[-1][-17:-3]).strftime('%Y-%m-%d')
    
    s_time = pd.to_datetime(ncFiles[0][-17:-3]).strftime('%Y-%m-%d')
    
    dir_p = outpath+'/'+f_name + '_FLEXPART_SRR_{}_{}'.format(s_time, e_time)
    try:
        os.mkdir(dir_p)
    except FileExistsError:

        shutil.rmtree(dir_p)
        os.mkdir(dir_p)
    dir_p = dir_p + '/'
    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()
    if locations == 'ALL':
        concat_output(ncFiles,outpath=dir_p,pointspec=i)
    else:
        for receptor in locations:
            if receptor in loc or receptor == str(i):
                concat_output(ncFiles, outpath=dir_p, pointspec=i, client=client)
            else:
                continue
