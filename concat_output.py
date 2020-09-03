import time
import xarray as xr

import pandas as pd
from functools import partial
import os
import shutil

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


def concat_output(ncfiles,locations = 'ALL',cluster=None, cluster_kwargs={}):
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
    if cluster == None:
        cluster = LocalCluster(**cluster_kwargs)
    
    
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

    

    d0 = xr.open_dataset(ncfiles[0])
    
    dend = xr.open_dataset(ncfiles[-1]).sel(pointspec=0, numpoint=0, numspec=0)

    pre = partial(pre_sum,pointspec=0)
    dsets = xr.open_mfdataset(ncfiles, preprocess=pre, decode_times=False,
                            combine='nested', concat_dim='time', parallel=True,data_vars='minimal')
    dsets = dsets.assign({
        'RELINT' : d0.RELEND - d0.RELSTART, 
        'START' : d0.RELSTART,
        'END' :  dend.RELSTART,
        'RELCOM' : d0.RELCOM,
        'RELLNG1' : d0.RELLNG1,
        'RELLNG2'  : d0.RELLNG2,
        'RELLAT1' : d0.RELLAT1,
        'RELZZ1' : d0.RELZZ1,
        'RELZZ2' : d0.RELZZ2,
        'RELKINDZ': d0.RELKINDZ
    })
    dsets = dsets.assign_attrs(d0.attrs)
    dsets = dsets.assign_attrs({'history': "Created " + time.ctime(time.time()),
                           'title' : 'FLEXPART backward simulation model output',
                           'iedate' : dend.iedate,
                           'ietime' : dend.ietime,
                           'varName' : 'spec001_mr'})

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Concat FLEXPART output from backward simulation along a single time dimmension')
    parser.add_argument('path', help='path to top directory containing flexpart output')
    parser.add_argument('--outpath', '--op', help='where the concatinated output should be stored', default='./out')
    parser.add_argument('--locations', '--loc', help='number of location to concatinate output for', default='ALL')
    args = parser.parse_args()
    outpath = 
    

    for i , com in enumerate(relCOMS):
    loc = str(com.values)[2:].strip().split()
    if locations == 'ALL':
        concat_output(ncFiles,outpath=outpath,pointspec=i)
    else:
        for receptor in locations:
            if receptor in loc or receptor == str(i):
                concat_output(outpath,ncFiles,flexdust,i, zlib=zlib)
            else:
                continue
