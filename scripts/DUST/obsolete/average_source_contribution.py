#!/usr/bin/env python

import xarray as xr
import argparse as ap
import dask
import pandas as pd
import time
import sys
import os 
from dask.distributed import LocalCluster, Client
def pre(ds):
    ds = ds.drop(['RELEND', 'RELSTART', 'RELCOM'])
    ds = ds.sum(dim='btime').mean(dim='time')
    return ds


if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Create monthly mean source contribution')
    parser.add_argument('paths_flexpart',nargs ='+', help='paths to flexpart output')
    parser.add_argument('--varName','--vn', help='name of variable to concatenate')
    parser.add_argument('--outpath', '--op',default='./', help='path to where output should be stored')
    
    args = parser.parse_args()

    paths = args.paths_flexpart
    outpath = args.outpath
    varName = args.varName
    #print(paths)
    paths.sort()
    cluster = LocalCluster(n_workers=4, threads_per_worker=1, memory_limit='16GB')
    client = Client(cluster)
   
    terminal_input = ' '.join(sys.argv)
    dsets = xr.open_mfdataset(paths,preprocess=pre,data_vars=[varName,'surface_sensitivity', 'RELPART'], 
                                concat_dim='time', combine='nested', parallel=True)
    times = pd.to_datetime([path.split('-')[-1][:-3] for path in paths])
    #times.freq = 'M'
    ds = xr.open_dataset(paths[0])
    relcom = ' '.join(str(ds['RELCOM'][0].values).strip().split(' ')[1:])
    dsets = dsets.assign_coords(time=times)
    dsets.attrs = ds.attrs
    dsets.attrs['RELCOM']=relcom
    dsets[varName].attrs = ds[varName].attrs
    dsets['surface_sensitivity'].attrs = ds['surface_sensitivity'].attrs
    dsets['lon'].attrs = ds['lon'].attrs
    dsets['lat'].attrs = ds['lat'].attrs
    #print(dsets.attrs)
    dsets.attrs['iedate'] = times[-1].strftime('%Y%m%d')  
    dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),terminal_input) + dsets.attrs['history']
    relcom = '_'.join(dsets.attrs['RELCOM'].split(' '))
    outFile_name = os.path.join(outpath, dsets.attrs['varName'] + '_' + relcom + '_' +  'monthly' + '_' + times[0].strftime('%Y%m%d') +'-'+ times[-1].strftime('%Y%m%d') + '.nc')
    print('writing to {}'.format(outFile_name))
    dsets.to_netcdf(outFile_name)


