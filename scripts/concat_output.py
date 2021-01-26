#!/usr/bin/env python
import os
from IPython import embed
from dask.distributed import Client, LocalCluster 
from DUST.read_data import read_multiple_flexpart_outputs
from DUST.utils.utils import arg_parser
from distributed import client
import pandas as pd
import shutil
import numpy as np
import sys
import argparse as ap
from dask.diagnostics import ProgressBar
import time
def slice_data(dset, x0, x1, y0, y1):
    dset = dset.sel(lon=slice(x0,x1),lat=slice(y0,y1))
    spatial_attrs_names = ['outlon0', 'outlon1', 'outlat0', 'outlat1']
    if y1 == None:
        y1 = len(dset.lat) * dset.dyout + dset.outlat0
    if x1 == None:
        x1 = len(dset.lon) * dset.dxout + dset.outlon0
    # Update attribute if output grid has been sliced
    for key, outloc in zip(spatial_attrs_names, [x0,x1, y0,y1]):
        if outloc != None:
            dset.attrs[key] = outloc
    return dset

def prepare_output_dataset(dset, height, point, fileName=None):

    dset = dset.squeeze()
    d0 = dset.time[0].dt.strftime('%Y%m%d%-%H')
    d1 = dset.time[-1].dt.strftime('%Y%m%d%-%H')
    Loc_name = str(dset.RELCOM.values).strip().split(' ')
    if fileName==None:
        file_name = "{}_{}_{}_{}.nc".format('_'.join(Loc_name),h, d0, d1) 
    else:
        file_name=fileName
    dset = dset.drop_vars('RELCOM')

    dset.attrs['relcom'] = Loc_name
    dset[dset.varName] = dset[dset.varName].astype(np.float32)
    dset.attrs['source'] = dset.attrs['source'] + ' '.join(nc_files)  
    dset.attrs['filename']=file_name 
    return dset 

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Concatenate FLEXPART simulations')
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--height', help='which levels to concatenate output for', type=int, default=None)
    parser.add_argument('paths', help='Paths to netcdf files to be concatenated',nargs='+')
    parser.add_argument('--outpath', '--op', help='Where the output should be stored', default='.')
    parser.add_argument('location', '--loc', help='Number corresponding to specific location in dataset', default=None, type=int)
    args = parser.parse_args()
    location = args.location
    nc_files = args.paths
    outpath = args.outpath
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0
    height = args.height
    
    # Load the netCDF files
    dset = read_multiple_flexpart_outputs(nc_files, height=height, location=location, parallel=False)
    dset = slice_data(dset, x0, x1, y0,y1)
    dset=prepare_output_dataset(dset, height, location)
    shape_dset = dset[dset.varName].shape
    encoding = {'zlib': True, 'complevel':9,
        'fletcher32': False,'contiguous': False, 'shuffle':False,
        'chunksizes':(1,shape_dset[1], shape_dset[2], shape_dset[3]),
    }
    if outpath==None:
        outfile_path = dset['filename']
    else:
        outfile_path = outpath

    dset = dset.drop_vars('RELCOM')

    dset[dset.varName] = dset[dset.varName].astype(np.float32)
    dset.attrs['source'] = dset.attrs['source'] + ' '.join(nc_files)    
    out_file = dset.to_netcdf(outfile_path,encoding={dset.varName:encoding}, 
                    unlimited_dims=['time'])

