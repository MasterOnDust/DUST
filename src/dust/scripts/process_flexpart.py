#!/usr/bin/env python
import numpy as np
import xarray as xr
import argparse as ap
import dust
import os
from dust.process_data_dust import process_per_pointspec, process_per_timestep, create_output
import pandas as pd
from IPython import embed
import dask
"""
AUTHOR
======
    Ove Haugvaldstad

    ovehaugv@outlook.com

"""

def main():
    parser = ap.ArgumentParser(description="""Multiply FLEXPART emissions sensitivities with modelled dust emissions from FLEXDUST,
        and save the output to a new NETCDF file file.""")
    parser.add_argument('path_flexpart', help='path to flexpart output')
    parser.add_argument('path_flexdust', help='path to flexdust output')
    parser.add_argument('--outpath', '--op', help='path to where output should be stored', default='./out')
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--height', help='height of lowest outgrid height', default=None, type=int)
    parser.add_argument('--use_dask', help='whether to use dask or not', action='store_true',
                        default=False)
    args = parser.parse_args()

    pathflexpart = args.path_flexpart
    pathflexdust = args.path_flexdust
    outpath = args.outpath
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0
    height = args.height
    use_dask = args.use_dask

    flexdust_ds = DUST.read_flexdust_output(pathflexdust)['dset']
    flexdust_ds = flexdust_ds.sel(lon=slice(x0,x1), lat=slice(y0,y1))
    # Check whether output is per time step or per release?

    if use_dask:
        from dask.distributed import Client, LocalCluster
        import dask
        cluster = LocalCluster(n_workers=8,
                       threads_per_worker=1,
                       memory_limit='64GB')
        client = Client(cluster)
        flexdust_ds = flexdust_ds.chunk({'time':124})
        ds = xr.open_dataset(pathflexpart,chunks={'pointspec':3})
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            if 'pointspec' in ds.dims:
                print('per receptor point')
                ds = xr.open_dataset(pathflexpart, chunks={'time':50, 'pointspec':20})
                ds, out_data, surface_sensitivity = process_per_pointspec(ds, flexdust_ds, x0, x1, y0, y1, height=height)
                relcom_str=str(ds.RELCOM[0].values.astype('U35')).strip().split(' ',2)[1:]
                ds.attrs['relcom']=[s.strip() for s in relcom_str]
            else:
                print('per timestep')
                ds = xr.open_dataset(pathflexpart, chunks={'time':50})
                ds, out_data, surface_sensitivity = process_per_timestep(ds, flexdust_ds, x0, x1, y0, y1, height=height)

            out_ds = create_output(out_data,surface_sensitivity,ds)

    else:
        ds = xr.open_dataset(pathflexpart)

        if 'pointspec' in ds.dims:
            print('per receptor point')
            #ds = xr.open_dataset(pathflexpart, chunks={'time':50, 'pointspec':20})
            ds = xr.open_dataset(pathflexpart)
            ds, out_data, surface_sensitivity = process_per_pointspec(ds, flexdust_ds, x0, x1, y0, y1, height=height)
            relcom_str=str(ds.RELCOM[0].values.astype('U35')).strip().split(' ',2)[1:]
            ds.attrs['relcom']=[s.strip() for s in relcom_str]
        else:
            print('per timestep')
            ds = xr.open_dataset(pathflexpart, chunks={'time':50})
            ds, out_data, surface_sensitivity = process_per_timestep(ds, flexdust_ds, x0, x1, y0, y1, height=height)

        out_ds = create_output(out_data,surface_sensitivity,ds)

    flexdust_ds.close()
    ds.close()
    spec_com = ds.spec001_mr.attrs['long_name']
    f_name = out_ds.attrs['varName']
    shape_dset = out_ds[f_name].shape
    encoding = {'zlib':True, 'complevel':9, 'chunksizes' : (1,10, shape_dset[2], shape_dset[3]),
    'fletcher32' : False,'contiguous': False, 'shuffle' : False}
    outFile_name = os.path.join(outpath,out_ds.attrs['filename'])
    print('writing to {}'.format(outFile_name))
    out_ds.to_netcdf(outFile_name, encoding={f_name:encoding, 'surface_sensitivity':encoding})
    if use_dask:
        client.close()
        cluster.close()

if __name__ == "__main__":
    main()
