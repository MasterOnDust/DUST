#!/usr/bin/env python
import os
from IPython import embed
from dask.distributed import Client, LocalCluster 
from DUST.read_data import read_multiple_flexpart_outputs
from DUST.utils.utils import arg_parser
import pandas as pd
import shutil
import numpy as np
import sys
import argparse as ap

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Concatenate FLEXPART simulations')
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--heights', help='which levels to concatenate output for', nargs='+', type=int, default=None)
    parser.add_argument('path', help='Path to top directory containing output')
    parser.add_argument('--outpath', '--op', help='Where the output should be stored', default='.')
    parser.add_argument('--locations', '--loc', help='Number of location to apply fuction to', default=None, nargs='+', type=int)
    parser.add_argument('--sdate', '--sd', help='Begining of time slice')
    parser.add_argument('--edate', '--ed', help='End of time slice')
    args = parser.parse_args()
    locations = args.locations
    path = args.path
    outpath = args.outpath
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0
    sdate = args.sdate
    edate = args.edate
    heights = args.heights
    try:
        df = pd.read_csv(os.path.join(path,'AVAILABLE_OUTPUT'), index_col=0).loc[slice(sdate, edate)]
        nc_files = [os.path.join(path,row['dir_paths'] + '/'+ row['ncfiles']) for index,row in df.iterrows()]

    except FileNotFoundError:
        raise FileNotFoundError('AVAIALBE output file not found, create AVAILABLE with list_available_output.py')
    date0 = df.index[0]
    date1 = df.index[-1]
    outdir = ''.join(df['ncfiles'][0].split('_')[:-1]) + '_' + date0 + '_' + date1
    path_to_dir = os.path.join(outpath, outdir)

    try: 
        os.mkdir(path_to_dir)
    except FileExistsError:
        askConfirmation = input("""This folder already exists,
        do you want to delete it? (y/n)""")
        if askConfirmation.strip() == 'y':
            shutil.rmtree(path_to_dir)
            os.mkdir(path_to_dir)
        else:
            sys.exit()

    cluster = LocalCluster(n_workers=8,threads_per_worker=1 ,memory_limit='32GB')
    client = Client(cluster)
    # Load the netCDF files
    dset = read_multiple_flexpart_outputs(nc_files, height=heights, location=locations)
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
    
    if 'height' not in dset.dims:
        heights = ['total']
        sel_height=False
    else:
        if dset.height.values.size == 1:
            heights = [int(dset.height.values)]
            sel_height=False
        else:
            heights = dset.height.values
            sel_height=True
    
    if dset.pointspec.values.size ==1:
        pointspecs = [int(dset.pointspec.values)]
        sel_point=False
    else:
        pointspecs = dset.pointspec.values
        sel_point=True
    for point in pointspecs:
        for h in heights:

            if sel_point==False and sel_height:
                temp_dset= dset.sel(height=h)
                h=int(h)
            elif sel_height==False and sel_point:
                temp_dset=dset.sel(numpoint=point,pointspec=point)
            elif sel_height and sel_point:
                temp_dset = dset.sel(numpoint=point, pointspec=point, height=h)
                h = int(h)
            else:
                temp_dset = dset
            temp_dset = temp_dset.squeeze()
            Loc_name = str(temp_dset.RELCOM.values).strip().split(' ')

            file_name = "{}_{}_{}_{}.nc".format('_'.join(Loc_name),h, date0, date1) 
            outfile_path = os.path.join(path_to_dir, file_name)
            temp_dset = temp_dset.drop_vars('RELCOM')

            temp_dset.attrs['relcom'] = ' '.join(Loc_name)
            temp_dset[temp_dset.varName] = temp_dset[temp_dset.varName].astype(np.float32)
            temp_dset.attrs['source'] = temp_dset.attrs['source'] + ' '.join(nc_files)    
            shape_dset = temp_dset[temp_dset.varName].shape
            encoding = {'zlib': True, 'complevel':9,
                'fletcher32': False,'contiguous': False, 'shuffle':False,
                'chunksizes':(1,shape_dset[1], shape_dset[2], shape_dset[3]),
            }
            #embed()
            #print(temp_dset)
            temp_dset.to_netcdf(outfile_path,encoding={temp_dset.varName:encoding}, 
                    unlimited_dims=['time'])
    #cluster.close()
    #client.close()

    #embed()

