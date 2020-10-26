#!/usr/bin/env python
import os 
from DUST.read_data import read_multiple_flexpart_outputs
from DUST.utils.utils import arg_parser, region_slice
import pandas as pd
import shutil

if __name__ == "__main__":
    parser = arg_parser('Concatenate FLEXPART simulations')
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--heights', help='which levels to concatenate output for', nargs='+', type=int, default=None)
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


    # Load the netCDF files
    dset = read_multiple_flexpart_outputs(nc_files, chunks={'pointspec':1})
    dset = region_slice(dset, x0, x1, y0, y1)
    spatial_attrs_names = ['outlon0', 'outlon1', 'outlat0', 'outlat1']

    # Update attribute if output grid has been sliced
    for key, outloc in zip(spatial_attrs_names, [x0,x1, y0,y1]):
        if outloc == None:
            dset.attrs[key] = outloc
    if heights == None:
        heights = dset.height.values
    if locations == 'ALL':
        locations = dset.point.values
    
    date0 = df.index[0]
    date1 = df.index[-1]
    outdir = ''.join(df['ncfiles'][0].split('_')[:-1]) + '_' + date0 + '_' + date1
    path_to_dir = os.path.join(outpath, outdir)
    try: 
        os.mkdir(path_to_dir)
    except:
        askConfirmation = input("""This folder already exists,
        do you want to delete it? (y/n)""")
        if askConfirmation.strip() == 'y':
            shutil.rmtree(path_to_dir)
            os.mkdir(path_to_dir)

    for location in locations:
        for height in heights:
            temp_dset = dset.sel(point=location, height=height)
            temp_dset = temp_dset.squeeze()
            # print(temp_dset.chunk())
            temp_dset = temp_dset.chunk({'time':10})
            # print(temp_dset.chunks)
            Loc_name = "_".join(location.split(' '))
            file_name = "{}_{}_{}_{}.nc".format(Loc_name, height, date0, date1) 
            outfile_path = os.path.join(path_to_dir, file_name)
            temp_dset['point'] = temp_dset.point.astype('S{}'.format(len(location)))
            temp_dset['RELCOM'] = temp_dset['RELCOM'].astype('S{}'.format(len(location)))
            temp_dset = temp_dset.compute()
            
            temp_dset.to_netcdf(outfile_path)
            temp_dset.close()

