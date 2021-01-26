#!/usr/bin/env python

import argparse as ap
import xarray as xr
import pandas as pd
import os 
from IPython import embed


def get_sensitvity_at_point(paths, fname,df,box_size_receptor=0, 
                    groupby_key='Mountain', outpath='./out', tag=None):
    """
    Need to figure out how it should store the data
    
    """
    
    out_df = pd.DataFrame(columns=df.index, index=range(1,13))
        
    for path in paths:

        ds = xr.open_dataset(path)
        if tag == None:
            out_tag=ds.height.values
        else:
            out_tag = tag
        groupby_loc = df.groupby('Mountain')
        month = int(ds.time.dt.month[0])
        #print(path,month)
        for index, location in groupby_loc:
            if index==fname:
                ems_sens = ds['spec001_mr'].sel(lon=slice(float(location['lon'].values-box_size_receptor),
                                       float(location['lon'].values+box_size_receptor)), 
                    lat=slice(float(location['lat'].values-box_size_receptor),
                              float(location['lat'].values+box_size_receptor))
                    ).sum(dim='btime').mean(dim=['lon','lat']).mean(dim='time').values
           
            else:
                ems_sens = ds['spec001_mr'].sel(lon=location['lon'], 
                    lat=location['lat'], method='nearest').sum(dim='btime').mean(dim='time').values
                
            out_df.loc[month,index] = float(ems_sens)
        full_outpath = os.path.join(outpath,'{}_{}_sensitivty_matrix.csv'.format(fname,out_tag))
        ds.close()
    out_df.to_csv(full_outpath)
            


if __name__=="__main__":
    parser = ap.ArgumentParser(description='Plot monthly mean maps')
    parser.add_argument('paths', nargs='+', help='Paths to merged output folder')
    parser.add_argument('--mountain_file', '--mf', help='path to where montain file is stored'
            ,default='~/DUST_scripts/scripts/Seeds_Africa/mountainlon_lat.csv')
    parser.add_argument('--box_size', '--bs', default=0.3, 
                            help='box size around the receptor location', type=float)
    parser.add_argument('--out_path', '--op', 
                        help='path to where the sensitivty matrix should be stored', 
                        default='./out')
    parser.add_argument('--tag', help='file name tag', default=None)
    args = parser.parse_args()
    paths = args.paths
    bs = args.box_size
    mf_path = args.mountain_file
    outpath=args.out_path
    tag = args.tag
        #paths = glob.glob(path, recursive=True)
    fname = paths[0].split('/')[-1].split('_')[0]
    paths.sort()
    df = pd.read_csv(mf_path, index_col=0, usecols=['Mountain','lon','lat'])
    print(fname)
    get_sensitvity_at_point(paths, fname, df, box_size_receptor=bs, outpath=outpath, tag=tag)
