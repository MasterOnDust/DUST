#!/usr/bin/env python

import xarray as xr
import argparse as ap
import pandas as pd
import time
import sys
import os 
from DUST.utils.resample import resample_monthly, concatenate_monthly




if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Create monthly mean source contribution')
    parser.add_argument('paths_flexpart',nargs ='+', help='paths to flexpart output')
    parser.add_argument('--outpath', '--op',default='./', help='path to where output should be stored')
    parser.add_argument('--filename', '--fn', 
                    help='Name of output files, default determine automatically', default=None)
    
    args = parser.parse_args()

    paths = args.paths_flexpart
    outpath = args.outpath
    file_name = args.filename 
    paths.sort()
    dsets_list=[resample_monthly(xr.open_dataset(path)) for path in paths]

    terminal_input = ' '.join(sys.argv)
    dsets=concatenate_monthly(dsets_list)


    dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),terminal_input) + dsets.attrs['history']
    if file_name:
        outFile_name = os.path.join(outpath, file_name)
    else:
        outFile_name = os.path.join(outpath, dsets.attrs['filename'])
    print('writing to {}'.format(outFile_name))
    dsets.to_netcdf(outFile_name)


