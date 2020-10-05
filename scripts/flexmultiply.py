#!/usr/bin/env python

import argparse as ap
from DUST.utils.multiply_emsfield import multi_flexpart_flexdust
import os
import glob
import xarray as xr
import pandas as pd

"""
AUTHOR
======
    Ove Haugvaldstad

    ovehaugv@outlook.com

"""
parser_description =("""
Multiply FLEXPART emissions sensitivities with modelled dust emissions from FLEXDUST,
and save the output to a new NETCDF file file.
""")

defaultOutPath = './out'
defaultLocations = 'ALL'

helpFlexpartPath = "Path to FLEXPART output top directory"
helpFlexdustPath = "Path to file/folder containing FLEXDUST emission fields"
helpLocation = "Name or integer correspoding to location defined in FLEXPART RELEASE file"
helpOutPath = "Path to directory where multiplied field should be stored"
helpZlib = "Whether to use zlib compression"

if __name__ == "__main__":

    parser = ap.ArgumentParser(description=parser_description)
    parser.add_argument("fpPath", help=helpFlexpartPath)
    parser.add_argument("fdPath", help=helpFlexdustPath)
    parser.add_argument("--zlib", "--zl", help=helpZlib, action='store_true')
    parser.add_argument("--locations", "--locs",default=defaultLocations, help=helpLocation, nargs='+')
    parser.add_argument("--out", "--o", default=defaultOutPath, help=helpOutPath)
    parser.add_argument("--nFiles", "--nF", default = 'ALL', help ='How many flexpart files to mulitply')
    args = parser.parse_args()
    path = args.fpPath
    flexdust = args.fdPath
    locations = args.locations
    outpath = args.out
    zlib = args.zlib
    nFiles = args.nFiles

    if os.path.isdir(outpath):
        pass
    else:
        os.mkdir(outpath)
    
    try:
        df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
        ncFiles = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
    except FileNotFoundError:
        ncFiles = glob.glob(path + "**/output/grid*.nc", recursive=True) #recursively find FLEXPART output files

    if nFiles == 'ALL':
        pass
    else:
        ncFiles = glob.glob(path + "/**/grid*.nc", recursive=True)[:int(nFiles)]

    ncFiles.sort()
    d = xr.open_dataset(ncFiles[0])
    relCOMS = d.RELCOM
    d.close()


    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()
        if locations == 'ALL':
            multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib)
        else:
            for receptor in locations:
                if receptor in loc or receptor == str(i):
                    multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib)
                else:
                    continue


