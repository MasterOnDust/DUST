import argparse as ap
from DUST.utils.multiply_emsfield import multi_flexpart_flexdust
import os
import glob
import xarray as xr
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

    args = parser.parse_args()
    path = args.fpPath
    flexdust = args.fdPath
    locations = args.locations
    outpath = args.out
    zlib = args.zlib


    if os.path.isdir(outpath):
        pass
    else:
        os.mkdir(outpath)
    ncFiles = glob.glob(path + "/**/grid*.nc", recursive=True) #recursively find FLEXPART output files
    
    d = xr.open_dataset(ncFiles[0])
    relCOMS = d.RELCOM
    d.close()
    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()[0]
        if locations == 'ALL':
            if loc or i in locations:
                multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib)    
        else:
            multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib)


