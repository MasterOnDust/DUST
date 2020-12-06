#!/usr/bin/env python

import xarray as xr
import sys
import argparse as ap
import DUST
import os
import time


"""
AUTHOR
======
    Ove Haugvaldstad

    ovehaugv@outlook.com

"""



if __name__ == "__main__":
    parser = ap.ArgumentParser(description="""Multiply FLEXPART emissions sensitivities with modelled dust emissions from FLEXDUST,
        and save the output to a new NETCDF file file.""")
    parser.add_argument('path_flexpart', help='path to flexpart output')
    parser.add_argument('path_flexdust', help='path to flexdust output')
    parser.add_argument('--outpath', '--op', help='path to where output should be stored', default='./out')
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    args = parser.parse_args()
    pathflexpart = args.path_flexpart
    pathflexdust = args.path_flexdust
    outpath = args.outpath
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0


    ds = xr.open_dataset(pathflexpart,chunks={'pointspec':1})

    ds = ds.drop_vars(['WD_spec001', 'DD_spec001'])
    ind_receptor = ds.ind_receptor

    if ind_receptor == 1:
        f_name = 'conc'
        field_unit = 'g/m^3'
        field_name = 'Concentration'
    elif ind_receptor == 4:
        f_name = 'drydep'
        field_unit = 'g/m^2 s'
        field_name = 'Dry depostion'
    elif ind_receptor == 3:
        f_name = 'wetdep'
        field_unit = 'g/m^2 s'
        field_name = 'Wet depostion'
    else:
        f_name = 'spec_mr'
        field_unit = 'g/m^3'
        field_name = 'Unknown'

    # only select surface sensitivity
    ds = ds.isel(height=0)
    height = ds.height.values
    ds = ds.sel(nageclass=0)
    # rename variables
    
    ds = ds.rename({'latitude':'lat','longitude':'lon'})

    # select output domain
    ds = ds.sel(lon=slice(x0,x1), lat=slice(y0,y1))

    # Read flexdust output

    flexdust_ds = DUST.read_flexdust_output(pathflexdust)
    flexdust_ds = flexdust_ds.sel(time=ds.time)
    flexdust_ds = flexdust_ds.sel(lon=slice(70,120), lat=slice(25,50))

    # Interpolate into same coordinates as FLEXPART
    flexdust_ds = flexdust_ds.interp_like(ds)

    # create output DataArray
    out_data = xr.DataArray(coords={'pointspec':ds.pointspec,'lon':ds['spec001_mr'].lon, 
                            'lat':ds['spec001_mr'].lat, 'time':ds.time}, 
                        dims={'pointspec':ds.pointspec,'time':ds.time,
                        'lat':ds.lat,'lon':ds.lon}, attrs=ds.spec001_mr.attrs)
    
    out_data.attrs['units'] = field_unit
    spec_com = ds.spec001_mr.attrs['long_name']
    out_data.attrs['spec_comment'] = spec_com
    out_data.attrs['long_name'] = field_name
    
    
    # Fill output Data arrary
    scale_factor = (1/height)*1000

    for i in range(len(ds.pointspec)):
        out_data[i] = ds.spec001_mr.sel(pointspec=i) * flexdust_ds.Emission * scale_factor
    
    terminal_input = ' '.join(sys.argv)

    out_ds = xr.Dataset({f_name : out_data, 'surface_sensitvity':ds.spec001_mr, 'RELEND':ds.RELEND, 'RELSTART':ds.RELSTART, 
                        'ORO':ds.ORO, 'RELPART':ds.RELPART.sum(), 'RELZZ1':ds.RELZZ1[0],
                        'RELZZ2': ds.RELZZ2[0], 'RELLAT':ds.RELLAT1[0], 'RELLNG':ds.RELLNG1[0]
                        ,'RELCOM':ds.RELCOM.astype('U35', copy=False)}, attrs=ds.attrs)
    ds.close()
    receptor_name = str(out_ds.RELCOM[0].values).strip().split(' ')[1]

    out_ds.attrs['title'] = 'FLEXPART/FLEXDUST model output'
    out_ds.attrs['references'] = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'
    out_ds.attrs['history'] = out_ds.attrs['history'] + '{} processed by {}'.format(time.ctime(time.time()),terminal_input)
    out_ds.attrs['varName'] = f_name
    shape_dset = out_ds[f_name].shape
    encoding = {'zlib':True, 'complevel':9, 'chunksizes' : (1,10, shape_dset[2], shape_dset[3]),
    'fletcher32' : False,'contiguous': False, 'shuffle' : False}
    outFile_name = os.path.join(outpath,f_name + '_' + receptor_name + '_' + spec_com + '_' + out_ds.iedate +'-'+ out_ds.ibdate + '.nc')
    
    out_ds.to_netcdf(outFile_name, encoding={f_name:encoding, 'surface_sensitvity':encoding})