#!/usr/bin/env python
import numpy as np
import xarray as xr
import sys
import argparse as ap
import DUST
import os
from netCDF4 import date2num, num2date
import pandas as pd
import time
#from dask.distributed import Client, LocalCluster
#import dask
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
    
    print(outpath)
    #cluster=LocalCluster(n_workers=4, threads_per_worker=1, memory_limit='8GB')
    #client=Client(cluster)
    
    ds = xr.open_dataset(pathflexpart,chunks={'pointspec':1})

    ds = ds.drop_vars(['WD_spec001', 'DD_spec001'])
    ind_receptor = ds.ind_receptor

    if ind_receptor == 1:
        f_name = 'conc'
        field_unit = 'g/m^3'
        sens_unit = 's'
        field_name = 'Concentration'
    elif ind_receptor == 4:
        f_name = 'drydep'
        field_unit = 'g/m^2 s'
        sens_unit = 'm'
        field_name = 'Dry depostion'
    elif ind_receptor == 3:
        f_name = 'wetdep'
        field_unit = 'g/m^2 s'
        sens_unit = 'm'
        field_name = 'Wet depostion'
    else:
        f_name = 'spec_mr'
        field_unit = 'g/m^3'
        sens_unit = 's'
        field_name = 'Unknown'

    spec_com = ds.spec001_mr.attrs['long_name']

    # only select surface sensitivity
    ds = ds.isel(height=0)
    height = ds.height.values
    ds = ds.sel(nageclass=0)
    # rename variables
    
    ds = ds.rename({'latitude':'lat','longitude':'lon'})

    # select output domain
    ds = ds.sel(lon=slice(x0,x1), lat=slice(y0,y1))

    # Determine simulation start
    sim_start = pd.to_datetime(ds.iedate+ds.ietime)

    # Determine size of backward time dimmension
    lout_step = abs(ds.attrs['loutstep'])
    btime_size = int(ds['LAGE']/lout_step*1e-9)
    lout_step_h = int(lout_step/(60*60))
    btime_array = -np.arange(lout_step_h,btime_size*lout_step_h+lout_step_h,lout_step_h)

    # Create new time forward time dimmension
    t0 = pd.to_datetime(ds.ibdate+ds.ibtime) + ds['LAGE'].values

    time_a = np.arange(0,len(ds['pointspec'])*3,3)
    time_var = xr.Variable('time',time_a, 
                    attrs=dict(units='hours since {}'.format(t0.strftime('%Y-%m-%d %H:%S')),
                    calendar="proleptic_gregorian"))

    # Read flexdust output
    flexdust_ds = DUST.read_flexdust_output(pathflexdust)['dset']
    flexdust_ds = flexdust_ds.sel(lon=slice(x0,x1), lat=slice(y0,y1))

    # Interpolate into same coordinates as FLEXPART
    flexdust_ds = flexdust_ds.interp_like(ds)

    # create output DataArray
    print('creating output array')
    out_data = xr.DataArray(np.zeros((len(ds['pointspec']),btime_size,len(ds['lat']),len(ds['lon'])),dtype=np.float32),
        dims=['time', 'btime', 'lat', 'lon'],
        coords={'time':time_var,
        'btime': ('btime',btime_array, dict(
            long_name='time along back trajectory',
            units='hours')), 
        'lon':('lon',ds['spec001_mr'].lon, ds.lon.attrs), 
        'lat':('lat',ds['spec001_mr'].lat, ds.lat.attrs)},
        attrs=dict(
            units=field_unit,
            spec_com=ds.spec001_mr.attrs['long_name'],
            long_name=field_name
        ) 
                        )
    surface_sensitvity = out_data.copy()
    surface_sensitvity.attrs['units'] = sens_unit
    print('created output array')

    
    
    # Fill output Data arrary
    scale_factor = (1/height)*1000
    #result=[]
    last_btime = out_data.btime[-1]
    first_btime =out_data.btime[0]
    time_units = out_data.time.units
    for i in range(len(out_data.time)):
        date0 = num2date(out_data[i].time + first_btime, time_units).strftime('%Y%m%d%H%M')
        date1 = num2date(out_data[i].time + last_btime, time_units).strftime('%Y%m%d%H%M')
        temp_data = ds['spec001_mr'].sel(time=slice(date0, date1), pointspec=i)
        emission_field = flexdust_ds['Emission'].sel(time=temp_data.time)
        #print(temp_data.time[0], emission_field.time[0])
        out_data[i] = temp_data.values*emission_field.values*scale_factor
        surface_sensitvity[i] = temp_data.values

    print('finish emsfield*sensitvity')

    terminal_input = ' '.join(sys.argv)
    
    out_ds = xr.Dataset({f_name : out_data, 'surface_sensitivity' : surface_sensitvity ,'RELEND':ds.RELEND, 'RELSTART':ds.RELSTART, 
                        'ORO':ds.ORO, 'RELPART':ds.RELPART.sum(), 'RELZZ1':ds.RELZZ1[0],
                        'RELZZ2': ds.RELZZ2[0], 'RELLAT':ds.RELLAT1[0], 'RELLNG':ds.RELLNG1[0]
                        ,'RELCOM':ds['RELCOM'].astype('U35', copy=False)}, attrs=ds.attrs)

    flexdust_ds.close()
    ds.close()
    receptor_name = str(out_ds.RELCOM[0].values).strip().split(' ')[1]

    out_ds.attrs['title'] = 'FLEXPART/FLEXDUST model output'
    out_ds.attrs['references'] = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'
    out_ds.attrs['history'] = '{} processed by {}, '.format(time.ctime(time.time()),terminal_input) + out_ds.attrs['history']
    out_ds.attrs['varName'] = f_name
    shape_dset = out_ds[f_name].shape
    encoding = {'zlib':True, 'complevel':9, 'chunksizes' : (1,10, shape_dset[2], shape_dset[3]),
    'fletcher32' : False,'contiguous': False, 'shuffle' : False}
    outFile_name = os.path.join(outpath,f_name + '_' + receptor_name + '_' + spec_com + '_' + out_ds.ibdate +'-'+ out_ds.iedate + '.nc')
    print('writing to {}'.format(outFile_name))
    out_ds.to_netcdf(outFile_name, encoding={f_name:encoding, 'surface_sensitivity':encoding})
