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
from IPython import embed
#from dask.distributed import Client, LocalCluster
#import dask
"""
AUTHOR
======
    Ove Haugvaldstad

    ovehaugv@outlook.com

"""
def process_per_pointspec(dset,flexdust_ds, x0,x1,y0,y1, height=None):
    dset = dset.drop_vars(['WD_spec001', 'DD_spec001'])
    dset = dset.chunk({'pointspec':1})
    # only select surface sensitivity
    if height == None:
        dset = dset.isel(height=0)
    else:
        height = dset.sel(height=height)
    height = dset.height.values
    dset = dset.sel(nageclass=0)
    # rename variables
    
    dset = dset.rename({'latitude':'lat','longitude':'lon'})

    # select output domain
    dset = dset.sel(lon=slice(x0,x1), lat=slice(y0,y1))

    #interpolate flexdust to match flexpart coordinates
    flexdust_ds = flexdust_ds.interp({'lon':dset.lon,'lat':dset.lat})    

    # Determine simulation start
    sim_start = pd.to_datetime(dset.iedate+dset.ietime)

    # Determine size of backward time dimmension
    lout_step = abs(dset.attrs['loutstep'])
    btime_size = int(dset['LAGE']/lout_step*1e-9)
    lout_step_h = int(lout_step/(60*60))
    btime_array = -np.arange(lout_step_h,btime_size*lout_step_h+lout_step_h,lout_step_h)
    # Create new time forward time dimmension
    t0 = pd.to_datetime(dset.ibdate+dset.ibtime) + dset['LAGE'].values

    time_a = np.arange(0,len(dset['pointspec'])*3,3)
    time_var = xr.Variable('time',time_a, 
                    attrs=dict(units='hours since {}'.format(t0.strftime('%Y-%m-%d %H:%S')),
                    calendar="proleptic_gregorian"))
    
    # create output DataArray
    print('creating output array')
    out_data = xr.DataArray(np.zeros((len(dset['pointspec']),btime_size,len(dset['lat']),len(dset['lon'])),
            dtype=np.float32),
        dims=['time', 'btime', 'lat', 'lon'],
        coords={'time':time_var,
        'btime': ('btime',btime_array, dict(
            long_name='time along back trajectory',
            units='hours')), 
        'lon':('lon',dset['spec001_mr'].lon, dset.lon.attrs), 
        'lat':('lat',dset['spec001_mr'].lat, dset.lat.attrs)},
        attrs=dict(
            spec_com=dset.spec001_mr.attrs['long_name'],
        ) 
                        )

    surface_sensitivity = out_data.copy()
    scale_factor = (1/height)*1000
    #result=[]
    last_btime = out_data.btime[-1]
    first_btime =out_data.btime[0]
    time_units = out_data.time.units

    for i in range(len(out_data.time)):
        date0 = num2date(out_data[i].time + first_btime, time_units).strftime('%Y%m%d%H%M')
        date1 = num2date(out_data[i].time + last_btime, time_units).strftime('%Y%m%d%H%M')
        temp_data = dset['spec001_mr'].sel(time=slice(date0, date1), pointspec=i)
        emission_field = flexdust_ds['Emission'].sel(time=temp_data.time)
        out_data[i] = temp_data.values*emission_field.values*scale_factor
        surface_sensitivity[i] = temp_data.values    
    
    print('finish emsfield*sensitvity')
    dset = dset.assign({
        'RELLAT1' : dset['RELLAT1'][0],
        'RELLNG1' : dset['RELLNG1'][0],
        'RELZ1' : dset['RELZZ1'][0],
        'RELZ2' : dset['RELZZ2'][0],
        'RELPART' : dset['RELPART'].sum(keep_attrs=True)
        })
    
    dset.attrs['ibdate'] = t0.strftime('%Y%m%d')
    dset.attrs['ibtime'] = t0.strftime('%H%M%S')
    receptor_name = str(dset.RELCOM[0].values).strip().split(' ')[1:]
    dset.attrs['relcom'] = receptor_name
    return dset,out_data, surface_sensitivity

def process_per_timestep(dset, flexdust_ds,x0,x1,y0,y1, height=None):
    dset = dset.sel(lon=slice(x0,x1), lat=slice(y0,y1))
    
    #interpolate flexdust to match flexpart coordinates
    flexdust_ds = flexdust_ds.interp({'lon':dset.lon,'lat':dset.lat})    
    
    if height == None:
        
        height = dset.height.values
    else:
        height = height
    scale_factor = (1/height)*1000
    print('creating output array')
    out_data = xr.zeros_like(dset['spec001_mr'])
    for i in range(len(out_data.time)):
        
        temp_data = dset['spec001_mr'].isel(time=i)
        time_steps = temp_data.time + temp_data.btime 
        emission_field = flexdust_ds['Emission'].sel(time=time_steps)
        out_data[i] = temp_data.values*emission_field.values*scale_factor
    print('finish emsfield*sensitvity')
    surface_sensitivity = dset['spec001_mr']
    iedate_stamp = dset.time[-1]
    ibdate_stamp = dset.time[0]
    dset.attrs['iedate'] = str(iedate_stamp.dt.strftime('%Y%m%d').values)
    dset.attrs['ietime'] = str(iedate_stamp.dt.strftime('%H%M%S').values)

    dset.attrs['ibdate'] = str(ibdate_stamp.dt.strftime('%Y%m%d').values)
    dset.attrs['ibtime'] = str(ibdate_stamp.dt.strftime('%H%M%S').values)
    return dset, out_data, surface_sensitivity


def create_output(out_data, surface_sensitivity, dset):
    ind_receptor = dset.ind_receptor
    f_name, field_unit, sens_unit, field_name = determine_units(ind_receptor)



    terminal_input = ' '.join(sys.argv)
    
    out_ds = xr.Dataset({f_name : out_data, 'surface_sensitivity' : surface_sensitivity ,'RELEND':dset.RELEND, 'RELSTART':dset.RELSTART, 
                        'RELPART':dset.RELPART, 'RELZZ1':ds.RELZZ1,
                        'RELZZ2': dset.RELZZ2, 'RELLAT':dset.RELLAT1, 'RELLNG':dset.RELLNG1
                        }, attrs=dset.attrs)

    out_ds[f_name].attrs['units'] = field_unit
    out_ds[f_name].attrs['long_name'] = field_name
    out_ds['surface_sensitivity'].attrs['units'] = sens_unit

    out_ds.attrs['title'] = 'FLEXPART/FLEXDUST model output'
    out_ds.attrs['references'] = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'
    out_ds.attrs['history'] = '{} processed by {}, '.format(time.ctime(time.time()),terminal_input) + out_ds.attrs['history']
    out_ds.attrs['varName'] = f_name
    receptor_name = out_ds.attrs['relcom'][0]
    spec_com = out_ds.attrs['relcom'][1]
    file_name = f_name + '_' + receptor_name + '_' + spec_com + '_' + out_ds.ibdate +'-'+ out_ds.iedate + '.nc'
 
    out_ds.attrs['filename'] = file_name

    
    
    return out_ds
    


def determine_units(ind_receptor):
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
    return f_name, field_unit, sens_unit, field_name

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
    parser.add_argument('--height', help='height of lowest outgrid height', default=None, type=int)
    args = parser.parse_args()

    pathflexpart = args.path_flexpart
    pathflexdust = args.path_flexdust
    outpath = args.outpath
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0
    height = args.height
    
    flexdust_ds = DUST.read_flexdust_output(pathflexdust)['dset']
    flexdust_ds = flexdust_ds.sel(lon=slice(x0,x1), lat=slice(y0,y1))
    # Check whether output is per time step or per release?
    ds = xr.open_dataset(pathflexpart)
    if 'pointspec' in ds.dims:
        print('per receptor point')
        ds, out_data, surface_sensitivity = process_per_pointspec(ds, flexdust_ds, x0, x1, y0, y1, height=height)
        ds.attrs['relcom'] = str(ds.RELCOM[0].values.astype('U35')).strip().split(' ')[1:]

    else:
        print('per timestep')
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
 
