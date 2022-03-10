import numpy as np
import xarray as xr
import sys
from netCDF4 import date2num, num2date
import pandas as pd
import time

"""
AUTHOR Ove Haugvaldstad

This file containts functions for processing data for working with FLEXDUST and FLEXPART model output

"""


def process_per_pointspec(dset,flexdust_ds, x0,x1,y0,y1, height=None):
    """
    DESCRIPTION:
    ===========
        Combines FLEXPART backward simulation output and with FLEXDUST emission field. 
        For FLEXPART output stored in a per pointspec format. Meaning that each FLEXPART 
        simulation has to only contain on recepotor point and the forward time dimension
        is defined by specifying a new release. 

    ARGUMENTS:
    ==========
        dset: xarray.dataset
            dataset containing FLEXPART model output
        flexdust_ds: xarray.dataset
            dataset containing FLEXDUST emission fields.
        x0: float
            longitude of lower left corner of cutout.
        x1: float
            longtude of upper right corner of cutout.
        y0: float
            latitude of lower left corner of cutout.
        y1: float 
            latidue of upper left corner of cutout. 
    KEYWORD ARGUMENTS:
    ==================
        height: int (default=None)
            Which height in the FLEXPART output to multiply with the FLEXDUST
            emissions fields. If None then the lowerst layer output layer is 
            selected. 
    
    RETURN:
    =======
        dset: xarray.Dataset
            Cleaned up input dataset.
        out_data: xarray.DataArray
            DataArray contatining FLEXPART and FLEXDUST multiplied.
        surface_sensitivity: xarray.DataArray
            DataArray conataining the senstivity of the lower FLEXPART 
            output layer.
            
    """

    
    dset = dset.drop_vars(['WD_spec001', 'DD_spec001'])
    dset = dset.chunk({'pointspec':1})
    # only select surface sensitivity
    if height == None:
        dset = dset.isel(height=0)
    else:
        dset = dset.sel(height=height)
    height = dset.height.values
    dset = dset.sel(nageclass=0)
    # rename variables
    
    dset = dset.rename({'latitude':'lat','longitude':'lon'})

    # select output domain
    dset = dset.sel(lon=slice(x0,x1), lat=slice(y0,y1))

    #interpolate flexdust to match flexpart coordinates
    flexdust_ds = flexdust_ds.interp({'lon':dset.lon,'lat':dset.lat})    


    # Assumes that the first part of the RELCOM string contains the date. 
    if isinstance(dset.RELCOM.values[0],np.bytes_):
        dset = dset.assign(RELCOM=dset.RELCOM.astype(str))

    # Determine size of backward time dimmension
    lout_step = abs(dset.attrs['loutstep'])
    btime_size = int(dset['LAGE']/lout_step*1e-9)
    lout_step_h = int(lout_step/(60*60))
    btime_array = -np.arange(lout_step_h,btime_size*lout_step_h+lout_step_h,lout_step_h)
    # Create new time forward time dimmension
    t0 = pd.to_datetime(dset.ibdate+dset.ibtime) + pd.to_timedelta(dset['LAGE'].values, unit='ns')
    # Assumes that the first part of the RELCOM string contains the date. 
    time_a = pd.to_datetime([date.split(' ')[0] for date in dset.RELCOM.values],format='%Y%m%d%H').to_pydatetime()
    time_a = date2num(time_a,units='hours since {}'.format(t0.strftime('%Y-%m-%d %H:%S')))
    time_var = xr.Variable('time',time_a, 
                    attrs=dict(units='hours since {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S')),
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
        'lon':('lon',dset['spec001_mr'].lon.data, dset.lon.attrs), 
        'lat':('lat',dset['spec001_mr'].lat.data, dset.lat.attrs)},
        attrs=dict(
            spec_com=dset.spec001_mr.attrs['long_name'],
        ) 
                        )

    surface_sensitivity = out_data.copy()
    if dset.ind_receptor == 3 or dset.ind_receptor == 4:
        scale_factor = (1/(height))*1000 # Deposition is accumulative 
    else:
        # Concentration is not  accumulative Units of FLEXDUST need to be g/m^3s
        scale_factor = (1/(height*lout_step))*1000 
    # print(scale_factor)
    last_btime = out_data.btime[-1].values
    first_btime =out_data.btime[0].values
    time_units = out_data.time.units
    
    dset = dset.sortby('time')

    for i in range(len(out_data.time)):
        date0 = out_data[i].time.values + first_btime
        date1 = out_data[i].time.values + last_btime
        date0 = num2date(date0, time_units).strftime('%Y%m%d%H%M%S')
        date1 = num2date(date1, time_units).strftime('%Y%m%d%H%M%S')
        temp_data = dset['spec001_mr'].sel(time=slice(date1, date0), pointspec=i)
        
        emission_field = flexdust_ds['Emission'].sel(time=temp_data.time)
        
        da = (temp_data*emission_field)*scale_factor
        da = da.rename(time='btime')
        da = da[::-1].assign_coords(btime=out_data.btime)
        out_data[i] = da
        surf_sens_da = temp_data.rename(time='btime')
        surf_sens_da = surf_sens_da[::-1].assign_coords(btime=out_data.btime)
        surface_sensitivity[i] = surf_sens_da    
    
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
    """
    DESCRIPTION:
    ============
        Setup the output dataset format
    
    ARGUMENTS:
    ==========
        out_data: xarray.DataArray
            DataArray containing the combined FLEXPART and FLEXDUST output
        surface_sensitivity: xarray.DataArray
            DataArray contatining the FLEXPART emission sensitivity of the lower most 
            output layer
        dset: xarray.Dataset
            Output dataset
    RETURN:
    =======

        output_ds:xarray.Dataset
            Complete output dataset with CF convensions applied. 

    """

    ind_receptor = dset.ind_receptor
    f_name, field_unit, sens_unit, field_name = determine_units(ind_receptor)



    terminal_input = ' '.join(sys.argv)
    
    out_ds = xr.Dataset({f_name : out_data, 'surface_sensitivity' : surface_sensitivity ,'RELEND':dset.RELEND, 'RELSTART':dset.RELSTART, 
                        'RELPART':dset.RELPART, 'RELZZ1':dset.RELZZ1,
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
    """ 
    Determine units of FLEXPART output
    """

    if ind_receptor == 1:
        f_name = 'conc'
        field_unit = 'g/m^3'
        sens_unit = 's'
        field_name = 'Concentration'
    elif ind_receptor == 4:
        f_name = 'drydep'
        field_unit = 'g/m^2'
        sens_unit = 'm'
        field_name = 'Dry depostion'
    elif ind_receptor == 3:
        f_name = 'wetdep'
        field_unit = 'g/m^2'
        sens_unit = 'm'
        field_name = 'Wet depostion'
    else:
        f_name = 'spec_mr'
        field_unit = 'g/m^3'
        sens_unit = 's'
        field_name = 'Unknown'
    return f_name, field_unit, sens_unit, field_name
