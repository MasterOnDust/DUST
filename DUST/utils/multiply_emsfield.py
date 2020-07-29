import glob
import xarray as xr
import pandas as pd
import DUST
from dask.delayed import delayed
import dask
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num, num2date
import os
import time
from IPython import embed
from numba import njit, prange



def multi_flexpart_flexdust(path, nc_files, flexdust_path, point_spec):
    """
    Description:
    ===========
        
        Multiply flexpart emissions sensitivities with modelled dust emissions from FLEXDUST, 
        and save the output to a new netcdf file.

    USAGE:
    =====
        
        path : path to where the netcdf file should be created
        nc_files : list containing the path to the FLEXPART netcdf files 
        flexdust_path : path to the directory containing the flexdust netcdf output
        point_spec : integer --> corrensponding to the release point in the flexpart netcdf file

        multi_flexpart_flexdust(path, nc_files, flexdust, point_spec)

    AUTHOR:
    =======
        Ove Haugvaldstad
        ovehaugv@outlook.com

"""

    # Read in the first flexpart netcdf file as reference for setting up the output file
    d0 = xr.open_dataset(ncfiles[0], decode_times=False)
    d0 = sel_location(d0,point_spec)
       
    try:
        ncfile = Dataset(path, 'w', format="NETCDF4")
    except PermissionError:
        # If netcdf file exist delete old one
        os.remove(path)
        ncfile = Dataset(path, 'w', format='NETCDF4') 
    
    #Setup attributes
    ncfile.title = 'Flexpart - Flexdust SSR'
    ncfile.history = "Created " + time.ctime(time.time())
    ncfile.flexpart_v = d0.source
    ncfile.receptor_name = str(d0.RELCOM.values)[2:].strip()
    ncfile.reference = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'

    #Copy variables 
    lats = d0.latitude.values
    lons = d0.longitude.values
    ts = d0.time.values
    ind_receptor = d0.ind_receptor
    dims = d0.dims
    

    #Spatial dims
    
    lat_dim = ncfile.createDimension('lat', dims['latitude'])
    lon_dim = ncfile.createDimension('lon', dims['longitude'])
    point_dim = ncfile.createDimension('npoint',1)
    #temporal dims
    
    time_dim = ncfile.createDimension('time', None)
    btime_dim = ncfile.createDimension('btime', dims['time'])
    
    
    #Setup lon/lat (spatial variables)
    
    lat = ncfile.createVariable('lat', 'f4', ('lat', ))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    
    lon = ncfile.createVariable('lon', 'f4', ('lon', ))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    
    #Set receptor location
    rellat = ncfile.createVariable('RELLAT', 'f4', ('npoint',))
    rellat.units = 'degrees_north'
    rellat.long_name = 'latitude_receptor'
    
    rellon = ncfile.createVariable('RELLON', 'f4', ('npoint',))
    rellon.units = 'degrees_east'
    rellon.long_name = 'longitude_receptor'
    
    lon[:] = lons
    lat[:] = lats
    
    rellat[:] = d0.RELLAT1.values
    rellon[:] = d0.RELLNG1.values
    
    
    # Setup temporal variable 
    btime = ncfile.createVariable('btime', 'i4', ('btime',))
    btime.units = 's'
    btime.long_name = 'seconds_since_release'
    
    btime[:] = d0.time.values
    
    time_var = ncfile.createVariable('time', 'f8', ('time',))
    time_var.units = "hours since 1980-01-01"
    time_var.long_name = 'time'
    
    # Determind which kind of output should be created
    if ind_receptor == 1:
        field = ncfile.createVariable('Concentration', 'f4',('time', 'btime', 'lat', 'lon'))
        field.units = 'kg/m^3'
        field.long_name = 'Concentration'
    elif ind_receptor == 4:
        field = ncfile.createVariable('DryDep', 'f4', ('time', 'btime', 'lat', 'lon'))
        field.units = 'kg/m^2s'
        field.long_name = 'Dry depostion'
    elif ind_receptor == 3:
        field = ncfile.createVariable('WetDep', 'f4', ('time', 'btime', 'lat', 'lon'))
        field.units = 'kg/m^2s '
        field.long_name = 'Wet depostion'
    else:
        field = ncfile.createVariable('Spec_mr', 'f4', ('time', 'btime', 'lat', 'lon'))
        field.units = 'kg/m^3'
    
    
    field.height = d0.height.values     #set the height of the lowest model output layer
    td = pd.to_timedelta(d0.time.values, unit='s')
    emsField = DUST.read_flexdust_output(flexdust_path)['dset'].Emission
    time_steps = []
    for n, nc_file in enumerate(nc_files):
        temp_ds = xr.open_dataset(nc_file, decode_times=False)
        temp_ds = _sel_location(temp_ds, point_spec)
        s_time = pd.to_datetime(temp_ds.iedate + temp_ds.ietime)
        time_steps.append(date2num(s_time, time_var.units))
        ems_sens = temp_ds.spec001_mr
        temp_ems = emsField.sel(time=s_time + td)
        #multiply emission and emission sensitivity
        temp_array = _multiply_emissions(temp_ems.values, ems_sens.values)  
        field[n,:,:,:] = temp_array
    time_var[:] = time_steps
    ncfile.close()

@njit(parallel=True)
def _multiply_emissions(ems, ems_sens):
    temp_array = np.zeros_like(ems_sens)
    for i in prange(ems.shape[0]):
        for j in range(ems.shape[1]):
            for k in range(ems.shape[2]):
                temp_array[i,j,k] = ems[i,j,k]*ems_sens[i,j,k]
    return temp_array


def _sel_location(ds,pointspec):
    ds = ds.sel(pointspec=pointspec)
    ds = ds.sel(numpoint=pointspec)
    ds = ds.sel(numspec=pointspec)
    ds = ds.sel(nageclass=0)
    ds = ds.isel(height=0)
    return ds


if __name__=="__main__":

    path = '/mnt/c/Users/oveha/Documents/FLEXPART_output'

    ncfiles = glob.glob(path + "/**/*.nc", recursive=True)

    _setup_netcdf4('/mnt/c/Users/oveha/Documents/test_parallel.nc', ncfiles, '/mnt/c/Users/oveha/Documents/dust_2019_ISRIC/',0)

