import glob
import xarray as xr
import xarray

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



def multi_flexpart_flexdust(path, nc_files, flexdust, point_spec, **kwargs):
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

        return: path to newly created netcdf4 file

    AUTHOR:
    =======
        Ove Haugvaldstad
        ovehaugv@outlook.com

    """

    # Read in the first flexpart netcdf file as reference for setting up the output file
    d0 = xr.open_dataset(nc_files[0], decode_times=False)
    d0 = _sel_location(d0,point_spec)


    #Copy variables
    lats = d0.latitude.values
    lons = d0.longitude.values
    ts = d0.time.values
    ind_receptor = d0.ind_receptor
    dims = d0.dims
    relcom = str(d0.RELCOM.values)[2:].strip()[:-1].split()
    sdate = num2date(0,d0.time.units).strftime('%Y%m%d%H%M')
    d0.close()
    if ind_receptor == 1:
        f_name = 'Conc'
        field_unit = 'kg/m^3'
        field_name = 'Concentration'
    elif ind_receptor == 4:
        f_name = 'DryDep'
        field_unit = 'kg/m^2 s'
        field_name = 'Dry depostion'
    elif ind_receptor == 3:
        f_name = 'WetDep'
        field_unit = 'kg/m^2 s'
        field_name = 'Wet depostion'
    else:
        f_name = 'spec_mr'
        field_units = 'kg/m^3'
        field_name = 'Unknown'
    # Read flexdust

    if isinstance(flexdust,xarray.Dataset):
        emsField=flexdust.Emission
        area = flexdust.area
    else:
        flexdust = DUST.read_flexdust_output(flexdust)['dset']
        emsField = flexdust.Emission
        area = flexdust.area
    t_d = pd.to_timedelta((emsField.time[1]-emsField.time[0]).values)
    time_int = t_d.total_seconds()

    # Create netcdf outfile

    outFileName = path + '/' + '_'.join(relcom) + f_name + sdate +'.nc'
    try:
        ncfile = Dataset(outFileName, 'w', format="NETCDF4")
    except PermissionError:
        # If netcdf file exist delete old one
        os.remove(outFileName)
        ncfile = Dataset(outFileName, 'w', format='NETCDF4')

    #Setup attributes
    ncfile.title = 'Flexpart - Flexdust SSR'
    ncfile.history = "Created " + time.ctime(time.time())
    ncfile.flexpart_v = d0.source
    ncfile.receptor_name =  ' '.join(relcom)
    ncfile.reference = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'
    ncfile.varName = f_name


    #Spatial dims

    lat_dim = ncfile.createDimension('lat', dims['latitude'])
    lon_dim = ncfile.createDimension('lon', dims['longitude'])
    point_dim = ncfile.createDimension('npoint',1)
    #temporal dims

    time_dim = ncfile.createDimension('time', None)
    btime_dim = ncfile.createDimension('btime', dims['time'])


    #Setup lon/lat (spatial variables)

    lat = ncfile.createVariable('lat', 'f4', ('lat', ),**kwargs)
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'

    lon = ncfile.createVariable('lon', 'f4', ('lon', ), **kwargs)
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'

    #Set receptor location
    rellat = ncfile.createVariable('RELLAT', 'f4', ('npoint',),**kwargs)
    rellat.units = 'degrees_north'
    rellat.long_name = 'latitude_receptor'

    rellon = ncfile.createVariable('RELLNG', 'f4', ('npoint',), **kwargs)
    rellon.units = 'degrees_east'
    rellon.long_name = 'longitude_receptor'

    lon[:] = lons
    lat[:] = lats

    rellat[:] = d0.RELLAT1.values
    rellon[:] = d0.RELLNG1.values

    # Setup area variable

    area_var = ncfile.createVariable('area', 'f4', ('lat','lon'),**kwargs)
    area_var.units = 'm^2'
    area_var.long_name = 'gridbox_area'

    area_var[:] = area.values
    # Setup temporal variable
    btime = ncfile.createVariable('btime', 'i4', ('btime',), **kwargs)
    btime.units = 'hours'
    btime.long_name = 'time along backtrajectory'

    btime[:] = ts/3600

    time_var = ncfile.createVariable('time', 'i4', ('time',), **kwargs)
    time_var.units = "hours since 1980-01-01"
    time_var.long_name = 'time'

    # Determind which kind of output should be created


    field = ncfile.createVariable(f_name, 'f4', ('time', 'btime', 'lat', 'lon'), **kwargs)
    field.units = field_unit
    field.long_name = field_name

    height = d0.height.values
    field.height = height     #set the height of the lowest model output layer
    field.lon0 = d0.RELLNG1.values
    field.lat0 = d0.RELLAT1.values
    field.name_location = '-'.join(relcom)
    field.time_int = time_int
    td = pd.to_timedelta(d0.time.values, unit='s')
    ncfile.close()
    #unit_correction = 1/height
    unit_correction = 1/(height*time_int)

    for n, nc_file in enumerate(nc_files):
        if n == 0:
            ncfile = Dataset(outFileName, 'a', format="NETCDF4")
            outfield = ncfile[f_name]
        temp_ds = xr.open_dataset(nc_file, decode_times=False)
        temp_ds = _sel_location(temp_ds, point_spec)
        s_time = pd.to_datetime(temp_ds.iedate + temp_ds.ietime)
        time_step = date2num(s_time, units = "hours since 1980-01-01")
        ems_sens = temp_ds.spec001_mr
        temp_ems = emsField.sel(time=s_time + td)

        #multiply emission and emission sensitivity
        temp_array = _multiply_emissions(temp_ems.values, ems_sens.values)
        #after multiplying the emissions and emission sensitivity has units kg.s/m^2 or kg/m
        temp_array = temp_array*unit_correction
        ncfile['time'][n] = time_step
        #write to file, concentration in kg/m^3 or total deposition in kg/m^2 s
        outfield[n,:,:,:] = temp_array

        if n % 5== 0:
            ncfile.close()
            print(n,'saved')
            ncfile = Dataset(outFileName, 'a', format="NETCDF4")
            outfield = ncfile[f_name]
    ncfile.close()
    return outFileName

@njit
def _multiply_emissions(ems, ems_sens):
    """python looping is slow..."""
    temp_array = np.zeros_like(ems_sens)
    for i in range(ems.shape[0]):
        for j in range(ems.shape[1]):
            for k in range(ems.shape[2]):
                # Emsission field has unit kg/m^2 ems_sens has unit m or s
                temp_array[i,j,k] = ems[i,j,k]*ems_sens[i,j,k]
    return temp_array


def _sel_location(ds,pointspec):
    ds = ds.sel(pointspec=pointspec)
    ds = ds.sel(numpoint=pointspec)
    ds = ds.sel(numspec=pointspec)
    ds = ds.sel(nageclass=0)
    ds = ds.isel(height=0)
    return ds

@dask.delayed
def multi_flexpart_flexdust_dask(xr_flexpart, xr_flexdust, point_spec):
    """
    NOT implemented yet!

    DESCRIPTION
    ===========

        Multiply flexpart and flexdust xarray dataset in paraell using dask, and returns
        dask delayed object.

    USAGE
    =====

        xr_multi = multi_flexpart_flexdust_dask(xr_flexpart, xr_flexdust, point_spec)

        xr_flexpart : xarray dataset containing flexpart model output(s)
        xr_flexdust : xarray dataset containing flexdust model output which is to be
                      matched and multiplied by the flexpart output


    """

    # Create output dataset
    attrs = dict(title ='FLEXPART * FLEXDUST multiplied SSR',
                history = "Created " + time.ctime(time.time()),
                reference = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482')


    xr_flexpart = _sel_location(point_spec)
    d_out = xr.Dataset(data_vars={'RELCOM' : xr_flexpart.RELCOM,
                                'RELLAT' : xr_flexpart.RELLAT1,
                                'RELLON' : xr_flexpart.RELLNG1,
                                'ORO' : xr_flexpart.ORO,
                                'RELPART' : xr.RELPART},
                    coords={'time':xr_flexpart.time,
                               'btime': xr_flexpart.btime,
                               'lon' : xr_flexpart.lon,
                               'lat' : xr_flexdust.lat
                            }, attrs=attrs)
    # for d_step in xr_flexpart:






    return xr_multi

if __name__=="__main__":

    path = '/mnt/c/Users/oveha/Documents/FLEXPART_output'

    ncfiles = glob.glob(path + "/**/*.nc", recursive=True)

    multi_flexpart_flexdust('/mnt/c/Users/oveha/Documents/test_parallel.nc', ncfiles, '/mnt/c/Users/oveha/Documents/dust_2019_ISRIC/',0, zlib=True)

