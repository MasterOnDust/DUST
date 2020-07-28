import glob
import xarray as xr
import pandas as pd
import DUST
from dask.delayed import delayed
import dask
import numpy as np
from netCDF4 import Dataset
import os
import time
from IPython import embed

"""
Description:
===========
    Create empty netcdf file for storing the merged flexpart/flexdust output.

USAGE:
=====
    path : path to where the netcdf file should be created
    locations: name of the locations contained in pointpsec (python list)
    d0 : xarray dataset necessary for setting up dimmensions and metadata

    _setup_netdcf4(path, locations, d0)


AUTHOR:
=======
    Ove Haugvaldstad
    ovehaugv@outlook.com

"""


def _setup_netcdf4(path,
                  locations,
                  d0):
    lat = d0.latitude.values
    lon = d0.longitude.values
    t = d0.time.values
    ind_receptor = d0.ind_receptor
    dims = d0.dims
    try:
        ncfile = Dataset(path, 'w', format="NETCDF4")
    except PermissionError:
        os.remove(path)
        ncfile = Dataset(path, 'w', format='NETCDF4')
    groups = []
    for location in locations:
        groups.append(ncfile.createGroup("/{}".format(location)))
    for group in groups:
        recordDim = group.createDimension('record', None)
        timeDim = group.createDimension('time', dims['time'])
        longitudeDim = group.createDimension('lon', dims['longitude'])
        latitudeDim = group.createDimension('lat', dims['latitude'])
        lats = group.createVariable('lat', 'f4', ('lat', ))
        lons = group.createVariable('lon', 'f4', ('lon', ))
        times = group.createVariable('time', 'i4', ('time',))
        lons[:] = lon
        lats[:] = lat
        times[:] = t
        s_time = pd.to_datetime(d0.iedate + d0.ietime)
        record[0] = s_timw

        if ind_receptor == 1:
            field = group.createVariable('Concentration', 'f4',('record', 'time', 'lat', 'lon'))
            field.units = 'kg/m^3'
        elif ind_receptor == 4:
            field = group.createVariable('DryDep', 'f4', ('record', 'time', 'lat', 'lon'))
            field.units = 'kg/m^2s '
        elif ind_receptor == 3:
            field = group.createVariable('WetDep', 'f4', ('record', 'time', 'lat', 'lon'))
            field.units = 'kg/m^2s '
        else:
            field = group.createVariable('Spec_mr', 'f4', ('record', 'time', 'lat', 'lon'))
            field.units = 'kg/m^3'
        group.title = 'FLEXPART/FLEXDUST SSR'
        group.history = "Created " + time.ctime(time.time())
        group.flexpart_v = d0.source
        group.reference = 'https://doi.org/10.5194/gmd-12-4955-2019, https://doi.org/10.1002/2016JD025482'
    ncfile.close()


def multiply_emissions(nc_files, flexdust_path, outfile, point_spec):
    emsField = DUST.read_flexdust_output(flexdust_path)['dset'].Emission
    for nc_file, n in enumerate(nc_files):
        ofile = Dataset(outfile, 'a')
        temp_loc = ofile.groups['SACOL']
        oVar = temp_loc.variables['DryDep']
        temp_ds = xr.open_dataset(nc_file, decode_times=False)
        s_time = pd.to_datetime(temp_ds.iedate + temp_ds.ietime)
        temp_ds = temp_ds.spec001_mr
        temp_ds = temp_ds.sel(height=100)
        temp_ds = temp_ds.sel(nageclass=0)
        temp_ds = temp_ds.sel(pointspec=point_spec)
        td = pd.to_timedelta(temp_ds.time.values, unit='s')
        temp_ems = emsField.sel(time=s_time + td)
        for i in range(temp_ds.time.shape[0]):
            for j in range(temp_ds.latitude.shape[0]):
                for k in range(temp_ds.longitude.shape[0]):
                    oVar[0,i,j,k] = temp_ems[i,j,k]*temp_ds[i,j,k]/temp_ds.height
        ofile.close()
        del oVar


if __name__=="__main__":

    path = '/cluster/work/users/ovewh/accumulated_avg_drydep_test'

    nc_files = glob.glob(path + "/**/*.nc", recursive=True)
    d0 = xr.open_dataset(nc_files[0])
    _setup_netcdf4('test.nc', ['SACOL', 'BAODE'], d0)
    fdp = '/cluster/projects/nn2806k/ovewh/flexdust/outputs/dust_2019_ISRIC/'
    emsField = DUST.read_flexdust_output(fdp)['dset']
    multiply_emissions(nc_files, fdp, 'test.nc',0)
