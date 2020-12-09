#!/usr/bin/env python
from netCDF4 import Dataset
from DUST import read_multiple_flexpart_output
import time
import xarray as xr

import pandas as pd
from functools import partial
import os
import shutil

import argparse as ap

from dask.distributed import Client, LocalCluster

def not_usefull(ds):
    essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

def pre_sum(ds, pointspec):
    ds = ds.rename(dict(longitude = 'lon', time= 'btime', latitude='lat'))
    ds = ds.sel(pointspec=pointspec, numpoint=pointspec, numspec=pointspec, nageclass=0)
    ds = ds.assign_coords(time=pd.to_datetime(ds.iedate + ds.ietime))
    ds = ds.drop(not_usefull(ds))
    ds = ds.sum(dim='btime')


    return ds



def avg_emission_senstivity(path, ncfiles, pointspec ,date_slice=None, data_var='spec001_mr',**kwargs):
    pre = partial(pre_sum,pointspec=pointspec)


    if date_slice == None:
        dsets = xr.open_mfdataset(ncfiles, preprocess=pre, decode_times=False,
                            combine='nested', concat_dim='time', parallel=True)
        s_time = pd.to_datetime(dsets.time[0].values).strftime('%Y%m%d')
        e_time = pd.to_datetime(dsets.time[-1].values).strftime('%Y%m%d')

    else:
        s_time = pd.to_datetime(date_slice.start).strftime('%Y%m%d')
        e_time = pd.to_datetime(date_slice.stop).strftime('%Y%m%d')
        ncdict = {'{}'.format(pd.to_datetime(path[-17:-3])):path for path in ncfiles}
        ncdict = pd.DataFrame.from_dict(ncdict, orient='index')
        ncdict = ncdict[date_slice]
        dsets = xr.open_mfdataset(ncdict[0].values, preprocess=pre, decode_times=False,
                            combine='nested', concat_dim='time', parallel=True)



    d0 = xr.open_dataset(ncfiles[0], decode_times=False)
    d0 = d0.sel(pointspec=pointspec, numpoint=pointspec, numspec=pointspec)
    lats = d0.latitude.values
    lons = d0.longitude.values
    dims = d0.dims
    height = d0.height
    relpart = d0.RELPART.values
    relcom = str(d0.RELCOM.values)[2:].strip()[:-1].split()
    d0.close()
    ind_receptor = d0.ind_receptor

    if ind_receptor == 1:
        f_name = 'Conc_emsens'
        field_unit = 's'
        field_name = 'Sensitivity to Concentration'
    elif ind_receptor == 4:
        f_name = 'DryDep_emsens'
        field_unit = 'm'
        field_name = 'Sensitivity to dry deposition'
    elif ind_receptor == 3:
        f_name = 'WetDep_emsens'
        field_unit = 'm'
        field_name = 'Sensitivity to wet deposition'
    else:
        field = 'Spec_mr'
        field.units = 's'
        field_name = 'Unknown'

    outFileName = path + '/' + '_'.join(relcom) +'_avg_{}_{}_{}'.format(f_name,s_time, e_time) + '.nc'
    try:
        ncfile = Dataset(outFileName, 'w', format="NETCDF4")
    except PermissionError:
        # If netcdf file exist delete old one
        os.remove(outFileName)
        ncfile = Dataset(outFileName, 'w', format='NETCDF4')

    ncfile.title = 'Flexpart emission sensitivity'
    ncfile.history = "Created " + time.ctime(time.time())
    # ncfile.flexpart_v = d0.source
    ncfile.receptor_name =  ' '.join(relcom)
    ncfile.reference = 'https://doi.org/10.5194/gmd-12-4955-2019'

    ncfile.avg_window_start = s_time
    ncfile.avg_window_end = e_time
    ncfile.particle_released = relpart
    ncfile.varName = f_name

    lat_dim = ncfile.createDimension('lat', dims['latitude'])
    lon_dim = ncfile.createDimension('lon', dims['longitude'])
    height_dim = ncfile.createDimension('height', dims['height'])
    point_dim = ncfile.createDimension('npoint',1)

    lat = ncfile.createVariable('lat', 'f4', ('lat', ), **kwargs)
    lon = ncfile.createVariable('lon', 'f4', ('lon',), **kwargs)
    height = ncfile.createVariable('height', 'i4', ('height', ), **kwargs)

    rellat = ncfile.createVariable('RELLAT', 'f4', ('npoint',),**kwargs)
    rellat.units = 'degrees_north'
    rellat.long_name = 'latitude_receptor'

    rellon = ncfile.createVariable('RELLON', 'f4', ('npoint',), **kwargs)
    rellon.units = 'degrees_east'
    rellon.long_name = 'longitude_receptor'

    lon[:] = lons
    lat[:] = lats


    rellat[:] = d0.RELLAT1.values
    rellon[:] = d0.RELLNG1.values

    height[:] = d0.height.values


    field = ncfile.createVariable(f_name, 'f4', ('height','lat', 'lon'), **kwargs)
    field.units = field_unit
    field.long_name = field_name

    field[:] = dsets[data_var].mean(dim='time').values

    ncfile.close()

if __name__ == "__main__":
    import glob
    defaultLocations = 'ALL'
    helpLocation = "Name or integer correspoding to location defined in FLEXPART RELEASE file"
    helpZlib = "Whether to use zlib compression"

    parser = ap.ArgumentParser(description='''Averages emission sensitivities
    from FLEXPART netCDF files''')

    parser.add_argument('--out_path', '--op', default='.', help='Path to where averaged output will be stored')
    parser.add_argument('flexpart_topdir', help='top directory containing the FLEXPART out')
    parser.add_argument("--zlib", "--zl", help=helpZlib, action='store_true')
    parser.add_argument("--locations", "--locs",default=defaultLocations, help=helpLocation, nargs='+')
    parser.add_argument('--etime', '--et', help='end of averaging window', default= None)
    parser.add_argument('--stime','--st', help = 'start of averaging window', default=None)
    parser.add_argument('--use_client', '--uc', help = 'create dask client', action='store_true')

    args = parser.parse_args()
    path = args.flexpart_topdir
    locations = args.locations
    outpath = args.out_path
    zlib = args.zlib
    e_time = args.etime
    s_time = args.stime
    create_client = args.use_client
    ncfiles = glob.glob(path + '/**/grid*.nc', recursive=True)
    ncfiles.sort()


    d = xr.open_dataset(ncfiles[0])
    relCOMS = d.RELCOM
    ind_receptor = d.ind_receptor
    d.close()

    if e_time == None:
        e_time = pd.to_datetime(ncfiles[-1][-17:-3]).strftime('%Y-%m-%d')
    if s_time == None:
        s_time = pd.to_datetime(ncfiles[0][-17:-3]).strftime('%Y-%m-%d')
    print(create_client)
    if create_client==True:
        cluster = LocalCluster(n_workers=32, threads_per_worker=1, memory_limit='16GB')
        client= Client(cluster)
        print(cluster)

    date_slice = slice(s_time, e_time)
    if ind_receptor == 1:
        f_name = 'Conc'
    elif ind_receptor == 3:
        f_name = 'WetDep'
    elif ind_receptor ==4:
        f_name = 'DryDep'
    else:
        f_name = 'Unknown'

    dir_p = outpath+'/'+f_name + '_mean_{}_{}'.format(s_time, e_time)
    try:
        os.mkdir(dir_p)
    except FileExistsError:

        shutil.rmtree(dir_p)
        os.mkdir(dir_p)

    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()
        if locations == 'ALL':
            avg_emission_senstivity(dir_p,ncfiles, i, date_slice=date_slice, zlib=zlib)
        else:
            for receptor in locations:
                if receptor in loc or receptor == str(i):
                    avg_emission_senstivity(dir_p,ncfiles, i, date_slice=date_slice,  zlib=zlib)
                else:
                    continue
