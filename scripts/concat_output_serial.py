import argparse as ap
from netCDF4 import Dataset, num2date, date2num
import xarray as xr
import os
import glob
import pandas as pd
import shutil

def _sel_location(ds,pointspec, height=None):
    ds = ds.sel(pointspec=pointspec)
    ds = ds.sel(numpoint=pointspec)
    ds = ds.sel(numspec=pointspec)
    ds = ds.sel(nageclass=0)
    if height != None:
        ds = ds.sel(height=height)
    return ds

def concat_output(ncfiles,pointspec,outpath, netCDF_kwargs={}):
    ncfiles.sort()
    
    d0 = xr.open_dataset(ncfiles[0], decode_times=False)
    d0 = _sel_location(d0,pointspec)


    #Copy variables
    dataVar = 'spec001_mr'
    lats = d0.latitude.values
    lons = d0.longitude.values
    heights = d0.height.values
    rel_lat = d0.RELLAT1.values
    rel_lon = d0.RELLNG1.values
    ts = d0.time.values
    ind_receptor = d0.ind_receptor
    dims = d0.dims
    relcom = str(d0.RELCOM.values)[2:].strip()[:-1].split()
    sdate = num2date(0,d0.time.units).strftime('%Y%m%d%H%M')
    f_unit = d0[dataVar].units
    f_longname = d0[dataVar].long_name
    version = d0.source
    out_lat0 = d0.outlat0
    out_lon0 = d0.outlon0
    ind_source = d0.ind_source
    ind_receptor = d0.ind_receptor
    lout_step = d0.loutstep
    lout_aver = d0.loutaver
    lsub_grid = d0.lsubgrid
    s_time = d0.ietime
    s_date = d0.iedate


    title = d0.title


    d0.close()
    name_str = '_'.join(ncfiles[0].split('/')[-1].split('_')[:2])
    outFileName = outpath + '/' + '_'.join(relcom) + name_str + sdate +'.nc'
    try:
        ncfile = Dataset(outFileName, 'w', format="NETCDF4")
    except PermissionError:
        # If netcdf file exist delete old one
        os.remove(outFileName)
        ncfile = Dataset(outFileName, 'w', format='NETCDF4')
    
    lat_dim = ncfile.createDimension('lat', dims['latitude'])
    lon_dim = ncfile.createDimension('lon', dims['longitude'])
    height_dim = ncfile.createDimension('height',dims['height'])
    point_dim = ncfile.createDimension('npoint',1)
    #temporal dims

    time_dim = ncfile.createDimension('time', None)
    btime_dim = ncfile.createDimension('btime', dims['time'])

   
    #Setup lon/lat (spatial variables)

    lat = ncfile.createVariable('lat', 'f4', ('lat', ),**netCDF_kwargs)
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'

    lon = ncfile.createVariable('lon', 'f4', ('lon', ), **netCDF_kwargs)
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'

    height = ncfile.createVariable('height', 'i4',('height',),**netCDF_kwargs)
    height.units = 'm'
    height.long_name = 'height above ground'

    #Set receptor location
    rellat = ncfile.createVariable('RELLAT', 'f4', ('npoint',),**netCDF_kwargs)
    rellat.units = 'degrees_north'
    rellat.long_name = 'latitude_receptor'

    rellon = ncfile.createVariable('RELLON', 'f4', ('npoint',), **netCDF_kwargs)
    rellon.units = 'degrees_east'
    rellon.long_name = 'longitude_receptor'

    btime = ncfile.createVariable('btime', 'i4', ('btime',), **netCDF_kwargs)
    btime.units = 's'
    btime.long_name = 'seconds_since_release'

    btime[:] = d0.time.values

    lon[:] = lons
    lat[:] = lats
    height[:] = height
    rellat[:] = rel_lat
    rellon[:] = rel_lon

    time_var = ncfile.createVariable('time', 'f8', ('time',), **netCDF_kwargs)
    time_var.units = "hours since 1980-01-01"
    time_var.long_name = 'time'

    field = ncfile.createVariable(name_str, 'f4', ('time', 'btime', 'height','lat', 'lon'), **netCDF_kwargs)
    field.units = f_unit
    field.spec_name = f_longname
    field.long_name = 'SRR {}'.format(name_str) 

    # ncfile attributes

    ncfile.title = title
    ncfile.info = 'Senstivity to emission from FLEXPART backwards simulation'
    ncfile.version = version
    ncfile.concatenated = ','.join(ncfiles)
    ncfile.dataVar = name_str
    ncfile.ind_receptor = ind_receptor
    ncfile.ind_source = ind_source
    ncfile.outlon0 = out_lon0
    ncfile.outlat0 = out_lat0
    ncfile.sdate = s_date
    ncfile.stime = s_time
    ncfile.lsubgrid = lsub_grid
    
    ncfile.close()


    for n,file_path in enumerate(ncfiles):
        if n == 0:
            ncfile = Dataset(outFileName, 'a', format="NETCDF4")
            outfield = ncfile[name_str]
        temp_ds = xr.open_dataset(file_path)[dataVar]
        temp_ds = _sel_location(temp_ds, pointspec)
        s_time = pd.to_datetime(temp_ds.iedate + temp_ds.ietime)
        time_step = date2num(s_time, units = "hours since 1980-01-01")
        temp_da = temp_ds[dataVar].values
        ncfile['time'][n] = time_step
        outfield[:] = temp_da
        if n % 5== 0:
            ncfile.close()
            print(n,'saved')
            ncfile = Dataset(outFileName, 'a', format="NETCDF4")
            outfield = ncfile[name_str]




if __name__ == "__main__":
    
    parser = ap.ArgumentParser(description='Concat FLEXPART output from backward simulation along a single time dimmension')
    parser.add_argument('path', help='path to top directory containing flexpart output')
    parser.add_argument('--outpath', '--op', help='where the concatinated output should be stored', default='.')
    parser.add_argument('--locations', '--loc', help='number of location to concatinate output for', default='ALL')

    args = parser.parse_args()
    path = args.path
    outpath = args.outpath
    locations = args.locations


    if path.endswith('/') == False:
        path = path +'/'
    #IF AVAILABLE_OUPUT file is created, use that before recursive search, slow on mounted system 
    try:
        df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
        ncFiles = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
    except FileNotFoundError:
        ncFiles = glob.glob(path + "**/output/grid*.nc", recursive=True) #recursively find FLEXPART output files

    d = xr.open_dataset(ncFiles[0])
    relCOMS = d.RELCOM
    ind_receptor = d.ind_receptor
    d.close()
    

    if ind_receptor == 1:
        f_name = 'Conc'
    elif ind_receptor == 3:
        f_name = 'WetDep'
    elif ind_receptor ==4:
        f_name = 'DryDep'
    else:
        f_name = 'Unknown'
    
    e_time = pd.to_datetime(ncFiles[-1][-17:-3]).strftime('%Y-%m-%d')
    
    s_time = pd.to_datetime(ncFiles[0][-17:-3]).strftime('%Y-%m-%d')
    
    dir_p = outpath+'/'+f_name + '_FLEXPART_SRR_{}_{}'.format(s_time, e_time)
    try:
        os.mkdir(dir_p)
    except FileExistsError:

        shutil.rmtree(dir_p)
        os.mkdir(dir_p)
    dir_p = dir_p + '/'
    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()
    if locations == 'ALL':
        concat_output(ncFiles,outpath=dir_p,pointspec=i)
    else:
        for receptor in locations:
            if receptor in loc or receptor == str(i):
                concat_output(ncFiles, outpath=dir_p, pointspec=i)
            else:
                continue
