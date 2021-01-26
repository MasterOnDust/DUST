#!/usr/bin/env python

from netCDF4 import Dataset, num2date, date2num
import xarray as xr
import os
import glob
import pandas as pd
import shutil
from subprocess import call


def setup_netcdf(p0, pointspec,outpath='.', height=None,**netCDF_kwargs):
    d0 = Dataset(p0, 'r')
    dataVar = 'spec001_mr'
    lats = d0.variables['latitude'][:]
    lons = d0.variables['longitude'][:]
    heights = d0.variables['height'][:]
    rel_lat = d0.variables['RELLAT1'][:][pointspec]
    rel_lon = d0.variables['RELLNG1'][:][pointspec]
    ts = d0.variables['time'][:]
    time_units = d0['time'].units
    relcom = d0.variables['RELCOM'][pointspec][:]

    dims = d0.dimensions
    rel_com_str = ''
    for char in relcom:
        rel_com_str += char.decode('UTF-8') 
    rel_com_str = rel_com_str.strip()
    sdate = num2date(0,time_units).strftime('%Y%m%d%H%M')
    
    f_unit = d0[dataVar].units
    f_longname = d0[dataVar].long_name

    ind_receptor = d0.ind_receptor
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

    if ind_receptor == 1:
        f_name = 'Conc'
        f_unit = 's'
    elif ind_receptor == 3:
        f_name = 'WetDep'
        f_unit = 'm'
    elif ind_receptor ==4:
        f_name = 'DryDep'
        f_unit = 'm'
    else:
        raise(ValueError('Model settings not recognized, ind_receptor {}'.format(ind_receptor)))
    
    name_str = '_'.join(p0.split('/')[-1].split('_')[:2])
    outFileName = outpath + '/' + '_'.join(rel_com_str.split(' ')) + name_str + sdate +'.nc'
    
    try:
        ncfile = Dataset(outFileName, 'w', format="NETCDF4")
    except PermissionError:
        # If netcdf file exist delete old one
        os.remove(outFileName)
        ncfile = Dataset(outFileName, 'w', format='NETCDF4')
    
    lat_dim = ncfile.createDimension('lat', dims['latitude'].size)
    lon_dim = ncfile.createDimension('lon', dims['longitude'].size)
    height_dim = ncfile.createDimension('height',dims['height'].size)
    point_dim = ncfile.createDimension('npoint',1)
    #temporal dims

    time_dim = ncfile.createDimension('time', None)
    btime_dim = ncfile.createDimension('btime', dims['time'].size)

   
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
    btime.units = 'hours'
    btime.long_name = 'time along backtrajectory'

    btime[:] = ts/3600

    lon[:] = lons
    lat[:] = lats
    height[:] = heights
    rellat[:] = rel_lat
    rellon[:] = rel_lon

    time_var = ncfile.createVariable('time', 'f8', ('time',), **netCDF_kwargs)
    time_var.units = ''
    time_var.long_name = 'time along back trajectory'

    field = ncfile.createVariable(f_name, 'f4', ('time', 'btime', 'height','lat', 'lon'), **netCDF_kwargs)
    field.units = f_unit
    field.spec_name = f_longname
    field.long_name = 'SRR {}'.format(name_str) 

    # ncfile attributes

    ncfile.title = title
    ncfile.info = 'Senstivity to emission from FLEXPART backwards simulation'
    ncfile.version = version
#     ncfile.concatenated = ','.join(ncfiles)
    ncfile.dataVar = f_name
    ncfile.ind_receptor = ind_receptor
    ncfile.ind_source = ind_source
    ncfile.outlon0 = out_lon0
    ncfile.outlat0 = out_lat0
    ncfile.sdate = s_date
    ncfile.stime = s_time
    ncfile.lsubgrid = lsub_grid
    
    ncfile.close()

    return {'path':outFileName, 'point': pointspec}



def concat_output(ncfiles,outpath, n_processes = None, netCDF_kwargs={}):
    ncfiles.sort()
    receptors = xr.open_dataset(ncfiles[0]).numpoint.values
    p0 = ncfiles[0]

    outpaths = [setup_netcdf(p0, receptor) for receptor in receptors]
    
    dt1 = pd.to_datetime(ncfiles[0][-17:-3])

    dt2 = pd.to_datetime(ncfiles[1][-17:-3])

    

    t_delta = int((dt2 - dt1).seconds/3600)
    for outfile_dict in outpaths:
        time_step = 0
        point_spec = outfile_dict['point']
        out_path = outfile_dict['path'] 
        print(out_path)
        for n, ncfile in enumerate(ncFiles):
            if n == 0:
                outfile = Dataset(out_path, 'a', format="NETCDF4")
                print(outfile)
                outfield = outfile[outfile.dataVar]
            print(ncfile)
            outfile['time'][n] = time_step
            # print(xr.open_dataset(ncfile))
            ds = Dataset(ncfile, 'r', format='NETCDF4')
            # print(ds)
            outfield[:] = ds.variables['spec001_mr'][0,point_spec,:,:,:,:]
            if n % 10== 0:
                outfile.close()
                print(n,'saved')
                outfile = Dataset(out_path, 'a', format="NETCDF4")
                outfield = outfile[outfile.dataVar]
            time_step = time_step + t_delta
        outfile.close()
        
    
    



if __name__ == "__main__":
    import argparse as ap
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
        call(['list_available_output.py', path])
        df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
        df.index = pd.to_datetime(df.index, format='%Y%m%d-%H')
        ncFiles = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
        # ncFiles = glob.glob(path + "**/output/grid*.nc", recursive=True) #recursively find FLEXPART output files

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
    concat_output(ncFiles,dir_p)