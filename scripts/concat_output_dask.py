import argparse as ap
from netCDF4 import Dataset, num2date, date2num
import xarray as xr
import os
import glob
import pandas as pd
import shutil
import DUST
from dask.distributed import Client, LocalCluster

def concat_output(ncfiles,outpath, locations='ALL', time_slice=None, netCDF_kwargs={}):
    ncfiles.sort()
    name_str = '_'.join(ncfiles[0].split('/')[-1].split('_')[:2])
    d = xr.open_dataset(ncfiles[0])
    relCOMS = d.RELCOM
    d.close()

    dsets = DUST.read_multiple_flexpart_output(path=ncfiles)

    if time_slice:
        dsets = dsets.sel(time=time_slice) 

    ind_receptor = dsets.ind_receptor
    if ind_receptor == 1:
        f_name = 'conc' 
        long_name = 'sensitivity to concentration'
        units = 's'
    elif ind_receptor == 3:
        f_name = 'wetdep'
        long_name = 'sensitvity to wet deposition'
        units = 'm'
    elif ind_receptor ==4:
        f_name = 'drydep'
        long_name = 'sensitvity to dry deposition'
        units = 'm'
    else:
        raise(ValueError('Value for ind_receptor not recognized {}'.format(ind_receptor)))   
    dsets = dsets.rename({'spec001_mr' : f_name})
    dsets = dsets.assign_attrs(dataVar=f_name)
    dsets[f_name].attrs['units'] = units
    dsets[f_name].attrs['long_name'] = long_name
    dsets = dsets.assign_attrs({'dataVar':f_name})

    outfileNames = []
    loc_data = []
    sdate = dsets.ibdate
    comp = dict(zlib=True, complevel=5)
    for i , com in enumerate(relCOMS):
        loc = str(com.values)[2:].strip().split()
        if locations == 'ALL':
            temp_dset = dsets.sel(pointspec=i, numpoint=i, nageclass=0)


            loc_data.append(temp_dset)
            outfileNames.append(outpath + '/' + '_'.join(loc[:2]) + name_str + sdate +'.nc')

        else:
            for receptor in locations:
                if receptor in loc or receptor == str(i):
                    temp_dset = dsets.sel(pointspec=i, numpoint=i, nageclass=0)

                    loc_data.append(temp_dset)

                    outfileNames.append(outpath + '/' + '_'.join(loc[:2]) + name_str + sdate +'.nc')
                    
                else:
                    continue
    for path, dset in zip(outfileNames, loc_data):
        dset.to_netcdf(path, encoding = {f_name:
                            {'dtype' :'f8' , 'zlib':True, 'complevel' : 6}}, unlimited_dims = 'time')

    

if __name__ == "__main__":
    
    parser = ap.ArgumentParser(description='Concat FLEXPART output from backward simulation along a single time dimmension')
    parser.add_argument('path',nargs='+', help='path to top directory containing flexpart output')
    parser.add_argument('--outpath', '--op', help='where the concatinated output should be stored', default='.')
    parser.add_argument('--locations', '--loc', help='number of location to concatinate output for', default='ALL')
    parser.add_argument('--bdate', '--bd', help='Beginning of time slice (YYYY-MM-DD)', default=None)
    parser.add_argument('--edate', '--ed', help='End of time slice (YYYY-MM-DD)', default=None)
    parser.add_argument('--memory_limit', default='4GB', help='memory limit local dask cluster')
    parser.add_argument('--n_worker', default=8, type=int, help='Number of dask workers')
    parser.add_argument('--use_cluster', '--uc', action='store_true', help='Weather to use cluster or not')
    parser.add_argument('--threads_per_worker', '--tpw', default=1, type=int)

    args = parser.parse_args()
    paths = args.path
    outpath = args.outpath
    locations = args.locations
    bdate = args.bdate
    edate = args.edate
    n_workers=args.n_worker
    memory_limit=args.memory_limit
    uc = args.use_cluster
    n_wthreads = args.threads_per_worker

    cluster = LocalCluster(n_workers=n_workers, memory_limit=memory_limit, threads_per_worker=n_wthreads)
    client = Client(cluster) 
    concat_output(ncFiles,dir_p, locations=locations,time_slice=time_slice)

    for path in paths:
        if bdate and edate == None:
            time_slice = None
        elif bdate and edate:
            e_time = edate
            s_time = bdate


        elif bdate:
            s_time = bdate
        else:
            e_time = edate

        if path.endswith('/') == False:
            path = path +'/'
        #IF AVAILABLE_OUPUT file is created, use that before recursive search, slow on mounted system 
        try:
            df = pd.read_csv(path+'AVAILABLE_OUTPUT', index_col=0)
            df.index = pd.to_datetime(df.index, format='%Y%m%d-%H')
            df = df[s_time:e_time]
            ncFiles = [path+row['dir_paths'] + '/'+ row['ncfiles'] for index,row in df.iterrows()]
            time_slice = None
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
