from dask.distributed import Client
import xarray as xr
import pandas as import pd
import dask
from dask.delayed import delayed
import numpy as np
from dask.diagnostics import ProgressBar
from dask_jobqueue import SLURMCluster

def not_usefull(ds):
    essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                  'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

def pre(ds):
    # add datetime to record dimension
    ds = ds.assign_coords(record=pd.to_datetime(ds.iedate + ds.ietime))
    ds=ds.sel(height=100)
    ds=ds.sel(pointspec=point)
    ds=ds.sel(nageclass=0)
    ds=ds.assign_coords(time=pd.to_timedelta(ds.time.values,unit='s'))
    return ds.drop(not_usefull(ds))

def calc_concentration(step, ems):
    return np.einsum('mnr,mnr->mnr', step, ems)

if if __name__ == "__main__":

    cluster =SLURMCluster(project='nn2806k', 
        local_directory='/cluster/work/users/ovewh',
        cores = 16,
        memory = '64 GB',
        walltime = '00:30:00'

    )

    """
    Shoud use argeparse at some point!
    """
    
    point = 0
    path_flexdust = 'path_to_flexdust'
    path_flexpart = '/opt/uio/flexpart/nc_files_apr2019/grid_time_20190403*.nc'
    outfile = 'some path '
    location = 'name of location'
    cluster.scale(1)
    client = Client(cluster)
    client.cluster
    dsets = xr.open_mfdataset(path_flexpart, parallel=True,combine='nested',
                          concat_dim='record', preprocess=pre, decode_times=False)
    ems = xr.open_dataarray(path_flexdust)

    spec_mr = dsets.spec001_mr
    time_steps = [step for step in spec_mr]
    ems_steps = [ems.sel(time=step.record + step.time) for step in spec_mr]
    conc = [delayed(calc_concentration)(sens_step.values, ems_step.values) for ems_step, sens_step in zip(ems_steps, time_steps)]
    conc_computed = client.compute(conc)
    da = xr.DataArray(data=conc_computed[0],coords=spec_mr.coords,dims=spec_mr.dims, attrs=spec_mr.attrs).chunk(spec_mr.chunks)
    dsets['Concentration'] = da
    dsets.encoding['unlimited_dims'] = 'record'
    dsets.Concentration.encoding = {
        'zlib' : True,
        'shuffle': False,
        'complevel': 9,
        'fletcher32': False,
        'contiguous': False,
        'chunksizes': (1, 1, 1, 4, 320, 680),
        'original_shape': (1, 7, 40, 4, 320, 680),
        'dtype': np.dtype('float32')
        }
    dsets.attrs['Location'] = location
    delayed_save = dsets.to_netcdf(outfile, compute=False)
    with ProgressBar():
        result = client.compute(delayed_save)