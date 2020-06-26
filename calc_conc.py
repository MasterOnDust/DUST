import numpy as np
import xarray as xr
import pandas as pd
from dask import delayed

"""
Both merging and writing to file want to be done in a parallell 
"""
def multiply_ems(ds, ems):
    ds = ds.spec001_mr
    conc = np.zeros_like(ds.values)
    for record in range(len(ds.record)):
        for time in range(len(ds.time)):
#                 
            conc[record,time,:,:] = (delayed((ds[record,time,:,:]*
                                            ems.sel(time=ds.record[record].values + ds.time[time])[:,:].values)/ds.height.values))

    return conc

def not_usefull(ds):
    essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                  'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

def pre(ds):
    ds = ds.assign_coords(record=pd.to_datetime(ds.iedate + ds.ietime))
    return ds.drop(not_usefull(ds))


if __name__ == "__main__":
            
    ds = xr.open_mfdataset('/cluster/projects/nn2806k/ovewh/flexpart/outputs/grid_time_2019*.nc',decode_times=False, 
                            combine='nested', concat_dim='record', preprocess=pre, parallel=True)
    #Subset dataset out of memory
    ds = ds.sel(height=100)
    ds = ds.isel(nageclass=0)
    ds = ds.sel(pointspec=0)

    ds = ds.assign_coords(time=pd.to_timedelta(ds.time.values,unit='s'))
    ems = xr.open_dataarray('/cluster/projects/nn2806k/ovewh/flexpart/dust_2019.nc', chunks={'time':40})

    ds['Concentration'] = multiply_ems(ds,ems)

    #Can I save dataset in parallel?
    ds.to_netcdf('/cluster/work/users/ovewh/merged.nc')