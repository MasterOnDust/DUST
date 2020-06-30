from numba import njit, prange
import xarray as xr
import numpy as np
import os

@njit(parallel=True)
def element_wise(rep, source, height):

    conc = np.zeros_like(rep)
    s = conc.shape

    for i in prange(s[0]):
        for j in prange(s[1]):
            for k in prange(s[2]):
                conc[i,j,k] = rep[i,j,k]*source[-i,j,k]
    return conc

def multiply_ems(ds, ems_f):
    data = ds.spec001_mr
    ems_f = ems_f.sel(time=ds.record+ ds.time)
#     ems_f = ems_f.reindex(time=ems.time[::-1])
    conc = element_wise(data.values,ems_f.values,data.height.values)
    print(conc.max())
#     conc = data.dot(ems_f)
#     for i in range(len(data.time)):
#         temp_ems = ems.isel(time=-i)
#         temp_data = data.isel(time=i)
#         print(temp_data.dot(tem))

#     print(ems.transpose('time', 'lon', 'lat', ))
#         conc =((data[time,:,:]*ems.sel(time=ds.record+ data.time[time])[:,:].values)/data.height.values)
#     ds = xr.Dataset(data_vars=[conc,ds.data_vars],coords=ds.coords,attrs=ds.attrs)
    ds = ds.assign(concentration = xr.DataArray(conc, coords=data.coords, 
                                                dims=data.dims,attrs={'units':'kg/m^3','long_name':data.long_name}))
    return ds

def not_usefull(ds):
    essentials = ['RELCOM','RELLNG1','RELLNG2','RELLAT1','RELLAT2','RELZZ1','RELZZ2',
                  'RELKINDZ','RELSTART','RELEND','RELPART','RELXMASS','LAGE','ORO', 'spec001_mr']
    return  [v for v in ds.data_vars if v not in essentials]

def pre(ds):
    # add datetime to record dimension
    ds = ds.assign_coords(record=pd.to_datetime(ds.iedate + ds.ietime))
    return ds.drop(not_usefull(ds))

def sel_data(ds,ems,pointspec, name_of_location, output):
    ds = pre(ds)
    ds = ds.sel(height=100)
    ds = ds.isel(nageclass=0)
    ds = ds.sel(pointspec=0)
    ds = ds.assign_coords(time=pd.to_timedelta(ds.time.values,unit='s'))
    ds = multiply_ems(ds,ems.Emission)
    ds.attrs['Location'] = name_of_location
    ds.attrs['Output'] = output
    return ds

if __name__ == "__main__":
    path = '/cluster/work/users/ovewh/'

    ems = xr.open_dataarray('/cluster/projects/nn2806k/ovewh/flexpart/dust_2019.nc')
    dsets = []
    for ncfile in os.listdir(path):
        if ncfile.endswith('.nc'):
            ds = xr.open_dataset(path+ncfile,decode_times=False)
            ds = sel_data(ds, ems, pointspec=0, 
                name_of_location='SACOL',output='Concentration')
            dsets.append(ds)
        else:
            continue
    dset = xr.concat(dsets,dim='record')
    dset.to_netcdf('/test.nc')

    