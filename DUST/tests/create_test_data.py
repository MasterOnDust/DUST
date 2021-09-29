from numpy.testing._private.utils import assert_approx_equal
import xarray as xr
import numpy as np
import pandas as pd
from DUST.process_data_dust import process_per_pointspec, create_output
import pytest

def create_flexdust_test_data(seed=None):
    time_vals = np.arange(10800, 10800*8*24,10800)
    time_var = xr.Variable('time',time_vals,attrs={'units':'seconds since 2000-03-01 00:00:00','calendar':'proleptic_gregorian'})
    longitude = xr.Variable('lon',np.arange(75,86,0.5))
    latitude = xr.Variable('lat', np.arange(30,35,0.5))
    da = xr.DataArray(data=np.random.rand(len(time_var),10,22),
        dims={'time' : len(time_var), 'lat':10, 'lon':22},
        coords={'time':time_var, 'lat':latitude, 'lon':longitude}
        )
    ds=xr.Dataset({'Emission':da},
                coords={'lon':longitude, 'time':time_var,
                    'lat':latitude})
    ds = xr.decode_cf(ds)
    return ds

def create_test_data(seed=None, ind_receptor=1):
    
    rs = np.random.RandomState(seed)
    time_var = pd.date_range("2000-03-10", periods=80, freq='3H')
    rel_com_time = time_var.strftime('%Y%m%d%H')[::-1][:9]
    time_var = xr.Variable('time',time_var[::-1][1:])
    
    longitude = xr.Variable('longitude',np.arange(75,86,0.5))
    latitude = xr.Variable('latitude', np.arange(30,35,0.5))
    height = xr.Variable('height',[100, 1000, 5000])
    relcom = xr.Variable('numpoint', ['{} TEST_RELEASE'.format(ts) for ts in rel_com_time])
    rellat = xr.Variable('numpoint', [32.3 for i in range(9)])
    rellng = xr.Variable('numpoint', [78.3 for i in range(9)])
    relz1 = xr.Variable('numpoint', [0 for i in range(9)])
    relz2 =xr.Variable('numpoint', [30 for i in range(9)])
    relpart =xr.Variable('numpoint', [3000 for i in range(9)])

    da = xr.DataArray(data=np.random.rand(1,9,79,3,10,22),
        dims={'nageclass':1, 'pointspec':9, 'time' : 79, 'height':3, 'latitude':10, 'longitude':22},
        coords={'time':time_var, 'latitude':latitude, 'longitude':longitude, 'height':height},
        attrs={'long_name':'Test-Spec'}
        )

    
    ds = xr.Dataset(data_vars={'spec001_mr':da, 'WD_spec001':da, 
                'DD_spec001':da,  
                'RELCOM':relcom, 'RELLNG1':rellng,
                'RELLAT1':rellat, 'RELZZ1':relz1, 'RELZZ2':relz2,'RELPART':relpart},
            coords={'longitude':longitude, 'time':time_var,
                    'latitude':latitude, 'height':height})

    ds = ds.assign(LAGE=int(pd.Timedelta(pd.to_timedelta('5d'),unit='D').value))
    ds.attrs['loutstep']=-10800
#     ds = ds.reindex(time=ds.time[::-1])
    t1 = ds.time[0]

    ds.attrs['iedate']=str(t1.dt.strftime('%Y%m%d').values)
    ds.attrs['ietime']=str(t1.dt.strftime('%H%M%S').values)
    t0=ds.time[-1]
    ds.attrs['ibdate']=str(t0.dt.strftime('%Y%m%d').values)
    ds.attrs['ibtime']=str(t0.dt.strftime('%H%M%S').values)
    ds.attrs['ind_receptor'] = ind_receptor

    return ds
def test_region_slicing():
    fds= create_flexdust_test_data(seed=None)
    fpds = create_test_data(seed=None)
    ds_orr,surface_sensitivity, pre_ds = process_per_pointspec(fpds,fds,x0=76, x1=80, y0=31, y1=34,height=100)
    ns_ds_orr,ns_surface_sensitivity, ns_pre_ds = process_per_pointspec(fpds,fds,
                        x0=None, x1=None, y0=None,y1=None,height=100)
    ns_pre_ds = ns_pre_ds.sel(lon=slice(x0=76,x1=80), lat=slice(y0=31,y1=34))
    assert assert_approx_equal(pre_ds.sum(), ns_pre_ds.sum())

if __name__ == "__main__":
    fds= create_flexdust_test_data(seed=None)
    fpds = create_test_data(seed=None)
