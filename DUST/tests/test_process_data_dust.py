from .create_test_data import create_flexdust_test_data, create_test_data
import xarray as xr
import numpy as np
import pandas as pd
from DUST.process_data_dust import process_per_pointspec, create_output
import pytest


def test_region_slicing():
    fds= create_flexdust_test_data(seed=None)
    fpds = create_test_data(seed=None)
    ds_orr, pre_ds,surface_sensitivity = process_per_pointspec(fpds,fds,x0=76, x1=80, y0=31, y1=34,height=100)
    ns_ds_orr, ns_pre_ds,ns_surface_sensitivity = process_per_pointspec(fpds,fds,
                        x0=None, x1=None, y0=None,y1=None,height=100)
    ns_pre_ds = ns_pre_ds.sel(lon=slice(76,80), lat=slice(31,34))
    assert pre_ds.sum().values == ns_pre_ds.sum().values

def test_multiplication():
    HEIGHT=100
    POINTSPEC=3
    BTIME=5
    scale_factor = (1/HEIGHT)*1000
    fds= create_flexdust_test_data(seed=None)
    fpds = create_test_data(seed=None)
    
    ds_orr,pre_ds, out_data = process_per_pointspec(fpds,fds,
                        x0=None, x1=None, y0=None,y1=None,height=HEIGHT)
    produced_ds = xr.decode_cf(pre_ds.to_dataset(name='spec001_mr'))
    fpds = fpds.rename({'latitude':'lat','longitude':'lon'})
    
    produced_ds=produced_ds.isel(time=POINTSPEC, btime=BTIME)
    test_time=produced_ds.time+produced_ds.btime
    
    fpds=fpds.sel(time=test_time.values, pointspec=POINTSPEC, height=HEIGHT, nageclass=0)['spec001_mr']
    fds=fds.sel(time=test_time.values)['Emission']
    fpds_times_fds = (fpds*fds).values*scale_factor
    assert pytest.approx(produced_ds['spec001_mr'].sum(dim=['lon','lat']).values,0.01) == fpds_times_fds.sum()
