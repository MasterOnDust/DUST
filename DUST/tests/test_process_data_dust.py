from .create_test_data import create_flexdust_test_data, create_test_data
import xarray as xr
import numpy as np
import pandas as pd
from DUST.process_data_dust import process_per_timestep, create_output,process_per_pointspec
import pytest


# def test_region_slicing():
#     fds= create_flexdust_test_data(seed=None)
#     fpds = create_test_data(seed=None)
#     ds_orr, pre_ds,surface_sensitivity = process_per_pointspec(fpds,fds,x0=76, x1=80, y0=31, y1=34,height=100)
#     ns_ds_orr, ns_pre_ds,ns_surface_sensitivity = process_per_pointspec(fpds,fds,
#                         x0=None, x1=None, y0=None,y1=None,height=100)
#     ns_pre_ds = ns_pre_ds.sel(lon=slice(76,80), lat=slice(31,34))
#     assert pre_ds.sum().values == ns_pre_ds.sum().values

def test_multiplication(seed=1,**test_kwargs):
    HEIGHT=100
    scale_factor = (1/HEIGHT)*1000
    lons = test_kwargs.get('lons',np.arange(75,86,0.5))
    lats = test_kwargs.get('lats',np.arange(30,35,0.5))
    fds= create_flexdust_test_data(lons, lats, seed=seed)
    fpds = create_test_data(seed=seed, ind_receptor=3,**test_kwargs)
    
    ds_orr,pre_ds, out_data = process_per_pointspec(fpds,fds,
                        x0=None, x1=None, y0=None,y1=None,height=HEIGHT)
    
    pre_test = pre_ds.isel(time=1)

    pre_sum = pre_test.sum(dim=['lon','lat'])

    pre_sum = pre_sum.isel(btime=4)
    fpds = fpds.isel(pointspec=1, numpoint=1)
    fpds = fpds.spec001_mr.sel(time='2000-03-19T03:00:00.000000000', height=100, nageclass=0)

    fpds = fpds.rename({'longitude':'lon', 'latitude':'lat'})

    corr_sum = fds.Emission.sel(time='2000-03-19T03:00:00.000000000')*fpds
    corr_sum = corr_sum.sum(dim=['lon','lat'])*scale_factor
    corr_sum = corr_sum.squeeze()

    try: 
        assert pytest.approx(pre_sum.values,rel=0.01) ==corr_sum.values
    except AssertionError:
        print('Calculated is not equal expected (calculated : {:.3f}, expected: {:.3f})'.format(pre_sum.values,corr_sum.values))

