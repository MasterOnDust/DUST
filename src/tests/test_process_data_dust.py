from .create_test_data import create_flexdust_test_data, create_test_data
import xarray as xr
import numpy as np
import pandas as pd
from DUST.process_data_dust import process_per_timestep, create_output,process_per_pointspec
import pytest



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
    pre_sum = pre_test
    pre_sum = pre_sum.isel(btime=14)
    test_ds = xr.decode_cf(pre_sum.to_dataset(name='drydep'))

    t_time = test_ds.time + test_ds.btime
    pre_sum = pre_sum.sel(lon=slice(76,77),lat=slice(30,31))

    fpds = fpds.isel(pointspec=1, numpoint=1)
    fpds = fpds.spec001_mr.sel(time=t_time, height=100, nageclass=0)

    fpds = fpds.rename({'longitude':'lon', 'latitude':'lat'})

    corr_sum = fds.Emission.sel(time=t_time)*fpds
    corr_sum = corr_sum.sel(lon=slice(76,77),lat=slice(30,31)).squeeze()*scale_factor
    np.testing.assert_allclose(pre_sum.values, corr_sum.values)


