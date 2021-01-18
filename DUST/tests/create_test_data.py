import xarray as xr
import numpy as np
import pandas as pd

def create_flexdust_test_data(seed=None):
    time_var = xr.Variable('time',pd.date_range("2000-03-01", periods=80, freq='6H'))
    longitude = xr.Variable('lon',np.arange(75,86,0.5))
    latitude = xr.Variable('lat', np.arange(30,35,0.5))
    da = xr.DataArray(data=np.random.rand(1,10,80,3,10,22),
        dims={'nageclass':1, 'pointspec':10, 'time' : 80, 'height':3, 'lat':10, 'lon':22},
        coords={'time':time_var, 'lat':latitude, 'lon':longitude}
        )
    ds=xr.Dataset({'Emission':da},
                coords={'lon':longitude, 'time':time_var,
                    'lat':latitude})
    
    return ds

def create_test_data(seed=None):
    
    rs = np.random.RandomState(seed)
    time_var = xr.Variable('time',pd.date_range("2000-03-10", periods=50, freq='6H'))
    longitude = xr.Variable('longitude',np.arange(75,86,0.5))
    latitude = xr.Variable('latitude', np.arange(30,35,0.5))
    height = xr.Variable('height',[100, 1000, 5000])
    relcom = xr.Variable('numpoint', ['Test Location' for i in range(10)])
    rellat = xr.Variable('numpoint', [32.3 for i in range(10)])
    rellng = xr.Variable('numpoint', [78.3 for i in range(10)])
    relz1 = xr.Variable('numpoint', [0 for i in range(10)])
    relz2 =xr.Variable('numpoint', [30 for i in range(10)])
    relpart =xr.Variable('numpoint', [3000 for i in range(10)])

    da = xr.DataArray(data=np.random.rand(1,10,50,3,10,22),
        dims={'nageclass':1, 'pointspec':10, 'time' : 50, 'height':3, 'latitude':10, 'longitude':22},
        coords={'time':time_var, 'latitude':latitude, 'longitude':longitude, 'height':height}
        )

    
    ds = xr.Dataset(data_vars={'spec001_mr':da, 'WD_spec001':da, 
                'DD_spec001':da,  
                'RELCOM':relcom, 'RELLNG1':rellng,
                'RELLAT1':rellat, 'RELZZ1':relz1, 'RELZZ2':relz2,'RELPART':relpart},
            coords={'longitude':longitude, 'time':time_var,
                    'latitude':latitude, 'height':height})

    ds = ds.assign(LAGE=pd.to_timedelta('5d','D').value)
    ds.attrs['loutstep']=21600
    ds = ds.reindex(time=ds.time[::-1])
    ds.attrs['iedate']=str(ds.time[0].dt.strftime('%Y%m%d').values)
    ds.attrs['ietime']=str(ds.time[0].dt.strftime('%H%M%S').values)
    #ds.attrs['iedate']=

    return ds

if __name__ == "__main__":
    create_flexdust_test_data(seed=None)
    create_test_data(seed=None)

