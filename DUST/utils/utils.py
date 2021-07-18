
import numpy as np
import xarray as xr
import pandas as pd
import argparse as ap

def arg_parser(description='FLEX parser'):
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('path', help='Path to top directory containing output')
    parser.add_argument('--outpath', '--op', help='Where the output should be stored', default='.')
    parser.add_argument('--locations', '--loc', help='Number of location to apply fuction to', default=None, nargs='+', type=int)
    parser.add_argument('--sdate', '--sd', help='Begining of time slice')
    parser.add_argument('--edate', '--ed', help='End of time slice')
    
    return parser




def set_varName(dset, varName):
    """
    DESCRIPTION
    ============
        Set the attribute which the processing and plotting methods
        will use to acess the data.

    """
    dset = dset.assign_attrs({'varName':varName})
    return dset
def multiply_area(dset, area=None):

    """
    DESCRIPTION:
    ===========
        Muliply gridded data by it's area, if the data set does not 
        contain an area variable, the area dataarray need to be provided. 

    """
    if area == None:
        try:
            area = dset.area.values
        except AttributeError:
            print('dataset does not contain a datavariable called area')

    with xr.set_options(keep_attributes=True):
        area_multiplied = dset[dset.varName]*area
    
    return area_multiplied



def region_slice(dset, x0=None, x1=None
                        , y0=None, y1=None):

    """
    DESCRIPTION:
    ===========
        Returns lon / lat rectangle slice of a xarray data set
        
    USAGE:
    =====
        sliced_dataset = region(dset, x0, x1, y0, y1)

        dset: xarray dataset
        x0: longitude of bottom corner, (default min longitude)
        x1: longitude of top corner, (default max longtitude)
        y0: latitude of bottom corner, (default min latitude)
        y1: latitude of top corner, (default max latitude)

        returns : xarray.dataset

    """
    
    if x0 == None: 
        x0 = dset.lon.min()
    if x1 == None:
        x1 = dset.lon.max()
    if y0 == None:
        y0 = dset.lat.min()
    if y1 == None:
        y1 = dset.lat.max()
    return dset.where((((dset.lon >= x0) & (dset.lon <= x1)) & 
                        ((dset.lat >= y0) & (dset.lat <= y1))), drop=True)




def _fix_time_flexdust(ncfile,**xarray_kwargs):
    """Fixes the time in FLEXDUST"""
    dset = xr.open_dataset(ncfile, decode_times=False,**xarray_kwargs)
    if dset.attrs.get('Date',None) == None:
        dset = xr.decode_cf(dset)
    else:
        s_date = dset.startdate.values 
        s_hour = dset.starthour.values
        s_dT =  pd.to_timedelta(s_hour,unit='h') 
        sTime = pd.to_datetime(s_date, format='%Y%m%d') + s_dT 
        time_index = np.unique(np.reshape(dset.Date.values,dset.Date.shape[0]*2))
        time_freq = int((time_index[1]- time_index[0])/60/60)
        nTimeSteps = len(time_index)-1
        time_index = pd.date_range(start='{}'.format(sTime.strftime('%Y%m%d %H:%M:%S').values[
            0]), periods=nTimeSteps, freq='{}h'.format(time_freq))
        dset['time'] = time_index
        dset = dset.drop_dims('time_s')
    return dset

