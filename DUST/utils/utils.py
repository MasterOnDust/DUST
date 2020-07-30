import numpy as np
import xarray as xr
import pandas as pd
# def select_data(ds):
#     """
#     DESCRIPTION
#     ===========

    
#     USAGE:
#     ======



#     """


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


# def merge_flexpart_flexdust(dset,ems, units='kg/m^3'):

#     """
#     DESCRIPTION:
#     ===========
        
#         Calculates surface concentrations / dry deposition / wet deposition 
#         depending on the kind of FLEXPART dataset

#         The FLEXDUST emissions has unit of kg/m^2 which is devided by the height of the lowest layer in the 
#         FLEXPART output

#         The FLEXPART emission sensitivity has unit of either 's', 'm' depending on where concentration or wet dep /
#         dry dep is specified. 

#     USAGE:
#     =====
        
#         da = marge_flexpart_flexdust(dset,ems,units)
#         dset: flexpart DataArray containing the emission sensitivity (xarray DataArray)
#         ems:  flexdust DataArray containing the emission flux output (xarray dataset)
#         units: the final units of the merged product


#         returns xarray.DataArray
#     """


#     da = xr.DataArray(coords=[dset.pointspec, dset.time, dset.latitude, dset.longitude],
#                       dims=['pointspec','time', 'latitude', 'longitude'], name='surf_combined')

#     for pointspec in range(len(dset.pointspec)):
#         for time in range(len(dset.time)):

#             da[pointspec,time,:,:] = (dset[0,pointspec,time,0,:,:]*(ems[time,:,:]/dset.height[0]).values)
                
    

#     return da

def _fix_time_flexdust(ncfile):
    """Fixes the time in FLEXDUST"""
    dset = xr.open_dataset(ncfile, decode_times=False)
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
    return dset

def _gen_log_clevs(dat_min, dat_max):
    """Creates a logarithmic color scale."""

    if dat_max > 0:
        dmx = int(np.round(np.log10(dat_max)))
    else:
        dmx = 1

    # TODO: What's the default value of dmn?
    if dat_min > 0:
        dmn = int(np.round(np.log10(dat_min)))
    elif dat_min == 0. or np.isnan(dat_min):
        dmn = dmx - 3

    # create equally spaced range
    # ERROR: dmn could be uninitialized
    if dmx == dmn:
        dmx = dmn + 1
    clevs = np.logspace(dmn, dmx, 100)

    return clevs


def _gen_flexpart_colormap(ctbfile=None, colors=None):
    """Generate the ast colormap for FLEXPART."""

    from matplotlib.colors import ListedColormap
    if ctbfile:
        try:
            colors = np.loadtxt(ctbfile)
        except:
            print("WARNING: cannot load ctbfile. using colors")
    if colors:
        name = 'user_colormap'
    else:
        # AST Colorset for FLEXPART
        colors = [
            1.0000000e+00, 1.0000000e+00, 1.0000000e+00,
            9.9607843e-01, 9.1372549e-01, 1.0000000e+00,
            9.8431373e-01, 8.2352941e-01, 1.0000000e+00,
            9.6470588e-01, 7.1764706e-01, 1.0000000e+00,
            9.3333333e-01, 6.0000000e-01, 1.0000000e+00,
            8.9019608e-01, 4.4705882e-01, 1.0000000e+00,
            8.3137255e-01, 2.0000000e-01, 1.0000000e+00,
            7.5686275e-01, 0.0000000e+00, 1.0000000e+00,
            6.6274510e-01, 0.0000000e+00, 1.0000000e+00,
            5.4901961e-01, 0.0000000e+00, 1.0000000e+00,
            4.0784314e-01, 0.0000000e+00, 1.0000000e+00,
            2.4705882e-01, 0.0000000e+00, 1.0000000e+00,
            7.4509804e-02, 0.0000000e+00, 1.0000000e+00,
            0.0000000e+00, 2.8235294e-01, 1.0000000e+00,
            0.0000000e+00, 4.8627451e-01, 1.0000000e+00,
            0.0000000e+00, 6.3137255e-01, 1.0000000e+00,
            0.0000000e+00, 7.4509804e-01, 1.0000000e+00,
            0.0000000e+00, 8.4705882e-01, 1.0000000e+00,
            0.0000000e+00, 9.3725490e-01, 1.0000000e+00,
            0.0000000e+00, 1.0000000e+00, 9.7647059e-01,
            0.0000000e+00, 1.0000000e+00, 8.9411765e-01,
            0.0000000e+00, 1.0000000e+00, 8.0000000e-01,
            0.0000000e+00, 1.0000000e+00, 6.9019608e-01,
            0.0000000e+00, 1.0000000e+00, 5.6470588e-01,
            0.0000000e+00, 1.0000000e+00, 4.0000000e-01,
            0.0000000e+00, 1.0000000e+00, 0.0000000e+00,
            3.9607843e-01, 1.0000000e+00, 0.0000000e+00,
            5.6470588e-01, 1.0000000e+00, 0.0000000e+00,
            6.9019608e-01, 1.0000000e+00, 0.0000000e+00,
            7.9607843e-01, 1.0000000e+00, 0.0000000e+00,
            8.9411765e-01, 1.0000000e+00, 0.0000000e+00,
            9.7647059e-01, 1.0000000e+00, 0.0000000e+00,
            1.0000000e+00, 9.4509804e-01, 0.0000000e+00,
            1.0000000e+00, 8.7450980e-01, 0.0000000e+00,
            1.0000000e+00, 7.9215686e-01, 0.0000000e+00,
            1.0000000e+00, 7.0588235e-01, 0.0000000e+00,
            1.0000000e+00, 6.0392157e-01, 0.0000000e+00,
            1.0000000e+00, 4.8235294e-01, 0.0000000e+00,
            1.0000000e+00, 3.1372549e-01, 0.0000000e+00,
            1.0000000e+00, 0.0000000e+00, 1.4901961e-01,
            1.0000000e+00, 0.0000000e+00, 3.3333333e-01,
            1.0000000e+00, 0.0000000e+00, 4.4705882e-01,
            1.0000000e+00, 0.0000000e+00, 5.3725490e-01,
            1.0000000e+00, 0.0000000e+00, 6.1176471e-01,
            9.7647059e-01, 0.0000000e+00, 6.6666667e-01,
            8.9411765e-01, 0.0000000e+00, 6.6666667e-01,
            7.9607843e-01, 0.0000000e+00, 6.3921569e-01,
            6.9019608e-01, 0.0000000e+00, 5.9215686e-01,
            5.6470588e-01, 0.0000000e+00, 5.0980392e-01,
            3.9607843e-01, 0.0000000e+00, 3.8039216e-01]
        colors = np.reshape(colors, (-1, 3))
        name = 'flexpart_cmap'
    cmap = ListedColormap(colors, name)
    return cmap