import xarray as xr

def resample_monthly(dset):
    """
    dset : List of xarray dataset

    """
    with xr.set_options(keep_attrs=True):
        dset = dset.sum(dim='time').mean(dim='btime') 

    return dset    

def concatenate_monthly(dsets):
    location_name = dsets[0].attrs['relcom']
    varName = dsets[0].attrs['varName']
    for dset in dsets:
        if dset.attrs['relcom']!=location_name:
            raise(ValueError('The dataset you want to concatenate have to be from the same location {} =! {}'.format(dset.attrs['relcom'],
                                                                                                        location_name)))
        if dset.attrs['varName']!=varName:
            raise(ValueError('You are trying to concatenate two different variables {} =! {}'.format(dset.attrs['varName'],
                                                                                                        varName)))

    out_dset = xr.concat(dsets, data_vars=[varName,'surface_sensitivity'])

    return out_dset    