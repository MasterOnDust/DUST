import pandas as pd

import numpy as np
import sys
import os
import glob

import xarray as xr
import xarray

from IPython import embed

from .utils.utils import multiply_area
from .utils.multiply_emsfield import multi_flexpart_flexdust


"""
This module contain all functionalilty working on the data.

"""

# Design interface first then, create the modeling functions!


def multiply_flexpart_flexdust(flexdust, outpath='./out', locations=None, ncFiles = None, path = None, zlib=True):
    """
    DESCRIPTION
    ===========

        Goes through all subdirectories in path looking for flexpart netcdf files with .nc
        file extension. Then multiply corresponding emission sensitivities with emission field

    USAGE:
    ======

        flexdust           : path to flexdust output folder/output or xarray.Dataset containing flexdust output
        outpath (optional) : directory for where the combined data are going to be stored
        locations(optional): receptor location in flexpart either interger or string containing the name,
                             correspoding to RELCOM in FLEXPART release file. If not provided all locations is used.
        ncFiles (optional) : List of paths to flexpart output files
        path (optional)    : path to top directory of flexpart output which is search recursively looking
                             for FLEXPART netcdf files, either ncFiles or path has to be provided!

        returns: python list containing paths to multiplied flexpart / flexdust output

    """

    if ncFiles:
        ncFiles = ncFiles
    elif  path:
        ncFiles = glob.glob(path + "/**/grid*.nc", recursive=True)
    else:
        raise(NameError("Both path and ncFiles is None"))

    files = []
    if os.path.isdir(outpath):
        pass
    else:
        os.mkdir(outpath)


    d = xr.open_dataset(ncFiles[0])
    for i , com in enumerate(d.RELCOM):
        loc = str(com.values)[2:].strip().split()[0]
        if locations:
            if loc or i in locations:
                files.append(multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib), )
        else:
            files.append(multi_flexpart_flexdust(outpath,ncFiles,flexdust,i, zlib=zlib))


    d.close()
    return files

def get_total(dset, unit='kg'):
    tot_emissions = dset.sum(dim='time')
    return tot_emissions

def time_series_to_df(dset):
    """
    DESCRIPTION
    ===========
        Sums lon lat 

    """

    # emissions_kg = integrate_area('kg')
    # emissions_kg_m2 = self._integrate_area('kg/m2')
    # df = pd.DataFrame(columns=['emissions(kg)','emissions(kg/m2)'], index=emissions_kg.time.values)
    # date0 = np.datetime_as_string(self._obj.time[0].values, unit='D')
    # date_end = np.datetime_as_string(self._obj.time[-1].values, unit='D')
    # df['emissions(kg)'] = emissions_kg.to_dataframe()
    # df['emissions(kg/m2)'] = emissions_kg_m2.to_dataframe()
    # return df
    pass

def emission_time_series_to_csv(self, filename=None):
    df = self.emission_time_series_to_df()

    if filename == None:
        filename = 'Emissions' + df.index[0].strftime(format='%b_%d_%Y') + df.index[-1].strftime(format='%b_%d_%Y') + '.csv'
    else:
        filename =filename
    df.to_csv(filename)



def resample_data(dset, freq, method='mean'):
    if method == 'mean':
        dset =  dset.resample(time=freq).mean()
    elif method =='sum':
        dset = dset.resample(time=freq).sum()
    else:
        raise ValueError("`method` param '%s' is not a valid one." % method)
    
    return dset


