import pandas as pd

import numpy as np
import sys
import os
import glob

import xarray as xr
import xarray

from .utils.utils import multiply_area


"""
This module contain all functionalilty working on the data.

"""

# Design interface first then, create the modeling functions!


def get_total(dset, unit="kg"):
    tot_emissions = dset.sum(dim="time")
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
        filename = (
            "Emissions"
            + df.index[0].strftime(format="%b_%d_%Y")
            + df.index[-1].strftime(format="%b_%d_%Y")
            + ".csv"
        )
    else:
        filename = filename
    df.to_csv(filename)


def resample_data(dset, freq, method="mean"):
    if method == "mean":
        dset = dset.resample(time=freq).mean()
    elif method == "sum":
        dset = dset.resample(time=freq).sum()
    else:
        raise ValueError("`method` param '%s' is not a valid one." % method)

    return dset
