import xarray as xr
import pandas as pd
import re
from netCDF4 import date2num, num2date
import numpy as np


def resample_monthly(dset, squeeze_numpoint=True):
    """
    dset : xarray dataset

    """
    with xr.set_options(keep_attrs=True):
        if dset.ind_receptor == 4 or dset.ind_receptor == 3:
            # Depostion is accumulative
            dset = dset.sum(dim="btime").sum(dim="time")
        else:
            # Concentration is not accumulative
            dset = dset.sum(dim="btime").mean(dim="time")
    if squeeze_numpoint:
        if "numpoint" in dset.dims:
            for data_var in dset.data_vars:
                if "numpoint" in dset[data_var].dims:
                    if data_var == "RELEND":
                        dset[data_var] = dset["RELSTART"][0]
                    elif data_var == "RELSTART":
                        dset[data_var] = dset["RELSTART"][-1]
                    else:
                        dset[data_var] = dset[data_var][0]
        dset = dset.squeeze()
    date_str = pd.to_datetime(re.search("-([0-9]{8})", dset.filename).group()[1:])
    num_date = int(
        date2num(date_str, "hours since 1900-01-01 00:00:00.0", calendar="gregorian")
    )
    dset = dset.assign_coords(time=np.array(num_date))
    dset["time"].attrs = dict(
        units="hours since 1900-01-01 00:00:00.0", calendar="gregorian"
    )
    dset[dset.varName] = dset[dset.varName].expand_dims("time")
    if "surface_sensitivity" in dset.data_vars:
        dset["surface_sensitivity"] = dset["surface_sensitivity"].expand_dims("time")

    return dset


def concatenate_monthly(dsets):
    """
    Takes as input a list of xarray dataset, which has
    already been averaged to monthly temporal resolution and concatenate them
    """

    location_name = dsets[0].attrs["relcom"]
    varName = dsets[0].attrs["varName"]

    for i, dset in enumerate(dsets):
        if dset.attrs["relcom"] != location_name:
            raise (
                ValueError(
                    "The dataset you want to concatenate have to be from the same location {} =! {}".format(
                        dset.attrs["relcom"], location_name
                    )
                )
            )
        if dset.attrs["varName"] != varName:
            raise (
                ValueError(
                    "You are trying to concatenate two different variables {} =! {}".format(
                        dset.attrs["varName"], varName
                    )
                )
            )
        dsets[i] = dset.drop(labels=["RELSTART", "RELEND", "RELPART"])
    out_dset = xr.concat(dsets, "time", data_vars=[varName, "surface_sensitivity"])

    out_dset.attrs["ibdate"] = dsets[0].attrs["ibdate"]
    out_dset.attrs["iedate"] = dsets[-1].attrs["iedate"]
    out_dset.attrs["ibtime"] = dsets[0].attrs["ibtime"]
    out_dset.attrs["ietime"] = dsets[-1].attrs["ietime"]
    out_dset["time"].attrs = dict(
        units="hours since 1900-01-01 00:00:00.0", calendar="gregorian"
    )
    fname = (
        out_dset.attrs["varName"]
        + "_monthly-mean_"
        + "_".join(out_dset.attrs["relcom"])
        + "_"
        + out_dset.attrs["ibdate"][:6]
        + "-"
        + out_dset.attrs["iedate"]
        + ".nc"
    )
    out_dset.attrs["filename"] = fname
    return out_dset
