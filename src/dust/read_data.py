import os
import pandas as pd
import numpy as np
import glob
import xarray as xr
import dask
from .utils.read_utils import (
    read_command_namelist,
    read_outGrid_namelist,
    read_release_namelist,
    read_flex_dust_summary,
)
from .utils.utils import _fix_time_flexdust
from functools import partial

"""
This file contain functions for preparing data for analysis

Standard conventions:
    "time" : allways forward time
    "btime" : time along backward trajectory


"""

xr.set_options(keep_attrs=True)


def read_multiple_flexpart_outputs(
    path,
    data_Vars="spec001_mr",
    time_step=None,
    ldirect=-1,
    height=None,
    location=None,
    parallel=False,
    **dset_kwargs
):
    """
    DESCRIPTION
    ===========
        Reads in FLEXPART output files in a parent directory,
        if there is an AVAILABLE_OUTPUT file it reads the path from there
        recursively looks for netcdf files with *.nc file extension

    USAGE
    =====
        dsets = read_multiple_flexpart_outputs(path, data_Vars, time_step, ldirect,
                                        height, location,**dset_kwargs)
        path : path to top directory containing the FLEXPART model output (simulation per time step)
        data_Vars : which variable to concatenate
        time_step : the timestep between each new FLEXPART simulation, default determined by the filename
        ldirect : wether the FLEXPART is run in forward or backward mode (default -1)
        height : which output height to concatenate, if None the output heights are summed, default None.
        location : the location of the receptor point, indexed by integer
        **dset_kwargs : **kwargs to xr.open_mfdataset()

        return : xarray dataset



        todo:: I don't know really whats the best way to concatenated netCDF files using xarray
                So, currently this fuction feels very slow and not particularly robust, but it might
                work fine for concatinating few files. At the moment it is better to first run concat_output.py
                , which create the a concatinated netCDF file of the output with the correct formatting and
                then to the analysis of the output.

    """

    if isinstance(path, list) == True:
        nc_files = path
    else:
        try:
            df = pd.read_csv(path + "AVAILABLE_OUTPUT", index_col=0)
            nc_files = [
                path + row["dir_paths"] + "/" + row["ncfiles"]
                for index, row in df.iterrows()
            ]
        except FileNotFoundError:
            nc_files = glob.glob(
                path + "**/output/grid*.nc", recursive=True
            )  # recursively find FLEXPART output files
    if time_step == None:
        # determine timestep of forward time dimmension based on file name
        time_step = int(
            int(nc_files[1].split("/")[-1].split("_")[-1][9])
            - int(nc_files[0].split("/")[-1].split("_")[-1][9])
        )
    else:
        time_step = time_step

    prep_func = partial(prepare_flexpart_dataset, dataVars=data_Vars, ldirect=ldirect)
    if isinstance(data_Vars, list) == False:
        data_Vars = [data_Vars]
    dsets = xr.open_mfdataset(
        nc_files,
        concat_dim="time",
        decode_times=False,
        data_vars=data_Vars,
        combine="nested",
        parallel=parallel,
        preprocess=prep_func,
        **dset_kwargs
    )
    if height != None:
        dsets = dsets.sel(height=height)
    else:
        dsets = dsets.sum(dim="height", keep_attrs=True)
    if location != None:
        dsets = dsets.sel(pointspec=location, numpoint=location)
    t_index = pd.Index(range(0, len(dsets.time) * time_step, time_step))
    dsets = dsets.assign_coords(time=t_index)
    dsets.attrs["source"] = (
        dsets.attrs["source"]
        + ", concatenated by DUST.read_data.read_multiple_flexpart_output"
    )
    dsets.time.attrs["units"] = "hours since {}".format(
        pd.to_datetime(dsets.iedate + dsets.ietime).strftime("%Y-%m-%d %H:%M")
    )
    dsets = xr.decode_cf(dsets, decode_times=True)
    if "RELCOM" in dsets.data_vars:
        dsets["RELCOM"] = dsets["RELCOM"].astype("U35")
    else:
        # Make sure RELCOM does not disappear!
        ds = xr.open_dataset(nc_files[0]).sel(numpoint=location, pointspec=location)
        relcom = ds.RELCOM.copy().astype("U35")
        relcom_dis = ds.RELCOM.attrs
        dsets = dsets.assign(RELCOM=relcom)
        dsets["RELCOM"].attrs = relcom_dis
    return dsets


def read_flexpart_metadata(path_output_folder):

    """
    DESCRIPTION
    ===========

        This should be rewritten so that it also can concatinated  or maybe I should
        just remove this since xarray.open_dataset does better? I could have one function
        called read metadata and one called read trajectories.

        Takes in the path to a flexpart output directory or a FLEXPART output file.
        If a path to a output directory is provided then it will read all the FLEXPART
        output files present in that directory and return it as a python directory. If a full path
        to a FLEXPART output file is provided then it will only read that specific output file.


    USAGE
    =====

        flexpart_output : str containing path to FLEXPART output directory or a FLEXPART
                          output file.

        fp = read_flexpart(flexpart_output)

        return : python directory / pandas.dataframe / xarray.dataset

    """
    outs = {}
    if os.path.isfile(path_output_folder):
        path = path_output_folder.split("/")
        files = [path[-1]]
        out_dir = str()
        for s in path[:-1]:
            out_dir += s + "/"
    elif os.path.isdirectory(path_output_folder):
        files = os.listdir(path_output_folder)
        if path_output_folder.endswith("/"):
            out_dir = path_output_folder
        else:
            out_dir = path_output_folder + "/"

    else:
        raise (
            ValueError(
                "flexpart_output has to be a string not {}".format(
                    type(path_output_folder)
                )
            )
        )

    for f in files:
        if "COMMAND.namelist" in f:
            com_dict = read_command_namelist(out_dir)
            outs["command"] = com_dict
        elif "RELEASES.namelist" in f:
            rel_df = read_release_namelist(out_dir)
            outs["releases"] = rel_df
        elif "OUTGRID.namelist" in f:
            out_dict = read_outGrid_namelist(out_dir)
            continue
    if len(outs.keys()) == 1:
        key = [key for key in outs.keys()][0]
        return outs[key]
    else:
        return outs


def load_trajectories(path, nclusters=5, skip_rows=24):
    cluster_list = []
    cluster_names = ["xcluster", "ycluster", "zcluster", "fcluster", "rmscluster"]
    for i in range(nclusters):
        for cn in cluster_names:
            cluster_list.append(cn + "(" + str(i) + ")")

    with open(path, "r") as trajecFile:
        header = trajecFile.readline().split(" ")
        trajecFile.readline()
        nlocs = int(trajecFile.readline().strip())
        lines = [next(trajecFile) for i in range(nlocs * 3)]
        locs = [line.strip().split(" ")[0] for i, line in enumerate(lines[2::3])]
        locs_dict = {i + 1: line.strip() for i, line in enumerate(lines[2::3])}

    intial_data = [
        line1.strip() + line2.strip() for line1, line2 in zip(lines[::3], lines[1::3])
    ]

    inital_data = [line.split() for line in intial_data]

    df_head = pd.DataFrame(np.array(inital_data)[:, 1:4])
    df_head = df_head.rename(columns={1: "time", 2: "lon", 3: "lat"})

    df_head["locations"] = locs

    s_time = pd.to_datetime(header[0] + header[1])

    cols = [
        "location",
        "time",
        "lon",
        "lat",
        "height",
        "mean topography",
        "mean mixing height",
        "mean tropopause height",
        "mean PV index",
        "rms distance",
        "rms",
        "zrms distance",
        "zrms",
        "fraction mixing layer",
        "fraction PV<2pvu",
        "fraction in troposphere",
    ] + cluster_list

    df = pd.read_csv(path, sep="\s+", skiprows=lambda x: x < skip_rows, names=cols)

    for key, location in locs_dict.items():
        df.loc[df.loc[:, "location"] == key, "location"] = location

    df.loc[:, ["location"]] = df.loc[:, ["location"]].astype("category")
    #     time_p_rel = s_time + pd.to_timedelta(sec_p_rel, unit = 's')
    df["start time"] = s_time
    return df, df_head


def read_flexpart_trajectories(path_to_top_directory, nclusters=5, skip_rows=24):
    """
    DESCRIPTION
    ===========
        Read all trajectories.txt files present the top folder and subfolders.
        The concatinate the data return return a dictionary containing the recptor
        location and pandas DataFrame containing the trajectory information.

    USAGE:
    =====
        trajectories, locations = read_flexpart_trajectories(path, nclusters)

        nclusters : how many clusters does the trajectory file contain? (default = 5)
        path_to_top_directory : path top directory containing FLEXPART model output

    """
    if isinstance(path_to_top_directory, list) == True:
        paths = path_to_top_directory
    else:
        try:
            df = pd.read_csv(path_to_top_directory + "AVAILABLE_OUTPUT", index_col=0)
            paths = [
                path_to_top_directory + row["dir_paths"] + "/trajectories.txt"
                for index, row in df.iterrows()
            ]
        except FileNotFoundError:
            paths = glob.glob(
                path_to_top_directory + "**/trajectories.txt", recursive=True
            )  # recursively find FLEXPART output files

    locations = load_trajectories(paths[0], skip_rows=skip_rows)[1]
    trajectories = [
        load_trajectories(path, nclusters, skip_rows=skip_rows)[0] for path in paths
    ]
    return pd.concat(trajectories, axis=0, ignore_index=True), locations


def read_flexdust_output(path_output, **xarray_kwargs):
    """
    DESCRIPTION:
    ===========
        Takes in a path to a FLEXDUST output folder/output file.
        If path directly to the FLEXDUST output file is provided then xarray.dataset is returned
        otherwise return a python dictionary, contaning the information in the summary file and
        a xarray.dataset

    ARGUMENTS:
    ==========
        path_output: str
            Path to where the FLEXDUST output is stored.d

    USAGE:
    =====
        fd = read_flexdust_output(path_output)

        returns : python dictionart/xarray.dataset

    """
    outdir = {}
    if path_output.endswith(".nc"):
        outdir = _fix_time_flexdust(path_output, **xarray_kwargs)
    else:
        for output_file in os.listdir(path_output):
            if output_file.endswith(".nc"):
                ncfile = path_output + output_file

                outdir["dset"] = _fix_time_flexdust(ncfile)
            elif output_file.endswith(".txt"):
                summary_dir = read_flex_dust_summary(path_output + output_file)
                outdir["Summary"] = summary_dir
            else:
                continue
    # if len(outdir.keys()) == 1:
    #    key = [key for key in outdir.keys()][0]
    #    return outdir[key]
    # else:
    return outdir


def read_flexpart_output(path_output, dataVars="spec001_mr", ldirect=-1, **dset_kwargs):
    """
    DESCRIPTION
    ===========

        Read flexpart netcdf output file and prepare the dataset for further analysis.

    ARGUMENTS:
    ==========
        path_output : str
            path to netCDF file containing flexpart output

    KEYWORD ARGUMENTS:
    =================
        dataVars: str,
            name of datavariable to read
        ldirect: str,
            which direction the FLEXPART is run
        dset_kwargs: dict
            dictionary of key word arguments for xr.open_dataset()



    USAGE
    =====

        dset = read_flexpart_output(path_output)

    """

    dset = xr.open_dataset(path_output, decode_times=False, **dset_kwargs)
    dset = prepare_flexpart_dataset(dset, dataVars, ldirect)
    return dset


def prepare_flexpart_dataset(dset, dataVars="spec001_mr", ldirect=-1):
    """
    DESCRIPTION
    ===========
        Base function cleaing the flexpart output and preparing the data for
        further processing.

    """

    if ldirect == -1:

        usefull = [
            "RELCOM",
            "RELLNG1",
            "RELLNG2",
            "RELLAT1",
            "RELLAT2",
            "RELZZ1",
            "RELZZ2",
            "RELKINDZ",
            "RELSTART",
            "RELEND",
            "RELPART",
        ]
        usefull.append(dataVars)
        not_usefull = [v for v in dset.data_vars if v not in usefull]
        dset = dset.drop(not_usefull)
        dset = dset.rename({"time": "btime"})
        dset = dset.assign_coords(btime=(dset.btime / 3600).astype(np.short))
        dset.btime.attrs["units"] = "hours"
        dset.btime.attrs["long_name"] = "time along back trajectory"
        dset = dset.assign_attrs({"varName": dataVars})
    dset = dset.rename({"longitude": "lon", "latitude": "lat"})

    dset = dset.squeeze()

    return dset


def read_data(path_output, **dset_kwargs):
    """
    DESCRIPTION
    ===========

        Function reading datafiles which have already been formated in the correct datafromat
        This is just a wrapper of xarray.open_dataset()
    """

    return xr.open_dataset(path_output, **dset_kwargs)


def read_grib_data(path, shortName=None, **dset_kwargs):
    """
    DESCRIPTION
    ===========
        Reading grib files for analysing model forcing.

    """
    if shortName == None:
        ds_grib = xr.open_dataset(
            path,
            engine="cfgrib",
            backend_kwargs={
                "filter_by_keys": {"typeOfLevel": "heightAboveGround", "edition": 1}
            },
        )
    else:
        ds_grib = xr.open_dataset(
            path,
            engine="cfgrib",
            backend_kwargs={"filter_by_keys": {"shortName": shortName}},
        )
        ds_grib = ds_grib.assign_attrs({"varName": shortName})

    return ds_grib


def read_wetprecip_data(path):
    """
    DESCRIPTION
    ===========
        Reads FLEXPART precip data wetscav_precip.txt
        returns  pandas dataframe

    """
    with open(path, 'r') as f:
        data = f.readlines()
        combined = []
        for p1, p2 in zip(data[::2], data[1::2]):
            combined.append(p1[:-1]+p2[:-1])
    combined = pd.DataFrame(combined)
    d = combined[0].str.extract(r'(?P<date>\d{8})(\s*\d{1}[1-9]{0,2})', expand=False)
    combined.index=pd.to_datetime(d['date'].str.cat(d[1].str.strip(), sep='-'), format='%Y%m%d-%H')
    combined = combined.drop_duplicates()
    combined = combined[0].str.extract(r'(?P<lsp>\s\d[.]\S\d*\S+\d{2})(?P<convprec>\s+\d[.]\S\d*\S+\d{2})')
    combined = combined.astype(float)
    combined = combined.sort_index()
    return combined