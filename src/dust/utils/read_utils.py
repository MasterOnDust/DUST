"""
Different function for reading flexpart/flexdust 
text files.

Works only for FLEXPART V10 or newer

AUTHOR: 
======

    Ove Haugvaldstad (ovehaugv@outlook.com)


"""
import pandas as pd


def read_command_namelist(output_dir):

    """
    DESCIRPTION:
    ===========

        Reads the Command namelist in the flexpart output directory
        and return the content as a dictionary

    USAGE:
    =====

        com_dict = read_command_namelist(output_dir)

            output_dir : directory containing FLEXPART output

            returns: python dictionary
    """
    comDict = {}
    lines = open(output_dir + "COMMAND.namelist", "r")
    lines.readline()
    for line in lines:
        var_val = line.split("=", 1)
        try:
            values = var_val[1].split(",", 1)[0].strip()
        except IndexError:
            continue
        var = var_val[0].strip()
        comDict[var] = values
    return comDict


def read_release_namelist(output_dir):
    """
    DESCRIPTION:
    ===========

        Reads the release namelist file in the output directory
        and return the content as a pandas dataframe

    USAGE:
    =====

        df = read_release_namelist(output_dir)

            output_dir: Path to directory with containing the output folder

            returns: pandas.Dataframe

    """
    relLocactions = pd.DataFrame()
    lines = open(output_dir + "RELEASES.namelist", "r")
    lines.readline()
    nspec = int(lines.readline().split("=", 1)[1].split(",", 1)[0].strip())
    lines.readline()
    lines.readline()
    relDict = {}
    for line in lines:
        if "&RELEASE" in line:
            relDict = {}
        elif "/" in line:
            df = pd.DataFrame(relDict, index=[0])
            relLocactions = relLocactions.append(df, ignore_index=True)
        else:
            var_val = line.split("=", 1)
            values = var_val[1].split(",", 1)[0].strip()
            var = var_val[0].strip()
            relDict[var] = values
    return relLocactions


def read_outGrid_namelist(output_dir):
    """
    DESCRIPTION:
    ===========

        Reads the OUTGRID.namelist file in the output directory and
        returns the data as python dictionary.

    USAGE:
    =====

        outGrid = read_outGrid_namelist(output_dir)

            output_dir: path to FLEXPART output directory.
            returns: python.dictonary

    """
    outGrDict = {}
    lines = open(output_dir + "OUTGRID.namelist")
    lines.readline()
    for line in lines:
        if "/" in line:
            continue
        elif "OUTHEIGHTS" in line:
            var_val = line.split("=", 1)
            values = var_val[1].split(",")[0].strip()
            var = var_val[0].strip()
            outGrDict[var] = values
        else:
            var_val = line.split("=", 1)
            values = var_val[1].split(",", 1)[0].strip()
            var = var_val[0].strip()
            outGrDict[var] = values

    return outGrDict


def read_trajectories(output_dir, nclusters=5):

    """
    DESCRIPTION
    ===========

        Read in the trajectories.txt and returns the data in a pandas dataframe

    USAGE
    =====
        df = read_trajectories(output_dir, nclusters=5)

            output_dir: Path to directory containing the FLEXPART output
            nclusters: The number of clusters used in the simulation. FLEXPART default is 5 (it might be
            written in the Trajectories.txt file)
    """
    cluster_list = []
    cluster_names = ["xcluster", "ycluster", "zcluster", "fcluster", "rmscluster"]
    for i in range(nclusters):
        for cn in cluster_names:
            cluster_list.append(cn + "(" + str(i) + ")")

    cols = [
        "release number",
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
    trajecFile = open(output_dir + "trajectories.txt", "r")
    header = trajecFile.readline().split(" ")
    s_time = pd.to_datetime(header[0] + header[1])
    print(s_time)

    df = pd.read_csv(
        output_dir + "trajectories.txt",
        sep="\s+",
        skiprows=lambda x: x < 24,
        names=cols,
    )

    sec_p_rel = df["time"]
    time_p_rel = s_time + pd.to_timedelta(sec_p_rel, unit="s")
    df.index = time_p_rel
    # df = df.drop(['time'], axis=1)
    return df


def read_flex_dust_summary(path_to_textfile):
    """Read FLEXDUST summary file and return metadata as a python dictionary"""
    text_file = open(path_to_textfile, "r")
    out_dict = {}
    for line in text_file.readlines():
        if len(line.split(":")) == 1:
            continue
        else:
            text_in_line = line.split(":")
            concat_str = ""
            for chunk in text_in_line[1].split(" "):

                if chunk != "":
                    chunk = chunk.strip()
                    concat_str += " " + chunk
            #                 chunks.append(chunk.strip())
            concat_str = concat_str.strip().split(" ")
            if len(concat_str) > 1:
                out_dict[text_in_line[0].strip()] = concat_str
            else:
                out_dict[text_in_line[0].strip()] = concat_str[0]
    return out_dict
