"""
functions for reading FLEXPART and FLEXDUST output

Works only for FLEXPART V10 or newer

AUTHOR: 
======

    Ove Haugvaldstad (ovehaugv@outlook.com)


"""
import pandas as pd
import xarray as xr
from IPython import embed 

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
def read_command_namelist(output_dir):
    comDict = {}
    lines = open(output_dir + '/COMMAND.namelist', 'r')
    lines.readline()
    for line in lines:
        var_val = line.split('=', 1)
        try:
            values = var_val[1].split(',',1)[0].strip()
        except IndexError:
            continue
        var = var_val[0].strip()
        comDict[var] = values
    return comDict
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
def read_release_namelist(output_dir):
    relLocactions = pd.DataFrame()
    lines = open(output_dir + 'RELEASES.namelist', 'r')
    lines.readline()
    nspec = int(lines.readline().split('=', 1)[1].split(',',1)[0].strip())
    lines.readline()
    lines.readline()
    for line in lines:
        if '&RELEASE' in line:
            relDict = {}
        elif '/' in line:
            df = pd.DataFrame(relDict, index=[0])
            relLocactions = relLocactions.append(df, ignore_index=True) 
        else:
            var_val = line.split('=', 1)
            values = var_val[1].split(',',1)[0].strip()
            var = var_val[0].strip()
            relDict[var] = values
    return relLocactions

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
def read_outGrid_namelist(output_dir):
    outGrDict = {}
    lines = open(output_dir + 'OUTGRID.namelist')
    lines.readline()
    for line in lines:
        if '/' in line:
            continue
        elif 'OUTHEIGHTS' in line:
            var_val = line.split('=', 1)
            values = var_val[1].split(',')[0].strip()
            var = var_val[0].strip()
            outGrDict[var] = values
        else:
            var_val = line.split('=', 1)
            values = var_val[1].split(',',1)[0].strip()
            var = var_val[0].strip()
            outGrDict[var] = values
        
    return outGrDict
        

"""
DESCRIPTION
===========

    Read in the trajectories.txt and returns the data in a pandas dataframe

USAGE
=====
    df = read_trajectories(output_dir, nclusters=5)

        output_dir: Path to directory containing the FLEXPART output
        nclusters: The number of cluster used in the simulation. (it might be 
        written in the Trajectories.txt file)
"""
def read_trajectories(output_dir, nclusters=5):
    cluster_list = []
    cluster_names = ['xcluster', 'ycluster', 'zcluster', 'fcluster',
    'rmscluster']
    for i in range(nclusters):
        for cn in cluster_names:
            cluster_list.append(cn + '(' +str(i)+ ')')


    cols = ['release number', 'time', 'lon', 'lat',
         'height', 'mean topography',
         'mean mixing height', 'mean tropopause height', 'mean PV index',
         'rms distance', 'rms', 'zrms distance', 'zrms',
         'fraction mixing layer', 'fraction PV<2pvu',
         'fraction in troposphere'] + cluster_list
    trajecFile = open(output_dir + 'trajectories.txt', 'r')
    header = trajecFile.readline().split(' ')
    s_time = pd.to_datetime(header[0] + header[1])

    df = pd.read_csv(output_dir + 'trajectories.txt', sep='\s+',
                    skiprows=lambda x: x <24, names=cols)
    
    sec_p_rel = df['time']
    time_p_rel = s_time + pd.to_timedelta(sec_p_rel, unit = 's')
    df.index = time_p_rel
    df = df.drop(['time'], axis=1)
    return df


if __name__ == "__main__":
    df = read_trajectories('/opt/uio/flexpart/Compleated_runs/20190306_15/output/')