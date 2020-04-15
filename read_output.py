"""
functions for reading FLEXPART and FLEXDUST output

Works only for FLEXPART V10 or newer

AUTHOR: 
======

    Ove Haugvaldstad (ovehaugv@outlook.com)


"""
import xarray as xr
import pandas as pd 

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
USAGE
=====
    
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
        
#Need to give Name to columns
"""
      write(unitouttraj,'(i5,i8,2f9.4,4f8.1,f8.2,4f8.1,3f6.1,&
           &5(2f8.3,f7.0,f6.1,f8.1))')&
           &j,itime-(ireleasestart(j)+ireleaseend(j))/2, &
           xcenter,ycenter,zcenter,topocenter,hmixcenter,tropocenter, &
           pvcenter,rmsdist,rms,zrmsdist,zrms,hmixfract,pvfract, &
           tropofract, &
           (xclust(k),yclust(k),zclust(k),fclust(k),rmsclust(k), &
           k=1,ncluster)
"""
def read_trajectories(output_dir, nclusters=5):
    cluster_list = []
    cluster_names = ['xcluster', 'ycluster', 'zcluster', 'fcluster',
    'rmscluster']
    for i in range(nclusters):
        for cn in cluster_names:
            cluster_list.append(cn + '(' +str(i)+ ')')


    cols = ['release number', 'seconds prior to release', 'lon', 'lat',
         'height', 'mean topography',
         'mean mixing height', 'mean tropopause height', 'mean PV index',
         'rms distance', 'rms', 'zrms distance', 'zrms',
         'fraction mixing layer', 'fraction PV<2pvu',
         'fraction in troposphere'] + cluster_list
    print(len(cols))
    print(cluster_list)
    
    print(len(cols))
    df = pd.read_csv(output_dir, sep='\s+',
                    skiprows=lambda x: x <24, names=cols)
    return df

        
