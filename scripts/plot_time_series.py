import xarray as xr
from dask.distributed import LocalCluster, Client
import DUST
import argparse as ap
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import glob

def create_time_series(path_2micron,path_20micron, outpath='.',seasonal=True,
                       e_time=None,s_time=None,to_netcdf=False,
                       source_strenght2m = 0.08, source_strenght20m = 0.03)):

    date_slices = []

    d0 = xr.open_dataset(path_2micron)
    relcom = d0.receptor_name.split()
    loc_name = relcom[0]
    data_var = d0.srr.var
    d0.close()
    if outpath.endswith('/'):
        pass
    else:
        outpath = outpath +'/'
    outpath = outpath + '{}_{}_2019'.format(data_var,loc_name)
    try:

        os.mkdir(outpath)
    except FileExistsError:
        shutil.rmtree(outpath)
        os.mkdir(outpath)

    source_strenghts=[source_strenght2m, source_strenght20m]



    p_sizes = ['2 $\mu m$', '20 $\mu m$']
    if seasonal == True:
        date_slices.append(slice('2019-03-01','2019-05-31'))
        date_slices.append(slice('2019-06-01', '2019-08-31'))
        date_slices.append(slice('2019-09-01','2019-10-31'))

    else:
        if e_time == None:
            e_time = pd.to_datetime(d0.time[-1].values).strftime('%Y-%m-%d')
        if s_time == None:
            s_time = pd.to_datetime(d0.time[0].values).strftime('%Y-%m-%d')

        date_slices.append(slice(e_time,s_time))

    dsets = []
    for path, source_strenght in zip([path_2micron, path_20micron], source_strenghts):
        dset = xr.open_dataset(path, chunks={'time': 20})
        dset = dset.srr.make_time_seires().persist()
        with xr.set_options(keep_attrs=True):
            dset = dset.assign({data_var:dset[data_var]*source_strenght})
        if to_netcdf:
            dset.to_netcdf(outpath +'/{}_{}_2019'.format(data_var, "_".join(dset.receptor_name.split())) + '.nc')
        dsets.append(dset)

    for date_slice in date_slices:
        fig, (ax,ax1) = plt.subplots(2,1,sharex=True, sharey=True, figsize=(16,5))

        for dset, ax_i, color, p_size in zip(dsets,(ax,ax1),['saddlebrown', 'sandybrown'], p_sizes):
            temp_dset = dset.sel(time=date_slice)
            ax_i = temp_dset.srr.plot_time_series(ax=ax_i, label = dset.receptor_name + ' ' + p_size, color=color)
            ax_i.legend(loc='upper right')
            ax_i.grid(linestyle='-')
        ax.axes.xaxis.label.set_visible(False)
        plt.savefig(outpath+'/{}_{}_{}_{}'.format(data_var,loc_name, date_slice.start, date_slice.stop) + '.png'
                                , dpi=300, bbox_inches='tight')
        plt.clf()

    for dset in dsets:
        dset.close()

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Plot timeseries of Concentration/Dry-/Wet deposition')
    parser.add_argument('path_2micron',  help='path to output 2micron')
    parser.add_argument('path_20micron',  help='path to output 20micron')
    parser.add_argument('--etime','--et', default=None, help='end date of timeseries plot')
    parser.add_argument('--stime', '--st', default=None, help='start data of timeseries plot')
    parser.add_argument('--seasonal','--s', action='store_true', help='''create seasonal time series
                        disregard etime and stime''')
    parser.add_argument('--out_dir', '--op' ,help='path to where output should be stored', default='.')
    parser.add_argument('--use_cluster', '--uc',action='store_true')
    parser.add_argument('--to_netCDF', '--to_nc', action='store_true', help = 'enable saving timeseries as netCDF')
    args = parser.parse_args()
    seasonal = args.seasonal
    e_time = args.etime
    s_time = args.stime
    path_2micron = args.path_2micron
    path_20micron = args.path_20micron
    outpath = args.out_dir
    use_cluster = args.use_cluster
    to_netcdf = args.to_netCDF


    if use_cluster == True:
        cluster = LocalCluster(n_workers=32, threads_per_worker=1, memory_limit='16GB')
        client= Client(cluster)
        print(cluster)

    files_2m = glob.glob(path_2micron +'**.nc', recursive=True)
    files_20m = glob.glob(path_20micron + '**.nc', recursive=True)

    path_dict_2m = {f_name.split('-')[0].split('/')[-1].split('_')[0] : f_name for f_name in files_2m}
    path_dict_20m = {f_name.split('-')[0].split('/')[-1].split('_')[0] : f_name for f_name in files_20m}

    for key, item in path_dict_2m.items():
        create_time_series(item,path_dict_20m[key], outpath=outpath,seasonal=seasonal,
                       e_time=e_time,s_time=s_time,to_netcdf=to_netcdf)


