import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import argparse as ap
import DUST.plot.plotting as dplot
import pandas as pd
from matplotlib import rcParams
import os

def process_data(ds, method='mean', area=None):
    time0 = ds.time[0].dt.strftime("%Y")
    time1 = ds.time[-1].dt.strftime("%Y")
    ds.attrs['x0'] = ds.lon.min(); ds.attrs['x1'] = ds.lon.max()
    ds.attrs['y0'] = ds.lat.min(); ds.attrs['y1'] = ds.lat.max()
    ds = ds.sum(dim=['lat','lon'],keep_attrs=True)
    ds.attrs['sdate'] = time0
    ds.attrs['edate'] = time1
    if method == 'mean':
        ds = ds.resample(time='Y').mean(dim='time',keep_attrs=True)
    elif method == 'sum':
        ds = ds.resample(time='Y').sum(dim='time', keep_attrs=True)
    else:
        raise(ValueError("method {} is invalid".format(method)))
    

def read_oni_dfj(path, sdate, edate):
    names=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep','Oct', 'Nov', 'Dec']
    df = pd.read_csv(path, sep=r"\s+", header=None, names=names, 
                    skiprows=[0,1], index_col=0, na_values=-99.99, skipfooter=8, engine='python')
    df_Jan = df['Jan'].loc[sdate:edate|]
    df_Feb = df['Feb'].loc[sdate:edate]
    df_Dec = df['Dec'].loc[sdate:edate]

    mean_oni = []
    for dec, jan, feb in zip(df_Dec, df_Jan, df_Feb):
        
        mean_oni.append((dec+jan+feb)/3)

    return mean_oni




if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Create yearly timeseries from monthly averaged data")
    parser.add_argument("paths",nargs="+", help="path to netcdf file to be processes")
    parser.add_argument("--outpath", "--op", help="path where figure should be stored", default="./")
    parser.add_argument("--title", "--t", default=None)
    parser.add_argument("--sdate", "--sd",default=None, help="start of time series")
    parser.add_argument("--edate", "--ed", default=None, help="end date of time series")
    parser.add_argument('--x0' ,help = 'longitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--y0', help='latitude of lower left corner of grid slice', default=None, type=int)
    parser.add_argument('--x1', help='longitude of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--y1', help='latidute of top right corner of grid slice', default=None, type=int)
    parser.add_argument('--psd', help='fraction of particle size distrubution corresponding to size of particle', 
                            default=1, type=int)
    parser.add_argument('--oni_index', '--oni', help='Include plot of the ONI index, if path to oni data file is given', 
                        default=None)
    parser.add_argument('--area_path', '--ap', help='path to netCDF file containing the area (stored in flexdust output)',
                            default=None)
    parser.add_argument('--method', '--m', help="wether to sum or average spring months", default='mean')
    parser.add_argument('--tag', '--tg', default='', 
                        help="tag to include in the beginning for file name")
    parser.add_argument('--file_extension', '--fe', default='png')
    source_strenght2m = 0.08, source_strenght20m = 0.03
    args = parser.parse_args()
    paths = args.paths
    outpath = args.outpath
    title = args.title
    sdate = args.sdate
    edate = args.edate  
    x0 = args.x0
    x1 = args.x1
    y1 = args.y1
    y0 = args.y0
    oniIndex = args.oni_index
    psd = args.psd
    area_path = args.area_path
    tag = args.tag
    file_extension=args.file_extension
    method = args.method
    area=None
    if area_path:
        area_ds = xr.open_dataarray(area_path, decode_times=False)
        area = area_ds['area']

    dsets = [process_data(xr.open_dataset(path).slice(lon=slice(x0,x1),lat=slice(y0,y1))
                            ,method=method, area=area) for path in paths]
    rcParams.update({'figure.autolayout': True})
    fig,ax =plt.subplot(figsize=(12,5))
    for dset in dsets:
        dplot.plot_emission_time_series(dset,ax=ax,linewidth=3)
    ax.grid()
    ax.set_xticks(dsets[0].time, rotation=40, ha='rights', fontsize=12)

    ax.set_title('lon0 {}, lat0 {}, lon1 {}, lat1 {}'.format(dsets[0].x0,dsets[0].y0, dsets[0].x1, dsets[0].y1))
    if oniIndex:
        ax2 = ax.twinx()
        sdate = dsets[0].attrs['sdate'] ; edate = dsets[0].attrs['edate']
        mean_oni = read_oni_dfj(oniIndex,sdate, edate)
        ax2.plot(dsets[0].time, mean_oni, linestyle = '--', color = 'blue', label='Averaged ONI index DJF')
        ax2.set_ylim(-3,3)
        ax2.set_ylabel('ONI Index SST anomalies $\degree C$')
        ax2.set_ylim(-3,3)
    
    fig.legend()
    outfilename = os.path.join(outpath,tag+"_time-series_{}_{}_{}.{}".format(dsets[0].varName, dsets[0].sdate,
                                                                            dsets[0].edate, file_extension))
    plt.savefig(outfilename, dpi=300, bbox_inches='tight')