import xarray as xr
import argparse as ap
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import DUST.plot.plotting as dplot
import DUST.plot.maps as dmap
import os 

if __name__ == "__main__":
    parser = ap.ArgumentParser(description='create averaged map of source contribution')
    parser.add_argument('path', help='Path to monthly FLEXPART source contribution netcdf file')
    parser.add_argument("--outpath", "--op", help="path where figure should be stored", default="./")
    parser.add_argument("--title", "--t", default=None)
    parser.add_argument("--sdate", "--sd",default=None, help="start of time series")
    parser.add_argument("--edate", "--ed", default=None, help="end date of time series")
    parser.add_argument('--psd', help='fraction of particle size distrubution corresponding to size of particle, 2micron 0.08, 20mciron 0.03', 
                            default=1, type=float)

    parser.add_argument('--tag', '--tg', default='', 
                        help="tag to include in the beginning for file name")
    parser.add_argument('--file_extension', '--fe', default='png')
    parser.add_argument('--varName', '--vn', help='Variable name in NETCDF4', default=None)
    args = parser.parse_args()
    path = args.path
    outpath = args.outpath
    title = args.title
    sdate = args.sdate
    edate = args.edate  
    psd = args.psd
    tag = args.tag
    file_extension=args.file_extension
    var_Name = args.varName

    ds = xr.open_dataset(path)
    ds = ds.sel(time=slice(sdate,edate))
    time0 = str(ds.time[0].dt.strftime("%Y").values)
    time1 = str(ds.time[-1].dt.strftime("%Y").values)
    ds.attrs['sdate'] = time0; ds.attrs['edate'] = time1
    with xr.set_options(keep_attrs=True):
        ds = ds.mean(dim="time")
        ds[ds.varName] = ds[ds.varName]*psd
    fig, ax = plt.subplots(figsize=(12,10), subplot_kw={'projection':ccrs.PlateCarree()})
    ax = dmap.map_terrain_china(ax)
    info_dict = {'Start Year': time0, 'End Year': time1}
    ax = dplot.plot_emission_sensitivity(ds, var_Name=var_Name, ax=ax, title=title, info_dict=info_dict)
    outfilename = os.path.join(outpath,tag+"_mean_source_contrib_{}_{}_{}.{}".format(ds.varName, ds.sdate,
                                                                            ds.edate, file_extension))
    plt.savefig(outfilename, dpi=300, bbox_inches='tight')


