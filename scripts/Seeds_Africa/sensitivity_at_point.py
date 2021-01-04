import argparse as ap
import xarray as xr
import pandas as pd

def get_sensitvity_at_point(paths, fname,df, groupby_key='Mountain', outpath='./out'):
    """
    Need to figure out how it should store the data
    
    """
    
    
    for path in paths:
        ds = xr.open_dataset(path)
        height=ds.height.values
        temp_df = df.drop(fname)
        groupby_loc = temp_df.groupby('Mountain')
        month = int(ds.time.dt.month[0].values)
        out_df = pd.DataFrame(columns=ds.index, index=range(13))
        for index, location in groupby_loc:
            ems_sens = ds.sel(lon=location['lon'], lat=location['lat'], method='nearest').sum(dim='btime').mean(dim='time').values
            out_df.loc[month,index] = ems_sens
        out_df.to_csv(outpath +'/{}_{}_sensitivty_matrix.csv')
            


if __name__=="__main__":
    parser = ap.ArgumentParser(description='Plot monthly mean maps')
    parser.add_argument('paths', nargs='+', help='Paths to merged output folder')
    
    args = parser.parse_args()
    paths = args.paths
    
        #paths = glob.glob(path, recursive=True)
    fname = paths[0].split('/')[-1].split('_')[0]
    paths.sort()
    get_sensitvity_at_point(paths, fname)