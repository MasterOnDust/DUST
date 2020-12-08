import xarray as xr
import argparse as ap
import dask
import pandas as pd
import time
import sys
import os 
def pre(ds):
    ds = ds.drop(['RELEND', 'RELSTART', 'RELCOM'])
    ds = ds.sum(dim='btime').mean(dim='time')
    ds = ds.assign_attrs(dict(RELCOM=' '.join(str(ds['RELCOM'][0].values).strip().split(' ')[1:])))
    return ds


if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Create monthly mean source contribution')
    parser.add_argument('paths_flexpart', help='paths to flexpart output')
    parser.add_argmuent('--outpath', '--op', help='path to where output should be stored')
    
    args = parser.parse_args()

    paths = args.paths_flexpart
    outpath = args.outpath
    paths.sort()
    terminal_input = ' '.join(sys.argv)
    dsets = xr.open_mfdataset(paths,preprocess=pre,data_vars=['drydep','surface_sensitivity', 'RELPART'], 
                                concat_dim='time', combine='nested', parallel=True)
    times = pd.to_datetime([path.split('-')[-1][:-3] for path in paths])
    times.freq = 'M'

    dsets = dsets.assign_coords(time=times)

    dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),terminal_input) + dsets.attrs['history']
    relcom = '_'.join(dsets['RELCOM'].splt(' '))
    outFile_name = os.path.join(outpath,relcom + '_' +  'monthly' + '_' + '_' + times[0].srtftime('%Y%m%d') +'-'+ times[-1].srtftime('%Y%m%d') + '.nc')
    print('writing to {}'.format(outFile_name))
    dsets.to_netcdf(outFile_name)


