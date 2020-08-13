import DUST
from DUST import read_flexdust_output, read_flexpart_output, read_multiple_flexpart_output
import xarray as xr
import pytest
from DUST.utils.multiply_emsfield import multi_flexpart_flexdust

def test_multiply_emsfield():
    flexdust = xr.open_dataset('test_files/FLEXDUST_emission_4th-11th_may_2019.nc')
    ref_emsens = read_flexpart_output('test_files/flexpart/20190510_00/output/grid_drydep_20190510000000.nc')

    ref_ems = flexdust.sel(time = '2019-05-05 12:00:00').Emission
    ref_emsens = ref_emsens.sel(btime=-388800)
    ref_emsens = ref_emsens.spec001_mr
    ref_emsens = ref_emsens.sel(pointspec=0, nageclass=0, height=100)
    multi_ems = ref_emsens.values * ref_ems.values
    ref_sum = multi_ems.sum()
    ncfiles = ['test_files/flexpart/20190510_00/output/grid_drydep_20190510000000.nc']
    fnames = multi_flexpart_flexdust('test_files',nc_files=ncfiles,
                                flexdust=flexdust, point_spec=0)
    multi_plied = xr.open_dataset(fnames)
    multi_plied = multi_plied.sel(btime=-388800)
    test_sum = multi_plied.DryDep.sum(dim=['lon','lat'])

    assert test_sum.values == ref_sum

