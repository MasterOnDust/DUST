import DUST
from DUST import read_flexdust_output, read_flexpart_output, read_multiple_flexpart_output
import xarray as xr
import unittest

class TEST_FP_METHODS(unittest.TestCase):

    def __init__(self):
        super().__init__()
        path = 'test_files/flexpart/20190510_00/output/grid_drydep_20190510000000.nc'
        ds_test = read_flexpart_output('test_files/flexpart/20190510_00/output/grid_drydep_20190510000000.nc')
        ds_test = ds_test.fp.select_receptor_point(0)
        ds_ref = xr.open_dataset(path)
        ds_ref = ds_ref.sel(pointspec=0, numspec=0, numpointspec=0, nageclass=0)

    a

