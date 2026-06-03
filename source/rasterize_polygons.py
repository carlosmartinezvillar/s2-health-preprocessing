'''
Script to burn census tracts to 10m-resolution rasters matching Sentinel-2 products.
Finds unique tiles in Sentinel-2. Iterates through each.
Overlapping polygons are converted to raster labels where channels correspond to
particular features. This is done via two files.
The first file is simply the labels, and corresponds to the disease prevalance.
This value is applied to each pixel in a chip covering a tract with that value.
The second file, "the mask", maps individual pixels to the remaining features of
RUCA code, population, area.
'''

import geopandas as gpd
import rasterio
from rasterio.features import rasterize

