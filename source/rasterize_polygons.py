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
############################################################
# LIBRARIES
############################################################
import glob
import os
import geopandas as gpd
from shapely.geometry import box
import rasterio
from rasterio.features import rasterize
import numpy as np

############################################################
# GLOBAL
############################################################
SCALING_FACTOR = 10

############################################################
# FUNCTIONS
############################################################

def find_unique_tiles():
	pass


def process_tile(s2_path,master_polygons):

	with rasterio.open(s2_path) as src:
		s2_crs       = src.crs
		s2_transform = src.transform
		s2_shape     = src.shape
		s2_bounds    = src.bounds
		s2_meta      = src.meta.copy()

	# Avoid 492px in border.


	# Make a bounding box polygon
	tile_box = box(s2_bounds.left, s2_bounds.bottom, s2_bounds.right, s2_bounds.top)

	# Reproject bounding box back to the master polygons' CRS
	tile_box_df = gpd.GeoDataFrame(geometry=[tile_box], crs=s2_crs).to_crs(master_polygons.crs)
	tile_geom   = tile_box_df.geometry.iloc[0]

	# Keep only the polygons intersecting this Sentinel-2 tile
	intersecting_polygons = master_polygons[master_polygons.intersects(tile_geom)]

	# Check not empty
	if intersecting_polygons.empty:
	    print(f"--> No matching polygons found for {filename}. Skipping.")
	    return    

	# Reproject to tile's CRS
	projected_polygons = intersecting_polygons.to_crs(s2_crs)
	projected_polygons['Data_Value'] = projected_polygons['Data_Value'] * SCALING_FACTOR

	# RASTERIZE LABEL
	shapes = [
		(g,v) for g,v in zip(projected_polygons.geometry,projected_polygons['Data_Value'])
	]

	rasterized_value = rasterize(
	    shapes,
	    out_shape=s2_shape,
	    transform=s2_transform,
	    fill=0,
	    all_touched=False,
	    dtype='uint16'
	)

	tile_str   = s2_path.split('/')[-1].split('_')[5]
	label_path = f"{tile_str}_diabetes.tif"
	s2_meta.update({
		"driver": "GTiff",
		"count": 1,
		"dtype": "uint16",
		"nodata": 0
	})

	os.makedirs("../masks", exist_ok=True)
	with rasterio.open(label_path, "w", **s2_meta) as dest:
	    dest.write(rasterized_array, 1)

	pass

############################################################
# MAIN
############################################################
if __name__ == '__main__':

	# FIND UNIQUE TILES
	with open('../other/search_results.tsv','r') as fp:
		s2_ids = [l.split('\t')[0] for l in fp.readlines()]
	mgrs = [s.split('_')[5] for s in s2_ids]
	unique_mgrs,first_index = np.unique(mgrs,return_index=True)
	unique_ids = np.array(s2_ids)[first_index]

	# SET PATHS FOR UNIQUE TILES

	# LOAD ALL POLYGONS
	master_polygons = gpd.read_file("../shapes/all_tracts/all_tracts.shp")
	#check crs

	# PROCESS -- BURN POLYGONS
	# for path in s2_files:
		# print(f"Processing {path}")
		# process_tile(path)

	pass