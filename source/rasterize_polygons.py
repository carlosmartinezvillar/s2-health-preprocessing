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
import argparse

############################################################
# GLOBAL
############################################################
SCALING_FACTOR = 10

############################################################
# FUNCTIONS
############################################################
def process_tile(s2_path:str,master_polygons:gpd.GeoDataFrame):

	# LOAD S2 RASTER
	with rasterio.open(s2_path) as src:
		s2_crs       = src.crs
		s2_transform = src.transform
		s2_shape     = src.shape
		s2_bounds    = src.bounds
		s2_meta      = src.meta.copy()


	# INTERSECT S2 GEOMETRY TO POLYGONS
	# Make a bounding box polygon
	tile_box = box(s2_bounds.left, s2_bounds.bottom, s2_bounds.right, s2_bounds.top)

	# Project bounding box back to the master polygons' CRS
	tile_box_df = gpd.GeoDataFrame(geometry=[tile_box], crs=s2_crs).to_crs(master_polygons.crs)
	tile_geom   = tile_box_df.geometry.iloc[0]

	# Keep only the polygons intersecting this Sentinel-2 tile
	intersecting_polygons = master_polygons[master_polygons.intersects(tile_geom)]

	# Check not empty
	if intersecting_polygons.empty:
	    print(f"--> No matching polygons found for {s2_path}. Skipping.")
	    return    

	# Reproject to tile's CRS
	projected_polygons = intersecting_polygons.to_crs(s2_crs)
	projected_polygons['Data_Value'] = projected_polygons['Data_Value'] * SCALING_FACTOR


	# RASTERIZE LABEL
	shapes = [
		(g,v) for g,v in zip(projected_polygons.geometry,projected_polygons['Data_Value'])
	]

	rasterized_shape = rasterize(
	    shapes,
	    out_shape=s2_shape,
	    transform=s2_transform,
	    fill=0,
	    all_touched=False,
	    dtype='uint16'
	)

	tile_str   = s2_path.split('/')[-1].split('_')[0]
	label_path = f"../masks/{tile_str}_diabetes.tif"
	s2_meta.update({
		"driver": "GTiff",
		"count": 1,
		"dtype": "uint16",
		"nodata": 0
	})

	os.makedirs("../masks", exist_ok=True)
	with rasterio.open(label_path, "w", **s2_meta) as dest:
	    dest.write(rasterized_shape, 1)
	print(f"Label file written to {label_path}.")



def get_band_path(s2_id:str,data_dir:str):
	date = s2_id.split('_')[2]
	y = date[0:4]
	m = date[4:6]
	d = date[6:8]
	band_regex = f"eodata/Sentinel-2/MSI/L2A/{y}/{m}/{d}/{s2_id}/GRANULE/*/IMG_DATA/R10m/*_B02_10m.jp2"
	path = glob.glob(band_regex,root_dir=data_dir)
	if len(path) == 0:
		print(f"File {band_regex} not found.")
		return None
	if len(path) > 1:
		print(f"Regex {band_regex} has multiple matches.")
		return None
	return path[0]

############################################################
# MAIN
############################################################
if __name__ == '__main__':

	# PARSE ARGV
	parser = argparse.ArgumentParser()
	parser.add_argument("--data-dir",required=True,default=None,help="Data directory.")
	args = parser.parse_args()

	if args.data_dir is None:
		print("Got None for data-dir argument.")
		sys.exit(1)

	if not os.path.isdir(args.data_dir):
		print(f"Data dir {args.data_dir} not found.")
		sys.exit(1)


	# FIND UNIQUE TILES
	# assume running inside s2-health-preprocessing/source/
	# with open('../other/search_results_2025.tsv','r') as fp:
	with open('../other/search_subset.tsv','r') as fp:
		s2_ids = [l.split('\t')[0] for l in fp.readlines()]
	mgrs_tiles = [s.split('_')[5] for s in s2_ids]
	unique_mgrs,first_index = np.unique(mgrs_tiles,return_index=True)
	unique_ids = np.array(s2_ids)[first_index]

	# 285 tiles,5322 products,03/01--10/31
	# 284 tiles,2344 products,05/01--08/31
	# 281 tiles,1637 products,05/15--08/15*


	# LOAD ALL POLYGONS
	master_polygons = gpd.read_file("../shapes/all_tracts/all_tracts.shp")
	if master_polygons.crs != "EPSG:4326":
		master_polygons.set_crs("EPSG:4326", inplace=True)


	# CHECK FILES IN DATA_DIR MATCH SEARCH RESULTS & DOWNLOADED QUEUE
	not_found_unique_tile_products = []
	for s2_id in unique_ids:
		band_path = get_band_path(s2_id,args.data_dir)
		if band_path is None:
			not_found_unique_tile_products.append(band_path)

	if len(not_found_unique_tile_products) > 0:
		for f in not_found_unique_tile_products:
			print(f"MISSING: {f}")
		sys.exit(1)
	else:
		print(f"All {len(unique_ids)} FOUND.")


	# PROCESS/BURN POLYGONS
	for s2_id in unique_ids:
		print(f"Processing {s2_id}")
		band_path = args.data_dir + '/' + get_band_path(s2_id,args.data_dir)
		print(band_path)
		process_tile(band_path,master_polygons)
