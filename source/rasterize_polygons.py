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
import subprocess as sp
import sys

############################################################
# GLOBAL
############################################################
SCALING_FACTOR = 10
REMOTE_PATH = "nrp:diabetes-s2"
# REMOTE_PATH = "/s2_volume"

############################################################
# FUNCTIONS
############################################################
def process_tile(s2_path:str,master_polygons:gpd.GeoDataFrame,data_dir:str) -> None:

	###################################
	# LOAD S2 RASTER
	###################################
	with rasterio.open(s2_path) as src:
		s2_crs       = src.crs
		s2_transform = src.transform
		s2_shape     = src.shape
		s2_bounds    = src.bounds
		s2_meta      = src.meta.copy()

	# Tile name. Base of large label raster names.
	tile_str   = s2_path.split('/')[-1].split('_')[0]


	###################################
	# INTERSECT S2 GEOMETRY TO POLYGONS
	###################################
	# Make a bounding box polygon
	tile_box = box(s2_bounds.left, s2_bounds.bottom, s2_bounds.right, s2_bounds.top)

	# Project tile box to the master polygons' CRS
	tile_box_df = gpd.GeoDataFrame(geometry=[tile_box], crs=s2_crs).to_crs(master_polygons.crs)
	tile_geom   = tile_box_df.geometry.iloc[0]

	# Keep only intersecting
	intersecting_polygons = master_polygons[master_polygons.intersects(tile_geom)]

	# Check not empty
	if intersecting_polygons.empty:
	    print(f"--> No matching polygons found for {tile_str}. Skipping.")
	    return    


	###################################
	# PROJECT POLYGONS TO TILE CRS
	###################################
	projected_polygons = intersecting_polygons.to_crs(s2_crs)


	###################################
	# ADJUST COLUMNS
	###################################
	# ADJUST DIABETES FROM PERCENT AND SINGLE DECIMAL TO 0-999 (FACTOR 10)
	projected_polygons['Data_Value'] = projected_polygons['Data_Value'] * SCALING_FACTOR
	# ADJUST LAND AREA TO 0-29674
	projected_polygons['LandArea'] = projected_polygons['LandArea'] * SCALING_FACTOR


	###################################
	# RASTERIZE LABEL (1 BAND)
	###################################
	label_shapes = [
		(g,v) for g,v in zip(projected_polygons.geometry,projected_polygons['Data_Value'])
	]

	rasterized_label = rasterize(
	    label_shapes,
	    out_shape=s2_shape,
	    transform=s2_transform,
	    fill=0,
	    all_touched=False,
	    dtype='uint16'
	)

	label_path = f"{data_dir}/masks/{tile_str}_diabetes.tif"
	s2_meta.update({
		"driver": "GTiff",
		"count": 1,
		"dtype": "uint16",
		"nodata": 0
	})

	os.makedirs(f"{data_dir}/masks", exist_ok=True)
	with rasterio.open(label_path, "w", **s2_meta) as dest:
	    dest.write(rasterized_label, 1)
	print(f"Label file written to {label_path}.")


	###################################
	# RASTERIZE FEATURES (4 BANDS)
	###################################
	cols_to_burn        = ["PrimaryRUC","Population","LandArea"]
	rasterized_features = np.zeros((4,rasterized_label.shape[0],rasterized_label.shape[1]),dtype=np.uint16)

	for i,col in enumerate(cols_to_burn):
		feature_shapes = [
			(g,v) for g,v in zip(projected_polygons.geometry,projected_polygons[col])
		]
		rasterized_features[i,:,:] = rasterize(
			feature_shapes,
			out_shape=s2_shape,
			transform=s2_transform,
			fill=0,
			all_touched=False,
			dtype='uint16'
		)

	index_shapes = [(g,v+1) for g,v in zip(projected_polygons.geometry,projected_polygons.index.tolist())]
	rasterized_features[3,:,:] = rasterize(
		index_shapes,
		out_shape=s2_shape,
		transform=s2_transform,
		fill=0,
		all_touched=False,
		dtype='uint16'
	)

	features_path = f"{data_dir}/masks/{tile_str}_features.tif"
	s2_meta.update({
		"count": 4
	})
	with rasterio.open(features_path, "w", **s2_meta) as dest:
	    dest.write(rasterized_features)
	print(f"Feature file written to {features_path}.")



def get_local_band_path(s2_id:str,data_dir:str) -> str:
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


def get_remote_band_path(s2_id:str) -> str:
	date = s2_id.split('_')[2]
	y = date[0:4]
	m = date[4:6]
	d = date[6:8]
	band_regex = f"eodata/Sentinel-2/MSI/L2A/{y}/{m}/{d}/{s2_id}/GRANULE/*/IMG_DATA/R10m/*_B02_10m.jp2"
	return band_regex


############################################################
# MAIN
############################################################
if __name__ == '__main__':

	# PARSE ARGV
	parser = argparse.ArgumentParser()
	parser.add_argument("--data-dir",required=True,default=None,
		help="Data directory.")
	parser.add_argument("--download",required=False,default=False,action=argparse.BooleanOptionalAction,
		help="Flag. If set, download Sentinel-2 products needed.")
	args = parser.parse_args()

	if args.data_dir is None:
		print("Got None for data-dir argument.")
		sys.exit(1)

	if not os.path.isdir(args.data_dir):
		print(f"Data dir {args.data_dir} not found.")
		sys.exit(1)


	# FIND UNIQUE TILES
	# assume running inside s2-health-preprocessing/source/
	with open('../other/search_results_2025.tsv','r') as fp:
	# with open('../other/search_subset.tsv') as fp: #--testing
		s2_ids = [l.split('\t')[0] for l in fp.readlines()]
	mgrs_tiles = [s.split('_')[5] for s in s2_ids]
	unique_mgrs,first_index = np.unique(mgrs_tiles,return_index=True)
	unique_ids = np.array(s2_ids)[first_index]
	# 285 tiles,5322 products,03/01--10/31
	# 284 tiles,2344 products,05/01--08/31
	# 281 tiles,1637 products,05/15--08/15*


	# IF TRUE, TRANSFER PRODUCTS CORRESPONDING TO UNIQUE TILE IDS
	if args.download:
		remote_paths = [get_remote_band_path(s) for s in unique_ids]

		# CP/UNIX
		# for rp in remote_paths:
		# 	path = REMOTE_PATH+'/'+rp
		# 	sp.run(["cp",path,args.data_dir]) #<--- need to add parent dir structure

		# S3/RCLONE
		with open("../other/temp_include.txt",'w') as fp:
			fp.write("\n".join(remote_paths))

		#run rclone download
		sp.run(["rclone","copy",REMOTE_PATH,args.data_dir,"--include-from","../other/temp_include.txt","--stats","10s","-v"])


	# LOAD ALL POLYGONS
	master_polygons = gpd.read_file("../shapes/all_tracts/all_tracts.shp")
	if master_polygons.crs != "EPSG:4326":
		master_polygons.set_crs("EPSG:4326", inplace=True)


	# CHECK FILES IN DATA_DIR MATCH SEARCH RESULTS & DOWNLOADED QUEUE
	not_found_unique_tile_products = []
	for s2_id in unique_ids:
		band_path = get_local_band_path(s2_id,args.data_dir)
		if band_path is None:
			not_found_unique_tile_products.append(s2_id)

	if len(not_found_unique_tile_products) > 0:
		# for f in not_found_unique_tile_products:
			# print(f"MISSING: {f}")
		print(f"missing {len(not_found_unique_tile_products)} products.")
		sys.exit(1)
	else:
		print(f"ALL .SAFE FOLDERS FOUND IN DATA DIR: {len(unique_ids)}\n")


	# PROCESS/BURN POLYGONS
	for i,s2_id in enumerate(unique_ids):
		print(f"Processing {s2_id} [{i+1}/{len(unique_ids)}]")
		band_path = args.data_dir + '/' + get_local_band_path(s2_id,args.data_dir)
		process_tile(band_path,master_polygons,args.data_dir)
