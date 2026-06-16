import rasterio
import matplotlib.pyplot as plt

# load POLYGON'S from TSV
import re
import pandas as pd
import geopandas as gpd
from shapely import wkt
import matplotlib.patches as mpatches
import numpy as np
 
# FILE PATHS
US_SHP_PATH      = "../figures/cb_2024_us_state_500k/cb_2024_us_state_500k.shp"
TRACTS_GEOM      = "../shapes/all_tracts/all_tracts.shp"
S2_PRODUCTS_GEOM = "../other/search_results_geometries_2025.tsv"

PLOT_CRS         = "EPSG:5070"
LABEL_MASK_LIST  = "../other/label_tiles.txt" #246 rasters w/ distinct tiles

def plot_tiles_and_tracts():
	'''
	Plot polygons of the United States, census tract polygons, and tile polygons 
	in dataset.
	'''
	# ---------------------------------------------------------------------------
	# 1. LOAD FILES
	# ---------------------------------------------------------------------------
	# 1.1 LOAD US STATES
	states      = gpd.read_file(US_SHP_PATH)
	territories = ['PR','AS','VI','MP','GU','AK','HI']
	contiguous  = states[~states['STUSPS'].isin(territories)]

	# 1.2 LOAD SENTINEL TILES USED --- load tsv file & clean.
	with open(S2_PRODUCTS_GEOM,'r') as fp:
		lines = fp.readlines()
	safe_ids  = [l.split('\t')[0] for l in lines]
	tile_geom = [l.split('\t')[1].split(';')[1].rstrip("'") for l in lines]

	# 1.3 LOAD CENSUS TRACTS
	all_tracts = gpd.read_file(TRACTS_GEOM)

	# ---------------------------------------------------------------------------
	# 2. Extract MGRS tile ID/Set to actual tiles used (in labels)
	# ---------------------------------------------------------------------------
	mgrs = [s.split("_")[5] for s in safe_ids]
	unique_mgrs,unique_mgrs_idx = np.unique(mgrs,return_index=True)
	unique_tile_geom = np.array(tile_geom)[unique_mgrs_idx]
	
	with open(LABEL_MASK_LIST,'r') as fp:
		label_tiles = [line.split('_')[0] for line in fp.readlines()]
 
	good_mgrs_mask = np.isin(unique_mgrs,label_tiles)
	good_mgrs      = unique_mgrs[good_mgrs_mask]
	good_tile_geom = unique_tile_geom[good_mgrs_mask]

	bad_mgrs      = unique_mgrs[~good_mgrs_mask]
	bad_tile_geom = unique_tile_geom[~good_mgrs_mask]

	# ---------------------------------------------------------------------------
	# 3. Parse WKT geometries
	#    Raw format: geography'SRID=4326;POLYGON ((...))' — strip the prefix.
	# ---------------------------------------------------------------------------
	good_tile_wkts = [wkt.loads(s) for s in good_tile_geom]
	tile_df = pd.DataFrame({"tile": good_mgrs})
	tile_gdf = gpd.GeoDataFrame(
	    tile_df,
	    geometry=good_tile_wkts,
	    crs="EPSG:4326"
	)

	bad_tile_wkts = [wkt.loads(s) for s in bad_tile_geom]
	bad_tile_df   = pd.DataFrame({"tile": bad_mgrs})
	bad_tile_gdf  = gpd.GeoDataFrame(
	    bad_tile_df,
	    geometry=bad_tile_wkts,
	    crs="EPSG:4326"
	)
	bad_tile_gdf['geometry'] = bad_tile_gdf['geometry'].make_valid()

	# ---------------------------------------------------------------------------
	# 4. PROJECT TO COMMON CRS
	# ---------------------------------------------------------------------------
	contiguous    = contiguous.to_crs(PLOT_CRS)
	tile_gdf      = tile_gdf.to_crs(PLOT_CRS)
	bad_tile_gdf  = bad_tile_gdf.to_crs(PLOT_CRS)
	all_tracts    = all_tracts.to_crs(PLOT_CRS)

	# ---------------------------------------------------------------------------
	# 5. PLOT LAYERS
	# ---------------------------------------------------------------------------
	fig, ax = plt.subplots(1,1,figsize=(24,20))

	contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.15)
	all_tracts.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.1)
	tile_gdf.plot(ax=ax,color='blue',alpha=0.1,edgecolor='blue',linewidth=1.0)
	bad_tile_gdf.plot(ax=ax,color='red',alpha=0.1,edgecolor='red',linewidth=1.0)

	# zoom in
	xmin, ymin, xmax, ymax = contiguous.total_bounds
	# print(f"X:{xmin}--{xmax} | Y: {ymin}--{ymax}")
	x_range = xmax - xmin
	y_range = ymax - ymin
	xmin = x_range*0.3 + xmin #last ~2/3 of contiguous US
	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)

	# title, layout, etc
	ax.set_title("Census Tracts & Intersecting Sentinel-2 Tiles",fontsize=24)
	ax.set_axis_off()
	plt.tight_layout()
	plt.savefig("../figures/tiles_tracts.png")


def plot_tiles():
	'''
	Plot polygons of the United States and S2 tile polygons dataset.
	'''
	# ---------------------------------------------------------------------------
	# 1. LOAD FILES
	# ---------------------------------------------------------------------------
	# 1.1 LOAD US STATES
	states      = gpd.read_file(US_SHP_PATH)
	territories = ['PR','AS','VI','MP','GU','AK','HI']
	contiguous  = states[~states['STUSPS'].isin(territories)]

	# 1.2 LOAD SENTINEL TILES USED --- load tsv file & clean.
	# df = pd.read_csv(S2_PRODUCTS_GEOM, sep="\t", header=None, names=["scene_id", "geometry_raw"])
	# print(f"Loaded {len(df)} rows.")
	with open(S2_PRODUCTS_GEOM,'r') as fp:
		lines = fp.readlines()
	safe_ids  = [l.split('\t')[0] for l in lines]
	tile_geom = [l.split('\t')[1].split(';')[1].rstrip("'") for l in lines]

	# ---------------------------------------------------------------------------
	# 2. Extract MGRS tile ID/Set to actual tiles in labels
	# ---------------------------------------------------------------------------
	mgrs = [s.split("_")[5] for s in safe_ids]
	unique_mgrs,unique_mgrs_idx = np.unique(mgrs,return_index=True)
	unique_tile_geom = np.array(tile_geom)[unique_mgrs_idx]
	
	with open(LABEL_MASK_LIST,'r') as fp:
		label_tiles = [line.split('_')[0] for line in fp.readlines()]
 
	good_mgrs_mask = np.isin(unique_mgrs,label_tiles)
	good_mgrs      = unique_mgrs[good_mgrs_mask]
	good_tile_geom = unique_tile_geom[good_mgrs_mask]

	bad_mgrs      = unique_mgrs[~good_mgrs_mask]
	bad_tile_geom = unique_tile_geom[~good_mgrs_mask]

	# ---------------------------------------------------------------------------
	# 3. Parse WKT geometries
	#    Raw format: geography'SRID=4326;POLYGON ((...))' — strip the prefix.
	# ---------------------------------------------------------------------------
	good_tile_wkts = [wkt.loads(s) for s in good_tile_geom]
	tile_df = pd.DataFrame({"tile": good_mgrs})
	tile_gdf = gpd.GeoDataFrame(
	    tile_df,
	    geometry=good_tile_wkts,
	    crs="EPSG:4326"
	)

	bad_tile_wkts = [wkt.loads(s) for s in bad_tile_geom]
	bad_tile_df   = pd.DataFrame({"tile": bad_mgrs})
	bad_tile_gdf  = gpd.GeoDataFrame(
	    bad_tile_df,
	    geometry=bad_tile_wkts,
	    crs="EPSG:4326"
	)
	bad_tile_gdf['geometry'] = bad_tile_gdf['geometry'].make_valid()

	# ---------------------------------------------------------------------------
	# 4. PROJECT TO COMMON CRS
	# ---------------------------------------------------------------------------
	contiguous    = contiguous.to_crs(PLOT_CRS)
	tile_gdf      = tile_gdf.to_crs(PLOT_CRS)
	bad_tile_gdf  = bad_tile_gdf.to_crs(PLOT_CRS)

	# ---------------------------------------------------------------------------
	# 5. PLOT LAYERS
	# ---------------------------------------------------------------------------
	#figure
	fig, ax = plt.subplots(1,1,figsize=(24,20))

	#plot
	contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.2)
	tile_gdf.plot(ax=ax,color='blue',alpha=0.1,edgecolor='blue',linewidth=1.0)
	bad_tile_gdf.plot(ax=ax,color='red',alpha=0.1,edgecolor='red',linewidth=1.0)

	# zoom in
	xmin, ymin, xmax, ymax = contiguous.total_bounds
	x_range = xmax - xmin
	y_range = ymax - ymin
	xmin = x_range*0.3 + xmin #last ~2/3 of contiguous US
	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)

	# title, layout, etc
	ax.set_title("Dataset Sentinel-2 Tiles",fontsize=24)
	ax.set_axis_off()
	plt.tight_layout()
	plt.savefig("../figures/tiles.png")


def plot_tracts():
	'''
	Plot US polygons and census tract polygons.
	'''
	# ---------------------------------------------------------------------------
	# 1. LOAD FILES
	# ---------------------------------------------------------------------------
	# 1.1 LOAD US STATES
	states      = gpd.read_file(US_SHP_PATH)
	territories = ['PR','AS','VI','MP','GU','AK','HI']
	contiguous  = states[~states['STUSPS'].isin(territories)]	
	all_tracts  = gpd.read_file(TRACTS_GEOM)

	# ---------------------------------------------------------------------------
	# 2. PROJECT TO COMMON CRS
	# ---------------------------------------------------------------------------
	contiguous    = contiguous.to_crs(PLOT_CRS)
	all_tracts    = all_tracts.to_crs(PLOT_CRS)

	# ---------------------------------------------------------------------------
	# 3. PLOT LAYERS
	# ---------------------------------------------------------------------------
	# ax figure
	fig, ax = plt.subplots(1,1,figsize=(24,20))

	# plot
	contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.2)
	all_tracts.plot(ax=ax,color='red',alpha=0.5,edgecolor='black',linewidth=0.1)

	# zoom in
	xmin, ymin, xmax, ymax = contiguous.total_bounds
	x_range = xmax - xmin
	y_range = ymax - ymin
	xmin = x_range*0.3 + xmin #last ~2/3 of contiguous US
	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)

	# title, layout, etc
	ax.set_title("Dataset Census Tracts",fontsize=24)
	ax.set_axis_off()
	plt.tight_layout()
	plt.savefig("../figures/tracts.png")


def plot_label(path):
	'''
	Plot a label tile.
	'''
	with rasterio.open(path) as src:
	    band_data = src.read(1, masked=False)

	band_data = band_data/10 #percentage stored as 0-999 to save memory with uint16
	unique_diabetes_values = np.unique(band_data)

	tile_str = path.split('/')[-1].split('_')[0]

	# 2. Plot the array with a specific data range and colormap
	plt.figure(figsize=(8, 6))
	plt.imshow(band_data, cmap='terrain', vmin=0, vmax=100)

	# 3. Add a colorbar and display the plot
	plt.colorbar(label='Prevalence (%/1000)')
	plt.title(f'Tile {tile_str} –– Diabetes Prevalence')
	plt.show()


def plot_features(path):
	'''
	Plot feature .tif for entire MGRS tile
	'''
	with rasterio.open(path) as src:
		band_data = src.read([1,2,3],masked=False)

	ruca_code   = band_data[0,:,:]
	population  = band_data[1,:,:]
	land_area   = band_data[2,:,:] / 10 #stored as uint16 by multiplying original single-decimal floats by 10
	# tract_index = band_data[3,:,:]

	unique_ruca       = np.unique(ruca_code)
	unique_population = np.unique(population)
	unique_land_area  = np.unique(land_area)
	# unique_indices    = np.unique(tract_index)

	tile_str = path.split('/')[-1].split('_')[0]

	plt.figure(figsize=(8, 6))
	plt.imshow(ruca_code, cmap='terrain', vmin=0, vmax=9) #0-9
	plt.colorbar(label='Primary RUCA Code')
	plt.title(f'Tile {tile_str}')
	plt.show()

	plt.figure(figsize=(8, 6))
	plt.imshow(population, cmap='terrain', vmin=0, vmax=8000) #0-9
	plt.colorbar(label='Population')
	plt.title(f'Tile {tile_str}')
	plt.show()

	plt.figure(figsize=(8, 6))
	plt.imshow(land_area, cmap='terrain', vmin=0, vmax=30000) #0-22k? or something like that, check filtered_diabetes.csv
	plt.colorbar(label='Area')
	plt.title(f'Tile {tile_str} -- Census Tract Area')
	plt.show()

	# plt.figure(figsize=(8, 6))
	# plt.imshow(tract_index, cmap='terrain', vmin=0, vmax=16473) #0-16473
	# plt.colorbar(label='Filtered CSV Tract Index')
	# plt.title(f'{tile_str}, Values: {str(list(unique_indices))}')
	# plt.show()


def plot_features_chip(path):
	'''
	Plot a feature .tif for a chip
	'''
	with rasterio.open(path) as src:
		band_data = src.read([1,2,3,4],masked=False)

	ruca_code   = band_data[0,:,:]
	population  = band_data[1,:,:]
	land_area   = band_data[2,:,:] / 10 #stored as uint16 by multiplying original single-decimal floats by 10
	tract_index = band_data[3,:,:]

	tile_str = path.split('/')[-1].split('_')[0]
	tile_row,tile_col = path.split('/')[-1].split('_')[3:5]

	unique_ruca       = np.unique(ruca_code)
	unique_population = np.unique(population)
	unique_land_area  = np.unique(land_area)
	unique_indices    = np.unique(tract_index)

	plt.figure(figsize=(8, 6))
	plt.imshow(ruca_code, cmap='terrain', vmin=0, vmax=9) #0-9
	plt.colorbar(label='Primary RUCA Code')
	plt.title(f'{tile_str}_{tile_row}_{tile_col}, Values: {str(list(unique_ruca))}')
	plt.show()

	plt.figure(figsize=(8, 6))
	plt.imshow(population, cmap='terrain', vmin=0, vmax=8000) #0-9
	plt.colorbar(label='Population')
	plt.title(f'{tile_str}_{tile_row}_{tile_col}, Values: {str(list(unique_population))}')
	plt.show()

	plt.figure(figsize=(8, 6))
	plt.imshow(land_area, cmap='terrain', vmin=0, vmax=30000) #0-22k? or something like that
	plt.colorbar(label='Area')
	plt.title(f'{tile_str}_{tile_row}_{tile_col}, Values: {str(list(unique_land_area))}')
	plt.show()

	# plt.figure(figsize=(8, 6))
	# plt.imshow(tract_index, cmap='terrain', vmin=0, vmax=16473) #0-16473
	# plt.colorbar(label='Filtered CSV Tract Index')
	# plt.title(f'{tile_str}, Values: {str(list(unique_indices))}')
	# plt.show()



if __name__ == "__main__":
	# plot_label('../masks/T13SGB_diabetes.tif')
	# plot_tiles_and_tracts()
	# plot_tiles()
	# plot_tracts()

	# plot_label('../../health_chips/chips/T15TUH_20250619T165849_R069_39_42_lbl.tif')
	# plot_features_chip('../../health_chips/features/T15TUH_20250619T165849_R069_39_42_ftr.tif')
	plot_features('../../health_chips/T15TUH_features.tif')
