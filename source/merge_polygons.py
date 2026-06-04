'''
Script to merge and filter all census tract polygons into a single
file.
'''
import geopandas as gpd
import pandas as pd
import glob
import os

# GET PATHS
shp_paths = glob.glob("../shapes/tl_2025_*_tract/tl_2025_*_tract.shp")

# LOAD SHAPE FILES
gdf_list = []
for path in shp_paths:
	gdf = gpd.read_file(path)
	gdf = gdf[['GEOID','ALAND','AWATER','geometry']]
	gdf_list.append(gdf)

# LOAD FEATURES
feature_cols = ['TractFIPS20','Data_Value','PrimaryRUCA','Population','LandArea']
features = pd.read_csv("../shapes/filtered_diabetes.csv",usecols=feature_cols)
features['TractFIPS20'] = features['TractFIPS20'].astype(str)

# CONCATENATE ALL GEOMETRIES
merged_polygons = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)
merged_polygons['GEOID'] = merged_polygons['GEOID'].astype(str)

# INNER JOIN TO DATA VALUES USING TRACT ID
merged = pd.merge(
	features,
	merged_polygons,
    left_on='TractFIPS20',
    right_on='GEOID',
    how='inner'
)
merged = gpd.GeoDataFrame(merged,geometry="geometry")
merged.crs = merged_polygons.crs
merged = merged.to_crs("EPSG:4326")

# SAVE
os.makedirs("../shapes/all_tracts", exist_ok=True)
out_path = "../shapes/all_tracts/all_tracts.shp"
merged.to_file(out_path)
print(f"Saved merged file to: {out_path}")