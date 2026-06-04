import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

# LOAD THE CENSUS BUREAU REGION SHAPEFILE
file_path = "../shapes/cb_2025_us_region_20m/cb_2025_us_region_20m.shp"
gdf = gpd.read_file(file_path)

# FILTER TABLE 
target_regions = ["Midwest"]
midwest_gdf    = gdf[gdf["NAME"].isin(target_regions)].copy()

# SET TO 4326 CRS
midwest_gdf    = midwest_gdf.to_crs("EPSG:4326")

# remove holes -- multipolygon to polygon
largest_polygon = Polygon(max(midwest_gdf["geometry"].iloc[0].geoms, key=lambda p: p.area).exterior.coords)
midwest_gdf["geometry"] = largest_polygon

# SIMPLIFY TO POLYGON
tolerance = 0.50 
midwest_gdf["geometry"] = midwest_gdf["geometry"].simplify(tolerance, preserve_topology=True)

# PRINT AS STRING
wkt_str = midwest_gdf["geometry"].to_wkt().iat[0]
print(wkt_str)
with open('../shapes/us_region_wkt.txt','w') as fp:
	fp.write(wkt_str)

# PLOT
midwest_gdf.plot(edgecolor='black', color='lightblue', figsize=(10, 6))
plt.title("Simplified Region Polygon")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.savefig("../figures/simplified_region.png")