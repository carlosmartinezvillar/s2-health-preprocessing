import rasterio
import matplotlib.pyplot as plt

# load POLYGON'S from TSV
import re
import pandas as pd
import geopandas as gpd
from shapely import wkt
import matplotlib.patches as mpatches
 

def plot_tiles():
	'''
	Plot a polygons of the United States along with the polygons of the 
	tiles in the dataset.
	'''

	# FILE PATHS
	US_SHP_PATH      = "../figures/cb_2024_us_state_500k/cb_2024_us_state_500k.shp"
	TRACTS_GEOM      = "../shapes/all_tracts/all_tracts.shp"
	S2_PRODUCTS_GEOM = "../other/search_results_geometries_2025.tsv"
	OUT_PATH         = "../figures/tiles.png"
	PLOT_CRS         = "EPSG:5070"

	# 1. LOAD FILES
	#---------------
	# 1.1 LOAD US STATES
	states      = gpd.read_file(US_SHP_PATH)
	territories = ['PR','AS','VI','MP','GU','AK','HI']
	contiguous  = states[~states['STUSPS'].isin(territories)]

	# 1.2 LOAD SENTINEL TILES USED --- load tsv file & clean.
	df = pd.read_csv(TSV_PATH, sep="\t", header=None, names=["scene_id", "geometry_raw"])
	print(f"Loaded {len(df)} rows.")

# ---------------------------------------------------------------------------
# 2. Extract MGRS tile ID from scene_id
#    Sentinel-2 naming: …_T14TNQ_… → tile = "T14TNQ"
# ---------------------------------------------------------------------------
MGRS_PATTERN = re.compile(r"(T\d{2}[A-Z]{3})")
 
df["mgrs_tile"] = df["scene_id"].str.extract(MGRS_PATTERN)
 
print(f"Unique MGRS tiles found: {df['mgrs_tile'].nunique()}")
print(df["mgrs_tile"].value_counts().head(10))
 
 
# ---------------------------------------------------------------------------
# 3. Define the MGRS tiles you want to keep
#    Edit this set to filter to whichever tiles you need.
# ---------------------------------------------------------------------------
TILES_OF_INTEREST = set(df["mgrs_tile"].dropna().unique())  # default: all tiles
# Example — uncomment and customise to restrict:
# TILES_OF_INTEREST = {"T14TNQ", "T14SMF", "T14TNM", "T15TUK"}
 
mask = df["mgrs_tile"].isin(TILES_OF_INTEREST)
df_filtered = df[mask].copy()
print(f"Rows after MGRS filter: {len(df_filtered)}")
 
 
# ---------------------------------------------------------------------------
# 4. Parse WKT geometries
#    Raw format: geography'SRID=4326;POLYGON ((...))' — strip the prefix.
# ---------------------------------------------------------------------------
def parse_geometry(raw: str):
    """Strip the geography'SRID=4326;' prefix and parse the WKT."""
    # Remove everything up to and including the first ';'
    wkt_str = re.sub(r"^geography'SRID=\d+;", "", raw).rstrip("'")
    try:
        return wkt.loads(wkt_str)
    except Exception:
        return None
 
 
df_filtered["geometry"] = df_filtered["geometry_raw"].apply(parse_geometry)
 
# Drop rows where geometry parsing failed
n_bad = df_filtered["geometry"].isna().sum()
if n_bad:
    print(f"Warning: {n_bad} geometries failed to parse and were dropped.")
df_filtered = df_filtered[df_filtered["geometry"].notna()]

 

	# tiles_gut = gpd.read_file(S2_) # read
	# tiles_gut['geometry'] = tiles_gut.geometry.apply(lambda x: x.geoms[0]) #POLYGON in GEOMETRYCOLLECTION
	# tiles_bad = gpd.read_file(S2_KML_PATH_2, driver='KML') # read
	# tiles_bad['geometry'] = tiles_bad.geometry.apply(lambda x: x.geoms[0]) #POLYGON in GEOMETRYCOLLECTION


	# 2. PROJECT TO COMMON CRS
	#-------------------------
	contiguous = contiguous.to_crs(common_crs)
	tiles_gut  = tiles_gut.to_crs(common_crs)
	tiles_bad  = tiles_bad.to_crs(common_crs)
	water      = water.to_crs(common_crs)

	# 3. PLOT LAYERS
	# --------------
	fig, ax = plt.subplots(1,1,figsize=(24,20))

	contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.2)
	tiles_gut.plot(ax=ax,color='red',alpha=0.3,edgecolor='red',linewidth=1.0)
	# tiles_bad.plot(ax=ax,color='blue',alpha=0.3,edgecolor='blue',linewidth=1.5)

	# zoom in
	xmin, ymin, xmax, ymax = contiguous.total_bounds
	print(f"X:{xmin}--{xmax} | Y: {ymin}--{ymax}")

	x_range = xmax - xmin
	y_range = ymax - ymin
	xmax = x_range*0.35 + xmin #~1/3 of US in plot
	xmin = xmin - x_range*0.02
	ymax = ymin + y_range*0.85
	ymin = ymin + y_range*0.30

	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)

	ax.set_title("Sentinel-2 (MGRS) Tiles",fontsize=24)
	ax.set_axis_off()
	plt.tight_layout()
	plt.savefig(OUT_PATH)


def plot_label(path):
	'''
	Plot a label tile.
	'''
	with rasterio.open(path) as src:
	    band_data = src.read(1, masked=False)

	tile_str = path.split('/')[-1].split('_')[0]

	# 2. Plot the array with a specific data range and colormap
	plt.figure(figsize=(8, 6))
	plt.imshow(band_data, cmap='terrain', vmin=0, vmax=1000)

	# 3. Add a colorbar and display the plot
	plt.colorbar(label='Diabetes Prevalence')
	plt.title(f'Tile {tile_str}')
	plt.show()	


if __name__ == "__main__":
	plot_label('../masks/T13SGB_diabetes.tif')
