# s2-health-preprocessing
Preprocessing of Sentinel-2 images and extraction of visual features for health-related analyses

## Source Files

* `simplify_us_region.py`: Takes the polygon shapes in .shp format for the 4 US regions (as defined by US census). Drops all polygons except the Midwest region. Simplifies the geometry by i) reducing the number of points in the shape, and ii) removing 'holes', i.e.: converting the multi-polygon shape to a single polygon shape. The resulting geometry is converted to WKT format and written to `../shapes/us_region_wkt.txt`. The number of points is arbitrarily reduced to an amount low enough to be used in a URL query.

* `merge_csv.py`: Merges the CSV tables in `../shapes/cdc_diabetes.csv` and `../shapes/RUCA-codes-2020-tract.csv`. The combining is done with an inner join on the census tract id's, `LocationID` and `TractFIPS20` respectively. Census tract rows with a `PrimaryRUCA` code of 10 are dropped. This results in a final set of 16472 census tracts (rows). The resulting table is written to `../shapes/filtered_diabetes.csv`. Columns are dropped. Final shape is 16472x9. For columns dropped see result table.

* `merge_polygons.py`: Loads all census tracts for each individual state, and merges them into a single dataframe. Census tract ids for these geometries are matched to those in `../shapes/filtered_diabetes.csv`. The combined 16472 geometries are written to `../shapes/all_tracts/all_tracts.shp`). The complete ESRI/shapefile output, including metadata, are all placed under `../shapes/all_tracts`.


* `rasterize_polygons.py`: Rasterizes/burns geometries outputted by `merge_polygons.py`. Census tract polygons are burned to raster images matching the 10-meter resolution of Sentinel-2 products. Temporal coverage of the Sentinel-2 products in our area-of-interest results in multiple Sentinel-2 images per MGRS tile. To avoid unnecessary overhead, a single raster is computed for each individual MGRS tile instead of individual S2 products. Pixel values inside each tract polygon are set to 'Data_Value', the diabetes crude prevalence value (estimated percentage over 1000 people) in `../shapes/filtered_diabetes.csv`.

* `create_chips.py`: Creates input-label chip pairs using Sentinel-2 images and rasterized census tract polygons created by `rasterize_polygons.py`. Output chip size is 224x224.