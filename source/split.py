import os
import argparse
import os
import numpy as np
import random
import glob
import shutil
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import geopandas as gpd
import rasterio as rio
from shapely.geometry import box

DEFAULT_SEED = 476
N_PROC = 8
EQUAL_AREA = "EPSG:5070"
PLOT_PATH  = "../figures/split.png"

RES = 20
CHIP_PX   = 256
BORDER_PX = 246
TILE_PX   = 5490

def parse_args():
	'''
	Get dir and check.
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument('--chip-dir',default=None,required=True,
		help='Dataset directory.')
	parser.add_argument('--feat-dir',default=None,required=True,
		help='Dataset directory.')	
	parser.add_argument('--out-dir',default='../other',
		help="Where to write split_chips.csv, split_tracts.csv and split_tiles.csv.")	
	parser.add_argument('--plot',action='store_true',
		help=f"Save a map of the split to {PLOT_PATH}.")
	args = parser.parse_args()

	assert args.chip_dir is not None, "split.py: Incorrect data dir argument."
	assert os.path.isdir(args.chip_dir), "split.py: Incorrect data dir argument."
	assert args.feat_dir is not None, "split.py: Incorrect feature dir argument."
	assert os.path.isdir(args.feat_dir), "split.py: Incorrect feature dir argument."
	assert args.out_dir is not None, "split.py: Incorrect out dir argument."
	assert os.path.isdir(args.out_dir), "split.py: Incorrect out dir argument."


	args.chip_dir = args.chip_dir.rstrip('/')
	args.feat_dir = args.feat_dir.rstrip('/')
	args.out_dir = args.out_dir.rstrip('/')

	return args


def build_tree(chip_names):
	'''
	Build a dict with the nested structure:

	{zone}
	|-> {tile}
		|-> {raster}
			|--> [chips]

	Passed string list follows:
		[T00ZBB_YYYYMMDDTHHMMSS_R000_00_00], which is:
		{tile}_{date}_{orbit}_{row}_{col}
	'''
	tree = {}

	# FOR EACH CHIP STR
	for chip in chip_names:

		# GET STRINGS
		raster = chip[0:27]
		tile   = chip[0:6]
		zone   = chip[1:4]

		# APPEND TO TREE -- I believe there's a better/pythonic way to do this
		if zone in tree:
			if tile in tree[zone]:
				if raster in tree[zone][tile]:
					tree[zone][tile][raster].append(chip)
				else:
					tree[zone][tile][raster] = []
					tree[zone][tile][raster].append(chip)
			else:
				tree[zone][tile] = {raster:[chip]}
		else:
			tree[zone] = {tile:{raster:[chip]}}

	return tree

# THIS LEAKS! -- Not used in final version.
def simple_split(args):

	# LIST OF CHIP IDs
	band_files  = sorted(glob.glob("*_rgb.tif",root_dir=args.chip_dir))
	chip_names  = [_[0:-8] for _ in band_files]

	#TREE -- UTM Zone -> MGRS Tile -> S2 Raster -> Chips from raster
	tree = build_tree(chip_names)

	# PRINT A BRIEF SUMMARY
	print(f"Zones: {len(tree.keys())}")
	print(f"Tiles per zone: {[len(tree[z]) for z in tree]}")
	chips_per_raster = []
	rasters_per_tile = []
	for z in tree:
		for t in tree[z]:
			for r in tree[z][t]:
				chips_per_raster.append(len(tree[z][t][r]))
			rasters_per_tile.append(len(tree[z][t].keys()))
	print(f"Avg. chips in raster: {sum(chips_per_raster)/len(chips_per_raster)}")
	print(f"Avg. rasters in tile: {sum(rasters_per_tile)/len(rasters_per_tile)}\n")



	#GET TEST TILES
	N_TEST_ZONES = 12
	N_TEST_TILES = 2
	test_zones = np.random.choice(sorted(tree),N_TEST_ZONES,replace=False)
	test_tiles = []
	print(f"\nsampling {N_TEST_TILES} from")
	for z in test_zones:

		sample = list(np.random.choice(sorted(tree[z]),N_TEST_TILES,replace=False))
		# test_tiles += list(np.random.choice(sorted(tree[z]),N_TEST_TILES,replace=False))
		test_tiles += sample
		print(f"Zone {z} of size {len(tree[z])} -- sample: {sample}")

	#GET VALIDATION
	N_VALID_ZONES = 12
	N_VALID_TILES = 2
	validate_zones = np.random.choice(sorted(tree),N_VALID_ZONES,replace=False)
	validate_tiles = []
	print(f"\nsampling {N_VALID_TILES} from")
	for z in validate_zones:
		not_in_test = set(tree[z].keys()).difference(test_tiles)
		if len(not_in_test) > 1:
			# print(f"Zone {z} of size {len(tree[z])}")
			# print(f" -- {len(not_in_test)} remaining.")
			sample = list(np.random.choice(sorted(not_in_test),N_VALID_TILES,replace=False))
			validate_tiles += sample
			print(f"Zone {z} of size {len(tree[z])} -- sample: {sample}")

	# GET TRAINING
	train_tiles = []
	for z in tree:
		for t in tree[z]:
			if t not in test_tiles and t not in validate_tiles:
				train_tiles.append(t)


	print(f"TEST TILES:\n{test_tiles}\n")
	print(f"VALIDATION:\n{validate_tiles}\n")
	print(f"TRAINING:  \n{train_tiles}\n")
	print(f"NR. TEST TILES: {len(test_tiles)}")
	print(f"NR. VALIDATION: {len(validate_tiles)}")	
	print(f"NR. TRAINING:   {len(train_tiles)}\n\n")

	# GET CHIPS PER SET
	v_bands_files = []
	for tile in validate_tiles:
		v_bands_files += glob.glob(f"{tile}_*_rgb.tif",root_dir=args.chip_dir)		
	v_label_files = [f.replace("rgb","lbl") for f in v_bands_files]
	v_masks_files = [f.replace("rgb","ftr") for f in v_bands_files]	

	t_bands_files = []
	for tile in test_tiles:
		t_bands_files += glob.glob(f"{tile}_*_rgb.tif",root_dir=args.chip_dir)		
	t_label_files = [f.replace("rgb","lbl") for f in t_bands_files]
	t_masks_files = [f.replace("rgb","ftr") for f in t_bands_files]

	train_bands_files = set(band_files).difference(v_bands_files + t_bands_files)
	train_label_files = [f.replace("rgb","lbl") for f in train_bands_files]
	train_masks_files = [f.replace("rgb","lbl") for f in train_bands_files]

	print(f"TEST CHIPS: {len(t_bands_files)} {len(t_bands_files)/len(band_files)*100:.2f}%")
	print(f"VALIDATION: {len(v_bands_files)} {len(v_bands_files)/len(band_files)*100:.2f}%")
	print(f"TRAINING:   {len(train_bands_files)} {len(train_bands_files)/len(band_files)*100:.2f}%")

	#copy
	pass


def copy_file(pair):
	src,dst = pair
	return shutil.copy2(src,dst)


def copy_threaded(file_queue,dest_dir):
	
	# List  matching file to destination dir
	file_pairs = list(zip(file_queue,[dest_dir]*len(file_queue)))

	# Copy
	with ThreadPoolExecutor(max_workers=N_PROC) as executor:
		return list(executor.map(copy_file,file_pairs))


def read_chip(feat_path):
	with rio.open(feat_path) as src:
		ftr = src.read(4)
		epsg = src.crs.to_epsg()
		left,top = src.bounds.left,src.bounds.top

	tract_ids = [int(t)-1 for t in np.unique(ftr) if t > 0] #initial index 1

	return os.path.basename(feat_path)[:-8],epsg,left,top,tract_ids


def get_tiles(chips):
	rc   = chips.chip.str.split('_')
	row  = rc.str[3].astype(int)
	col  = rc.str[4].astype(int)
	ul_x = chips.left - (BORDER_PX + col*CHIP_PX)*RES
	ul_y = chips.top  + (BORDER_PX + row*CHIP_PX)*RES

	corners = pd.DataFrame({'tile':chips.tile,'epsg':chips.epsg,'ul_x':ul_x,'ul_y':ul_y})
	corners = corners.groupby('tile').agg({'epsg':'first','ul_x':'median','ul_y':'median'})	

	b,t     = BORDER_PX*RES,TILE_PX*RES
	parts   = []
	for epsg,sub in corners.groupby('epsg'):
		area = [box(x+b, y-t+b, x+t-b, y-b) for x,y in zip(sub.ul_x,sub.ul_y)]
		g    = gpd.GeoDataFrame({'tile':sub.index,'epsg':epsg}, geometry=area, crs=f"EPSG:{epsg}")
		g['center'] = gpd.GeoSeries(g.geometry.centroid).to_crs(EQUAL_AREA).values
		parts.append(g.to_crs(EQUAL_AREA))
	return gpd.GeoDataFrame(pd.concat(parts), crs=EQUAL_AREA).set_index('tile').sort_index()	


def assign_tiles(tiles):
	rng = np.random.default_rng(DEFAULT_SEED)
	order = rng.permutation(tiles.index.to_numpy())
	n_val  = round(len(order)*0.10)
	n_test = round(len(order)*0.10)
	split = pd.Series('train',index=tiles.index)
	split[order[:n_val]] = 'val'
	split[order[n_val:n_val+n_test]] = 'test'
	return split


def assign_tracts(tracts,tiles):
	'''
	Home tile per tract: tiles whose chipped area contains the tract first, then val/test over
	train, then the largest share of the tract, then the nearest tile center.
	'''
	pairs = gpd.sjoin(tracts[['geometry']], tiles[['geometry']], predicate='intersects')
	pairs = pairs.reset_index(names='tract').rename(columns={'index_right':'tile'})
	tile  = tiles.loc[pairs.tile]
	area  = gpd.GeoSeries(tile.geometry.to_numpy(), index=pairs.index, crs=EQUAL_AREA)
	ctr   = gpd.GeoSeries(tile.center.to_numpy(), index=pairs.index, crs=EQUAL_AREA)

	# CONDITIONS PER (TRACT,TILE) PAIR
	pairs['is_in']   = pairs.within(area)
	# pairs['holdout'] = pairs.split != 'train'
	pairs['frac']    = pairs.intersection(area).area / pairs.area
	pairs['dist']    = pairs.centroid.distance(ctr)

	# SELECTION RULE: INSIDE -> FRACTION IN -> CENTROID DISTANCE
	pairs = pairs.sort_values(['tract','is_in','frac','dist'], ascending=[True,False,False,True])
	best  = pairs.drop_duplicates('tract').set_index('tract').reindex(tracts.index)
	
	# CHECK OK
	assert best.tile.notna().all(), "split_tiles.py: tract outside every chipped area."
	
	return best.tile, best.is_in


def chip_split(tract_ids,tract_split):
	'''
	Assign chips to unique train/val/test split or marked as 'dropped'
	'''
	splits = {tract_split[t] for t in tract_ids}
	return splits.pop() if len(splits) == 1 else 'dropped'


def plot_split(chips,tracts,tiles,seed,out_path):
	import matplotlib.pyplot as plt
	from matplotlib.patches import Patch, Rectangle
	import plots

	t     = CHIP_PX*RES
	parts = [gpd.GeoDataFrame(sub, geometry=[box(l, tp-t, l+t, tp) for l,tp in zip(sub.left,sub.top)], crs=f"EPSG:{e}").to_crs(EQUAL_AREA)
	         for e,sub in chips.groupby('epsg')]
	chips = gpd.GeoDataFrame(pd.concat(parts), crs=EQUAL_AREA)

	col    = {'train':'#2a78d6','val':'#eb6834','test':'#1baf7a','dropped':'#8a8983'}
	states = gpd.read_file(plots.US_SHP_PATH)
	states = states[~states.STUSPS.isin(['PR','AS','VI','MP','GU','AK','HI'])].to_crs(EQUAL_AREA)
	N,NT   = len(chips),len(tracts)
	nd     = (chips.split == 'dropped').sum()

	# zoom on the val/test tile with most dropped chips
	dz = chips[chips.split == 'dropped'].tile.value_counts()
	zt = next((t for t in dz.index if tiles.split[t] != 'train'), tiles.index[0])
	zc = tiles.center[zt]; Z = 260_000
	zoom = (zc.x-Z/2, zc.y-Z/2, zc.x+Z/2, zc.y+Z/2)

	fig,(ax,az) = plt.subplots(1, 2, figsize=(26,12), gridspec_kw={'width_ratios':[1.55,1]})
	fig.patch.set_facecolor('white')
	for a,edge,tlw in [(ax,'none',0.6),(az,'white',1.8)]:
		states.plot(ax=a, color='#f7f7f5', edgecolor='#8a8983', linewidth=0.4)
		for k in ['train','val','test']:
			tracts[tracts.split == k].plot(ax=a, color=col[k], edgecolor='none', alpha=0.2)
			chips[chips.split == k].plot(ax=a, color=col[k], edgecolor=edge, linewidth=0.5)
		chips[chips.split == 'dropped'].plot(ax=a, color=col['dropped'], edgecolor=edge, linewidth=0.5)
		for k in ['val','test']:
			tiles[tiles.split == k].boundary.plot(ax=a, color=col[k], linewidth=tlw, zorder=5)
		a.set_axis_off()
	tiles[tiles.split == 'train'].boundary.plot(ax=az, color=col['train'], linewidth=1.0, linestyle='--', zorder=5)
	tracts.boundary.plot(ax=az, color='#2b2b29', linewidth=0.3, alpha=0.4)
	for t,r in tiles.iterrows():
		if zoom[0] < r.center.x < zoom[2] and zoom[1] < r.center.y < zoom[3]:
			az.text(r.center.x, r.center.y, f"{t}\n{r.split}", ha='center', va='center', fontsize=10, fontweight='bold',
			        color='#2b2b29', bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=col[r.split], lw=1.5, alpha=0.85), zorder=6)
	xmin,ymin,xmax,ymax = tracts.total_bounds; pad = 60_000
	ax.set_xlim(xmin-pad, xmax+pad); ax.set_ylim(ymin-pad, ymax+pad)
	ax.add_patch(Rectangle(zoom[:2], Z, Z, fill=False, edgecolor='#2b2b29', linewidth=2))
	az.set_xlim(zoom[0], zoom[2]); az.set_ylim(zoom[1], zoom[3])
	az.set_title(f"Zoom on {zt} · outlines = chipped area of tiles (dashed = train) · label at tile center", fontsize=15, loc='left')

	def lab(k):
		nt,nc = (tracts.split == k).sum(),(chips.split == k).sum()
		return f"{k.capitalize()}: {(tiles.split==k).sum()} tiles, {nt:,} tracts ({nt/NT*100:.0f}%), {nc:,} chips ({nc/N*100:.0f}%)"
	h = [Patch(color=col[k], label=lab(k)) for k in ['train','val','test']]
	h.append(Patch(color=col['dropped'], label=f"Dropped chips: {nd:,} ({nd/N*100:.0f}%) — tracts of 2+ splits"))
	ax.legend(handles=h, loc='lower left', fontsize=13, framealpha=0.95,
	          title="solid = chip · pale = rest of tract · outlined = val/test tiles", title_fontsize=12)
	fig.suptitle(f"Tile-block split · random MGRS tiles · chipped area ({BORDER_PX*RES/1000:g} km border removed) (seed {seed})",
	             fontsize=20, x=0.01, ha='left')
	plt.tight_layout(rect=(0,0,1,0.96))
	plt.savefig(out_path, dpi=100)
	print(f"Plot saved to {out_path}")


def boundary_aware_split(args):

	# PATHS
	chip_ids = sorted(f[:-8] for f in glob.glob("*_rgb.tif",root_dir=args.chip_dir))
	feat_paths = [f"{args.feat_dir}/{c}_ftr.tif" for c in chip_ids]

	# READ CHIP GEOMETRIES
	print(f"Reading chip geometries...")
	with ThreadPoolExecutor(N_PROC) as ex:
		rows = list(ex.map(read_chip,feat_paths))
	chips = pd.DataFrame(rows,columns=['chip','epsg','left','top','tracts'])
	chips['tile'] = chips.chip.str.split('_').str[0]

	# LOAD TRACTS
	print(f"Loading tracts...")
	all_tracts = gpd.read_file("../shapes/all_tracts/all_tracts.shp")

	# TRAIN/VAL/SPLIT OF TILES
	print(f"Loading tiles and assigning tiles to data split...")
	tiles          = get_tiles(chips)
	tiles['split'] = assign_tiles(tiles)

	# SELECT CORRECT TRACTS IN CHIPS
	print(f"Selecting subset of tracts...")
	used   = sorted({t for ids in chips.tracts for t in ids})
	tracts = all_tracts.loc[used,['GEOID','geometry']].to_crs(EQUAL_AREA)

	# MAIN FUNCTION -- SPLIT TRACTS BY TILE AND SUBSEQUENT RULES
	print(f"Assigning tracts to data split...")
	tracts['home'],tracts['inside'] = assign_tracts(tracts,tiles)
	tracts['split'] = tracts.home.map(tiles.split) # id -> tile -> split

	# APPLY SPLIT TO CHIPS
	print(f"Assigning chips to data split...")
	tract_split    = tracts.split.to_dict()
	chips['split'] = chips.tracts.map(lambda ids: chip_split(ids,tract_split))

	# SAVE CHIP SPLIT
	out = chips[['chip','tile','split']].copy()
	out['tracts'] = chips.tracts.map(lambda ids: ' '.join(map(str,ids)))
	out.to_csv(f"{args.out_dir}/split_chips.csv",index=False)

	# SAVE TRACT SPLIT
	tracts.drop(columns='geometry').rename_axis('tract_index').to_csv(f"{args.out_dir}/split_tracts.csv")

	# SAVE TILE SPLIT
	tiles[['epsg','split']].to_csv(f"{args.out_dir}/split_tiles.csv")


	# LOG
	kept   = {t for ids,s in zip(chips.tracts,chips.split) if s != 'dropped' for t in ids}
	N,NT   = len(chips),len(tracts)
	print(f"\n{N} chips, {NT} tracts, {len(tiles)} tiles (seed {DEFAULT_SEED}). "
		f"Tracts entirely within a chipped area: {tracts.inside.mean()*100:.1f}%")
	for s in ['train','val','test']:
		nt = (tracts.split == s).sum()
		nc = (chips.split == s).sum()
		no_chip = ((tracts.split == s) & ~tracts.index.isin(kept)).sum()
		print(f"{s.upper():7s} tiles {(tiles.split == s).sum():4d}  tracts {nt:6d} ({nt/NT*100:4.1f}%)  "
			f"chips {nc:6d} ({nc/N*100:4.1f}%)  tracts w/o chip {no_chip}")
	nd = (chips.split == 'dropped').sum()
	print(f"DROPPED chips {nd:6d} ({nd/N*100:4.1f}%) -- cover tracts of 2+ splits")
	print(f"\nWritten {args.out_dir}/split_chips.csv, split_tracts.csv and split_tiles.csv")	

	# PLOT
	if args.plot:
		plot_split(chips,tracts,tiles,DEFAULT_SEED,PLOT_PATH)

if __name__ == '__main__':

	# ARGS
	args = parse_args()

	# SPLIT
	boundary_aware_split(args)




