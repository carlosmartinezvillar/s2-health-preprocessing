'''
A script to produce 224x224 chips from a directory of Sentinel-2 images and 
matching labels.
'''
import os
import rasterio as rio
from rasterio.windows import Window
import numpy as np
import matplotlib.pyplot as plt
import glob
import math
import multiprocessing as mp
import time
from PIL import Image
import sys
import argparse
import subprocess as sp

import shutil
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import tifffile as tiff

__spec__ = None

#DIRS SET HERE BECAUSE THREAD ACCESS
WORK_DIR  = None #FAST VOLUME 100 S2 TIFFs: ~32GB, 100 MASK TIFFs: ~100GB
LABEL_DIR = None #SLOW VOLUME ~277GB
CHIP_DIR  = None #FAST VOLUME (inside working dir)
S2_DIR    = None #SLOW VOLUME ~338GB
# CHIP_REMOTE = "nrp:diabetes-chips"

# PIXEL LIMITS
CHIP_SIZE = 256
STRIDE    = 256

# NR OF PROCESSES PER RASTER
N_PROC = 8

# BANDS
BANDS = ["B01","B02","B03","B04","B05","B06","B07","B11","B12","B8A"]

# PRE-TRAIN LIMITS? NOT ATM
URBAN_THRESHOLD = CHIP_SIZE*CHIP_SIZE/2


####################################################################################################
# STRINGS+PARSING
####################################################################################################
def get_datastrip_id(str):
	pass


def get_granule_id(str):
	pass


def get_dynamicworld_id(s2_id: str) -> str:
	datastrip = None
	date,tile = s2_id.split('_')[2:6:3]
	gee_id    = '_'.join([date,datastrip,tile])
	return gee_id


def get_local_band_path(s2_id:str,data_dir:str) -> str:
	'''
	B02 (R20m) path of a .SAFE product, relative to data_dir. None if not found.
	Same as in rasterize_polygons.py.
	'''
	date = s2_id.split('_')[2]
	y = date[0:4]
	m = date[4:6]
	d = date[6:8]

	if y == '2023':
		prod_series = 'L2A_N0500'
	else:
		prod_series = 'L2A'

	band_regex = f"eodata/Sentinel-2/MSI/{prod_series}/{y}/{m}/{d}/{s2_id}/GRANULE/*/IMG_DATA/R20m/*_B02_20m.jp2"

	path = glob.glob(band_regex,root_dir=data_dir)
	if len(path) == 0:
		print(f"File {band_regex} not found.")
		return None
	if len(path) > 1:
		print(f"Regex {band_regex} has multiple matches.")
		return None
	return path[0]


####################################################################################################
# RASTER PROCESSING
####################################################################################################
def clean_dynamicworld_borders(src: rio.DatasetReader) -> dict:
	'''
	Take a rasterio DatasetReader for a dynamicworld image and get the indices 
	where non-zeros begin at the top, bottom, left, and right.

	Parameters
	----------
	src: rasterio.DatasetReader
		Dataset reader for a dynamic world array (which has zeroes where S2
		still has data, making it redundant to check for zeroes in the S2 array).

	Returns
	-------
	dict
		dictionary with indices of first non-zero values at top, left, right, 
		bottom

	'''
	top    = 0
	bottom = src.height-1
	left   = 0
	right  = src.width-1

	while(True):
		row = src.read(1,window=rio.windows.Window(0,top,src.width,1))
		if row.sum() == 0:
			top += 1
		else:
			break

	while(True):
		row = src.read(1,window=rio.windows.Window(0,bottom,src.width,1))
		if row.sum() == 0:
			bottom -= 1
		else:
			break

	while(True):
		col = src.read(1,window=rio.windows.Window(left,0,1,src.height))
		if col.sum() == 0:
			left += 1
		else:
			break

	while(True):
		col = src.read(1,window=rio.windows.Window(right,0,1,src.height))
		if col.sum() == 0:
			right -= 1
		else:
			break

	return {'top':top, 'bottom':bottom, 'left':left, 'right':right}


def align_dynamicworld(s2_src: rio.DatasetReader,dw_src: rio.DatasetReader) -> tuple:
	'''
	Do everything: match indices and remove borders.
	'''
	# 1. REMOVE DW NO-DATA BORDERS(~1-2px each side)
	dw_ij = remove_dynamicworld_borders(dw_src) # <---- THIS CAN BE COMBINED

	# 2. MATCH DW to S2 (DW has ~20px less on each side) 
	# DW ij's (px index) -> DW xy's (coords)
	dw_xy_ul = dw_src.xy(dw_ij['top'],dw_ij['left'],offset='center')
	dw_xy_lr = dw_src.xy(dw_ij['bottom'],dw_ij['right'],offset='center')
	# DW xy's (coords) -> S2 ij's (px index)
	s2_ij = {}
	s2_ij['top'],s2_ij['left']     = s2_src.index(dw_xy_ul[0],dw_xy_ul[1],op=math.floor)
	s2_ij['bottom'],s2_ij['right'] = s2_src.index(dw_xy_lr[0],dw_xy_lr[1],op=math.floor)

	# 3. TRIM S2 -- REMOVE S2 TILE OVERLAP & ADJUST DW
	if s2_ij['top'] < 492: #shift top down
		delta        = 492 - s2_ij['top']
		s2_ij['top'] = 492
		dw_ij['top'] = dw_ij['top'] + delta

	if s2_ij['bottom'] > 10487: #shift bottom up
		delta           = s2_ij['bottom'] - 10487
		s2_ij['bottom'] = 10487	
		dw_ij['bottom'] = dw_ij['bottom'] - delta

	if s2_ij['left'] < 492: #shift left right
		delta         = 492 - s2_ij['left']
		s2_ij['left'] = 492	
		dw_ij['left'] = dw_ij['left'] + delta

	if s2_ij['right'] > 10487: #shift right left
		delta          = s2_ij['right'] - 10487
		s2_ij['right'] = 10487		
		dw_ij['right'] = dw_ij['right'] - delta

	return s2_ij,dw_ij	


def get_strided_windows(borders):
	'''
	Given a dicts of boundaries, returns an array list with tuples (i,j) for block indices i,j and 
	window objects corresponding to the block i,j while considering only the area of the raster
	within the boundaries defined by the indices in the dict. For example, if the array had two rows
	and a column of no data (top and left) the blocks are offseted and defined as:

			    left   stride   stride*1
				| 0 0 ..  	      |
				| 0 0... |		  | 
	stride   ---+--------+--------+----
		    0 0 |        |        |
		    0 0 | (0, 0) | (0, 1) |
		     .  |        |        |
		     .  +--------+--------+
		     .  |        |        |
		        | (1, 0) | (1, 1) |
		        |        |        |
	stride*1 ---+--------+--------+---
				|                 |


	Parameters
	----------
	borders: dict
		The dictionary containing the first and last indices of usable data in
		both directions.

	Returns
	-------
	List of shape [(str,str),Window]. Contains Window objects to be read by
	rasterio.DatasetReaders and indices for the position of these objects in 
	original size raster.

	'''	
	# number of pixel rows and cols accounting for boundaries
	n_px_rows = borders['bottom'] + 1 - borders['top']
	n_px_cols = borders['right'] + 1 - borders['left']

	#nr of blocks in each direction
	block_rows = (n_px_rows - CHIP_SIZE) // STRIDE + 1
	block_cols = (n_px_cols - CHIP_SIZE) // STRIDE + 1

	#total blocks
	N = block_rows * block_cols

	windows = []

	for k in range(N):
		i = k // block_cols
		j = k % block_cols
		row_start = i * STRIDE + borders['top']
		col_start = j * STRIDE + borders['left']
		W = Window(col_start,row_start,CHIP_SIZE,CHIP_SIZE)
		windows += [[(str(i),str(j)),W]]

	return windows


def get_windows(borders):
	'''
	Given a dicts of boundaries, returns an array list with tuples (i,j) for block indices i,j and 
	window objects corresponding to the block i,j while considering only the area of the raster
	within the boundaries defined by the indices in the dict. For example, if the array had two rows
	and a column of no data (top and left) the blocks are offseted and defined as:

			    left    224      448
				| 0 0 ..  	      |
				| 0 0... |		  | 
	    top ----+--------+--------+----
		    0 0 |        |        |
		    0 0 | (0, 0) | (0, 1) |
		     .  |        |        |
		     .  +--------+--------+
		     .  |        |        |
		        | (1, 0) | (1, 1) |
		        |        |        |
		448 ----+--------+--------+---
				|                 |


	Parameters
	----------
	borders: dict
		The dictionary containing the first and last indices of usable data in
		both directions.

	Returns
	-------
	List of shape [(str,str),Window]. Contains Window objects to be read by
	rasterio.DatasetReaders and indices for the position of these objects in 
	original size raster.

	'''

	#number of px row columns accounting for outer bounds
	n_px_rows = borders['bottom'] + 1 - borders['top']
	n_px_cols = borders['right'] + 1 - borders['left']

	#nr of blocks in each direction
	block_rows = n_px_rows // CHIP_SIZE
	block_cols = n_px_cols // CHIP_SIZE

	#total blocks
	N = block_rows * block_cols

	#Set return list and append Window objects
	windows = []
	for k in range(N):
		i = k // block_cols
		j = k % block_cols
		row_start = i * CHIP_SIZE + borders['top']
		col_start = j * CHIP_SIZE + borders['left']
		W = Window(col_start,row_start,CHIP_SIZE,CHIP_SIZE)
		windows += [[(str(i),str(j)),W]]

	return windows


def copy_single_file(pair):
	src,dst = pair
	return shutil.copy2(src,dst)


def copy_threaded(file_queue,dest_dir):

	# List of tuples containing (source_path, destination_path)
	file_pairs = list(zip(file_queue,[dest_dir]*len(file_queue)))

	# Copy files simultaneously using N_PROC worker threads
	# list() consumes results so worker exceptions are raised here
	with ThreadPoolExecutor(max_workers=N_PROC) as executor:
		return list(executor.map(copy_single_file, file_pairs))


def chip_image(s2_readers,label_path,feature_path,base_id,index,N):

	# STDOUT
	print(f'[{index+1}/{N}] PROCESSING {base_id}')
	start_time = time.time()

	# LOAD BAND ARRAYS, CLIP, & NORMALIZE
	bands = []
	for reader in s2_readers:

		# LOAD
		band_array  = reader.read(1)

		# IF ONLY NO DATA, SKIP PRODUCT -- SOME EMPTY ARRAYS!?
		if int(band_array.sum()) == 0:
			print(f"EMPTY BAND ARRAY in {reader.files[0]} -- SKIPPING.")
			return

		# CLIP & NORMALIZE
		zero_mask   = band_array == 0
		high_cutoff = int(np.percentile(band_array[~zero_mask],99))
		low_cutoff  = int(np.percentile(band_array[~zero_mask],1)) #This might have to be lower?
		if high_cutoff == low_cutoff:
			print(f"CONSTANT BAND ARRAY in {reader.files[0]} -- SKIPPING.")
			return
		band_array  = np.clip(band_array,low_cutoff,high_cutoff)
		band_array  = np.round((band_array-low_cutoff)/(high_cutoff-low_cutoff)*254+1).astype(np.uint8)
		band_array  = np.where(zero_mask,0,band_array)
		bands.append(band_array)
	bands = np.array(bands)

	# SET WINDOWS
	# range 10m = 10980
	# range 20m = 5490 --> lastindex = N - 1 - BORDER = 5489-246 = 5243 = 5489 - 246
	# s2_borders = {'top': 492, 'bottom': 10487, 'left': 492, 'right': 10487}
	s2_borders = {'top': 246, 'bottom': 5243, 'left': 246, 'right': 5243}
	s2_windows = get_strided_windows(s2_borders)

	# SPLIT WINDOWS INTO WORKER SECTIONS
	process_share = len(s2_windows) // N_PROC
	leftover      = len(s2_windows) % N_PROC
	start         = [i*process_share for i in range(N_PROC)]
	stop          = [i*process_share+process_share for i in range(N_PROC)]
	stop[-1]      += leftover
	s2_window_chunks = [s2_windows[s0:s1] for s0,s1 in zip(start,stop)]

	# THROW WORKERS AT WINDOW SECTIONS
	# lock = mp.Lock() #lock to log stuff
	processes = []
	for i in range(N_PROC):
		p = mp.Process(
			target=chip_image_worker,
			args=(bands,label_path,feature_path,s2_window_chunks[i],base_id)
		)
		p.start()
		processes.append(p)

	for p in processes:
		p.join(timeout=60)

	# STDOUT	
	exec_time = time.time() - start_time
	print(f"All workers done ({exec_time:.3f} secs). ")


def chip_image_worker(band_arrays,label_path,feature_path,windows,base_id):

	# Distinct rio.DatasetReader for thread/race conditions
	lbl_rdr = rio.open(label_path,'r',tiled=True) #1 band, uint16
	ftr_rdr = rio.open(feature_path,'r',tiled=True) #3 bands, uint16

	ftr_dir = base_id.replace("chips","features")
	# Log chip info?
	# stats = []

	for k,(rowcol,w) in enumerate(windows):

		# LOAD (ONLY WINDOW SECTION) LABEL & FEATURES
		lbl_array = lbl_rdr.read(1,window=w)
		ftr_array = ftr_rdr.read([1,2,3,4],window=w) #this flushes with proc exit I suppose..

		# IF LABEL NO DATA -- SKIP CHIP 
		if (lbl_array == 0).any():
			continue

		# CHECK PROPORTION OF RUCA CODE==10 -- No
		# n_10 = (ftr_array[0]==10).sum()
		# if n_10 > URBAN_THRESHOLD:
			# continue

		# LOAD BANDS/IF NO DATA IN BAND SKIP CHIP
		b_array = band_arrays[0,w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]
		if (b_array == 0).any():
			continue

		# SLICE 
		arr = band_arrays[:,w.row_off:w.row_off+CHIP_SIZE,w.col_off:w.col_off+CHIP_SIZE]

		# GOOD -- SAVE BANDS
		row,col = rowcol
		outfile = f"{base_id}_{row.zfill(2)}_{col.zfill(2)}_rgb.tif"
		tiff.imwrite(
			outfile,
			arr,
			photometric="minisblack"
		)

		# SAVE LABEL
		outfile = f'{base_id}_{row.zfill(2)}_{col.zfill(2)}_lbl.tif'
		# img = Image.fromarray(lbl_array)
		# img.save(outfile)
		tiff.imwrite(
			outfile,
			lbl_array,
			photometric="minisblack"
		)

		# SAVE FEATURES
		ftr_meta = ftr_rdr.meta.copy()
		ftr_meta.update({
			"height": CHIP_SIZE,
			"width": CHIP_SIZE,
			"transform": rio.windows.transform(w,ftr_rdr.transform) #chip's own origin, not the tile's
		})
		outfile = f'{ftr_dir}_{row.zfill(2)}_{col.zfill(2)}_ftr.tif'
		with rio.open(outfile,"w",**ftr_meta) as dst:
			dst.write(ftr_array)

		# STATS/LOG?
		# diabetes = lbl_array.mean() #or weighted mean.. something like that.
		# stats.append(f'{outfile.split('/')[-1][:-8]}\t{diabetes}')

	# LOG?
	# lock.acquire()
	# with open(f'{CHIP_DIR}/stats.txt','a') as fp:
	# 	fp.write('\n'.join(stats))
	# lock.release()


if __name__ == '__main__':

	########## ARGV CONFIG ##########
	parser = argparse.ArgumentParser(
		prog="chips.py",
		description="Large Sentinel-2 and labels to 256x256 images.")

	# PATHS
	parser.add_argument('--work-dir',default=None,
		help="Temporary directory to load/offload data.")
	parser.add_argument('--chip-dir',default=None,
		help="Output directory for resulting chips")
	parser.add_argument('--s2-dir',default=None,
		help="Source directory for raw Sentinel-2 products.")
	parser.add_argument('--label-dir',default=None,
		help="Source directory for mask rasters.")
	# parser.add_argument('--selected-products',default='../other/selected_products.txt',
		# help="List of S2 product ids (one per tile) written by rasterize_polygons.py.")


	########## SET ARGS ##########
	args = parser.parse_args()
	WORK_DIR  = args.work_dir
	CHIP_DIR  = args.chip_dir
	S2_DIR    = args.s2_dir 
	LABEL_DIR = args.label_dir
	# SELECTED_PRODUCTS = args.selected_products

	if not os.path.isdir(WORK_DIR):
		print(f"WORK_DIR {WORK_DIR} not found. EXIT(1).")
		sys.exit(1)
	if WORK_DIR[-1] == '/':
		WORK_DIR = WORK_DIR.rstrip('/')

	if CHIP_DIR is None:
		os.makedirs(WORK_DIR + '/chips',exist_ok=True)
		os.makedirs(WORK_DIR + '/features',exist_ok=True)
		CHIP_DIR = WORK_DIR + '/chips'
	if not os.path.isdir(CHIP_DIR):
		print(f"CHIP_DIR in {CHIP_DIR} not found. EXIT(1).")
		sys.exit(1)

	if not os.path.isdir(S2_DIR):
		print("S2_DIR not found. EXITING.")
		sys.exit(1)
	if S2_DIR[-1] == '/':
		S2_DIR = S2_DIR.rstrip('/')

	if not os.path.isdir(LABEL_DIR):
		print("LABEL_DIR not found. EXITING.")
		sys.exit(1)
	if LABEL_DIR[-1] == '/':
		LABEL_DIR = LABEL_DIR.rstrip('/')

	print(f"WORK_DIR set to:  {WORK_DIR}")
	print(f"CHIP_DIR set to:  {CHIP_DIR}")
	print(f"S2_DIR set to:    {S2_DIR}")
	print(f"LABEL_DIR set to: {LABEL_DIR}")

	# if not os.path.isfile(SELECTED_PRODUCTS):
	# 	print(f"SELECTED_PRODUCTS {SELECTED_PRODUCTS} not found. Run rasterize_polygons.py first. EXIT(1).")
	# 	sys.exit(1)


	########## GET UNIQUE TILES FROM LABEL DIR ###############
	# label_tiffs  = glob.glob('*.tif',root_dir=LABEL_DIR) #arg/masks
	# label_tiles = [s.split('_')[0] for s in label_tiffs]

	########## GET PRODUCT INTERSECTION ##########
	# one product per tile (largest footprint), as selected by rasterize_polygons.py
	# with open(SELECTED_PRODUCTS,'r') as fp:
		# selected_ids = [l.strip() for l in fp.readlines() if l.strip()]
	# selected_ids = [s for s in selected_ids if s.split('_')[5] in label_tiles]

	# s2_good_products = selected_ids
	# missing_products = []
	# for s2_id in selected_ids:
	# 	b2_path = get_local_band_path(s2_id,S2_DIR)
	# 	if b2_path is None:
	# 		missing_products.append(s2_id)
	# 	else:
	# 		s2_good_products.append(b2_path)

	# if len(missing_products) > 0:
	# 	print(f"MISSING {len(missing_products)} SELECTED PRODUCTS IN S2_DIR. EXIT(1).")
	# 	sys.exit(1)
	# print(f"PRODUCTS MATCHING LABELS: {len(s2_good_products)}.")

	########## GO THRU ALL PRODUCTS ##########
	with open('../other/search_results_2023.tsv','r') as fp:
		lines = fp.readlines()
	safe_folder_ids = [l.split('\t')[0] for l in lines]
	band2_paths     = [get_local_band_path(s,CHIP_DIR) for s in safe_folder_ids] 

	sys.exit(1)

	########## SPLIT AND QUEUE ################
	chunk_size  = 50
	N_chunks    = len(band2_paths) // chunk_size
	remainder   = len(band2_paths) % chunk_size
	chunk_queue = []
	for i in range(N_chunks):
		chunk_queue.append(band2_paths[i*chunk_size:i*chunk_size+chunk_size])
	if remainder != 0:
		chunk_queue.append(band2_paths[N_chunks*chunk_size:])

	########## PROCESS  #######################
	for i,chunk in enumerate(chunk_queue):

		print(f"Chunk {i+1}/{len(chunk_queue)}")

		chip_base_paths = []
		tiles_in_chunk  = []
		copy_band_queue = []

		########## DOWNLOAD/COPY ####################
		for b2_path in chunk:

			# GET SOME STRINGS
			for b in BANDS:
				band_path = b2_path.replace("_B02_",f"_{b}_")
				copy_band_queue.append(f"{S2_DIR}/{band_path}")
			
			tile  = b2_path.split('/')[-1].split('_')[0]
			date  = b2_path.split('/')[-1].split('_')[1]
			orbit = b2_path.split('/')[7].split('_')[4]

			chip_base_paths.append(f"{CHIP_DIR}/{tile}_{date}_{orbit}")
			tiles_in_chunk.append(tile)

		print(f"Copying {len(chunk)} .jp2 band files in chunk...")
		copy_threaded(copy_band_queue,WORK_DIR)

		# COPY ONLY NECESSARY LABELS
		copy_mask_queue = []
		for t in list(np.unique(tiles_in_chunk)):
			copy_mask_queue.append(f"{LABEL_DIR}/{t}_diabetes.tif")
			copy_mask_queue.append(f"{LABEL_DIR}/{t}_features.tif")

		print(f"Copying {len(copy_mask_queue)} label+features in chunk.")
		copy_threaded(copy_mask_queue,WORK_DIR)

		########## CHIP ####################
		for i,product in enumerate(chunk):

			# PATHS & READERS
			local_b2_path = product.split('/')[-1]
			band_readers  = []

			for b in BANDS:
				local_band_path = local_b2_path.replace("_B02_",f"_{b}_")
				band_readers.append(rio.open(f"{WORK_DIR}/{local_band_path}",'r',tiled=True))

			label_path   = f"{WORK_DIR}/{tiles_in_chunk[i]}_diabetes.tif"
			feature_path = f"{WORK_DIR}/{tiles_in_chunk[i]}_features.tif"

			# CHIP
			try:
				chip_image(band_readers,label_path,feature_path,chip_base_paths[i],i,len(chunk))
			finally:
				for reader in band_readers:
					reader.close()


		########## DELETE FILES ############
		print("Deleting.tif files...")
		for file_path in glob.glob(os.path.join(WORK_DIR, '*.tif')):
		    os.remove(file_path)

		print("Deleting .jp2 files...")
		for file_path in glob.glob(os.path.join(WORK_DIR, '*.jp2')):
			os.remove(file_path)


	print("DONE.")