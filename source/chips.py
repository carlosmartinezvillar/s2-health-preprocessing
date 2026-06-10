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

# Typing
# ndarray = np.ndarray

#DIRS SET HERE BECAUSE THREAD ACCESS
WORK_DIR  = None #FAST VOLUME 100 S2 TIFFs: ~32GB, 100 MASK TIFFs: ~100GB
LABEL_DIR = None #SLOW VOLUME ~277GB
CHIP_DIR  = None #FAST VOLUME (inside working dir)
S2_DIR    = None #SLOW VOLUME ~338GB
# CHIP_REMOTE = "nrp:diabetes-chips"

# PIXEL LIMITS
CHIP_SIZE = 224
STRIDE    = 224

# NR OF PROCESSES PER RASTER
N_PROC = 8

####################################################################################################
# CLASSES
####################################################################################################
class EmptyLabelError(Exception):
	pass

class IncompleteDirError(Exception):
	pass

# class Product():
# 	'''
# 	An object referencing a single Sentinel-2 product in the ESA database.

# 	Parameters
# 	----------
# 	id:
# 	tile:
# 	date:
# 	orbit:
# 	s2_fnames:
# 	s2_readers:
# 	gee_id:
# 	dw_path:
# 	dw_reader:
# 	s2_borders:
# 	dw_borders:
# 	base_chip_id:

# 	Methods
# 	-------
# 	get_band_filenames()
# 	get_gee_id()

# 	'''
# 	def __init__(self,safe_id):
# 		self.id    = safe_id
# 		self.tile  = self.id[38:44]
# 		self.date  = self.id[11:26]
# 		self.orbit = self.id[33:37]

# 		#1.1 ID -> BAND READERS
# 		self.s2_fnames  = self.get_band_filenames() #sorted
# 		self.s2_readers = []
# 		for f in self.s2_fnames:
# 			band_path = f'{WORK_DIR}/{safe_id}/{f}'
# 			if not os.path.isfile(band_path):
# 				raise IncompleteDirError(f"Missing band file {f}")
# 			self.s2_readers += [rio.open(band_path,'r',tiled=True)]

# 		#1.2 ID -> XML PATH
# 		#2.XML -> DW PATH
# 		#3.DW PATH -> DW READER
# 		self.gee_id    = self.get_gee_id()		
# 		self.dw_reader = rio.open(self.dw_path,'r',tiled=True)

# 		#Check label
# 		if self.dw_reader.statistics(1).max == 0:
# 			raise EmptyLabelError("Label is zero everywhere.")

# 		#4.DW READER -> BOUNDS DW
# 		#5.DW READER+BAND2 READER -> BOUNDS S2 & BOUNDS DW
# 		self.s2_borders,self.dw_borders = align(self.s2_readers[0],self.dw_reader)
	
# 		#format: DATE_DSTRIP_TILE_ROTATION_WINROW_WINCOL_B0*.tif
# 		#format: DATE_DSTRIP_TILE_ROTATION_WINROW_WINCOL_LBL.tif	
# 		self.base_chip_id = self.gee_id + '_' + self.orbit

# 	def get_band_filenames(self):
# 		return [f'{self.tile}_{self.date}_{b}_10m.jp2' for b in ['B02','B03','B04','B08']]

# 	def get_gee_id(self):
# 		datastrip = None
# 		return '_'.join([self.date,datastrip,self.tile])


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


def chip_image(s2_readers,label_path,feature_path,base_id,index,N):

	# STDOUT
	print(f'[{index+1}/{N}] PROCESSING {base_id}')
	start_time = time.time()

	# LOAD BAND ARRAYS, CLIP, & NORMALIZE
	rgb = []
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
		band_array  = np.clip(band_array,low_cutoff,high_cutoff)
		band_array  = (band_array/(high_cutoff-low_cutoff)*255).astype(np.uint8)
		band_array  = np.where(zero_mask,0,band_array)
		rgb.append(band_array)

	# SET WINDOWS
	s2_borders = {'top': 492, 'bottom': 10487, 'left': 492, 'right': 10487}
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
			args=(rgb,label_path,feature_path,s2_window_chunks[i],base_id)
		)
		p.start()
		processes.append(p)

	for p in processes:
		p.join(timeout=60)

	# STDOUT	
	exec_time = time.time() - start_time
	print(f"All workers done ({exec_time:.3f} secs). ")


def chip_image_worker(rgb,label_path,feature_path,windows,base_id):

	# Distinct rio.DatasetReader for thread/race conditions
	lbl_rdr = rio.open(label_path,'r',tiled=True) #1 band, uint16
	ftr_rdr = rio.open(feature_path,'r',tiled=True) #3 bands, uint16

	# Log chip info?
	# stats = []

	for k,(rowcol,w) in enumerate(windows):

		# LOAD (ONLY WINDOW SECTION) LABEL & FEATURES
		lbl_array = lbl_rdr.read(1,window=w)
		ftr_array = ftr_rdr.read(1,window=w)

		# IF LABEL NO DATA -- SKIP CHIP
		if (lbl_array == 0).any():
			continue

		# LOAD RGB/IF NO DATA IN RGB SKIP CHIP
		r_array = rgb[0][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]
		if (r_array == 0).any():
			continue
		g_array = rgb[1][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]
		b_array = rgb[2][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]

		# GOOD -- SAVE BANDS
		row,col = rowcol
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_rgb.tif'
		r = Image.fromarray(r_array)
		g = Image.fromarray(g_array)
		b = Image.fromarray(b_array)
		img = Image.merge('RGB',(r,g,b))
		img.save(outfile)

		# SAVE LABEL
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_lbl.tif'
		img = Image.fromarray(lbl_array)
		img.save(outfile)

		# SAVE FEATURES <<<<----- FIX
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_ftr.tif'
		img = Image.fromarray(ftr_array)
		img.save(outfile)		

		# STATS/LOG?
		# diabetes = lbl_array.mean() #or weighted mean.. something like that.
		# stats.append(f'{outfile.split('/')[-1][:-8]}\t{diabetes}')

	# LOG?
	# lock.acquire()
	# # print(f'Worker {mp.current_process()} done.')	
	# with open(f'{CHIP_DIR}/stats.txt','a') as fp:
	# 	fp.write('\n'.join(stats))
	# lock.release()


if __name__ == '__main__':

	########## ARGV CONFIG ##########
	parser = argparse.ArgumentParser(
		prog="chips.py",
		description="Large Sentinel-2 and labels to 224x224 images.")

	# PATHS
	parser.add_argument('--work-dir',default='/cache',
		help="Temporary directory to load/offload data.")
	parser.add_argument('--chip-dir',default=None,
		help="Output directory for resulting chips")
	parser.add_argument('--s2-dir',default=None,
		help="Source directory for raw Sentinel-2 products.")
	parser.add_argument('--label-dir',default=None,
		help="Source directory for 10980x10980 mask rasters.")


	########## SET ARGS ##########
	args = parser.parse_args()
	WORK_DIR  = args.work_dir
	CHIP_DIR  = args.chip_dir
	S2_DIR    = args.s2_dir 
	LABEL_DIR = args.label_dir

	if not os.path.isdir(WORK_DIR):
		print(f"WORK_DIR {WORK_DIR} not found. EXIT(1).")
		sys.exit(1)
	if WORK_DIR[-1] == '/':
		WORK_DIR = WORK_DIR.rstrip('/')

	if CHIP_DIR is None:
		os.makedirs(WORK_DIR + '/chips',exist_ok=True)
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


	########## GET UNIQUE TILES FROM LABEL DIR ###############
	label_tiffs  = glob.glob('*.tif',root_dir=LABEL_DIR) #arg/masks
	unique_tiles = [s.split('_')[0] for s in label_tiffs]

	########## GET PRODUCT INTERSECTION ##########
	band2_regex = "eodata/Sentinel-2/MSI/L2A/*/*/*/*.SAFE/GRANULE/*/IMG_DATA/R10m/*_B02_10m.jp2" #1637
	s2_tiffs         = glob.glob(band2_regex,root_dir=S2_DIR)
	s2_tiles         = [s.split('/')[-1].split('_')[0] for s in s2_tiffs]
	intersection     = np.isin(s2_tiles,unique_tiles)
	s2_good_products = np.array(s2_tiffs)[intersection]
	print(f"PRODUCTS MATCHING LABELS: {len(s2_good_products)}.")

	########## SPLIT AND QUEUE ################
	chunk_size  = 100
	N_chunks    = len(s2_good_products) // chunk_size
	remainder   = len(s2_good_products) % chunk_size
	chunk_queue = []
	for i in range(N_chunks):
		chunk_queue.append(s2_good_products[i*chunk_size:i*chunk_size+chunk_size])
	if remainder != 0:
		chunk_queue.append(s2_good_products[N_chunks*100:])

	########## PROCESS  #######################
	for chunk in chunk_queue:

		base_ids       = []
		tiles_in_chunk = []

		########## DOWNLOAD/COPY ####################
		for b2_path in chunk:

			# GET SOME STRINGS
			b3_path = b2_path.replace("_B02_","_B03_")
			b4_path = b2_path.replace("_B02_","_B04_")
			# bands_regex = '/'.join(b2_path.split('/')[0:-1]) + '/*.jp2' #or [:-34]?
			tile  = b2_path.split('/')[-1].split('_')[0]
			date  = b2_path.split('/')[-1].split('_')[1]
			orbit = b2_path.split('/')[7].split('_')[4]
			base_ids.append(f"{tile}_{date}_{orbit}")
			tiles_in_chunk.append(tile)

			# COPY 3 BANDS
			# sp.run(["cp",f"{S2_DIR}/{bands_regex}",WORK_DIR,"-v"])
			sp.run(["cp",f"{S2_DIR}/{b2_path}",WORK_DIR,"-v","-n"])
			sp.run(["cp",f"{S2_DIR}/{b3_path}",WORK_DIR,"-v","-n"])
			sp.run(["cp",f"{S2_DIR}/{b4_path}",WORK_DIR,"-v","-n"])

		# COPY NECESSARY LABELS
		for t in list(np.unique(tiles_in_chunk)):
			sp.run(["cp",f"{LABEL_DIR}/{t}_diabetes.tif",WORK_DIR,"-v","-n"])
			sp.run(["cp",f"{LABEL_DIR}/{t}_features.tif",WORK_DIR,"-v","-n"])


		########## CHIP ####################
		for i,product in enumerate(chunk):

			# PATHS & READERS
			local_b2_path = product.split('/')[-1]
			local_b3_path = local_b2_path.replace("_B02_","_B03_")
			local_b4_path = local_b2_path.replace("_B02_","_B03_")
			b2_reader = rio.open(f"{WORK_DIR}/{local_b2_path}",'r',tiled=True)
			b3_reader = rio.open(f"{WORK_DIR}/{local_b3_path}",'r',tiled=True)
			b4_reader = rio.open(f"{WORK_DIR}/{local_b4_path}",'r',tiled=True)
			rgb_readers = [b4_reader,b3_reader,b2_reader]
			label_path   = f"{WORK_DIR}/{tiles_in_chunk[i]}_diabetes.tif"
			feature_path = f"{WORK_DIR}/{tiles_in_chunk[i]}_features.tif"

			# CHIP
			chip_image(rgb_readers,label_path,feature_path,base_ids[i],i,len(chunk))


		########## DELETE FILES ############
		sp.run(["rm",f"{WORK_DIR}/*.tif","-v"])
		sp.run(["rm",f"{WORK_DIR}/*.jp2","-v"])


	print("DONE.")