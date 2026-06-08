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

# Typing
# from typing import Tuple, List
# ndarray = np.ndarray

DATA_DIR  = None
LABEL_DIR = None
CHIP_DIR  = None
#SET DIRS HERE BECAUSE THREAD ACCESS

# PIXEL LIMITS
CHIP_SIZE = 224
STRIDE    = 112

# NR OF PROCESSES PER RASTER
N_PROC    = 16

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
# 			band_path = f'{DATA_DIR}/{safe_id}/{f}'
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
def remove_dynamicworld_borders(src: rio.DatasetReader) -> dict:
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


def align_dynamicworld(s2_src: rio.DatasetReader,dw_src: rio.DatasetReader) -> Tuple:
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
	'''

	# number of rows and cols takin' the boundaries into acct
	n_px_rows = borders['bottom'] + 1 - borders['top']
	n_px_cols = borders['right'] + 1 - borders['left']

	#nr of blocks in each direction
	block_rows = n_px_rows // CHIP_SIZE
	block_cols = n_px_cols // CHIP_SIZE

	#total blocks
	N = block_rows * block_cols

	windows = []

	for k in range(N):
		i = k // block_cols
		j = k % block_cols
		row_start = i * CHIP_SIZE + borders['top']
		col_start = j * CHIP_SIZE + borders['left']
		W = Window(col_start,row_start,CHIP_SIZE,CHIP_SIZE)
		windows += [[(str(i),str(j)),W]]

	return windows


def chip_image(s2_readers,label_reader,features_reader,index,N):
	print(f'[{index}/{N-1}] PROCESSING {product.id} ')
	start_time = time.time()

	# LOAD ARRAYS AND NORMALIZE BANDS
	rgb = []
	for reader in s2_readers:
		band_array  = reader.read(1)
		zero_mask   = band_array == 0
		high_cutoff = int(np.percentile(band_array[~zero_mask],99)) # DO NOT PASS FLOAT TO CLIP HERE!!!
		low_cutoff  = int(np.percentile(band_array[~zero_mask],1)) #This might have to be lower?
		band_array  = np.clip(band_array,low_cutoff,high_cutoff)
		band_array  = (band_array/(high_cutoff-low_cutoff)*255).astype(np.uint8)
		band_array  = np.where(zero_mask,0,band_array)
		rgb.append(band_array)

	# SET WINDOWS
	s2_borders = {'top': 492, 'bottom': 10487, 'left': 492, 'right': 10487}
	s2_windows = get_windows_strided(s2_borders)

	# SPLIT WINDOWS INTO WORKER SECTIONS
	process_share = len(s2_windows) // N_PROC
	leftover      = len(s2_windows) % N_PROC
	start         = [i*process_share for i in range(N_PROC)]
	stop          = [i*process_share+process_share for i in range(N_PROC)]
	stop[-1]      += leftover
	s2_window_chunks = [s2_windows[s0:s1] for s0,s1 in zip(start,stop)]

	#THROW WORKERS AT ARRAYS
	lock = mp.Lock()
	processes = []
	for i in range(N_PROC):
		p = mp.Process(
			target=chip_image_worker,
			args=(rgb,label_reader,feature_reader,s2_window_chunks[i],base_id,lock)
			)
		p.start()
		processes.append(p)

	for p in processes:
		p.join(timeout=60)
	print("All workers done. ",end='')
	exec_time = time.time() - start_time
	print(f"({exec_time:.3f} seconds).")


def chip_image_worker(rgbn,label_reader,feature_reader,windows,base_id,lock):

	stats = []
	# lbl_rdr = rio.open(dw_path,'r',tiled=True)

	for k,(rowcol,w) in enumerate(s2_windows):

		lbl_arr = lbl_rdr.read(1,window=windows[k][1])

		# CHECK LABEL NO DATA
		if (lbl_arr == 0).any():
			continue

		r_array = rgb[0][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]

		if (r_array == 0).any():
			continue

		g_array = rgb[1][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]
		b_array = rgb[2][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE]

		diabetes_prevalence = lbl_arr.mean()

		# ALL GOOD -- SAVE BANDS IN SINGLE [R,G,B,NIR] FILE (NIR stored in alpha)
		row = rowcol[0]
		col = rowcol[1]
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_rgb.tif'
		r = Image.fromarray(r_array)
		g = Image.fromarray(g_array)
		b = Image.fromarray(b_array)
		img = Image.merge('RGB',(r,g,b))
		img.save(outfile)

		# SAVE LABEL
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_lbl.tif'
		img = Image.fromarray(lbl_arr)
		img.save(outfile)


		# SAVE FEATURES
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_fea.tif'
		img = Image.fromarray(fea_arr)
		img.save(outfile)		

		stats.append(f'{outfile.split("/")[-1]}\t{diabetes_prevalence}')

	# LOG
	lock.acquire()
	# print(f'Worker {mp.current_process()} done.')	
	with open(f'{CHIP_DIR}/stats.txt','a') as fp:
		fp.write('\n'.join(stats))
	lock.release()


if __name__ == '__main__':

	########## ARGV CONFIG ##########
	parser = argparse.ArgumentParser(
		prog="chips.py",
		description="Chip Sentinel-2 and labels to 224x224 images.")

	# PATHS
	parser.add_argument('--data-dir',default='./dat',
		help="Dataset directory")
	parser.add_argument('--chip-dir',
		help="Chip directory")

	# LOAD
	args = parser.parse_args()

	########## SET ARGS ##########
	DATA_DIR  = args.data_dir
	LABEL_DIR = DATA_DIR+'/dynamicworld' #<-- fix this at some point...
	CHIP_DIR  = args.chip_dir

	if not os.path.isdir(DATA_DIR):
		print("DATA_DIR not found. EXITING.")
		sys.exit()
	print(f"DATA_DIR:  {DATA_DIR}")	

	if len(glob.glob('*.SAFE',root_dir=DATA_DIR)) == 0 :
		print("EMPTY DATA_DIR")
		sys.exit()

	print(f"LABEL_DIR set to: {LABEL_DIR}")
	print(f"CHIP_DIR set to:  {CHIP_DIR}")


	#.SAFE folders in data directory
	folders = glob.glob('*.SAFE',root_dir=DATA_DIR)
	paths   = glob.glob(DATA_DIR+'/*.SAFE')

	# Check everything is there
	if not os.path.isdir(LABEL_DIR):
		print("LABEL_DIR not found. EXITING.")
		sys.exit()

	#make chip dir if not already there
	if not os.path.isdir(CHIP_DIR):
		os.mkdir(CHIP_DIR) 

	# clean log file
	if os.path.isfile(f"{CHIP_DIR}/stats.txt"):
		os.remove(CHIP_DIR+'/stats.txt')


	########## PROCESS .SAFE FOLDERS ##########
	N = len(folders)
	for i,f in enumerate(folders):
		try:
			product = Product(f) #load metadata
		except (EmptyLabelError,IncompleteDirError) as e:
			print(f'ERROR: {e}')
			print(f'---> SKIPPING {f}')
			with open(f'{CHIP_DIR}/errored.txt','a') as fp:
				fp.write(f'{f}\n')
			continue

		# <----- CHIP ----->
		chip_image(product,i,N)

	print("DONE.")