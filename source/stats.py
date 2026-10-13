import os
import numpy as np
import tifffile as tiff
import argparse
import glob
from concurrent.futures import ThreadPoolExecutor
import csv
import itertools

N_PROC = 8
BANDS = ["B01","B02","B03","B04","B05","B06","B07","B8A","B11","B12"]
N_BANDS = len(BANDS)


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--chip-dir',default=None,required=True,
		help='Dataset directory.')
	parser.add_argument('--split-csv',default=None,required=True,
		help='CSV with train/val/test split details.')	
	args = parser.parse_args()

	assert os.path.isdir(args.chip_dir), f"Dir path not found in {args.chip_dir}"
	args.chip_dir = args.chip_dir.rstrip('/')

	assert os.path.isfile(args.split_csv), f"CSV not found in {args.split_csv}" 

	return args


def sum_band_chip(chip_path):

	bands = tiff.imread(chip_path).astype(np.float64)
	pixels = bands.reshape(N_BANDS,-1)

	sum_   = pixels.sum(axis=1)
	sum_sq = (pixels**2).sum(axis=1)
	count  = pixels.shape[1]
	return sum_,sum_sq,count	


def sum_tract_px(chip_path):
	features = tiff.imread(chip_path).astype(np.uint16)
	tracts,counts = np.unique(features[...,3],return_counts=True)
	return zip((tracts.astype(int)-1).tolist(),counts.tolist())


def calculate_band_mean_std(chip_dir,split_csv):

	# LOAD TRAIN/VAL/TEST SPLIT
	with open(split_csv,'r',encoding='utf-8') as fp:
		reader = csv.reader(fp)
		header = next(reader)
		rows = [_ for _ in reader]

	# ADJUST ROWS
	table = np.array(rows)
	split = table[:,2]
	full_paths = [f"{chip_dir}/{r[0]}_rgb.tif" for r in table]
	split_paths = np.array(full_paths)[split == 'train'].tolist()

	# READ CHIPS IN SPLIT
	with ThreadPoolExecutor(N_PROC) as ex:
		rows = list(ex.map(sum_band_chip,split_paths))

	# UNWRAP
	sums, sq_sums, px_counts = zip(*rows)

	# ADD INDIVIDUAL COUNTERS
	overall_sum =  np.stack(sums).sum(axis=0)
	overall_sq_sum = np.stack(sq_sums).sum(axis=0)
	overall_px_count = sum(px_counts)

	# GET
	means   = overall_sum / overall_px_count
	stddevs = np.sqrt(overall_sq_sum/overall_px_count - means ** 2)

	print("MEANS: ")
	print(f"[{','.join([str(b) for b in means])}]")
	print("STD DEVS: ")
	print(f"[{','.join([str(b) for b in stddevs])}]")


def calculate_px_per_tract(chip_dir,split_csv):
	
	# LOAD SPLIT
	with open(split_csv,'r',encoding='utf-8') as fp:
		reader = csv.reader(fp)
		header = next(reader)
		rows = [_ for _ in reader]
	table = np.array(rows)

	# ADJUST PATHS
	feat_dir = os.path.join(os.path.dirname(chip_dir),'features')
	feat_ids = table[:,0][table[:,2]!='dropped'].to_list()
	chip_paths = [f"{feat_dir}/{i}_ftr.tif" for i in feat_ids]

	# COUNT PER CHIP
	with ThreadPoolExecutor(N_PROC) as ex:
		nested_list = list(ex.map(sum_tract_px,chip_paths))

	# MERGE AND ADD
	n_tracts = 18407 #fix to read from filtered_diabetes.csv
	single_list = list(itertools.chain.from_iterable(nested_list))
	tract_counts = np.zeros(n_tracts).astype(np.float64) 

	for tract,count in single_list:
		tract_counts[tract] += count


if __name__ == '__main__':
	args = parse_args()
	calculate_band_mean_std(args.chip_dir,args.split_csv)