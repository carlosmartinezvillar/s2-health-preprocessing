#!/bin/bash
#THIS SCRIPT CONVERTS SENTINEL-2 IMAGES INTO
#A DATASET OF INPUT-LABEL CHIP PAIRS.
#IT DOES SO BY REPEATEDLY TRANSFERING A PORTION OF THE ORIGINAL 
#SIZE BANDS FROM IN AN S3 BUCKET INTO A TEMPORARY DIRECTORY
#(OFTEN ATTACHED TO A CONTAINER), SPLITTING THEM, AND PUSHING
#THE CHIPS TO A DIFFERENT S3 BUCKET.

DATA_DIR=/cache
N_TRANSFERS=16
DEST_IMAGE_BUCKET=
ORIG_IMAGE_BUCKET=nrp:diabetes-s2

#====================================================================================================
# SET CHUNKS OF .SAFE FOLDERS TO DOWNLOAD
#====================================================================================================
# Rclone dir to list
string_list=$(rclone lsf ${ORIG_BUCKET} | grep .SAFE | awk '{print substr($0,1,length($0)-1)}')
array=(${string_list})

# Chunk math
chunk_size=50
n_chunks=$(((${#array[@]}) / chunk_size))
remainder=$((${#array[@]} - chunk_size * n_chunks))

#====================================================================================================
# DOWNLOAD WHOLE LABEL SET IN REMOTE 
#====================================================================================================
rclone copy ${ORIG_BUCKET}/dynamicworld ${DATA_DIR}/dynamicworld -P --transfers ${N_TRANSFERS}

#====================================================================================================
# ITERATE--DOWNLOAD, CHIP & PUSH
#====================================================================================================

# EVEN CHUNKS OF BAND FOLDERS
#----------------------------- 
for (( i=0; i<n_chunks; i++ )); do
	echo "BATCH: ${i}"

	#0. SAVE CHUNK LIST TO TEMP FILE
	start=$((i * chunk_size))
	chunk=("${array[@]:$start:$chunk_size}")
	printf "%s/**\n" "${chunk[@]}" > "${DATA_DIR}/chunk.txt"

	#1. DOWLOAD/READ CHUNK
	rclone copy --include-from ${DATA_DIR}/chunk.txt ${ORIG_BUCKET} ${DATA_DIR} -P --transfers ${N_TRANSFERS}

	#2. CHIP IMAGES
	python3 chipping.py --data-dir ${DATA_DIR} --chip-dir ${DATA_DIR}/chips --chips

	#3. CLEAN UP .SAFE
	rm -r ${DATA_DIR}/*.SAFE

	#4. PUSH CHIPS (SO FAR) TO S3
	rclone copy ${DATA_DIR}/chips ${DEST_BUCKET} -P --transfers ${N_TRANSFERS}

	#CLEAN UP TIFFs
	rm  ${DATA_DIR}/chips/*.tif
done

# REMAINING BATCH
#--------------------- 
echo "BATCH: ${n_chunks}"

#0. CHUNK LIST TO TEMP FILE
start=$((n_chunks * chunk_size))
chunk=("${array[@]:$start:$remainder}")
printf "%s/**\n" "${chunk[@]}" > "${DATA_DIR}/chunk.txt"

#1. .SAFE TO LOCAL CACHE
rclone copy --include-from ${DATA_DIR}/chunk.txt ${ORIG_BUCKET} ${DATA_DIR} -P	--transfers ${N_TRANSFERS}

#2. CHIP IMAGES
python3 chipping.py --data-dir ${DATA_DIR} --chip-dir ${DATA_DIR}/chips --chips

#3. CLEAN UP -- REMOVE .SAFE
rm -r ${DATA_DIR}/*.SAFE

#4. PUSH CHIPS (SO FAR) TO S3
rclone copy ${DATA_DIR}/chips ${DEST_BUCKET} -P --transfers ${N_TRANSFERS}

#CLEAN UP CHIP TIFFs AND LIST IN WORKING DIR
rm  ${DATA_DIR}/chips/*.tif
rm ${DATA_DIR}/chunk.txt
