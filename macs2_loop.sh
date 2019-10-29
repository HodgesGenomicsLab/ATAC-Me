#!/bin/bash

# grab the files, and export it so the 'child' sbatch jobs can access it

# Directory location of bam files output from walt after conversion from sam file type
BAMDIR=/data/hodges_lab/EH4889/stand_atac/trimmed_reads/mapped_reads
export FILESBAM=($(ls -1 ${BAMDIR}/*.hum))

# Output directory of macs2 peak calls
export MACS_DIR=/data/hodges_lab/EH4889/revisions/macs
if test -d ${MACS_DIR}; then  echo "exist"; else mkdir ${MACS_DIR} && echo created; fi

# get size of array
NUMBAM=${#FILESBAM[@]}
# now subtract 1 as we have to use zero-based indexing (first cell is 0)
ZBNUMBAM=$(($NUMBAM - 1))
for i in `seq 0 $ZBNUMBAM`;
do
echo ${FILESBAM[$i]}
export BAM=${FILESBAM[$i]}
sbatch macs2_comparison_stand.slrm
done
