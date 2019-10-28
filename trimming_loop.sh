#!/bin/bash

# Line 1 of this script must not be altered so that the shell interpreter knows what to do
# grab the files, and export it so the 'child' sbatch jobs can access it
# Set directory to where fastq files are located
FQDIR=/data/hodges_lab/EH4889/fastqs

# Export read1 files
export FILESR1=($(ls -1 ${FQDIR}/R1_WGBS_THP1_PMA_Rep1A_0.fastq.gz))

# Export read2 files
export FILESR2=($(ls -1 ${FQDIR}/R2_WGBS_THP1_PMA_Rep1A_0.fastq.gz))

# Set directory where trimmed fastq files will be output
export TRIM_DIR=/data/hodges_lab/EH4889/fastqs/trimmed_reads

# Set directory where fastqc files will be output
export FASTQC_DIR=/data/hodges_lab/EH4889/fastqs/trimmed_reads_fastqc

# Check to see if the trimming and fastqc output directories have been created
# If not, create the directories
if test -d ${TRIM_DIR}; then  echo "exist"; else mkdir ${TRIM_DIR} && echo created; fi
if test -d ${FASTQC_DIR}; then  echo "exist"; else mkdir ${FASTQC_DIR} && echo created; fi

# get size of array
NUMFASTQ=${#FILESR1[@]}
# now subtract 1 as we have to use zero-based indexing (first cell is 0)

ZBNUMFASTQ=$(($NUMFASTQ - 1))
for i in `seq 0 $ZBNUMFASTQ`;
do
echo ${FILESR1[$i]}
export R1=${FILESR1[$i]}
export R2=${FILESR2[$i]}

# sbatch submits the SLURM job scheduler jobs
# trim.slrm is the job script which indicates the actual commands you want to run
sbatch trim.slrm
done
