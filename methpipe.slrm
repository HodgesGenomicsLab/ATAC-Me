#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --output=MethPipe.out
#SBATCH --error=MethPipe.error
#SBATCH --mail-user=yourEmail@university.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="MethPipe"

# Necessary to set this for proper sorting
export LC_ALL=C

# Sets up naming scheme for output files
MAPPED=`basename $SAM _val_1.fq.sam`

# Loads samtools via LMOD from cluster computing environment
module load Intel/2016.3.210-GCC-5.4.0-2.26 SAMtools/1.6

# Convert sam output from walt to bam filetype
samtools view -Sb  ${SAM}  >  ${PROCESSED_DIR}/${MAPPED}.bam

# Convert bam output to mapped read format used by methpipe tools
to-mr -o ${PROCESSED_DIR}/${MAPPED}.mr -m general ${PROCESSED_DIR}/${MAPPED}.bam

# Sort mapped read format file
sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -T ${PROCESSED_DIR} ${PROCESSED_DIR}/${MAPPED}.mr -o ${PROCESSED_DIR}/${MAPPED}.sort.mr

# Remove duplicate reads 
duplicate-remover -s -A -S ${PROCESSED_DIR}/${MAPPED}.dupstats -o ${PROCESSED_DIR}/${MAPPED}.mr.dremove ${PROCESSED_DIR}/${MAPPED}.sort.mr

# Estimate bisulfite conversion rate for sequencing library 
bsrate  -c ${CHROM_DIR} ${PROCESSED_DIR}/${MAPPED}.mr.dremove -o ${PROCESSED_DIR}/${MAPPED}.bsrate

# Generate epiread formate files for use in other downstream methpipe tools such as finding DMRs
methstates -c ${CHROM_DIR} ${PROCESSED_DIR}/${MAPPED}.mr.dremove -o ${PROCESSED_DIR}/${MAPPED}.epiread

# Calculate methylation rate for all individual cytosines in the genome
methcounts -c ${CHROM_DIR} -o ${PROCESSED_DIR}/${MAPPED}.all.meth ${PROCESSED_DIR}/${MAPPED}.mr.dremove

# Calculate statistics on the coverage and methylation levels for different context cytosines
levels -o ${PROCESSED_DIR}/${MAPPED}.levels  ${PROCESSED_DIR}/${MAPPED}.all.meth

# Get methylation rate for only symmetric CpG context cytosines
symmetric-cpgs -m -o ${PROCESSED_DIR}/${MAPPED}.meth ${PROCESSED_DIR}/${MAPPED}.all.meth

# Scan for hypomethylated regions
# This tool is geared primarily for WGBS data
# Still in development for use with ATACMe data
hmr -o ${PROCESSED_DIR}/${MAPPED}.hmr  -p ${PROCESSED_DIR}/${MAPPED}.hmrparams  ${PROCESSED_DIR}/${MAPPED}.meth

# Generate bigwig of symmetric CpG covering reads
awk '{OFS="\t"; print $1,$2,$2+1,$6}' < ${PROCESSED_DIR}/${MAPPED}.meth | wigToBigWig /dev/stdin ${SIZES_DIR} ${TRACK_DIR}/${MAPPED}.read.bw

# Generate bigwig of methylation levels for symmetric CpGs
awk '{OFS="\t"; print $1,$2,$2+1,$5}' < ${PROCESSED_DIR}/${MAPPED}.meth | wigToBigWig /dev/stdin ${SIZES_DIR} ${TRACK_DIR}/${MAPPED}.meth.bw

# Generate bigBed file of genomic intervals called as hypomethylated regions
cut -f 1-3 ${PROCESSED_DIR}/${MAPPED}.hmr > ${PROCESSED_DIR}/${MAPPED}.hmr.tmp
bedToBigBed ${PROCESSED_DIR}/${MAPPED}.hmr.tmp ${SIZES_DIR} ${TRACK_DIR}/${MAPPED}.hmr.bb && rm ${PROCESSED_DIR}/${MAPPED}.hmr.tmp
