FQDIR=/data/hodges_lab/EH4889/stand_atac/trimmed_reads
export FILESR1=($(ls -1 ${FQDIR}/*_val_1.fq))
export FILESR2=($(ls -1 ${FQDIR}/*_val_2.fq))
export INDEX=/home/hodgese/data/drosophila_genome/bt2/hg19+dm3
export MAPPED_DIR=/data/hodges_lab/EH4889/stand_atac/trimmed_reads/mapped_reads
if test -d ${MAPPED_DIR}; then  echo "exist"; else mkdir ${MAPPED_DIR} && echo created; fi

# get size of array
NUMFASTQ=${#FILESR1[@]}
# now subtract 1 as we have to use zero-based indexing (first cell is 0)
ZBNUMFASTQ=$(($NUMFASTQ - 1))
for i in `seq 0 $ZBNUMFASTQ`;
do
echo ${FILESR1[$i]}
export R1=${FILESR1[$i]}
export R2=${FILESR2[$i]}
sbatch stATAC_bowtie2.slrm
done
