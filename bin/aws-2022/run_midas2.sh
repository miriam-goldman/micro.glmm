#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 SAMPLE"
    exit 1
fi

sample_name="$1"

basedir="/mnt/cz/ibd/midas2_output"
mkdir -p ${basedir}

mkdir -p /mnt/cz/ibd/midas2_log

~/bin/time-1.9/time -v python -m iggtools midas_run_snps \
 --sample_name ${sample_name} \
 -1 /mnt/cz/ibd/qc/${sample_name}_1.fastq.gz -2 /mnt/cz/ibd/qc/${sample_name}_2.fastq.gz \
 --midasdb_name uhgg --midasdb_dir /mnt/cz/ibd/midasdb_v1.2 \
 --num_cores 32 --select_threshold=-0.5 --aln_mapq 10 --fragment_length 1000 --paired_only \
 --prebuilt_bowtie2_indexes /mnt/cz/ibd/bowtie2_indexes/repgenomes \
 --prebuilt_bowtie2_species /mnt/cz/ibd/bowtie2_indexes/repgenomes.species \
 --site_depth 5 --chunk_size 1000000 ${basedir} &> /mnt/cz/ibd/midas2_log/${sample_name}_32.log
