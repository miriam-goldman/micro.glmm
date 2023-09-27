#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 SAMPLE"
    exit 1
fi

sample_name="$1"

basedir="/mnt/cz/ibd/instrain_output"
mkdir -p ${basedir}

midasdir="/mnt/cz/ibd/midas2_output"

~/bin/time-1.9/time -v inStrain profile -o $basedir/$sample_name -p 32 \
  -l 0.94 --min_mapq 10 --pairing_filter paired_only -c 5 -f 0.05 \
  --skip_plot_generation --window_length 1000000 --skip_mm_profiling \
  -s /mnt/cz/ibd/bowtie2_indexes/repgenomes.stb --database_mode  \
  /mnt/cz/ibd/raw_bamfiles/${sample_name}.bam \
  /mnt/cz/ibd/bowtie2_indexes/repgenomes.fa &> /mnt/cz/ibd/instrain_log/${sample_name}_32.log
