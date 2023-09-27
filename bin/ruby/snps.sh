#!/bin/bash
# MIDAS2's SNPs Workflow on Ruby
# Chunyu Zhao

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 PROJECT_NAME SAMPLE_NAME"
    exit 1
fi

proj_name="$1"
sample_name="$2"

base_dir="/pollard/data/projects/mwas-projects"
midasdb_dir="/pollard/data/midas2-db/midas2db-uhgg-v1.5"

proj_dir="${base_dir}/${proj_name}"
data_dir="${proj_dir}/qc/03_decontam"
midas2_dir="${proj_dir}/midas2_output"

r1="${data_dir}/${sample_name}_1.fastq.gz"
r2="${data_dir}/${sample_name}_2.fastq.gz"

scratch_dir="/pollard/scratch/czhao/midas2_output"
OUTDIR="${scratch_dir}/${proj_name}"
mkdir -p ${OUTDIR}

lr1="${OUTDIR}/${sample_name}_1.fastq.gz"
lr2="${OUTDIR}/${sample_name}_2.fastq.gz"

cp -r ${r1} ${lr1}
cp -r ${r2} ${lr2}

cp -r ${midas2_dir}/${sample_name} ${OUTDIR}

midas2 run_snps \
  --sample_name ${sample_name} -1 ${lr1} -2 ${lr2} \
  --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
  --num_cores 8 \
  --select_by median_marker_coverage,unique_fraction_covered \
  --select_threshold=2,0.5 \
  --fragment_length 1400 --fragment_ratio 1 \
  --paired_only --ignore_ambiguous --advanced \
  ${OUTDIR}

cp -r ${OUTDIR}/${sample_name} ${midas2_dir}
