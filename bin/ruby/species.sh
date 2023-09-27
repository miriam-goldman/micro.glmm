#!/bin/bash
# Illumina MIDAS2's Species Workflow
# Chunyu Zhao

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 PROJECT_NAME SAMPLE_NAME"
    exit 1
fi

#proj_name="he2017_PMID28655159"

proj_name="$1"
sample_name="$2"

base_dir="/pollard/data/projects/mwas-projects"
midasdb_dir="/pollard/data/midas2-db/midas2db-uhgg-v1.5"

proj_dir="${base_dir}/${proj_name}"
data_dir="${proj_dir}/qc/03_decontam"
out_dir="${proj_dir}/midas2_output"

r1="${data_dir}/${sample_name}_1.fastq.gz"
r2="${data_dir}/${sample_name}_2.fastq.gz"

midas2 run_species --sample_name ${sample_name} \
  -1 ${r1} -2 ${r2} --num_cores 8 \
  --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
  ${out_dir}

echo "${sample_name} DONE"
