#!/bin/bash
# MIDAS2's SNPs Workflow on Ruby
# Chunyu Zhao

set -e

if [ $# -ne 2  ]; then
    echo "Usage: $0 SPECIES_ID cluster_pid"
    exit 1
fi

species_id="$1"
cluster_pid="$2"

base_dir="/pollard/data/projects/mwas-projects"
midasdb_dir="/pollard/data/midas2-db/midas2db-uhgg-v1.5"

out_dir="${base_dir}/merge_genes_cid.${cluster_pid}/${species_id}"
mkdir -p $out_dir

los="${base_dir}/samples_list_genes/${species_id}.txt"
cp $los $out_dir

midas2 merge_genes \
  --samples_list $los \
  --species_list ${species_id} \
  --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
  --num_cores 16 \
  --cluster_pid $cluster_pid \
  --min_copy 0.35 \
  --sample_counts 40 \
  --genome_depth 2 \
  ${out_dir}
