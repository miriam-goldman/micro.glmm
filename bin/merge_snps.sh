#!/bin/bash
# MIDAS2's SNPs Workflow on AWS
# Chunyu Zhao

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 SPECIES_ID"
    exit 1
fi

species_id="$1"

base_dir="/mnt/cz/IBD"
midasdb_dir="/mnt/cz/midas2db-uhgg-v1.5"
GNUTIME="/home/ubuntu/bin/time-1.9/bin/time"

data_dir="/mnt/cz/samples_list_snps"

OUTDIR="/mnt/cz/merge_snps/${species_id}"
mkdir -p ${OUTDIR}

cp $data_dir/${species_id}.txt ${OUTDIR}/samples.list

## GOOD, double check the sample list, because I don't want those samples with genome_over_marker_coverage
${GNUTIME} -v python -m midas2 merge_snps \
  --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
  --samples_list $data_dir/${species_id}.txt \
  --species_list ${species_id} \
  --genome_depth 5.0 --genome_coverage 0.4 --sample_counts 40 \
  --site_prev 0.7 --snp_pooled_method prevalence \
  --site_depth 5 --site_ratio 5 --snp_maf 0.05 \
  --chunk_size 100000 --advanced \
  --num_cores 32 \
  ${OUTDIR} 2> ${OUTDIR}/merge_snps.log
