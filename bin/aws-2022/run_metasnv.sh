#!/bin/bash
set -e

outdir="/mnt/cz/ibd/metasnv_output"
repfa="/mnt/cz/ibd/bowtie2_indexes/repgenomes.fa"
lob="/mnt/cz/ibd/lob"

#--min_pos_cov: minimum coverage (mapped reads) per position for snpCall. (4) (sum across all samples)
#--min_pos_snvs: minimum number of non-reference nucleotides per position for snpCall. (4)
#rm -rf ${outdir}

mkdir -p /mnt/cz/ibd/metasnv_log

#~/bin/time-1.9/time -v metaSNV.py --min_pos_snvs 5 --min_pos_snv 2 --threads 32 ${outdir} ${lob} ${repfa} &> /mnt/cz/ibd/metasnv_log/snv.log

~/bin/time-1.9/time -v metaSNV_Filtering.py -b 0.4 -d 5.0 -m 44 -c 5.0 --ind -p 0.8 --n_threads 32 ${outdir} &> /mnt/cz/ibd/metasnv_log/filter.log

#-b FLOAT: Coverage breadth: minimal horizontal genome coverage percentage per sample per species (default: 40.0)
#-d FLOAT: Coverage depth: minimal average vertical genome coverage per sample per species (default: 0.0)
#-m INT: Minimum number of samples per species (default: 1)
#-c FLOAT: FILTERING STEP II:minimum coverage per position per sample per species (default: 5.0)
#-p FLOAT: FILTERING STEP II:required proportion of informative samples (coverage non-zero) perposition (default: 0.5)


