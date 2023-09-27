#!/bin/bash                         #-- what is the language of this shell
#                                   #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/home/pollard/czhao/scratch/snps_out                       #-- output directory (fill in)
#$ -e /wynton/home/pollard/czhao/scratch/snps_err                       #-- error directory (fill in)
#$ -cwd                             #-- tell the job that it should start in your working directory
##$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -pe smp 16
#$ -l mem_free=12G                   #-- submits on nodes with enough free memory (required)
#$ -l arch=lx-amd64                 #-- SGE resources (CPU type)
#$ -l scratch=50G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-3                        #-- Array job: submit XX jobs to the job queues at once.
#$ -tc 3                         #-- specify the maximum number of concurrent tasks

module load CBI miniconda3/23.3.1-0-py39
conda activate midas2v1.1.0

JB_LAB=species_${SGE_TASK_ID}
TREADS=${NSLOTS:-1}

echo "Hello world, I’m running on node ${HOSTNAME} for species ${SGE_TASK_ID}"

base_dir="/pollard/data/projects/mwas-projects"
midasdb_dir="/pollard/data/midas2-db/midas2db-uhgg-v1.5"

GNUTIME="/wynton/home/pollard/czhao/local/bin/time-1.9/bin/time"
CWDIR="/wynton/home/pollard/czhao/mwas-project/bin"
GB_IN="/wynton/home/pollard/czhao/mwas-project/sge_jobs/${JB_LAB}"

read -r SPECIES_ID < ${GB_IN}/species

OUTDIR="/scratch/czhao/midas2_output_${SPECIES_ID}"
if [ ! -d $OUTDIR ]; then
  mkdir -p ${OUTDIR}
fi

cd "$OUTDIR" # Use a temporary working directory

cp ${GB_IN}/* ${OUTDIR}

list_of_samples="${OUTDIR}/los"
cat ${list_of_samples} | xargs -Ixx -P 4 bash -c "mkdir -p ${OUTDIR}/run/xx/snps"

list_of_files="${OUTDIR}/lof"
cat ${list_of_files} | xargs -l -P 4 bash -c 'cp $0 $1'

${GNUTIME} -v python -m midas2 merge_snps \
  --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
  --samples_list ${OUTDIR}/samples_list \
  --species_list ${SPECIES_ID} \
  --genome_depth 5.0 --genome_coverage 0.4 --sample_counts 40 \
  --site_prev 0.7 --snp_pooled_method prevalence \
  --site_depth 5 --site_ratio 5 --snp_maf 0.05 \
  --chunk_size 10000 --advanced \
  --num_cores ${TREADS} \
  ${OUTDIR} 2> ${OUTDIR}/merge_snps_${SPECIES_ID}.log

DEST="/pollard/home/czhao/mwas/merge_snps"
mkdir -p $DEST
cp -r ${OUTDIR}/merge/snps/${SPECIES_ID} ${DEST}
cp ${OUTDIR}/merge_snps_${SPECIES_ID}.log ${DEST}/${SPECIES_ID}
cp -r ${GB_IN} ${DEST}/${SPECIES_ID}


SUB_DIR=${CWDIR}
if ls -A ${SUB_DIR}/ | grep core; then
  rm -v ${SUB_DIR}/core.*
else
  echo “no core dumping detected”
fi

conda deactivate
