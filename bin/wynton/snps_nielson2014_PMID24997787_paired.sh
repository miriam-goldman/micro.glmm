#!/bin/bash                         #-- what is the language of this shell
#                                   #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/home/pollard/czhao/scratch/snps_out                       #-- output directory (fill in)
#$ -e /wynton/home/pollard/czhao/scratch/snps_err                       #-- error directory (fill in)
#$ -cwd                             #-- tell the job that it should start in your working directory
##$ -r y                            #-- tell the system that if a job crashes, it should be restarted
##$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -pe smp 8
#$ -l mem_free=10G                   #-- submits on nodes with enough free memory (required)
#$ -l arch=lx-amd64                 #-- SGE resources (CPU type)
#$ -l scratch=40G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=16:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-383                        #-- Array job: submit XX jobs to the job queues at once.
#$ -tc 383                         #-- specify the maximum number of concurrent tasks

module load CBI miniconda3/23.3.1-0-py39
conda activate midas2v1.1.0

JB_LAB=sample_${SGE_TASK_ID}
TREADS=${NSLOTS:-1}

echo "Hello world, I’m running on node ${HOSTNAME} for sample ${SGE_TASK_ID}"

proj_name="nielson2014_PMID24997787_paired"

base_dir="/pollard/data/projects/mwas-projects"
midasdb_dir="/pollard/data/midas2-db/midas2db-uhgg-v1.5"

GNUTIME="/wynton/home/pollard/czhao/local/bin/time-1.9/bin/time"
CWDIR="/wynton/home/pollard/czhao/mwas-project/bin"
GB_IN="/wynton/home/pollard/czhao/mwas-project/sge_jobs/${proj_name}/${JB_LAB}"

# Pollard server
proj_dir="${base_dir}/${proj_name}"
data_dir="${proj_dir}/qc/05_combined"
midas2_dir="${base_dir}/${proj_name}/midas2_output"

# local scratch disk as output directory
IN_DIR="/scratch/czhao/${proj_name}/05_combined"
if [ ! -d $IN_DIR ]; then
  mkdir -p $IN_DIR
fi

OUTDIR="/scratch/czhao/${proj_name}/midas2_output"
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi

cd "$OUTDIR" # Use a temporary working directory

while IFS= read -r sample_name
do

  r1="${data_dir}/${sample_name}_1.fastq.gz"
  r2="${data_dir}/${sample_name}_2.fastq.gz"

  lr1="${IN_DIR}/${sample_name}_1.fastq.gz"
  lr2="${IN_DIR}/${sample_name}_2.fastq.gz"

  cp -r ${r1} ${lr1}
  cp -r ${r2} ${lr2}

  cp -r ${midas2_dir}/${sample_name} ${OUTDIR}

  ${GNUTIME} -v python -m midas2 run_snps \
    --sample_name ${sample_name} -1 ${lr1} -2 ${lr2} \
    --midasdb_name newdb --midasdb_dir ${midasdb_dir} \
    --num_cores ${TREADS} \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=2,0.5 \
    --fragment_length 1400 --fragment_ratio 1 \
    --paired_only --ignore_ambiguous --advanced \
    ${OUTDIR} &> ${OUTDIR}/${sample_name}/snps_time.log

  cp -r ${OUTDIR}/${sample_name} ${midas2_dir}
  rm -f ${lr1} ${lr2}

done < ${GB_IN}

SUB_DIR=${CWDIR}
if ls -A ${SUB_DIR}/ | grep core; then
  rm -v ${SUB_DIR}/core.*
else
  echo “no core dumping detected”
fi

conda deactivate
