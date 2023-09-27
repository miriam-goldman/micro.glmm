#!/bin/sh
# properties = {properties}

# source ./env
#
# mkdir -p /scratch/bsmith
#

{exec_job}
_status=$?

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2
exit $_status
