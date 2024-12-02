#!/bin/bash
# run: export EXPERIMENT='...'; export CONTIGMAP='...'; sbatch /opt/xchem-fragalysis-2/mferla/library_making/third_pass.slurm.sh

#SBATCH --job-name=rescore
#SBATCH --chdir=/opt/xchem-fragalysis-2/mferla
#SBATCH --output=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
#SBATCH --error=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
#SBATCH --partition=main
#SBATCH --cpus-per-task=6
##SBATCH --nodes=1
##SBATCH --exclusive
##SBATCH --priority=0
#SBATCH --nice=10000
#SBATCH --export=NONE,EXPERIMENT,SLACK_WEBHOOK

# -------------------------------------------------------

export SUBMITTER_HOST=$HOST
export HOST=$( hostname )
export USER=${USER:-$(users)}
export DATA=/opt/xchem-fragalysis-2
export HOME2=$DATA/mferla;
export HOME=$HOME2;
source /etc/os-release;

echo "Running $SLURM_JOB_NAME ($SLURM_JOB_ID) as $USER in $HOST which runs $PRETTY_NAME submitted from $SUBMITTER_HOST"
echo "Request had cpus=$SLURM_JOB_CPUS_PER_NODE mem=$SLURM_MEM_PER_NODE tasks=$SLURM_NTASKS jobID=$SLURM_JOB_ID partition=$SLURM_JOB_PARTITION jobName=$SLURM_JOB_NAME"
echo "Started at $SLURM_JOB_START_TIME"
echo "job_pid=$SLURM_TASK_PID job_gid=$SLURM_JOB_GID topology_addr=$SLURM_TOPOLOGY_ADDR home=$HOME cwd=$PWD"

export CONDA_PREFIX=$DATA/mferla/waconda-slurm
export CONDA_QUIET=true
export CONDA_YES=true

source $CONDA_PREFIX/etc/profile.d/conda.sh
#export CONDA_ENVS_PATH="$CONDA_PREFIX/envs:$DATA/mferla/waconda/envs"

conda activate base;  # not RFdiffusion2

echo "************************"
echo "HELLO $SLURM_JOB_NAME!"
echo "************************"
echo "Greetings from $SLURM_JOB_NAME script ${0} as $USER in $HOST which runs on $PRETTY_NAME with $CONDA_PREFIX"

# END OF FLUFF
# -------------------------------------------------------
# ACTUAL COMMAND

echo $EXPERIMENT;
if [ -z "$EXPERIMENT" ]; then
  echo "Error: EXPERIMENT is not set."
  curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' missing EXPERIMENT var"}' $SLACK_WEBHOOK
  exit 1
fi
export OUTPUT_FOLDER=$HOME2/crysalin-redux/output/$EXPERIMENT;
cd $OUTPUT_FOLDER;
export PATH="$HOME2/crysalin-redux/repo/code:$PATH";
$CONDA_PREFIX/bin/python $HOME2/crysalin-redux/rescore_v2.py $OUTPUT_FOLDER $HOME2/crysalin-redux/pentakaihemimer_renumbered2.pdb

# END OF ACTUAL COMMAND
# -------------------------------------------------------
# FINAL FLUFF

#curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' complete"}' $SLACK_WEBHOOK
wget --quiet \
     --method=POST \
     --header="Content-Type: application/json" \
     --body-data='{"text":"'"$SLURM_JOB_NAME"' complete"}' \
     "$SLACK_WEBHOOK"