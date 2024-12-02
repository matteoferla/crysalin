#!/bin/bash
# run: export EXPERIMENT='...'; export CONTIGMAP='...'; sbatch /opt/xchem-fragalysis-2/mferla/library_making/third_pass.slurm.sh

#SBATCH --job-name=RFdiffusion
#SBATCH --chdir=/opt/xchem-fragalysis-2/mferla
#SBATCH --output=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
#SBATCH --error=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
#SBATCH --partition=gpu
##SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --priority=0
#SBATCH --export=NONE,EXPERIMENT,CONTIGMAP,NUMDESIGNS,HOTSPOTS,SLACK_WEBHOOK

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

conda activate RFdiffusion2;

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
echo $CONTIGMAP;
if [ -z "$CONTIGMAP" ]; then
  echo "Error: CONTIGMAP is not set."
  curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' missing CONTIGMAP var"}' $SLACK_WEBHOOK
  exit 1
fi
export NUMDESIGNS=${NUMDESIGNS:-1};

export RFDIFFUSION_CONFIG=$HOME2/.cache/RFdiffusion_config/inference;
export OUTPUT_FOLDER=$HOME2/crysalin-redux/output/$EXPERIMENT;
mkdir -p $OUTPUT_FOLDER;

run_inference.py \
--config-path=$RFDIFFUSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin-redux/output \
inference.output_prefix=$OUTPUT_FOLDER/$EXPERIMENT \
inference.input_pdb=$HOME2/crysalin-redux/pentakaihemimer_renumbered2.pdb \
inference.write_trajectory=false \
ppi.hotspot_res="$HOTSPOTS" \
contigmap.contigs="$CONTIGMAP" \
inference.num_designs=$NUMDESIGNS;


# END OF ACTUAL COMMAND
# -------------------------------------------------------
# FINAL FLUFF

#curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' complete"}' $SLACK_WEBHOOK
wget --quiet \
     --method=POST \
     --header="Content-Type: application/json" \
     --body-data='{"text":"'"$SLURM_JOB_NAME"' complete"}' \
     "$SLACK_WEBHOOK"