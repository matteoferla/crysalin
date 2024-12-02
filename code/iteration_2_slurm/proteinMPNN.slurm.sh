#!/bin/bash
# run: proteinMPNN.slurm.sh

#SBATCH --job-name=MPNN
#SBATCH --chdir=/opt/xchem-fragalysis-2/mferla
#SBATCH --output=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
#SBATCH --error=/opt/xchem-fragalysis-2/mferla/crysalin-redux/logs/slurm-error_%x_%j.log
##SBATCH --partition=gpu
#SBATCH --partition=main
##SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --cpus-per-task=6
##SBATCH --gres=gpu:1
#SBATCH --priority=0
#SBATCH --export=WORKPATH
##SBATCH --nodelist=host-192-168-222-202,host-192-168-222-211,host-192-168-222-213,host-192-168-222-214,host-192-168-222-215,
##SBATCH --mem=74962M

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

nvidia-smi;

printenv;

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
#export WORKPATH=${WORKPATH:-$HOME2/crysalin-redux/output}
export PROTEINMPNN_WEIGHTS=$HOME2/.cache/ProteinMPNN_weights/soluble_model_weights
cd $HOME2/crysalin-redux;
echo $WORKPATH
export PATH="$HOME2/crysalin-redux/repo/code:$PATH";
export SERIALNUMBER=$RANDOM;
#[ ! -f "$WORKPATH/chains_definitions.jsonl" ] && python $HOME2/crysalin-redux/repo/code/prep_MPNN_v3.py "$WORKPATH"
python $HOME2/crysalin-redux/repo/code/prep_MPNN_v3.py "$WORKPATH" "$SERIALNUMBER" 50;
python $CONDA_PREFIX/bin/protein_mpnn_run.py \
        --jsonl_path "$WORKPATH/chains_definitions$SERIALNUMBER.jsonl" \
        --chain_id_jsonl "$WORKPATH/fixed_chains$SERIALNUMBER.json" \
        --fixed_positions_jsonl "$WORKPATH/fixed_positions$SERIALNUMBER.json" \
        --out_folder $WORKPATH \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --seed 888 \
        --batch_size 1 \
        --path_to_model_weights $PROTEINMPNN_WEIGHTS

# END OF ACTUAL COMMAND
# -------------------------------------------------------
# FINAL FLUFF

#curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' complete"}' $SLACK_WEBHOOK
wget --quiet \
     --method=POST \
     --header="Content-Type: application/json" \
     --body-data='{"text":"'"$SLURM_JOB_NAME"' complete"}' \
     "$SLACK_WEBHOOK"
