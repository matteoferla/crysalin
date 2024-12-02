# job_arthor.sh
: << USAGE
export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh;
export APPTAINER_CONTAINER=/data/xchem-fragalysis/mferla/singularity/rocky3.sif;
export APPTAINERENV_CONDA_PREFIX=/data/xchem-fragalysis/mferla/waconda;
export JOB_INNER_SCRIPT=/data/xchem-fragalysis/mferla/crysalin/job_thread.sh;

condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-107.novalocal")';
USAGE

# nohup sleep 7000 && let JOB_ITERATION=30 && condor_submit /data/xchem-fragalysis/shared/target_script.condor;

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME"
# ---------------------------------------------------------------

source /data/xchem-fragalysis/mferla/.bashrc;
cd $HOME2/crysalin;
conda activate RFdiffusion;

echo '\n\nFILTER\n\n'
python $HOME2/crysalin/filter.py;

echo '\n\nGEN SEQS\n\n'
python $HOME2/crysalin/prep_MPNN.py
export PROTEINMPNN_WEIGHTS=$HOME2/.cache/ProteinMPNN_weights/vanilla_model_weights
python $CONDA_PREFIX/bin/protein_mpnn_run.py \
        --jsonl_path 'output_MPNN2/chains_definitions.jsonl' \
        --chain_id_jsonl 'output_MPNN2/fixed_chains.json' \
        --fixed_positions_jsonl 'output_MPNN2/fixed_positions.json' \
        --out_folder output_MPNN \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1 \
        --path_to_model_weights $PROTEINMPNN_WEIGHTS

echo '\n\nTHREADING & RELAX\n\n'
python $HOME2/crysalin/thread.py;

curl -X POST -H 'Content-type: application/json' --data '{"text":"Done threading"}' $SLACK_WEBHOOK