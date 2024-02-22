: << USAGE
export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh;
export APPTAINER_CONTAINER=/data/xchem-fragalysis/mferla/singularity/rocky3.sif;
export APPTAINERENV_CONDA_PREFIX=/data/xchem-fragalysis/mferla/waconda;
export JOB_INNER_SCRIPT=/data/xchem-fragalysis/mferla/crysalin/job_remodel.sh;

export APPTAINERENV_EXPERIMENT=ins139_141
condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-38.novalocal")';
export APPTAINERENV_EXPERIMENT=ins59_61
condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-39.novalocal")';
export APPTAINERENV_EXPERIMENT=ins36_41
condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-40.novalocal")';
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
conda activate compchem;
python $HOME2/crysalin/output/remodel_job.py;

curl -X POST -H 'Content-type: application/json' --data '{"text":"remodel done'$EXPERIMENT'"}' $SLACK_WEBHOOK