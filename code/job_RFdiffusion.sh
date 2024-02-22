: << USAGE
export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh;
export APPTAINER_CONTAINER=/data/xchem-fragalysis/mferla/singularity/rocky3.sif;
export APPTAINERENV_CONDA_PREFIX=/data/xchem-fragalysis/mferla/waconda;
export JOB_INNER_SCRIPT=/data/xchem-fragalysis/mferla/crysalin/job_rf.sh;
export APPTAINERENV_DGLBACKEND=pytorch;

# export APPTAINERENV_EXPERIMENT=partial
# export APPTAINERENV_EXPERIMENT=Aless
# export APPTAINERENV_EXPERIMENT=mini
# export APPTAINERENV_EXPERIMENT=mega
# export APPTAINERENV_EXPERIMENT=noise
#export APPTAINERENV_EXPERIMENT=epsilon
#export APPTAINERENV_EXPERIMENT=theta
export APPTAINERENV_EXPERIMENT=eta
export APPTAINERENV_EXPERIMENT=iota
export APPTAINERENV_EXPERIMENT=testpolish
export APPTAINERENV_EXPERIMENT=proteinpmnn
condor_submit /data/xchem-fragalysis/shared/target_script.condor -a 'Requirements=(machine == "orpheus-worker-gpu-18.novalocal")';
USAGE

# nohup sleep 7000 && let JOB_ITERATION=30 && condor_submit /data/xchem-fragalysis/shared/target_script.condor;

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME"
# ---------------------------------------------------------------

source /data/xchem-fragalysis/mferla/.bashrc;
conda activate RFdiffusion;
export HOME2=/data/xchem-fragalysis/mferla;
export RFDIFFUSSION_CONFIG=$HOME2/.cache/RFdiffusion_config/inference
cd $HOME2/crysalin;
pwd;

if [ "$EXPERIMENT" == "noise" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_midnoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &

CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_highnoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=1.0 denoiser.noise_scale_frame=1.0;

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_midhighnoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.8 denoiser.noise_scale_frame=0.8 &

CUDA_VISIBLE_DEVICES=3 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_midlownoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.4 denoiser.noise_scale_frame=0.4;

curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_highnoise done"}' $SLACK_WEBHOOK

elif [ "$EXPERIMENT" == "mega" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_mega \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-500/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_megadone"}' $SLACK_WEBHOOK

elif [ "$EXPERIMENT" == "mini" ]; then

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_mini \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-100/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=1.0 denoiser.noise_scale_frame=1.0
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_mini done"}' $SLACK_WEBHOOK

elif [ "$EXPERIMENT" == "Aless" ]; then

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_Aless \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_Aless done"}' $SLACK_WEBHOOK

elif [ "$EXPERIMENT" == "partial" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_partial1 \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 diffuser.partial_T=20 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_partial2 \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 diffuser.partial_T=50;

curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_partial done"}' $SLACK_WEBHOOK

elif [ "$EXPERIMENT" == "epsilon" ]; then
# -----------------
CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/epsilon_ins36_41 \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[A1-35 6-20 A42-367/0 B1-121/0 C1-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/epsilon_ins59_61 \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[A1-58 6-20 A62-367/0 B1-121/0 C1-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &

CUDA_VISIBLE_DEVICES=3 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/epsilon_ins139_141 \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[A1-138 6-15 A142-367/0 B1-121/0 C1-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/epsilon_combo \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[A1-35 6-20 A42-58 6-20 A62-138 6-20 A142-367/0 B1-121/0 C1-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5;

curl -X POST -H 'Content-type: application/json' --data '{"text":"epsilon done"}' $SLACK_WEBHOOK

# ------------------
elif [ "$EXPERIMENT" == "eta" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/eta/eta_fore_hot \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-200/A142-367/0 B13-135/0 C13-135/0 a202-331/0]' \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98]' \
inference.num_designs=1000 &


CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/eta/eta_fore \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-200/A147-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/eta/eta_comboplus \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[A5-40/6-20/A47-63/6-20/A67-143/6-25/A149-157/9-20/A166-189/9-20/A198-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=3 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/eta/eta_combo \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[A5-40/6-20/A47-63/6-20/A67-143/6-20/A147-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000;

curl -X POST -H 'Content-type: application/json' --data '{"text":"eta done"}' $SLACK_WEBHOOK


elif [ "$EXPERIMENT" == "theta" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98]' \
inference.output_prefix=$HOME2/crysalin/output/theta/theta_hot \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/theta/theta_cold \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000;

curl -X POST -H 'Content-type: application/json' --data '{"text":"theta done"}' $SLACK_WEBHOOK

# ------
elif [ "$EXPERIMENT" == "iota" ]; then

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/iota/iota_hot \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/iota/iota_cold \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000;

elif [ "$EXPERIMENT" == "testpolish" ]; then

export PROTEINMPNN_WEIGHTS=$HOME2/.cache/ProteinMPNN_weights/vanilla_model_weights
python $CONDA_PREFIX/bin/protein_mpnn_run.py \
        --jsonl_path 'tests/chains_definitions.jsonl' \
        --chain_id_jsonl 'tests/fixed_chains.json' \
        --fixed_positions_jsonl 'tests/fixed_positions.json' \
        --out_folder tests \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1 \
        --path_to_model_weights $PROTEINMPNN_WEIGHTS

elif [ "$EXPERIMENT" == "proteinpmnn" ]; then

export PROTEINMPNN_WEIGHTS=$HOME2/.cache/ProteinMPNN_weights/vanilla_model_weights
python $CONDA_PREFIX/bin/protein_mpnn_run.py \
        --jsonl_path 'output_MPNN/chains_definitions.jsonl' \
        --chain_id_jsonl 'output_MPNN/fixed_chains.json' \
        --fixed_positions_jsonl 'output_MPNN/fixed_positions.json' \
        --out_folder output_MPNN \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1 \
        --path_to_model_weights $PROTEINMPNN_WEIGHTS

else
    echo "EXPERIMENT value is not recognized"
    # Add default commands or error handling here

fi
echo 'COMPLETE'