## Combined dumpyard of RFdiffusion commands

This is for my personal reference. 

```bash

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/rfdiff \
inference.input_pdb=$HOME2/crysalin/dimer.pdb \
'contigmap.contigs=[1-150/A202-371/0 B13-135]' \
'ppi.hotspot_res=[A222,A223,A224,A227,A228,B35,B37,B55,B59,B60,B62,B63,B64,B78,B80,B81,B82,B83,B85,B87]' \
inference.num_designs=10 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/beta \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-150/A198-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=10 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_lownoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_lownoise done"}' $SLACK_WEBHOOK

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_midnoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_midnoise done"}' $SLACK_WEBHOOK


run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_highnoise \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=1.0 denoiser.noise_scale_frame=1.0
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_highnoise done"}' $SLACK_WEBHOOK


run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_mega \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-500/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_megadone"}' $SLACK_WEBHOOK


run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_mini \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-100/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=1.0 denoiser.noise_scale_frame=1.0
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_mini done"}' $SLACK_WEBHOOK


run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_Aless \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_Aless done"}' $SLACK_WEBHOOK


run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/gamma_partial \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[1-200/A196-367/0 B1-123/0 C1-123/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,B43,B68,B69,B70,B71,B73,B75,C23,C25,C47,C48,C50,C51,C52,C66,C53,C54,C55]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 diffuser.partial_T=20
curl -X POST -H 'Content-type: application/json' --data '{"text":"gamma_partial done"}' $SLACK_WEBHOOK

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/delta_full \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[150-400/A196-367/0 B3-121/0 C3-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=10 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &



CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/delta_full \
inference.input_pdb=$HOME2/crysalin/trimer_tweaked.pdb \
'contigmap.contigs=[150-400/A196-367/0 B3-121/0 C3-121/0]' \
'ppi.hotspot_res=[A222,A223,A224,A227,A228,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=10 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 diffuser.partial_T=20 &

# APO

CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/trimer-renumbered.pdb \
inference.input_pdb=$HOME2/crysalin/trimer_tweaked.pdb \
'contigmap.contigs=[150-400/0 B3-121/0 C3-121/0]' \
'ppi.hotspot_res=[C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=10 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 diffuser.partial_T=20 &

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/epsilon_combo \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[A1-35 6-20 A42-58 6-20 A62-138 6-20 A142-367/0 B1-121/0 C1-121/0]' \
'ppi.hotspot_res=[A197,A198,A199,A202,A203,C1,C2,C3,C23,C25,C47,C48,C50,C51,C52,C53,C54,C55,C66,C84,C86,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100]' \
inference.num_designs=100 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 &

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
inference.output_prefix=$HOME2/crysalin/output/zeta_no_fluff \
inference.input_pdb=$HOME2/crysalin/trimer-renumbered.pdb \
'contigmap.contigs=[150-400/A196-367/0 B3-121/0 C3-121/0]' \
inference.num_designs=1000

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/eta/eta_test2 \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[A5-40/6-20/A47-63/6-20/A67-143/6-20/A147-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1 &

# ----

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

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98]' \
inference.output_prefix=$HOME2/crysalin/output/theta/theta_hot_ \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/theta/theta_cold_ \
inference.input_pdb=$HOME2/crysalin/trikaihemimer.pdb \
'contigmap.contigs=[150-300/A202-331/0 B13-135/0 C13-135/0 a202-331/0]' \
inference.num_designs=1000;

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

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/kappa/kappa \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[50-150/A150-157/70-150/A202-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/lambda/lambda_fore_hotcorrect \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
'contigmap.contigs=[150-200/A142-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 &


CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/lambda/lambda_fore \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[150-200/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/lambda/lambda_comboplus \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[A5-40/6-20/A47-63/6-20/A67-143/6-25/A149-157/9-20/A166-189/9-20/A198-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=3 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/lambda/lambda_combo \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[A5-40/6-20/A47-63/6-20/A67-143/6-20/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000;

CUDA_VISIBLE_DEVICES=1 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/mu/mu_partial \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[10-75/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 & 

CUDA_VISIBLE_DEVICES=2 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/mu/mu_micro \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[10-20/A150-157/10-75/A202-331/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 & 

CUDA_VISIBLE_DEVICES=3 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/mu/mu_full1 \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[10-75/A150-157/10-75/A202-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000 &

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
inference.output_prefix=$HOME2/crysalin/output/mu/mu_full2 \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'contigmap.contigs=[10-75/A150-157/10-75/A202-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000;

CUDA_VISIBLE_DEVICES=0 run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/nu/nu \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'ppi.hotspot_res=[A203,A204,A207,a271,a205]' \
'contigmap.contigs=[A5-24/6-12/A32-80/6-12/A89-107/6-12/A114-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/omicron/omicron \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
'contigmap.contigs=[A5-9/6-20/A24-40/6-20/A47-63/6-20/A67-143/6-20/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000;

run_inference.py \
--config-path=$RFDIFFUSSION_CONFIG \
hydra.output_subdir=$HOME2/crysalin/output \
inference.output_prefix=$HOME2/crysalin/output/omicron/omicron2 \
inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb \
'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B80,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' \
'contigmap.contigs=[A5-9/10-20/A24-40/6-20/A47-63/6-20/A67-143/6-20/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' \
inference.num_designs=1000;

CUDA_VISIBLE_DEVICES=0 run_inference.py --config-path=$RFDIFFUSSION_CONFIG hydra.output_subdir=$HOME2/crysalin/output inference.output_prefix=$HOME2/crysalin/output/omicron/omicron3 inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb 'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B30,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' 'contigmap.contigs=[A5-9/10-20/A24-40/6-20/A47-63/6-20/A67-143/3-5/A147-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' inference.num_designs=1000

#??? Where is pi?

CUDA_VISIBLE_DEVICES=1 run_inference.py --config-path=$RFDIFFUSSION_CONFIG hydra.output_subdir=$HOME2/crysalin/output inference.output_prefix=$HOME2/crysalin/output/rho/rho inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb 'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B30,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' 'contigmap.contigs=[A5-9/10-20/A24-40/6-20/A47-51/17-30/A69-139/7-9/A149-158/8-12/A167-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' inference.num_designs=1000

CUDA_VISIBLE_DEVICES=1 run_inference.py --config-path=$RFDIFFUSSION_CONFIG hydra.output_subdir=$HOME2/crysalin/output inference.output_prefix=$HOME2/crysalin/output/rho/rho inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb 'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B30,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' 'contigmap.contigs=[A5-9/10-20/A24-40/6-20/A47-51/17-30/A69-139/7-9/A149-158/8-12/A167-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' inference.num_designs=1000

CUDA_VISIBLE_DEVICES=2 run_inference.py --config-path=$RFDIFFUSSION_CONFIG hydra.output_subdir=$HOME2/crysalin/output inference.output_prefix=$HOME2/crysalin/output/sigma/sigma inference.input_pdb=$HOME2/crysalin/pentakaihemimer.relax.pdb 'ppi.hotspot_res=[A203,A204,A207,a271,a205,B25,B27,B43,B45,B46,B47,B49,B52,B54,B55,B79,B30,B81,B82,B83,B84,B85,B86,B87,B88,B90,B92,B108,B110,B112,C13,C14,C15,C35,C37,C59,C60,C62,C63,C64,C65,C66,C67,C78,C96,C98,L24,L25]' 'contigmap.contigs=[A5-9/10-20/A24-40/6-20/A47-51/20-40/A69-139/7-9/A149-158/8-12/A167-331/0 B13-135/0 C13-135/0 K13-135/0 L13-135/0 F304-332/0 a202-331/0]' inference.num_designs=1000

```