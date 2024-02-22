#
export NEW_CONDA_PREFIX=RFdiffusion2
export RFDIFFUSSION_CONFIG=$HOME2/.cache/RFdiffusion_config/inference
# use scratch tmp:
export REPO_ROOT=/data/outerhome/tmp/RFdiffusion

# common fluff
unset LD_LIBRARY_PATH
unset CUDA_HOME
unset CUDA_DIR
unset XLA_FLAGS
unset CONDA_CHANNELS

# install
conda create --name $NEW_CONDA_PREFIX -y python=3.9 ipykernel
conda activate $NEW_CONDA_PREFIX

# common fluff:
conda env config vars set LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX:/.singularity.d/libs
conda env config vars set PYTHONUSERBASE=$CONDA_PREFIX
conda deactivate
conda activate $NEW_CONDA_PREFIX

# install pyrosetta (for my stuff)
pip install -q https://levinthal:paradox@graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python39.ubuntu.wheel/pyrosetta-2024.4+release.a56f3bb973-cp39-cp39-linux_x86_64.whl
pip install -q rdkit #py3dmol

# install stuff
pip install -q --no-cache-dir \
                      dgl==1.0.2+cu116 -f https://data.dgl.ai/wheels/cu116/repo.html \
                      torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116 \
                      e3nn==0.3.3 \
                      wandb==0.12.0 \
                      pynvml==11.0.0 \
                      git+https://github.com/NVIDIA/dllogger#egg=dllogger \
                      decorator==5.1.0 \
                      hydra-core==1.3.2 \
                      pyrsistent==0.19.3

# install RFdiffusion
conda install -y anaconda::git
# /usr/libexec/git-core/git-remote-https: symbol lookup error: /lib64/libldap.so.2: undefined symbol: EVP_md2, version OPENSSL_3.0.0
conda install -y anaconda::openldap
git clone https://github.com/RosettaCommons/RFdiffusion.git $REPO_ROOT
cd $REPO_ROOT
pip install -q --no-cache-dir  $REPO_ROOT/env/SE3Transformer
pip install -q --no-cache-dir $REPO_ROOT --no-deps

# install conda
export CONDA_OVERRIDE_CUDA="11.6.2";
conda install -y nvidia/label/cuda-11.6.2::cuda-toolkit
conda install -y nvidia/label/cuda-11.6.2::cuda-tools
conda install -y nvidia/label/cuda-11.6.2::cuda-nvrtc
conda install -y nvidia/label/cuda-11.6.2::libcufile
conda install -y nvidia/label/cuda-11.6.2::libcusparse
conda install -y nvidia/label/cuda-11.6.2::cuda-cudart
conda install -y nvidia/label/cuda-11.6.2::cuda-cudart-dev

conda env config vars set RFDIFFUSSION_CONFIG=$RFDIFFUSSION_CONFIG $HOME2/.cache/RFdiffusion_config/inference
cp


????

!mv models $HOME2/.cache/RFDiffusion_models

!mkdir $HOME2/.cache/ProteinMPNN_weights
!mv ../ProteinMPNN/vanilla_model_weights $HOME2/.cache/ProteinMPNN_weights/vanilla_model_weights
!mv ../ProteinMPNN/soluble_model_weights $HOME2/.cache/ProteinMPNN_weights/soluble_model_weights
!mv ../ProteinMPNN/ca_model_weights $HOME2/.cache/ProteinMPNN_weights/ca_model_weights
!cp ../ProteinMPNN/*.py $CONDA_PREFIX/bin/
cp /scripts/* $CONDA_PREFIX/bin/
export PROTEINMPNN_WEIGHTS=$HOME2/.cache/ProteinMPNN_weights/vanilla_model_weights