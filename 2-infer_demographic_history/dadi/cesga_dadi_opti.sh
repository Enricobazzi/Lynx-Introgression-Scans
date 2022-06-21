#!/bin/bash

module load cesga/2020 && module load miniconda3/4.11.0

source activate dadi

# 1st arg = model
model=($(echo ${1}))

# 2nd arg = pop-pair
pop_pair=($(echo ${2}))

# 3rd arg = rep number
n=($(echo ${3}))

# define working folder
wfold=dadi_${model}_${pop_pair}_${n}

# create new working folder and go there
cp -r dadi_proto_folder ${wfold}
cd ${wfold}

# run one round of optimization
python dadi_run_models_optimize.py ${model} ${pop_pair}
