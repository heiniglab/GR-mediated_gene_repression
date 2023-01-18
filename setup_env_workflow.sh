#!/bin/bash

set -eo pipefail

eval "$(conda shell.bash hook)" 

# create env that starts the workflow and runs most of the shell commands
#-----------------------------------------------------------------------
#mamba env create --file src/conda/gr_env.yml

# create memes env which is used to later call the command line meme
#-----------------------------------------------------------------------
# conda version of memes is too old, so we install it through BiocManager instead and change LIB_PATH for that
#mamba env create --file src/conda/memes.yml

conda activate memes

# copy scripts into the new conda env that make sure the R lib paths are set correctly
cp src/conda/env_vars_activate.sh $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh 
cp src/conda/env_vars_deactivate.sh $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

# create site-library directory to avoid breaking r-base when installing packages through R (instead conda)
# see https://githubmemory.com/repo/conda-forge/r-base-feedstock/issues/169 for explanation
mkdir -p $CONDA_PREFIX/lib/R/site-library

# de- and reactivate the env for changes to R path to take effect
conda deactivate
conda activate memes

# let's see if we even need ChIPUtils in the tidied up workflow
R -e "BiocManager::install('memes')"

