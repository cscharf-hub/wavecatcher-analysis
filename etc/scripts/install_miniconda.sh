#!/bin/bash

# Define the version of Miniconda to install
MINICONDA_VERSION="latest"

# Determine the platform (linux or macOS)
PLATFORM=$(uname)

# Determine the processor architecture
ARCH=$(uname -m)

if [[ $ARCH!="x86_64" ]] || [[ $ARCH!="x86" ]]; then
    printf "INFO: Your system architecture (%s) might not be fully supported\n" $(uname -m)
    printf "In case of problems please follow the instructions at https://github.com/conda-forge/miniforge" 
fi    

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$PLATFORM-$ARCH.sh
source $HOME/.bashrc
rm Miniforge3-$PLATFORM-$ARCH.sh

# previous intructions
# wget https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-$PLATFORM-$ARCH.sh -O miniconda.sh
# bash miniconda.sh -b -p $HOME/miniconda
# echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc
# source $HOME/.bashrc
# rm miniconda.sh

# Initialize the base environment
$HOME/miniconda/bin/conda init bash
$HOME/miniconda/bin/conda init zsh

# Please restart shell after this step