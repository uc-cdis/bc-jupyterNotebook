#!/bin/bash

export FREESURFER_HOME=/opt/freesurfer/
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/mnt/CVB/ds000030_R1.0.4

# Get subject
sub=$1
input=$2

# Run FreeSurfer
recon-all -i $input -s $sub -sd freesurfer/ -all -openmp 8
