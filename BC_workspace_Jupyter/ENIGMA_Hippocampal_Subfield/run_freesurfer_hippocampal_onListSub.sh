#!/bin/bash

export FREESURFER_HOME=/mnt/freesurfer/
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/home/ubuntu/enigma/output/
sub_list=$1
# Get subject
#exec < /home/ubuntu/enigma/input/list_subjects.txt
exec < sub_list
while read x; do
# Run FreeSurfer
recon-all -s $x -hippocampal-subfields-T1 
done
