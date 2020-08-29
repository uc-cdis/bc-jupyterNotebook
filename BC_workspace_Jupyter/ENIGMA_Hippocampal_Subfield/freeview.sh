#!/bin/bash
export FREESURFER_HOME='/mnt/freesurfer/'
source $FREESURFER_HOME/SetUpFreeSurfer.sh
sub=$1
path="/home/ubuntu/enigma/output/$sub/mri"
cd $path
#echo $path
#freeview -v /Users/xingyankuang/Downloads/output/subj10159/tmp/hippoSF_T1_v10_right/hippoAmygBinaryMask_autoCropped.mgz imageDump.mgz
#freeview -v /Users/xingyankuang/Downloads/output/subj10159/tmp/hippoSF_T1_v10_right/hippoAmygBinaryMask_autoCropped.mgz
freeview -v nu.mgz -v lh.hippoSfLabels-T1.v10.mgz:colormap=lut -v rh.hippoSfLabels-T1.v10.mgz:colormap=lut
