#!/bin/bash
export FREESURFER_HOME=/home/ubuntu/tools/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
#./mkIQIrange.sh > jnk.txt
cat jnk.txt | grep "has" |  awk -F/ ' { print $NF } ' | awk ' { print $1 } '| sort | uniq > jnk2.txt
cat jnk2.txt | wc -l
hd=`pwd`

for img in `cat jnk2.txt`
do echo $img 
cd /home/ubuntu/enigma/output/$img/mri/
mri_convert --out_orientation RAS --in_type mgz --out_type nii T1.mgz T1.nii ; 
mri_convert --out_orientation RAS --in_type mgz --out_type nii aseg.mgz aseg.nii ; 
more $hd/jnk.txt | grep $img
fslview T1.nii aseg.nii -t 0.2 -l "MGH-Subcortical"; 
rm *.nii
done
