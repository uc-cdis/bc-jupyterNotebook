#!/bin/bash
export SUBJECTS_DIR=/home/jovyan/pd/freesurfer
sub=$1
out_path=$2

echo "labl_import_annotation \"aparc.annot\"" > tmp.tcl
echo "scale_brain 1.35" >>tmp.tcl
echo "redraw" >>tmp.tcl
echo "save_tiff "$out_path"/"$sub".lh.lat.hr.tif" >>tmp.tcl
echo "rotate_brain_y 180.0"  >>tmp.tcl
echo "redraw" >>tmp.tcl
echo "save_tiff "$out_path"/"$sub".lh.med.hr.tif" >>tmp.tcl
echo "exit 0" >>tmp.tcl
tksurfer $sub lh pial -tcl tmp.tcl

echo "labl_import_annotation \"aparc.annot\"" > tmp.tcl
echo "scale_brain 1.35" >>tmp.tcl
echo "redraw"  >>tmp.tcl
echo "save_tiff "$out_path"/"$sub".rh.lat.hr.tif" >>tmp.tcl
echo "rotate_brain_y 180.0" >>tmp.tcl
echo "redraw" >>tmp.tcl
echo "save_tiff "$out_path"/"$sub".rh.med.hr.tif" >>tmp.tcl
echo "exit 0" >>tmp.tcl
tksurfer $sub rh pial -tcl tmp.tcl

rm tmp.tcl

