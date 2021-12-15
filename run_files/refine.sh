#!/bin/bash

dir=$(pwd)
loc=$(echo ${dir%/*})
DATABASE=rosetta_database

#phenix.rosetta_refine ../input_files/lilly.mtz ../input_files/lilly.pdb refinement.input.xray_data.labels=IMEAN_1,SIGIMEAN_1 input.xray_data.r_free_flags.generate=True script=../input_files/spades_rosetta_refine.xml nproc=1 post_refine=True  database=${DATABASE} verbose=True 

cd rosetta_1
mkdir hydrate_structure
cd hydrate_structure

cp ../lilly_rosetta_phenix_001.mtz .
cp ../lilly_rosetta_phenix_001.pdb .

# run rehydration
${loc}/run_files/rescore.sh lilly_rosetta_phenix_001.pdb ${loc}

mv reu_lilly_rosetta_phenix_001_copy_0001.pdb lilly_refined_final.pdb
