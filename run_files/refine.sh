#!/bin/bash

dir=$(pwd)
loc=$(echo ${dir%/*})
DATABASE=/PATH/TO/ROSETTA/DATABASE
ROSETTA_LOC=/PATH/TO/ROSETTA
INPUT_FILES=${loc}/input_files

phenix.rosetta_refine ${INPUT_FILE}/lilly.mtz ${INPUT_FILE}/lilly.pdb refinement.input.xray_data.labels=IMEAN_1,SIGIMEAN_1 input.xray_data.r_free_flags.generate=True script=${INPUT_FILE}/spades_rosetta_refine.xml nproc=2 post_refine=True  database=${DATABASE} verbose=True spades_args=${INPUT_FILE}/spades_args.txt

cd rosetta_1
mkdir hydrate_structure
cd hydrate_structure

cp ../lilly_rosetta_phenix_001.mtz .
cp ../lilly_rosetta_phenix_001.pdb .

# run rehydration
${loc}/run_files/rescore.sh lilly_rosetta_phenix_001.pdb ${loc} ${ROSETTA_LOC}
${loc}/run_files/pdb_remove_far_away_waters.sh reu_lilly_rosetta_phenix_001_copy_0001.pdb lilly_refined_final.pdb

