#!/bin/bash

loc=${2}
input=${loc}/input_files
ROSETTA=${3}

pdb=${1}
pdbname=$(echo "${pdb}" | cut -d'.' -f1)
copy_pdb=${pdbname}_copy.pdb

cp ${pdb} ${copy_pdb}
${loc}/run_files/pdb_only_atom.sh ${copy_pdb} ${copy_pdb}2
mv ${copy_pdb}2 ${copy_pdb}
    
${ROSETTA}/main/source/bin/hydrate.python.linuxgccrelease \
        @${input}/flags \
	-database ${ROSETTA}/main/database \
        -s ${copy_pdb} \
        -hydrate:just_score \
        -overwrite \
        -out:prefix tmp_ 
mv tmp_*.pdb ${copy_pdb}
touch tmp.tmp.tmp
rm *tmp*
    
${loc}/run_files/pdb_only_atom.sh ${copy_pdb} ${copy_pdb}2
mv ${copy_pdb}2 ${copy_pdb} 

phenix.get_cc_mtz_pdb mtz_in=${pdbname}.mtz pdb_in=${pdb} > $(basename ${pdb} .pdb).get_cc_mtz_pdb.log
mv cc.log $(basename ${pdb} .pdb).get_cc_mtz_pdb.txt
phenix.model_vs_data ${pdbname}.mtz ${pdbname}_offset.pdb > $(basename ${pdb} .pdb).model_vs_data.txt
phenix.molprobity ${copy_pdb} outliers_only=False keep_hydrogens=True > $(basename ${pdb} .pdb).molprobity.txt

${ROSETTA}/main/source/bin/hydrate.python.linuxgccrelease \
        @${input}/flags \
        -s ${copy_pdb} \
	-database ${ROSETTA}/main/database \
        -overwrite \
        -out:prefix reu_ >> reu_$(basename ${pdb} .pdb).log 2>&1


