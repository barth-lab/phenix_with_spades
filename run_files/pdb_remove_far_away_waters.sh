#!/bin/bash

if [[ -z $2 ]]
then
    echo "Usage: $0 input.pdb output.pdb"
    echo "  task: removes hetatm that are too far away"
    exit
fi

touch "${2}"
rm "${2}"

farAwayWaters=0
while read -r line
do
    if [[ "${line:0:6}" == 'HETATM' && $(echo "${line:29:9} > 100000.0" | bc) == 1 ]]
    then
        farAwayWaters=$(echo "${farAwayWaters} + 1" | bc)
    else
        echo "${line}" >> ${2}
    fi
done < ${1}

if [[ $farAwayWaters > 1 ]]
then
    echo "Found ${farAwayWaters} atoms that are far away" 
fi

