#!/bin/bash

if [[ -z $2 ]]
then
    echo "Usage: $0 input.pdb output.pdb"
    echo "  task: removes all lines except those that start with ATOM"
    exit
fi

touch "${2}"
rm "${2}"
cat "${1}" | grep -e '^ATOM' > "${2}"
echo "TER " >> "${2}"
