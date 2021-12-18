#!/bin/bash

dir=$(pwd)
loc=$(echo ${dir%/*})
input=${loc}/input_files
lilly_filename=${input}/lilly
database_loc=${loc}/rosetta/main/database/

sed "s|lilly|${lilly_filename}|g" ${input}/spades_args_tmp.txt  > ${input}/spades_args.txt
sed "s|rosetta_database|${database_loc}|g" -i refine.sh
# N.B. These are genearl hydrate weights, however to get sed working cleanly we call them lilly.wts
sed "s|lilly|${lilly_filename}|g" ${input}/flags_tmp > ${input}/flags
