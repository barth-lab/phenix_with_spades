#!/bin/bash

dir=$(pwd)
loc=$(echo ${dir%/*})
input=${loc}/input_files
lilly_filename=${input}/lilly

sed "s|lilly|${lilly_filename}|g" ${input}/spades_args_tmp.txt  > ${input}/spades_args.txt
# N.B. These are general hydrate weights, however to get sed working cleanly we call them lilly.wts
sed "s|lilly|${lilly_filename}|g" ${input}/flags_tmp > ${input}/flags
