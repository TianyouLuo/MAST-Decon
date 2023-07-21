#!/bin/bash

for slide in HCC-5A HCC-5B HCC-5C HCC-5D
do
for c in med
do
for ref in GuilliamsMa
do
name="${slide}_${ref}"
sbatch -t 128:00:00 --mem=40g -p general --job-name=${name} --wrap="Rscript /proj/yunligrp/users/tianyou/MASTDecon/HCCDB/codes/run_smoothing.R ${slide} ${c} ${ref}"
done
done
done

