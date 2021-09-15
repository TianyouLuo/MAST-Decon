#!/bin/bash

for slice in ST8059048
do
for sub in 1000 1500 2500
do
for c in low med high
do
name="${slice}_sub${sub}_${c}C"
sbatch -t 8:00:00 --mem=4g -p general --job-name=${name} --wrap="Rscript /pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/MAST-Decon/run_smoothing_max4.R ${slice} ${sub} ${c}"
done
done
done

