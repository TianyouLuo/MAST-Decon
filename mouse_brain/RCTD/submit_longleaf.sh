#!/bin/bash

for slice in ST8059049 ST8059051 ST8059052
do
for sub in 1000 1500 2500 3500 5000 7500 10000
do
name1="${slice}_sub${sub}_full"
name2="${slice}_sub${sub}_multi"
sbatch -t 4:00:00 --mem=10g -p general --job-name=${name1} --wrap="Rscript /pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/codes/RCTD.R ${slice} ${sub}"
sbatch -t 4:00:00 --mem=10g -p general --job-name=${name2} --wrap="Rscript /pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/RCTD/codes/RCTD_multi.R ${slice} ${sub}"
done
done

