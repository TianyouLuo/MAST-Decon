# ST-deconvolution
This is the project for Spatial Transcriptomics data deconvolution.


--------------------------------------------------
# Documentation of results
For smoothing, the easiest way to keep only the top celltypes above certain threshold is just to sort the result and only keep top ones and set the others to 0, and then also eliminate the small ones below threshold after normalization and renormalize them. We tried to do this step in every iteration and set max_CT=4 and threshold=0.02. These results are in the '/old_max_threshold/' folders, but they are actually much worse than expected.

