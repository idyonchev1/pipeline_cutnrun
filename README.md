Snakemake pipeline for the normalization of CUT&RUN data, followed by track generation and MACS2 peak-calling.

Normalization is performed against regions of high background using the hg38 CUT&RUN greenlist developed by Fabio N de Mello, et al. 2024 https://doi.org/10.1093/bib/bbad538

Requires a comma separated samples.txt file with columns Sample,condition,line,replicate, where at least one sample for each line+replicate is an IgG sample, and a conda environment containing bedtools deeptools base-r tidyr and dplyr. Each sample is normalized by dividing by the sum of background reads. I find samples normalized using this method to give more convincing sets of peaks during MACS2-peakcalling with disabled local lambda.   
