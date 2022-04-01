#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
source ~/.bashrc;
source activate ATACseq_env;
while getopts o:i flag
do
    case "${flag}" in
        o) output=${OPTARG};;
        i) input=${OPTARG};;
    esac
done
awk \
"{
	if (NR>1)
		{print $2 "\t" "merged_peaks" "\t" "Peak" "\t" $3 "\t" $4 "\t" $6 "\t" $5 "\t" "." "\t" "peak_name=" $1}
}" /gpfs/home/snagaraja/ATACseqPipeline/data/merged/merged_peaks_HOMER.bed > /gpfs/home/snagaraja/ATACseqPipeline/data/merged/merged_peaks_HOMER.gff;