#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-38:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



alkesgrp_2024_trait_pair_file="$1"
alkesgrp_2024_sumstats_dir="$2"
ldsc_code_dir="$3"
ldscore_hg19_dir="$4"
sldsc_h19_weights_dir="$5"
ldsc_rg_results_dir="$6"
job_number="$7"
num_jobs="$8"

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12


############################################
# Do some stuff to calculate starting and ending line number for parallelization purposes
############################################
# Count up number of lines in trait pair file
num_lines=-1  # Start at -1 to avoid header
# Loop through traits
while IFS=$'\t' read -r -a myArray
do
    ((num_lines++))
done < $alkesgrp_2024_trait_pair_file
lines_per_job=$(( num_lines / num_jobs +1))
start_line=$(( lines_per_job * job_number +1 ))
end_line=$(( lines_per_job * (job_number+1) ))


echo $job_number
echo $num_jobs

############################################
# Run LDSC genetic correlation analysis for all pairs of traits in this parrallel run
############################################
# Read and discard the first line
line_count=-1

# Loop through traits
while IFS=$'\t' read -r -a myArray
do

    # Increment the line count
    ((line_count++))
    # Check if it's the first line
    if [ $line_count -eq 0 ]; then
        continue  # Skip processing the first line
    fi
    # Skip line if not in this parallel run
    if (( line_count < start_line )); then
    	continue
    fi
    if (( line_count > end_line )); then
    	continue
    fi

    # Get names of pair of traits
	trait1_identifier="${myArray[0]}"
	trait2_identifier="${myArray[1]}"

	echo $trait1_identifier":"$trait2_identifier

	# Get gwas sumstat filename corresponding to trait identifiers
	gwas_trait1_file_name=${alkesgrp_2024_sumstats_dir}${trait1_identifier}".sumstats"
	gwas_trait2_file_name=${alkesgrp_2024_sumstats_dir}${trait2_identifier}".sumstats"

	# LDSC genetic correlation (rg) command
	python ${ldsc_code_dir}ldsc.py --rg ${gwas_trait1_file_name}","${gwas_trait2_file_name} --ref-ld-chr ${ldscore_hg19_dir} --w-ld-chr ${sldsc_h19_weights_dir} --out ${ldsc_rg_results_dir}${trait1_identifier}":"${trait2_identifier}"_ldsc_rg"

done < $alkesgrp_2024_trait_pair_file

