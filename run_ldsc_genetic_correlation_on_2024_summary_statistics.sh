#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



alkesgrp_2024_trait_pair_file="$1"
alkesgrp_2024_sumstats_dir="$2"
ldsc_code_dir="$3"
ldscore_hg19_dir="$4"
sldsc_h19_weights_dir="$5"
ldsc_rg_results_dir="$6"

source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12



# Read and discard the first line
line_count=0

# Loop through traits
while IFS=$'\t' read -r -a myArray
do

    # Increment the line count
    ((line_count++))
    # Check if it's the first line
    if [ $line_count -eq 1 ]; then
        continue  # Skip processing the first line
    fi

    # Get names of pair of traits
	trait1_identifier="${myArray[0]}"
	trait2_identifier="${myArray[1]}"


	# Get gwas sumstat filename corresponding to trait identifiers
	gwas_trait1_file_name=${alkesgrp_2024_sumstats_dir}${trait1_identifier}".sumstats"
	gwas_trait2_file_name=${alkesgrp_2024_sumstats_dir}${trait2_identifier}".sumstats"

	# LDSC genetic correlation (rg) command
	python ${ldsc_code_dir}ldsc.py --rg ${gwas_trait1_file_name}","${gwas_trait2_file_name} --ref-ld-chr ${ldscore_hg19_dir} --w-ld-chr ${sldsc_h19_weights_dir} --out ${ldsc_rg_results_dir}${trait1_identifier}":"${trait2_identifier}"_ldsc_rg"


done < $alkesgrp_2024_trait_pair_file
