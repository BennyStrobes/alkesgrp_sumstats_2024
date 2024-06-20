#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-45:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)



alkesgrp_2024_sumstats_summary_file="$1"
alkesgrp_2024_sumstats_dir="$2"
ldsc_code_dir="$3"
ldsc_baseline_ld_hg19_annotation_dir="$4"
ref_1kg_hg19_genotype_dir="$5"
sldsc_h19_weights_dir="$6"
sldsc_h2_results_dir="$7"


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

	# Run S-LDSC in each trait independently
	trait_identifier="${myArray[1]}"
	echo $trait_identifier

	# Get gwas sumstat filename corresponding to trait identifier
	gwas_file_name=${alkesgrp_2024_sumstats_dir}${trait_identifier}".sumstats"

	# S-LDSC command
	python ${ldsc_code_dir}ldsc.py --h2 ${gwas_file_name} --ref-ld-chr ${ldsc_baseline_ld_hg19_annotation_dir} --w-ld-chr ${sldsc_h19_weights_dir} --overlap-annot --print-coefficients --frqfile-chr ${ref_1kg_hg19_genotype_dir} --out ${sldsc_h2_results_dir}${trait_identifier}"_sldsc_baselineLD_v2.2"


done < $alkesgrp_2024_sumstats_summary_file


