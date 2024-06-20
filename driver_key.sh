#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)


#########################################################
# Ben Strober
#########################################################
# Scripts to generate Alkesgrp 2024 formatted summary statistics
#########################################################







####################
# Input data
####################
# Directory containing AlkesGrp 2021 summary statistics
alkesgrp_2021_sumstats_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"

# File describing AlkesGrp 2021 summary statistics
# This file was made by saving the first page of /n/groups/price/ldsc/sumstats_formatted_2021/Description_070121.xlsx as a csv
alkesgrp_2021_sumstats_summary_file="/n/groups/price/ldsc/sumstats_formatted_2024/reference_files/summary_statistics_2021_summary.csv"

# Directory containing Gazal-lab 2024 set of 231 summary statistics 
gazal_lab_sumstats_dir="/n/groups/price/ldsc/sumstats_formatted_2024/reference_files/gazal_lab_sumstats_231/"

# File describing Gazal lab 2021 summary statistics
gazal_lab_sumstats_summary_file="/n/groups/price/ldsc/sumstats_formatted_2024/reference_files/gazal_lab_sumstats_231/gazal_lab_sumstats_231_file_description.csv"

# File containing Hapmap3 rsids
# Needed because we are going to limit to HapMap3 variants
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/1000G_hm3_noMHC.snplist"

# Ldscore regression code
ldsc_code_dir="/n/groups/price/ben/tools/ldsc/"

# LDSC baseline LD Dir
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD."
ldsc_eas_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EAS_Phase3/baselineLD_v2.2/baselineLD."

# ldscore direcotry (used for genetic correlation analysis)
ldscore_hg19_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/LDscore/LDscore."
ldscore_eas_hg19_dir="/n/groups/price/ldsc/reference_files/1000G_EAS_Phase3/LDscore/LDscore."

# LDSC 1KG genotype files (hg19)
ref_1kg_hg19_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC."
ref_1kg_eas_hg19_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EAS_Phase3/plink_files/1000G.EAS.QC."


# hg19 sldsc weights
sldsc_h19_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC."
sldsc_eas_h19_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EAS_Phase3/weights/weights.EAS.hm3_noMHC."



####################
# Output data
####################
# Output root directory containing AlkesGrp 2024 formatted summary statistics
alkesgrp_2024_sumstats_root_dir="/n/groups/price/ldsc/sumstats_formatted_2024/"

# Output directory containing summary stat files
alkesgrp_2024_sumstats_dir=${alkesgrp_2024_sumstats_root_dir}"sumstats/"

# Output directory containing S-LDSC h2 estimates on 2024 summary statistics
sldsc_h2_results_dir=${alkesgrp_2024_sumstats_root_dir}"sldsc_h2/"

# Output directory containing LDSC genetic correlation (rg) estimates on 2024 summary statistics
ldsc_rg_results_dir=${alkesgrp_2024_sumstats_root_dir}"ldsc_rg/"

# Output directory containing QC analysis comparing h2 parameter to 2021 updates amongst overlapping traits
qc_comparison_to_2021_update_dir=${alkesgrp_2024_sumstats_root_dir}"qc_comparison_to_2021_update/"



if false; then
source /home/bes710/.bash_profile
fi

####################
# 1. Generate Alkesgrp 2024 formatted summary statistics
# This involves concatenating data from (a) alkesgrp 2021 formatted summary statistics and (b) Gazal lab summary statistics
####################
# File to be created describing AlkesGrp 2024 summary statistics (to be made in same format as $alkesgrp_2021_sumstats_summary_file)
alkesgrp_2024_sumstats_summary_file=${alkesgrp_2024_sumstats_root_dir}"summary_statistics_2024.txt"
alkesgrp_2024_x_trait_ldsc_sumstats_summary_file=${alkesgrp_2024_sumstats_root_dir}"signed_summary_statistics_2024.txt"
if false; then
python3 generate_alkesgrp_2024_formatted_sumstats.py $alkesgrp_2021_sumstats_dir $alkesgrp_2021_sumstats_summary_file $gazal_lab_sumstats_dir $gazal_lab_sumstats_summary_file $hapmap3_rsid_file $ldsc_code_dir $alkesgrp_2024_sumstats_dir $alkesgrp_2024_sumstats_summary_file $alkesgrp_2024_x_trait_ldsc_sumstats_summary_file
fi

# Split summary statistics into EUR and EAS groups
alkesgrp_2024_EUR_sumstats_summary_file=${alkesgrp_2024_sumstats_root_dir}"reference_files/summary_statistics_EUR_2024.txt"
alkesgrp_2024_EAS_sumstats_summary_file=${alkesgrp_2024_sumstats_root_dir}"reference_files/summary_statistics_EAS_2024.txt"
if false; then
python3 split_gwas_studies_into_eur_and_eas.py $alkesgrp_2024_sumstats_summary_file $alkesgrp_2024_EUR_sumstats_summary_file $alkesgrp_2024_EAS_sumstats_summary_file
fi

# Extract all pairs of traits within EUR and EAS groups
# EUR
alkesgrp_2024_EUR_trait_pair_file=${alkesgrp_2024_sumstats_root_dir}"reference_files/summary_statistics_trait_pairs_EUR_2024.txt"
if false; then
python3 extract_all_pairs_of_traits.py $alkesgrp_2024_EUR_sumstats_summary_file $alkesgrp_2024_EUR_trait_pair_file
fi
# EAS
alkesgrp_2024_EAS_trait_pair_file=${alkesgrp_2024_sumstats_root_dir}"reference_files/summary_statistics_trait_pairs_EAS_2024.txt"
if false; then
python3 extract_all_pairs_of_traits.py $alkesgrp_2024_EAS_sumstats_summary_file $alkesgrp_2024_EAS_trait_pair_file
fi



####################
# 2. Run s-LDSC on each trait to compute heritabilities
####################
if false; then
# EUR
sbatch run_sldsc_on_all_2024_summary_statistics.sh $alkesgrp_2024_EUR_sumstats_summary_file $alkesgrp_2024_sumstats_dir $ldsc_code_dir $ldsc_baseline_ld_hg19_annotation_dir $ref_1kg_hg19_genotype_dir $sldsc_h19_weights_dir $sldsc_h2_results_dir
# EAS
sbatch run_sldsc_on_all_2024_summary_statistics.sh $alkesgrp_2024_EAS_sumstats_summary_file $alkesgrp_2024_sumstats_dir $ldsc_code_dir $ldsc_eas_baseline_ld_hg19_annotation_dir $ref_1kg_eas_hg19_genotype_dir $sldsc_eas_h19_weights_dir $sldsc_h2_results_dir
fi


####################
# 3. Run LDSC-genetic correlation on all pairs of traits
####################
if false; then
# EAS
sbatch run_ldsc_genetic_correlation_on_2024_summary_statistics.sh $alkesgrp_2024_EAS_trait_pair_file $alkesgrp_2024_sumstats_dir $ldsc_code_dir $ldscore_eas_hg19_dir $sldsc_eas_h19_weights_dir $ldsc_rg_results_dir
# EUR
num_jobs="30"
for job_number in $(seq 0 29); do 
	sbatch run_ldsc_genetic_correlation_on_2024_summary_statistics_parrallelized.sh $alkesgrp_2024_EUR_trait_pair_file $alkesgrp_2024_sumstats_dir $ldsc_code_dir $ldscore_hg19_dir $sldsc_h19_weights_dir $ldsc_rg_results_dir $job_number $num_jobs 
done
fi

####################
# 4. Organize S-LDSC heritability results
####################
# EUR
alkesgrp_2024_EUR_h2_summary_file=${alkesgrp_2024_sumstats_root_dir}"sldsc_heritabilities_EUR_2024.txt"
if false; then
python3 organize_sldsc_h2_results.py $alkesgrp_2024_EUR_sumstats_summary_file $sldsc_h2_results_dir $alkesgrp_2024_EUR_h2_summary_file
fi
# EAS
alkesgrp_2024_EAS_h2_summary_file=${alkesgrp_2024_sumstats_root_dir}"sldsc_heritabilities_EAS_2024.txt"
if false; then
python3 organize_sldsc_h2_results.py $alkesgrp_2024_EAS_sumstats_summary_file $sldsc_h2_results_dir $alkesgrp_2024_EAS_h2_summary_file
fi

####################
# 5. Organize genetic correlation results
####################
# EUR
alkesgrp_2024_EUR_genetic_corr_summary_file=${alkesgrp_2024_sumstats_root_dir}"ldsc_genetic_correlations_trait_pairs_EUR_2024.txt"
if false; then
python3 organize_ldsc_genetic_correlation_results.py $alkesgrp_2024_EUR_trait_pair_file $ldsc_rg_results_dir $alkesgrp_2024_EUR_genetic_corr_summary_file
fi
# EAS
alkesgrp_2024_EAS_genetic_corr_summary_file=${alkesgrp_2024_sumstats_root_dir}"ldsc_genetic_correlations_trait_pairs_EAS_2024.txt"
if false; then
python3 organize_ldsc_genetic_correlation_results.py $alkesgrp_2024_EAS_trait_pair_file $ldsc_rg_results_dir $alkesgrp_2024_EAS_genetic_corr_summary_file
fi


####################
# 6. QC analysis comparing h2 parameter to 2021 updates amongst overlapping traits
####################
if false; then
sh qc_comparison_to_2021_sumstat_update.sh $alkesgrp_2024_EUR_h2_summary_file $alkesgrp_2024_EAS_h2_summary_file $alkesgrp_2024_EUR_genetic_corr_summary_file $alkesgrp_2024_EAS_genetic_corr_summary_file $alkesgrp_2024_sumstats_summary_file $qc_comparison_to_2021_update_dir
fi



####################
# 7. Generate non-redundent trait list
####################
# EUR
alkesgrp_2024_EUR_non_redundent_trait_file=${alkesgrp_2024_sumstats_root_dir}"non_redundent_traits_EUR_2024.txt"
if false; then
python3 generate_non_redundent_trait_list.py $alkesgrp_2024_sumstats_summary_file $alkesgrp_2024_EUR_h2_summary_file $alkesgrp_2024_EUR_genetic_corr_summary_file $alkesgrp_2024_EUR_non_redundent_trait_file
fi
# EAS
alkesgrp_2024_EAS_non_redundent_trait_file=${alkesgrp_2024_sumstats_root_dir}"non_redundent_traits_EAS_2024.txt"
if false; then
python3 generate_non_redundent_trait_list.py $alkesgrp_2024_sumstats_summary_file $alkesgrp_2024_EAS_h2_summary_file $alkesgrp_2024_EAS_genetic_corr_summary_file $alkesgrp_2024_EAS_non_redundent_trait_file
fi

