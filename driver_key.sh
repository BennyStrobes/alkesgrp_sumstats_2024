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
gazal_lab_sumstats_dir="/n/groups/price/ldsc/sumstats_formatted_2024/gazal_lab_sumstats_231/"

# File describing Gazal lab 2021 summary statistics
gazal_lab_sumstats_summary_file="/n/groups/price/ldsc/sumstats_formatted_2024/gazal_lab_sumstats_231/gazal_lab_sumstats_231_file_description.csv"

# File containing Hapmap3 rsids
# Needed because we are going to limit to HapMap3 variants
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/list.txt"

# File containing MHC rsids
# Needed because we exclude MHC rsids
mhc_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/MHC.rsid"



####################
# Output data
####################
# Directory containing AlkesGrp 2024 formatted summary statistics
alkesgrp_2024_sumstats_dir="/n/groups/price/ldsc/sumstats_formatted_2024/"




source /home/bes710/.bash_profile


####################
# 1. Generate Alkesgrp 2024 formatted summary statistics
# This involves concatenating data from (a) alkesgrp 2021 formatted summary statistics and (b) Gazal lab summary statistics
####################
# File to be created describing AlkesGrp 2024 summary statistics (to be made in same format as $alkesgrp_2021_sumstats_summary_file)
alkesgrp_2024_sumstats_summary_file=${alkesgrp_2024_sumstats_dir}"summary_statistics_2024.txt"

python3 generate_alkesgrp_2024_formatted_sumstats.py $alkesgrp_2021_sumstats_dir $alkesgrp_2021_sumstats_summary_file $gazal_lab_sumstats_dir $gazal_lab_sumstats_summary_file $hapmap3_rsid_file $mhc_rsid_file $alkesgrp_2024_sumstats_dir $alkesgrp_2024_sumstats_summary_file

