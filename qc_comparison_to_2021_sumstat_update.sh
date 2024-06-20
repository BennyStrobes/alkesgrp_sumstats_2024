#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)






alkesgrp_2024_EUR_h2_summary_file="$1"
alkesgrp_2024_EAS_h2_summary_file="$2"
alkesgrp_2024_EUR_genetic_corr_summary_file="$3"
alkesgrp_2024_EAS_genetic_corr_summary_file="$4"
alkesgrp_2024_sumstats_summary_file="$5"
qc_comparison_to_2021_update_dir="$6"

if false; then
python3 organize_data_for_qc_comparison_to_2021_sumstat_update.py $alkesgrp_2024_EUR_h2_summary_file $alkesgrp_2024_EAS_h2_summary_file $alkesgrp_2024_EUR_genetic_corr_summary_file $alkesgrp_2024_EAS_genetic_corr_summary_file $alkesgrp_2024_sumstats_summary_file $qc_comparison_to_2021_update_dir
fi


Rscript visualize_qc_analyses.R $alkesgrp_2024_EUR_h2_summary_file $alkesgrp_2024_EAS_h2_summary_file $alkesgrp_2024_EUR_genetic_corr_summary_file $alkesgrp_2024_EAS_genetic_corr_summary_file $qc_comparison_to_2021_update_dir
