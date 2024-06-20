import numpy as np
import os
import sys
import pdb


def extract_genetic_correlation_and_covariance_from_ldsc_rg_file(ldsc_rg_file):
	# File exists
	f = open(ldsc_rg_file)

	# Booleans to keep track that file parsing works
	genetic_cov_bool = False
	genetic_corr_bool = False

	for line in f:
		line = line.rstrip()
		if line.startswith('Genetic Covariance'):
			tmp = f.readline()
			tmp = f.readline()
			genetic_cov = tmp.split(': ')[1].split(' ')[0]
			genetic_cov_se = tmp.split(': ')[1].split('(')[1].split(')')[0]
			genetic_cov_bool = True
		if line.startswith('Genetic Correlation: '):
			genetic_corr = line.split(': ')[1].split(' ')[0]
			genetic_corr_se = line.split(': ')[1].split('(')[1].split(')')[0]
			tmp = f.readline()
			genetic_corr_z = tmp.rstrip().split(': ')[1]
			tmp = f.readline()
			genetic_corr_p = tmp.rstrip().split(': ')[1]
			genetic_corr_bool = True
	f.close()

	# Make sure parsing was successful
	# 16 of the pairs had error: "FloatingPointError: invalid value encountered in sqrt"
	# This was consistent with Martin's 2021 sumstat analysis
	if genetic_cov_bool == False or genetic_corr_bool == False:
		return 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'


	return genetic_corr, genetic_corr_se, genetic_corr_z, genetic_corr_p, genetic_cov, genetic_cov_se










############################
# Command line args
############################
trait_pairs_input_file = sys.argv[1]  # File containing all pairs of traits
ldsc_rg_results_dir = sys.argv[2]  # Directory containing ldsc genetic correlation results
alkesgrp_2024_EUR_genetic_corr_summary_file = sys.argv[3]  # Ouput file to keep summarized, organized results


# Open output file handle
t = open(alkesgrp_2024_EUR_genetic_corr_summary_file,'w')
# Write header
t.write('trait1\ttrait2\tgenetic_correlation\tgenetic_correlation_se\tgenetic_correlation_z\tgenetic_correlation_p\tgenetic_covariance\tgenetic_covariance_se\n')



# Loop through all pairs of traits
f = open(trait_pairs_input_file)
# Used to skip header
head_count = 0
# Stream file
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract trait pair for this line
	trait_1 = data[0]
	trait_2 = data[1]

	# LDSC rg file corresponding to the trait pair
	ldsc_rg_file = ldsc_rg_results_dir + trait_1 + ':' + trait_2 + '_ldsc_rg.log'


	# Extract genetic covariance and correlation estimates from LDSC
	genetic_correlation, genetic_correlation_se, genetic_correlation_z, genetic_correlation_p, genetic_covariance, genetic_covariance_se = extract_genetic_correlation_and_covariance_from_ldsc_rg_file(ldsc_rg_file)


	# Print result to output file
	t.write(trait_1 + '\t' + trait_2 + '\t' + genetic_correlation + '\t' + genetic_correlation_se + '\t' + genetic_correlation_z + '\t' + genetic_correlation_p + '\t' + genetic_covariance + '\t' + genetic_covariance_se + '\n')


f.close()