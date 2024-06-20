import numpy as np
import os
import sys
import pdb





def extract_h2_estimates_from_sldsc_log_file(sldsc_log_file):
	f = open(sldsc_log_file)

	# Boolean variables to make sure file has corrrect lines
	h2_bool = False
	intercept_bool = False

	for line in f:
		line = line.rstrip()

		if line.startswith('Total Observed scale h2'):
			h2 = line.split(': ')[1].split(' (')[0]
			h2_se = line.split(': ')[1].split(' (')[1].split(')')[0]
			h2_z = str(round(float(h2)/float(h2_se), 5))
			h2_bool = True

		if line.startswith('Intercept'):
			intercept = line.split(': ')[1].split(' (')[0]
			intercept_bool = True

	f.close()

	if h2_bool == False or intercept_bool == False:
		print('error in extracting h2 from sldsc log file')
		pdb.set_trace()

	return h2, h2_se, h2_z, intercept




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
		return 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'

	return genetic_corr, genetic_corr_se, genetic_corr_z, genetic_corr_p, genetic_cov, genetic_cov_se





def compare_rgs(alkesgrp_2024_EUR_genetic_corr_summary_file, alkesgrp_2024_EAS_genetic_corr_summary_file, rg_2021_eur_dir,rg_2021_eas_dir, pre_2024_traits, output_file):
	t = open(output_file,'w')
	t.write('trait_name1\ttrait_name2\trg_2024\trg_2021\n')


	#EUR
	f = open(alkesgrp_2024_EUR_genetic_corr_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait1 = data[0]
		trait2 = data[1]
		genetic_correlation = data[2]
		genetic_correlation_se = data[3]
		if trait1 not in pre_2024_traits or trait2 not in pre_2024_traits:
			continue

		if genetic_correlation == 'nan':
			continue

		ldsc_rg_2021_log_file1 = rg_2021_eur_dir + 'cor.' + trait1 + '.' + trait2 + '.log'
		ldsc_rg_2021_log_file2 = rg_2021_eur_dir + 'cor.' + trait2 + '.' + trait1 + '.log'

		if os.path.isfile(ldsc_rg_2021_log_file1):
			ldsc_rg_2021_log_file = ldsc_rg_2021_log_file1
		elif os.path.isfile(ldsc_rg_2021_log_file2):
			ldsc_rg_2021_log_file = ldsc_rg_2021_log_file2
		else:
			print('assumption error: 2021 log file doesnt exist')
			pdb.set_trace()


		genetic_corr_2021, genetic_corr_se_2021, genetic_corr_z_2021, genetic_corr_p_2021, genetic_cov_2021, genetic_cov_se_2021 = extract_genetic_correlation_and_covariance_from_ldsc_rg_file(ldsc_rg_2021_log_file)

		t.write(trait1 + '\t' + trait2 + '\t' + genetic_correlation + '\t' + genetic_corr_2021 + '\n')

	f.close()

	#EAS
	f = open(alkesgrp_2024_EAS_genetic_corr_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait1 = data[0]
		trait2 = data[1]
		genetic_correlation = data[2]
		genetic_correlation_se = data[3]
		if trait1 not in pre_2024_traits or trait2 not in pre_2024_traits:
			continue

		if genetic_correlation == 'nan':
			continue

		ldsc_rg_2021_log_file1 = rg_2021_eas_dir + 'cor.' + trait1 + '.' + trait2 + '.log'
		ldsc_rg_2021_log_file2 = rg_2021_eas_dir + 'cor.' + trait2 + '.' + trait1 + '.log'

		if os.path.isfile(ldsc_rg_2021_log_file1):
			ldsc_rg_2021_log_file = ldsc_rg_2021_log_file1
		elif os.path.isfile(ldsc_rg_2021_log_file2):
			ldsc_rg_2021_log_file = ldsc_rg_2021_log_file2
		else:
			print('assumption error: 2021 log file doesnt exist')
			pdb.set_trace()


		genetic_corr_2021, genetic_corr_se_2021, genetic_corr_z_2021, genetic_corr_p_2021, genetic_cov_2021, genetic_cov_se_2021 = extract_genetic_correlation_and_covariance_from_ldsc_rg_file(ldsc_rg_2021_log_file)

		t.write(trait1 + '\t' + trait2 + '\t' + genetic_correlation + '\t' + genetic_corr_2021 + '\n')

	f.close()


	t.close()
	return


def compare_h2s(alkesgrp_2024_EUR_h2_summary_file, alkesgrp_2024_EAS_h2_summary_file, h2_2021_eur_dir, h2_2021_eas_dir, pre_2024_traits, output_file):
	t = open(output_file,'w')
	t.write('trait_name\th2_2024\th2_2021\th2_se_2024\th2_se_2021\th2_z_2024\th2_z_2021\tintercept_2024\tintercept_2021\n')

	# EUR
	f = open(alkesgrp_2024_EUR_h2_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_identifier = data[0]
		h2 = data[1]
		h2_se = data[2]
		h2_z = data[3]
		intercept = data[4]
		if trait_identifier not in pre_2024_traits:
			continue

		sldsc_2021_log_file = h2_2021_eur_dir + trait_identifier + '.baselineLD_1000G.log'

		h2_2021, h2_se_2021, h2_z_2021, intercept_2021 = extract_h2_estimates_from_sldsc_log_file(sldsc_2021_log_file)

		t.write(trait_identifier + '\t' + h2 + '\t' + h2_2021 + '\t' + h2_se + '\t' + h2_se_2021 + '\t' + h2_z + '\t' + h2_z_2021 + '\t' + intercept + '\t' + intercept_2021 + '\n')
	f.close()

	# EAS
	f = open(alkesgrp_2024_EAS_h2_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_identifier = data[0]
		h2 = data[1]
		h2_se = data[2]
		h2_z = data[3]
		intercept = data[4]
		if trait_identifier not in pre_2024_traits:
			continue

		sldsc_2021_log_file = h2_2021_eas_dir + trait_identifier + '.baselineLD_1000G.log'

		h2_2021, h2_se_2021, h2_z_2021, intercept_2021 = extract_h2_estimates_from_sldsc_log_file(sldsc_2021_log_file)

		t.write(trait_identifier + '\t' + h2 + '\t' + h2_2021 + '\t' + h2_se + '\t' + h2_se_2021 + '\t' + h2_z + '\t' + h2_z_2021 + '\t' + intercept + '\t' + intercept_2021 + '\n')
	f.close()


	t.close()
	return



def extract_pre_2024_traits(alkesgrp_2024_sumstats_summary_file):
	f = open(alkesgrp_2024_sumstats_summary_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 9:
			print('assumption eroorro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count +1
			continue
		if data[2] == 'New 2024':
			continue
		dicti[data[1]] = 1
	f.close()

	return dicti





#####################
# Command line args
#####################
alkesgrp_2024_EUR_h2_summary_file = sys.argv[1]
alkesgrp_2024_EAS_h2_summary_file = sys.argv[2]
alkesgrp_2024_EUR_genetic_corr_summary_file = sys.argv[3] 
alkesgrp_2024_EAS_genetic_corr_summary_file = sys.argv[4]
alkesgrp_2024_sumstats_summary_file = sys.argv[5]
qc_comparison_to_2021_update_dir = sys.argv[6]


# File directories corresponding to 2021 update
h2_2021_eur_dir='/n/groups/price/ldsc/sumstats_formatted_2021/baselineLD_v2.2/res/'
h2_2021_eas_dir='/n/groups/price/ldsc/sumstats_formatted_2021/baselineLD_v2.2_EAS/res/'
rg_2021_eur_dir='/n/groups/price/ldsc/sumstats_formatted_2021/genet_cor/cor/'
rg_2021_eas_dir='/n/groups/price/ldsc/sumstats_formatted_2021/genet_cor_EAS/cor/'


# Extract pre-2024 trait identifiers
pre_2024_traits = extract_pre_2024_traits(alkesgrp_2024_sumstats_summary_file)


######################################
# Compare heritabilities
######################################
h2_comparison_file=qc_comparison_to_2021_update_dir + 'h2_comparison.txt'
compare_h2s(alkesgrp_2024_EUR_h2_summary_file, alkesgrp_2024_EAS_h2_summary_file, h2_2021_eur_dir, h2_2021_eas_dir, pre_2024_traits, h2_comparison_file)


######################################
# Compare genetic correlations
######################################
rg_comparison_file=qc_comparison_to_2021_update_dir + 'rg_comparison.txt'
compare_rgs(alkesgrp_2024_EUR_genetic_corr_summary_file, alkesgrp_2024_EAS_genetic_corr_summary_file, rg_2021_eur_dir,rg_2021_eas_dir, pre_2024_traits, rg_comparison_file)



