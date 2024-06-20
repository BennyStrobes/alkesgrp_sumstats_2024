import numpy as np
import os
import sys
import pdb






def extract_valid_traits(alkesgrp_2024_EUR_genetic_corr_summary_file, alkesgrp_2024_sumstats_summary_file):
	head_count = 0
	nan_traits = {}
	valid_traits = {}
	f = open(alkesgrp_2024_EUR_genetic_corr_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait1 = data[0]
		trait2 = data[1]
		corr = data[2]
		valid_traits[trait1] =1
		valid_traits[trait2] = 1

		if corr == 'nan':
			nan_traits[trait1] = 1
			nan_traits[trait2] = 1
	f.close()

	final_trait_dicti = {}
	final_trait_arr = []
	for trait_name in [*valid_traits]:
		if trait_name in nan_traits:
			continue
		final_trait_dicti[trait_name] = 1
		final_trait_arr.append(trait_name)
	return np.asarray(final_trait_arr), final_trait_dicti

def create_mapping_from_trait_pair_to_squared_rg(alkesgrp_2024_EUR_genetic_corr_summary_file, valid_trait_dictionary):
	f = open(alkesgrp_2024_EUR_genetic_corr_summary_file)
	head_count = 0
	trait_pair_to_squared_rg = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait1 = data[0]
		trait2 = data[1]
		corr = data[2]
		if trait1 not in valid_trait_dictionary or trait2 not in valid_trait_dictionary:
			continue
		if corr == 'nan':
			print('assumption eroror')
			pdb.set_trace()

		squared_rg = np.square(float(corr))

		if trait1 + ':' + trait2 in trait_pair_to_squared_rg or trait2 + ':' + trait1 in trait_pair_to_squared_rg:
			print('assumption eroor')
			pdb.set_trace()

		trait_pair_to_squared_rg[trait1 + ':' + trait2] = squared_rg
		trait_pair_to_squared_rg[trait2 + ':' + trait1] = squared_rg

	f.close()

	return trait_pair_to_squared_rg


def create_mapping_from_trait_name_to_h2_info(alkesgrp_2024_EUR_h2_summary_file, valid_trait_dictionary):
	f = open(alkesgrp_2024_EUR_h2_summary_file)
	head_count = 0
	dicti = {}

	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		trait_name = data[0]
		h2 = data[1]
		h2_se = data[2]
		h2_z = data[3]
		h2_intercept = data[4]

		if trait_name in dicti:
			print('assumption erororo')
			pdb.set_trace()

		dicti[trait_name] = (h2, h2_se, h2_z, h2_intercept)

	f.close()
	return dicti


def make_non_redundent_trait_file(alkesgrp_2024_EUR_non_redundent_trait_file, valid_trait_list, trait_pair_to_squared_rg, trait_name_to_h2_info):
	# Open output file handle
	t = open(alkesgrp_2024_EUR_non_redundent_trait_file,'w')
	t.write('trait_identifier\th2\th2_se\th2_z\tintercept\n')

	# Initialize list of traits and z-scores
	current_trait_list = []
	current_z_score_list = []
	for trait_name in valid_trait_list:
		sldsc_h2_z = float(trait_name_to_h2_info[trait_name][2])
		if sldsc_h2_z > 6:
			current_trait_list.append(trait_name)
			current_z_score_list.append(sldsc_h2_z)
	current_trait_list = np.asarray(current_trait_list)
	current_z_score_list = np.asarray(current_z_score_list)


	# Keep looping until we run out of traits
	while len(current_trait_list) > 0:
		# Sort traits by sldsc-h2 z-score
		ordered_current_z_score_list = current_z_score_list[np.argsort(-current_z_score_list)]
		ordered_current_trait_list = current_trait_list[np.argsort(-current_z_score_list)]

		# Add best trait to output file
		best_trait = ordered_current_trait_list[0]
		best_trait_h2_info = trait_name_to_h2_info[best_trait]
		t.write(best_trait + '\t' + best_trait_h2_info[0] + '\t' + best_trait_h2_info[1] + '\t' + best_trait_h2_info[2] + '\t' + best_trait_h2_info[3] + '\n')

		# Now make dictionary list of which traits to remove this round
		traits_to_remove = {}
		for trait_name in current_trait_list:
			if trait_name == best_trait:
				traits_to_remove[trait_name]=1
			else:
				squared_correlation = trait_pair_to_squared_rg[trait_name + ':' + best_trait]
				if squared_correlation > .1:
					traits_to_remove[trait_name] =1

		# Now create list of traits for next round
		new_traits = []
		new_trait_z_scores = []
		for ii, trait_name in enumerate(current_trait_list):
			trait_z_score = current_z_score_list[ii]
			if trait_name not in traits_to_remove:
				new_traits.append(trait_name)
				new_trait_z_scores.append(trait_z_score)

		# Update current trait list and current z-score list
		current_trait_list = np.asarray(new_traits)
		current_z_score_list = np.asarray(new_trait_z_scores)
	t.close()
	return



########################
# Command line args
########################
alkesgrp_2024_sumstats_summary_file = sys.argv[1]  # Full list of summary stats
alkesgrp_2024_EUR_h2_summary_file = sys.argv[2]  # Heritability file
alkesgrp_2024_EUR_genetic_corr_summary_file = sys.argv[3]  # genetic correlation file
alkesgrp_2024_EUR_non_redundent_trait_file = sys.argv[4]  # Output file



# First get list of valid traits to consider
# We are going to limit to those that:
## 1. have signed summary statistics
## 2. Have allele info
## 3. Had successful convergence of cross-trait LDSC with all other traits
valid_trait_list, valid_trait_dictionary = extract_valid_traits(alkesgrp_2024_EUR_genetic_corr_summary_file, alkesgrp_2024_sumstats_summary_file)


# Second, create mapping from trait pair name to squared genetic correlation
trait_pair_to_squared_rg = create_mapping_from_trait_pair_to_squared_rg(alkesgrp_2024_EUR_genetic_corr_summary_file, valid_trait_dictionary)

# Third, create mapping from trait name to h2 info
trait_name_to_h2_info = create_mapping_from_trait_name_to_h2_info(alkesgrp_2024_EUR_h2_summary_file, valid_trait_dictionary)


# Make non-redudent traits file
make_non_redundent_trait_file(alkesgrp_2024_EUR_non_redundent_trait_file, valid_trait_list, trait_pair_to_squared_rg, trait_name_to_h2_info)
