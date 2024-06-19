import numpy as np
import os
import sys
import pdb




def print_all_pairs_of_traits_to_trait_pair_output_file(alkesgrp_2024_sumstats_summary_file, trait_pair_output_file):
	# First extract list of all traits
	traits = []
	f = open(alkesgrp_2024_sumstats_summary_file)
	head_count = 0
	# Stream summary stat summary file
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Error checking
		if len(data) != 9:
			print('assumption error')
		# Only consider traits with signed z-scores for genetic correlation analysis
		if data[7] != 'YES':
			continue
		# Only consider traits with allele info for genetic correlation analysis
		if data[8] != 'YES':
			continue
		# Use trait-identifier as trait_name
		traits.append(data[1])
	f.close()
	traits = np.asarray(traits)


	# Print all pairs of traits to output file
	t = open(trait_pair_output_file,'w')
	# header
	t.write('Trait1\tTrait2\n')
	# Loop through all pairs
	for trait1_index, trait1 in enumerate(traits):
		for trait2 in traits[trait1_index + 1:]:
			# Print to output
			t.write(trait1 + '\t' + trait2 + '\n')
			# Quick error check
			if trait1 == trait2:
				print('trait pair assumption error')
				pdb.set_trace()
	t.close()
	return








####################
# Command line args
####################
trait_summary_file = sys.argv[1]  # Input file
trait_pair_file = sys.argv[2]  # output file



print_all_pairs_of_traits_to_trait_pair_output_file(trait_summary_file, trait_pair_file)