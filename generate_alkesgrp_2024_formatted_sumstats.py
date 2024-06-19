import numpy as np
import os
import sys
import pdb
import gzip
import scipy.stats


def extract_dictionary_list_of_hapmap3_rsids(hapmap3_rsid_file):
	# Now extract dictionary list of hapmap3 rsids
	hm3_rsids = {}
	f = open(hapmap3_rsid_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[0]
		a1 = data[1]
		a2 = data[2]
		if rsid.startswith('rs') == False:
			print('rsid assumption eroror')
			pdb.set_trace()
		if rsid in hm3_rsids:
			print('repeat hm3 assumption error')
			pdb.set_trace()
		if a1 == a2:
			print('assumption erooror')
			pdb.set_trace()
		hm3_rsids[rsid] = a1 + '_' + a2
	f.close()
	return hm3_rsids

def check_if_sumstat_file_has_a1_a2(sumstat_2021_file):
	if sumstat_2021_file.endswith('.gz'):
		f = gzip.open(sumstat_2021_file)
		zipped = True
	else:
		f = open(sumstat_2021_file)
		zipped = False
	sumstat_has_a1_a2 = False
	for line in f:
		if zipped:
			line = line.rstrip().decode("utf-8")
		else:
			line = line.rstrip()
		data = line.split()
		data = np.asarray(data)
		if 'A1' in data:
			sumstat_has_a1_a2 = True

		break

	f.close()

	return sumstat_has_a1_a2


def print_gazal_lab_summary_statistics_into_2024_sumstats_dir(gazal_lab_sumstats_dir, gazal_lab_sumstats_summary_file, hapmap3_snplist_file, hapmap3_rsids, ldsc_code_dir, alkesgrp_2024_sumstats_dir, t_summary_2024):
	# First extract names of traits according to gazal lab sumstats summary file
	trait_names = []
	trait_identifiers = []
	gazal_trait_identifiers = []
	trait_origins = []
	trait_populations = []
	trait_references = []	
	f = open(gazal_lab_sumstats_summary_file)
	head_count = 0
	for line in f:
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		line = line.rstrip()
		data = line.split(',')
		trait_name = data[2]
		trait_identifier = data[1]
		gazal_trait_identifier = data[0]
		trait_reference = data[3]
		# Skip this file as there is clearly something wrong with the z-scores
		if trait_identifier == 'PASS_VenousThromboembolism_Ghouse2023':
			continue
		# Add to global arrays
		trait_names.append(trait_name)
		trait_identifiers.append(trait_identifier)
		gazal_trait_identifiers.append(gazal_trait_identifier)
		trait_origins.append('New 2024')
		trait_populations.append('EUR')
		trait_references.append(trait_reference)
	f.close()

	# Convert to numpy arrays
	trait_names = np.asarray(trait_names)
	trait_identifiers = np.asarray(trait_identifiers)
	gazal_trait_identifiers = np.asarray(gazal_trait_identifiers)
	trait_origins = np.asarray(trait_origins)
	trait_populations = np.asarray(trait_populations)
	trait_references = np.asarray(trait_references)



	# Now loop through traits, and for each trait:
	# 1. Print trait summary statistic file to alkesgrp_2024_sumstats_dir in summary stat format
	# 2. Add relevent fields to t_summary_2024 file handle

	for trait_iter, trait_name in enumerate(trait_names):
		trait_identifier = trait_identifiers[trait_iter]
		trait_origin = trait_origins[trait_iter]
		trait_population = trait_populations[trait_iter]
		trait_reference = trait_references[trait_iter]
		gazal_trait_identifier = gazal_trait_identifiers[trait_iter]
		print(trait_identifier)


		# 1. Print trait summary statistic file to alkesgrp_2024_sumstats_dir in summary stat format
		sumstat_gazal_file = gazal_lab_sumstats_dir + gazal_trait_identifier + '.sumstats.gz'  # input file
		sumstat_2024_file = alkesgrp_2024_sumstats_dir + trait_identifier + '.sumstats'  # output file

		# Check if 2021 file has A1 and A2 columns
		sumstat_has_a1_a2 = check_if_sumstat_file_has_a1_a2(sumstat_gazal_file)

		if sumstat_has_a1_a2:
			MM, NN, signed_sumstat = print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir(sumstat_gazal_file, sumstat_2024_file, hapmap3_rsids)
			allele_info = 'YES'
		else:
			MM, NN, signed_sumstat = print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir_w_o_a1_a2(sumstat_gazal_file, sumstat_2024_file, hapmap3_rsids)
			allele_info = 'NO'

		t_summary_2024.write(trait_name + '\t' + trait_identifier + '\t' +  trait_origin + '\t' + trait_reference + '\t' + str(NN) + '\t' + str(MM) + '\t' + trait_population + '\t' + signed_sumstat + '\t' + allele_info + '\n')
		t_summary_2024.flush()

	return t_summary_2024



def print_alkesgrp_2021_summary_statistics_into_2024_sumstats_dir(alkesgrp_2021_sumstats_dir, alkesgrp_2021_sumstats_summary_file, hapmap3_rsid_file, hapmap3_rsids, alkesgrp_2024_sumstats_dir, t_summary_2024):
	# First extract names of traits according to 2021 sumstats summary file
	trait_names = []
	trait_identifiers = []
	trait_origins = []
	trait_populations = []
	trait_references = []
	f = open(alkesgrp_2021_sumstats_summary_file)
	head_count = 0
	for line in f:
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		line = line.rstrip()
		# Extract relevent fields
		data = line.split(',')
		trait_identifier = data[1]
		trait_identifiers.append(trait_identifier)
		trait_name = data[0]
		# Fix trat name issue in two traits
		if trait_name.startswith('"'):
			trait_name = trait_name.split('"')[0]
		# Name of trait
		trait_names.append(trait_name)
		# Origin (2019, 2021, 2024, etc)
		trait_origin = data[2]
		trait_origins.append(trait_origin)
		# Gwas population
		trait_population= data[8]
		if trait_population != 'EUR' and trait_population != 'EAS':
			trait_population = data[7]
			if trait_population != 'EUR' and trait_population != 'EAS':
				print('assumption eroror')
		trait_populations.append(trait_population)

		# Get reference (ie. citation) for trait
		if trait_origin.startswith('UKBio'):
			trait_reference = 'Bycroft et al., 2018 Nature'
		else:
			if len(data[3].split('"')) != 2:
				trait_reference = data[3]
			else:
				trait_reference = data[3].split('"')[1] + ',' + data[4].split('"')[0]
		trait_references.append(trait_reference)

	f.close()
	# Convert to numpy arrays
	trait_names = np.asarray(trait_names)
	trait_identifiers = np.asarray(trait_identifiers)
	trait_origins = np.asarray(trait_origins)
	trait_populations = np.asarray(trait_populations)
	trait_references = np.asarray(trait_references)


	# Now loop through traits, and for each trait:
	# 1. Print trait summary statistic file to alkesgrp_2024_sumstats_dir in summary stat format
	# 2. Add relevent fields to t_summary_2024 file handle

	for trait_iter, trait_name in enumerate(trait_names):
		trait_identifier = trait_identifiers[trait_iter]
		trait_origin = trait_origins[trait_iter]
		trait_population = trait_populations[trait_iter]
		trait_reference = trait_references[trait_iter]

		print(trait_identifier)

		# 1. Print trait summary statistic file to alkesgrp_2024_sumstats_dir in summary stat format
		sumstat_2021_file = alkesgrp_2021_sumstats_dir + trait_identifier + '.sumstats'  # input file
		sumstat_2024_file = alkesgrp_2024_sumstats_dir + trait_identifier + '.sumstats'  # output file
		# Check if 2021 file has A1 and A2 columns
		sumstat_has_a1_a2 = check_if_sumstat_file_has_a1_a2(sumstat_2021_file)
		if sumstat_has_a1_a2:
			MM, NN, signed_sumstat = print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir(sumstat_2021_file, sumstat_2024_file, hapmap3_rsids)
			allele_info = 'YES'
		else:
			MM, NN, signed_sumstat = print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir_w_o_a1_a2(sumstat_2021_file, sumstat_2024_file, hapmap3_rsids)
			allele_info = 'NO'

		#os.system('python ' + ldsc_code_dir + 'munge_sumstats.py --sumstats ' + sumstat_2024_file + ' --out ' + alkesgrp_2024_sumstats_dir + trait_identifier + '.sumstats_tmp ' + '--merge-alleles ' + hapmap3_rsid_file)

		t_summary_2024.write(trait_name + '\t' + trait_identifier + '\t' +  trait_origin + '\t' + trait_reference + '\t' + str(NN) + '\t' + str(MM) + '\t' + trait_population + '\t' + signed_sumstat + '\t' + allele_info + '\n')
		t_summary_2024.flush()
	return t_summary_2024


def print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir_w_o_a1_a2(sumstat_2021_file, sumstat_2024_file, hapmap3_rsids):
		# Open output file handle
		t = open(sumstat_2024_file,'w')
		t.write('SNP\tN\tZ\n')

		# Keep track of various parameters describing sumstat file
		MM = 0 # number of snps
		n_vec = []  # vector to keep track of sample size at each snp
		pos_signed = False  # snp has positive signed z-score for some snp
		neg_signed = False  # snp has negative signed z-score for some snp

		# Open input file handle
		if sumstat_2021_file.endswith('.gz'):
			f = gzip.open(sumstat_2021_file)
			zipped = True
		else:
			f = open(sumstat_2021_file)
			zipped = False

		# Stream 2021 sum stat file
		head_count = 0
		for line in f:
			if zipped:
				line = line.rstrip().decode('utf-8')
			else:
				line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get column index of various identifiers for this sumstat file
				snp_index = np.where(np.asarray(data)=='SNP')[0][0]
				N_index = np.where(np.asarray(data)=='N')[0][0]
				Z_index = np.where(np.asarray(data)=='Z')[0][0]
				continue
			# Standard line of summary stat file
			# Each line corresponds to a snp

			# Extract relevent fields for this line
			snp = data[snp_index]
			snp_NN = data[N_index]
			Z = data[Z_index]

			# Quick error check to make sure snp is a hapmap3 snp
			if snp not in hapmap3_rsids:
				continue


			# Print to output file handle
			t.write(snp + '\t' + snp_NN + '\t' + Z + '\n')

			# Update other quantities
			MM = MM + 1
			n_vec.append(float(snp_NN))
			Z_float = float(Z)
			if Z_float > 0:
				pos_signed = True
			if Z_float < 0:
				neg_signed = True
		t.close()
		t.close()

		# Convert to neat numpy array
		n_vec = np.asarray(n_vec)
		if pos_signed and neg_signed:
			signed_sumstat = 'YES'
		else:
			signed_sumstat = 'NO'

		# Get average sample size
		avg_NN = np.round(np.mean(n_vec)).astype(int)

		return MM, avg_NN, signed_sumstat

def alleles_match(a1, a2, hapmap3_a1, hapmap3_a2):
	COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	booler = False

	# Stand and ref match
	if a1 == hapmap3_a1 and a2 == hapmap3_a2:
		booler = True
	# Ref flip and strand match
	if a1 == hapmap3_a2 and a2 == hapmap3_a1:
		booler = True
	# Ref match, strand flip
	if a1 == COMPLEMENT[hapmap3_a1] and a2 == COMPLEMENT[hapmap3_a2]:
		booler = True
	# Strand and ref flip
	if a1 == COMPLEMENT[hapmap3_a2] and a2 == COMPLEMENT[hapmap3_a1]:
		booler = True

	return booler


def print_single_alkesgrp_2021_summarystat_file_to_sumstat_2024_dir(sumstat_2021_file, sumstat_2024_file, hapmap3_rsids):
		# Open output file handle
		t = open(sumstat_2024_file,'w')
		t.write('SNP\tA1\tA2\tN\tZ\n')

		# Keep track of various parameters describing sumstat file
		MM = 0 # number of snps
		n_vec = []  # vector to keep track of sample size at each snp
		pos_signed = False  # snp has positive signed z-score for some snp
		neg_signed = False  # snp has negative signed z-score for some snp

		# Open input file handle
		if sumstat_2021_file.endswith('gz'):
			f = gzip.open(sumstat_2021_file)
			zipped = True
		else:
			f = open(sumstat_2021_file)
			zipped = False

		# Stream 2021 sum stat file
		head_count = 0
		allele_mismatches = 0
		for line in f:
			if zipped:
				line = line.rstrip().decode('utf-8')
			else:
				line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get column index of various identifiers for this sumstat file
				snp_index = np.where(np.asarray(data)=='SNP')[0][0]
				a1_index = np.where(np.asarray(data)=='A1')[0][0]
				a2_index = np.where(np.asarray(data)=='A2')[0][0]
				N_index = np.where(np.asarray(data)=='N')[0][0]
				Z_index = np.where(np.asarray(data)=='Z')[0][0]
				continue
			# Standard line of summary stat file
			# Each line corresponds to a snp

			# Some stat files currently have some lines with only a snp id and no information. Presumably for missing snps. filter those lines out
			if len(data) == 1:
				continue

			# Extract relevent fields for this line
			snp = data[snp_index]
			a1 = data[a1_index]
			a2 = data[a2_index]
			snp_NN = data[N_index]
			Z = data[Z_index]

			# Quick error check to make sure snp is a hapmap3 snp
			if snp not in hapmap3_rsids:
				continue
			# also make sure alleles match
			hapmap3_alleles = hapmap3_rsids[snp]
			hapmap3_a1 = hapmap3_alleles.split('_')[0]
			hapmap3_a2 = hapmap3_alleles.split('_')[1]

			if alleles_match(a1, a2, hapmap3_a1, hapmap3_a2) == False:
				allele_mismatches = allele_mismatches + 1
				continue


			# Print to output file handle
			t.write(snp + '\t' + a1 + '\t' + a2 + '\t' + snp_NN + '\t' + Z + '\n')

			# Update other quantities
			MM = MM + 1
			n_vec.append(float(snp_NN))
			Z_float = float(Z)
			if Z_float > 0:
				pos_signed = True
			if Z_float < 0:
				neg_signed = True
		t.close()
		t.close()

		# Convert to neat numpy array
		n_vec = np.asarray(n_vec)
		if pos_signed and neg_signed:
			signed_sumstat = 'YES'
		else:
			signed_sumstat = 'NO'

		# Get average sample size
		avg_NN = np.round(np.mean(n_vec)).astype(int)

		print(str(MM) + '\t' + str(allele_mismatches))

		return MM, avg_NN, signed_sumstat




############################
# Command line arguments
#############################
alkesgrp_2021_sumstats_dir = sys.argv[1]
alkesgrp_2021_sumstats_summary_file = sys.argv[2]
gazal_lab_sumstats_dir = sys.argv[3]
gazal_lab_sumstats_summary_file = sys.argv[4]
hapmap3_rsid_file = sys.argv[5]  # limit to these rsids
ldsc_code_dir = sys.argv[6]
alkesgrp_2024_sumstats_dir = sys.argv[7] # Output dir
alkesgrp_2024_sumstats_summary_file = sys.argv[8]  # Output summary file


# Extract dictionary list of hapmap3 rsids and MHC rsids
hapmap3_rsids = extract_dictionary_list_of_hapmap3_rsids(hapmap3_rsid_file)


# Open file-handle for writing to alkesgrp 2024 sumstats summary file
t_summary_2024 = open(alkesgrp_2024_sumstats_summary_file,'w')
t_summary_2024.write('Trait Name\tTrait Identifier\tOrigin\tReference\tN\tM\tPopulation\tSigned SumStat\tAllele Info\n')

# Print alkesgrp 2021 summary statistics into alkesgrp 2024 sumstats dir
t_summary_2024 = print_alkesgrp_2021_summary_statistics_into_2024_sumstats_dir(alkesgrp_2021_sumstats_dir, alkesgrp_2021_sumstats_summary_file,hapmap3_rsid_file, hapmap3_rsids, alkesgrp_2024_sumstats_dir, t_summary_2024)


# Print gazal lab summary statistics into alkesgrp 2024 sumstats dir
t_summary_2024 = print_gazal_lab_summary_statistics_into_2024_sumstats_dir(gazal_lab_sumstats_dir, gazal_lab_sumstats_summary_file, hapmap3_rsid_file, hapmap3_rsids, ldsc_code_dir, alkesgrp_2024_sumstats_dir, t_summary_2024)

# Close file handle
t_summary_2024.close()



