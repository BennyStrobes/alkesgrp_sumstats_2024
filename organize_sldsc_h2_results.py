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











##############################
# Command line args
##############################
trait_list_file = sys.argv[1]  # File containing a list of traits analyzed
sldsc_results_dir = sys.argv[2]  # Directory containing sldsc resutls
output_sldsc_results_summary_file = sys.argv[3]  # Output file



# Open output file handle
t = open(output_sldsc_results_summary_file,'w')
# Write header
t.write('trait_identifier\th2\th2_se\th2_z\tintercept\n')


# Open filehandle containing trait names
f = open(trait_list_file)
# Variable used to skip header
head_count = 0
# Stream trait names file
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue
	# Extract trait identifier from line
	trait_identifier = data[1]

	# sldsc results filename corresponding to trait identifier
	sldsc_log_file = sldsc_results_dir + trait_identifier + '_sldsc_baselineLD_v2.2.log'

	# Extract h2 estimates from sldsc log file
	h2, h2_se, h2_z, intercept = extract_h2_estimates_from_sldsc_log_file(sldsc_log_file)

	# Print to output file
	t.write(trait_identifier + '\t' + h2 + '\t' + h2_se + '\t' + h2_z + '\t' + intercept + '\n')


f.close()
t.close()

