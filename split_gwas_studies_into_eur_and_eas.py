import numpy as np
import os
import sys
import pdb













######################
# Command line args
######################
alkesgrp_2024_sumstats_summary_file = sys.argv[1]
alkesgrp_2024_EUR_sumstats_summary_file = sys.argv[2]
alkesgrp_2024_EAS_sumstats_summary_file = sys.argv[3]


# Split alkesgrp_2024_sumstats_summary_file into alkesgrp_2024_EUR_sumstats_summary_file and alkesgrp_2024_EAS_sumstats_summary_file
f = open(alkesgrp_2024_sumstats_summary_file)
t_eur = open(alkesgrp_2024_EUR_sumstats_summary_file,'w')
t_eas = open(alkesgrp_2024_EAS_sumstats_summary_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		# Copy header for both output files
		t_eur.write(line + '\n')
		t_eas.write(line + '\n')
		continue
	if data[6] == 'EUR':
		t_eur.write(line + '\n')
	elif data[6] == 'EAS':
		t_eas.write(line + '\n')
	else:
		print('assumption error')
		pdb.set_trace()

f.close()
t_eur.close()
t_eas.close()