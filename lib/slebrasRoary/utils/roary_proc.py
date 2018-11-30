import sys
import os
import json

def run_roary(scratch, data_path, params):
	"""
	"""
	data_path_files = os.path.join(data_path,'*.gff')
	out_dir = os.path.join(scratch,'outputs')

	# arguments we probably want to include in UI
	blastp_percentage_identity = str(params['blast_p_percentage'])
	no_split_paralogs = False
	# default number
	max_num_clusters = str(params['max_num_clusters'])

	# internal arguments
	graphing = False
	num_threads = str(4)

	# find out what other arguments we want to use for this from Paramvir
	args = ['roary','-f', out_dir, '-p', num_threads, '-i', \
			blastp_percentage_identity,'-g', max_num_clusters] 
	if no_split_paralogs:
		args.append('-s')

	if graphing:
		args.append('-r')
	# finally add files on which to run Roary
	args.append(data_path_files)

	# Roary seems to return code 255 regardless of if it worked or not so we'll use os.system instead
	# of subprocess
	os.system(' '.join(args))

	return out_dir
