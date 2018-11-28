import sys
import os
import subprocess
import json

def run_roary(scratch, data_path):
	"""
	"""
	data_path_files = data_path + '/*.gff'
	out_dir = scratch + '/outputs/'

	# arguments we probably want to include in UI
	blastp_percentage_identity = str(0.98)
	no_split_paralogs = False
	# default number
	max_num_clusters = str(50000)

	# internal arguments
	graphing = False
	num_threads = str(4)

	# find out what other arguments we want to use for this from Paramvir
	args = ['roary','-f', out_dir, '-p', num_threads, '-i', blastp_percentage_identity,\
			'-g', max_num_clusters] 
	if no_split_paralogs:
		args.append('-s')

	if graphing:
		args.append('-r')
	# finally add files on which to run Roary
	args.append(data_path_files)
	# not aware of the errors for Roary yet, so we'll just use .call for now
	proc_status = subprocess.call(args)
	if proc_status == 0:
		continue
	else:
		raise RuntimeError('Roary or subprocess ran into an error')

	return out_dir
