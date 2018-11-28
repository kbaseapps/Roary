'''
Download a GenomeSet
'''
import os
import sys
import subprocess

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil

def download_gffs(genome_set_ref, cb_url, scratch):
	"""
	Args:
		refs - workspace references in the form 'workspace_id/object_id/obj_version'
    	cb_url - callback server URL
	Returns the path to the folder containing .gff files


	we want to first handle a GenomeSet Object "KBaseSearch.GenomeSet" or "KBaseSets.GenomeSet"
	"""

	# Get our utilities
	dfu = DataFileUtil(cb_url)
	au  = AssemblyUtil(cb_url)
	gfu = GenoeFileUtil(cb_url)


	obj_data = dfu.get_objects({'object_refs':[genome_set_ref]})
	gs_obj = obj_data['data'][0]
	obj_type = obj_data['info'][2]

	if 'KBaseSets.GenomeSet' in obj_type:
		refs = [gsi['ref'] for gsi in gs_obj['items']]
	elif 'KBaseSearch.GenomeSet' in obj_type:
		refs = [gse['ref'] for gse in gs_obj['elements']]
	else:
		raise TypeError('provided input must of type KBaseSets.GenomeSet or KBaseSearch.GenomeSet not '+ str*(obj_type))

	# name the output directory
	temp_dir = scratch + '/temp'
	final_dir = scratch + '/gff'

	os.mkdir(final_dir)
	os.mkdir(temp_dir)

	# write file that will help us cat the gff and fasta files
	cat_path = scratch + '/fast_cat.txt'

	cat_file = open(cat_path,'w')	
	cat_file.write("##FASTA")
	cat_file.close()

	path_to_ref = {}
	for ref in refs:
		gen_obj = dfu.get_objects({'object_refs':[ref]})['data'][0]

		fasta_path = temp_dir + "/" + gen_obj['id'] + ".fa"
		fasta_file = au.get_assembly_as_fasta({'ref':gen_obj['assembly_ref'], 'filename': fasta_path})
		gff_file = gfu.genome_to_gff({'genome_ref':ref, 'target_dir',temp_dir})

		new_file_path = final_dir + "/" + gen_obj['id'] + ".gff"
		# next we make a new 'gff' file that contains both the gff and fasta information
		args = ['cat',fasta_file['path'], cat_path, gff_file['path'], '>', new_file_path]
		subprocess.call(args)

		path_to_ref[new_file_path] = ref

	return final_dir, path_to_ref