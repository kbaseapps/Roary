'''
Download a GenomeSet
'''
import os
import sys
import subprocess

from multi_key_dict import multi_key_dict

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil


def download_gffs(cb_url, scratch, genome_set_ref):
	"""
	Args:
	cb_url - callback server URL
	scratch - scratch work folder
	genome_set_ref - reference to genome_set object in workspace
	Returns the path to the folder containing .gff files


	we want to first handle a GenomeSet Object "KBaseSearch.GenomeSet" or "KBaseSets.GenomeSet"
	"""

	# Get our utilities
	dfu = DataFileUtil(cb_url)
	au = AssemblyUtil(cb_url)
	gfu = GenomeFileUtil(cb_url)

	obj_data = dfu.get_objects({'object_refs': [genome_set_ref]})['data'][0]
	gs_obj = obj_data['data']
	obj_type = obj_data['info'][2]

	if 'KBaseSets.GenomeSet' in obj_type:
		refs = [gsi['ref'] for gsi in gs_obj['items']]
	elif 'KBaseSearch.GenomeSet' in obj_type:
		refs = [gse['ref'] for gse in gs_obj['elements'].values()]
	else:
		raise TypeError(
			'provided input must of type KBaseSets.GenomeSet or KBaseSearch.GenomeSet not ' +
			str *
			(obj_type))

	if len(refs) < 2:
		raise ValueError("Must provide GenomeSet with at least 2 Genomes.")

	# name the output directory
	temp_dir = scratch + '/temp'
	final_dir = scratch + '/gff'

	os.mkdir(final_dir)
	os.mkdir(temp_dir)

	# write file that will help us cat the gff and fasta files
	cat_path = scratch + '/fast_cat.txt'

	cat_file = open(cat_path, 'w')
	cat_file.write("##FASTA")
	cat_file.close()

	path_to_ref_and_ID_pos_dict = {}
	for ref in refs:
		gen_obj = dfu.get_objects({'object_refs': [ref]})['data'][0]['data']

		# figure out if genome is of acceptable type.
		# NO Eukaryotes, NO Fungi,
		# yes bacateria, yes archaea, yes(?) virus
		if gen_obj['domain'] not in ['Bacteria', 'Archaea']:
			raise TypeError(
				'Provided Genomes are not labeled as Bacteria or Archaea. Roary is only equipped to handle Archaea or Bacteria')

		fasta_path = temp_dir + "/" + gen_obj['id'] + ".fa"
		gff_file = gfu.genome_to_gff({'genome_ref': ref, 'target_dir': temp_dir})
		if 'assembly_ref' not in gen_obj.keys():
			raise TypeError("All genomes must contain an 'assembly_ref'")
		else:
			fasta_file = au.get_assembly_as_fasta(
				{'ref': gen_obj['assembly_ref'], 'filename': fasta_path})
			# check that fasta_file exists
			if not os.path.isfile(fasta_file['path']):
				raise ValueError('An input Genome does not have an associated FASTA file.')

		# need to figure out if FASTA is already in gff file
		# not sure if we need to do this step.
		if 'path' in gff_file:
			gff_file_path = gff_file['path']
		elif 'file_path' in gff_file:
			gff_file_path = gff_file['file_path']
		elif 'gff_file' in gff_file:
			gff_file_path = gff_file['gff_file']['path']
		else:
			raise ValueError("No GFF File Path found.")
		gff_file_path, ID_to_pos, contains_fasta = filter_gff(gff_file_path)

		new_file_path = final_dir + "/" + gen_obj['id'] + ".gff"

		if contains_fasta:
			args = ['mv', gff_file_path, new_file_path]
			subprocess.call(args)
		else:
			# NOTE: We have to pipe output of cat call to the new_file_path
			# next we make a new 'gff' file that contains both the gff and fasta information
			args = ['cat', fasta_file['path'], cat_path, gff_file_path]
			catted_files = subprocess.check_output(args, stdout=subprocess.PIPE)
			f = open(new_file_path, 'w')
			f.write(catted_files.decode('utf-8'))
			f.close()

		path_to_ref_and_ID_pos_dict[new_file_path] = (ref, ID_to_pos)

	return final_dir, path_to_ref_and_ID_pos_dict


def make_some_trash(gff_file):
	f = open(gff_file)
	output = []
	IDs = set([])
	cds_to_pos = {}
	gene_pos = 0
	for l in f:
		if '##FASTA' in l:
			output.append(l)
			output += [j for j in f]
			break
		elif l[:2] == '##':
			output.append(l)
			continue
		ID = l.split(';')[0].split('=')[-1]
		feat_type = l.split()[2]
		if feat_type == 'CDS':
			cds_to_pos[ID] = gene_pos
			gene_pos += 1
		output.append(l)
	f.close()

	# write output to file
	f = open(gff_file, 'w')
	for l in output:
		f.write(l)
	f.close()

	return gff_file, cds_to_pos, True


def filter_gff(gff_file):
	# if there is a duplicate ID's toss one of them out (if it is programmed frameshift toss one)
	# either the smaller or the second one.
	# for feature pos, use the order of the gff 'gene' items. create a mapping
	# from ID to order for feature pos
	f = open(gff_file)
	output = []
	IDs = set([])
	# definitely a better way of doing this

	# ideally we would do the feature mapping o th genes, but there are nt always 'gene' features in the gff files
	# so to get around this we base the gene position in the gff off the CDS position because they are most often
	# consecutively ordered and paired. For the organisms that Roary is supposed to service, there should only be
	# a one to one pairing of 'gene' to 'CDS'

	# gene_to_pos = {}
	# cds_to_gene = {}

	cds_to_pos = {}
	gene_pos = 0
	contains_fasta = False
	for l in f:
		if "##FASTA" in l:
			# add FASTA sequence to output if it exists and then exit
			contains_fasta = True
			output.append(l)
			output += [j for j in f]
			break
		# if we find a comment add it to output
		if l[0:2] == '##':
			output.append(l)
			continue

		# get ID of feature
		ID = l.split(';')[0].split('=')[-1]
		# determine type of feature
		feat_type = l.split()[2]
		# if feat_type == 'gene':
		# 	gene_to_pos[ID] = gene_pos
		# 	gene_pos += 1
		if feat_type == 'CDS':
			cds_to_pos[ID] = gene_pos
			gene_pos += 1
			# # find parent
			# parent_split = l.split('Parent=')
			# if len(parent_split) == 2:
			# 	parent_gene_ID = parent_split[1].split(';')[0]
			# 	cds_to_gene[ID] = parent_gene_ID
			# elif len(parent_split) == 1:
			# 	# there is no parent to this feature.

		ids_len = len(IDs)
		IDs.add(ID)
		if len(IDs) == ids_len:
			# we found a duplicate and its the second one. get rid of it
			continue
		else:
			# make sure that the feature we get is a 'CDS' object
			if feat_type == 'CDS':
				output.append(l)
	f.close()

	# write output to file
	f = open(gff_file, 'w')
	for l in output:
		f.write(l)
	f.close()

	return gff_file, cds_to_pos, contains_fasta
