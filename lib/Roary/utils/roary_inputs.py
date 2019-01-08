'''
Download a GenomeSet
'''
import os
import sys
import subprocess
import random

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
	cat_file.write("##FASTA\n")
	cat_file.close()

	path_to_ref_and_ID_pos_dict = {}
	for ref in refs:
		gen_obj = dfu.get_objects({'object_refs': [ref]})['data'][0]['data']

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

		assert(os.path.isfile(gff_file_path)), "Could not find input GFF file for object with workspace reference: %s"%ref

		# oki doki, here we wanna make sure that the ID's in the genome object match up with
		# ID's in the gff file. This is importatnt because the pangenome object uses the genome
		# objects (in the pangenomeviewer).

		gff_file_path, ID_to_pos, gffid_to_genid, contains_fasta = filter_gff(gff_file_path, gen_obj)

		new_file_path = final_dir + "/" + gen_obj['id'] + '.gff'

		if contains_fasta:
			args = ['mv', gff_file_path, new_file_path]
			subprocess.call(args)
		else:
			# NOTE: We have to pipe output of cat call to the new_file_path
			# next we make a new 'gff' file that contains both the gff and fasta information

			args = ['cat', gff_file_path, cat_path, fasta_file['path']]
			catted_files = subprocess.check_output(args)
			f = open(new_file_path, 'w')
			f.write(catted_files.decode('utf-8'))
			f.close()

		path_to_ref_and_ID_pos_dict[new_file_path] = (ref, ID_to_pos, gffid_to_genid)

	return final_dir, path_to_ref_and_ID_pos_dict

def mapping_func(gff_id, gen_id):
	'''
	function to map gff IDs to genome IDs.
	'''
	if gff_id in gen_id:
		return True
	return False


def filter_gff(gff_file, genome_obj, overwrite=True):
	# if there is a duplicate ID's toss one of them out (if it is programmed frameshift toss one)
	# Here we throw out the second duplicate ID.


	# ideally we would do the feature mapping to the genes, but there aren't always 'gene' features in the gff files
	# so to get around this we base the gene position in the gff off the CDS position because they are most often
	# consecutively ordered and paired. For the organisms that Roary is supposed to service, there should only be
	# a one to one pairing of 'gene' to 'CDS', so this shouldn't be much of an issue

	features = genome_obj['features']

	# we want to make sure that the feature ID's lineup with the gff ID's
	# and we'll change the feature_pos argument to be the position of the ID
	# in the feature array.
	gen_ids = set([])
	ID_to_pos = {}
	for feat_pos, feature in enumerate(features):
		id_ = feature["id"]
		gen_ids.add(id_)
		ID_to_pos[id_] = feat_pos

	f = open(gff_file)
	output = []
	gff_ids = set([])
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

		# get ID of feat
		ID = l.split('ID=')[-1].split(';')[0]

		# determine type of feature
		feat_type = l.split()[2]

		ids_len = len(gff_ids)
		gff_ids.add(ID)
		if len(gff_ids) == ids_len:
			# we found a duplicate and its the second one. don't include it
			continue
		else:
			# make sure that the feature we get is a 'CDS' object
			if feat_type == 'CDS':
				output.append(l)

	f.close()

	# we want to make sure that the gen_ids and gff_ids are the same.
	# diff = gen_ids.symmetric_difference(gff_ids)

	# instead of the symmetric difference, I think its better to make sure the genome has all the gff ids
	# the reason we do this is to 
	overlap = gen_ids.intersection(gff_ids)
	gffid_to_genid = {}
	if len(overlap) == len(gen_ids) or len(overlap) == len(gff_ids):
		# all the ids overlap
		for o in overlap: 
			gffid_to_genid[o] = o
	else:
		# we should make a mapping from the gff ID's to the genome ID's
		# what is a good way to do that.... hmmm...
		for o in overlap: 
			gffid_to_genid[o] = o


		# get non overlapping ones
		diff = gen_ids - gff_ids
		gff_diff = gff_ids - gen_ids


		mapping = defaultdict(lambda:list)
		for gen_id in diff:
			# find gff_id that matches with the associated genome id
			# Im 85% certain that the gff_id will be a substring of
			# the genome id
			for gff_id in gff_diff:
				if mapping_func(gff_id, gen_id):					
					# this is where they overlap
					mapping[gen_id].append(gff_id)

		prob_map = {key:mapping[key] for key in mapping if len(mapping[key]) > 1}
		used_set = set([mapping[key][0] for key in mapping if len(mapping[key]) == 1])

		# iteratively find pairings for the key values.
		prob_map_copy = dict(prob_map)
		iters = 0
		while len(prob_map) > 1:
			# 1.) filtering step
			for key, v_list in prob_map.items():
				for index, item in enumerate(v_list):
					if item in used_set:
						v_list.pop(index)
				# 2.) add items that have 1 to 1 mapping to used list.
				#     and remove them from copied dict.
				if len(v_list) == 1:
					mapping[key] = v_list
					used_set.add(v_list[0])
					prob_map_copy.pop(key, None)
				else:
					prob_map_copy[key] = v_list
			prob_map = dict(prob_map_copy)
			iters+=1
			if iters > 20:
				raise ValueError("Could not resolve mapping of KBaseGenomes.Genome object IDs to GFF file IDs.")
		
		# now we should have a complete 1 to 1 mapping.
		for key in mapping:
			gffid_to_genid[mapping[key][0]] = key

	# now we can use gffid_to_genid for when we construct the pangenome object.
	if len(gffid_to_genid) != len(gen_ids) or len(gffid_to_genid) != len(gff_ids):
		gff_diffs = set(gff_ids) - set(gffid_to_genid.keys())
		gen_diffs = set(gen_ids) - set(gffid_to_genid.values())
		n = 100
		samp_keys = random.sample(list(gffid_to_genid), n)
		sampled = {k:gffid_to_genid[k] for k in samp_keys}
		raise ValueError("Genome object with id %s cannot match all of its IDs to an ID in its GFF File. "%genome_obj['id'],
			"Heres the whole dang mapping: ",sampled,"[[[[THE MAPPING FROM GFF TO GENOME ID ENDS HERE]]]]]",# "GFF ids not in mapping: ",
			#gff_diffs, "Genome IDs not in mapping: ", gen_diffs, 
			"length of genome features:", len(genome_obj['features']), 'length of mapping:', len(gffid_to_genid),
			"example of genome id:",genome_obj['features'][0]['id'], "number of gff ids:", len(gff_ids), "number of gen ids",len(gen_ids))


	if len(diff) != 0:
		# here is where we see they have different ids.
		raise ValueError("Genome object with id %s does not having matching ID's to gff file, output difference: "%genome_obj['id'], diff,
						 "gff ids length %i, genome ids length %i, difference length %i"%(len(gff_ids), len(gen_ids), len(diff)))


	assert(len(output) > 1), "Could not succesfully filter %f. It may be empty or contain no CDS information."%gff_file.split('/')[-1]

	if overwrite:
		# write output to file
		f = open(gff_file, 'w')
		for l in output:
			f.write(l)
		f.close()

	return gff_file, ID_to_pos, gffid_to_genid, contains_fasta
