import os
import sys

def filter_gff(gff_file):
	# if there is a duplicate ID's toss one of them out (if it is programmed frameshift toss one)
	# either the smaller or the second one.
	# for feature pos, use the order of the gff 'gene' items. create a mapping from ID to order for feature pos
	f = open(gff_file)
	output = []
	IDs = set([])
	# definitely a better way of doing this

	# ideally we would do the feature mapping o th genes, but there are nt always 'gene' features in the gff files
	# so to get around this we base the gene position in the gff off the CDS position because they are most often
	# consecutively ordered and paired. For the organisms that Roary is supposed to service, there should only be
	# a one to one pairing of 'gene' to 'CDS' 

	gene_to_pos = {}
	cds_to_gene = {}

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
	for l in output: f.write(l)
	f.close()
	# get mapping from CDS to gene_position

	return gff_file, cds_to_pos, contains_fasta

files = ['000008285.gff','000026705.gff','000168815.gff', \
		 '000021185.gff','000168635.gff','000196035.gff']
i = 0
for f in files:
	filter_gff(f)
	i+=1
	print('finished file #; ', i, f)