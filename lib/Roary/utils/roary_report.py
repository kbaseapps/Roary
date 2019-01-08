import os
import uuid
import pandas as pd
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WSLargeDataIOClient import WsLargeDataIO

# utils
from .roary_output import format_output_html

def get_col_name_from_path(path):
	return os.path.splitext(path.split('/')[-1])[0]

def generate_pangenome(gene_pres_abs, path_to_ref_and_ID_pos_dict, pangenome_id, pangenome_name):
	'''
		params:
			gene_pres_abs               : file path to gene_presence_absence.csv output from Roary
			path_to_ref_and_ID_pos_dict : dictionary mapping gff file path to a tuple of (workspace ref, {gene ID :-> file position}, {genome object id -> gff id})
			pangenome_id				: pangenome identifier
			pangenome_name				: pangenome display name
		Returns:
			Pangenome 					: KBaseGenomes.Pangenome like object (see type spec) https://narrative.kbase.us/#spec/type/KBaseGenomes.Pangenome
	'''
	Pangenome = {}

	Pangenome['genome_refs'] = [tup[0] for tup in path_to_ref_and_ID_pos_dict.values()]
	Pangenome['id'] = pangenome_id
	Pangenome['name'] =  pangenome_name
	Pangenome['type'] = None

	OrthologFamilyList = []

	consistent_cols  = ['Gene','Non-unique Gene name','Annotation','No. isolates','No. sequences',\
						'Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment',\
						'Accessory Order with Fragment','QC','Min group size nuc','Max group size nuc','Avg group size nuc']
	# load gene_pres_abs data and format it into Pangenome dictionary
	df = pd.read_csv(gene_pres_abs)
	# ignore errors to only drop columns that are in the dataframe

	cols = set(df.columns.values) - set(consistent_cols)

	col_to_ref = {}
	for col in cols:
		start_len = len(col_to_ref)
		for path in path_to_ref_and_ID_pos_dict:
			if col == get_col_name_from_path(path):
				col_to_ref[col] = path_to_ref_and_ID_pos_dict[path]
				break
		if len(col_to_ref) == start_len:
			raise ValueError("could not find file name match for " + col + " column")
	for i, row in df.iterrows():
		OrthologFamily = {}

		# put in standard arguments for an OrhtologFamily as found in Pangenome Spec file.
		OrthologFamily['id'] = row['Gene']
		OrthologFamily['type'] = None
		OrthologFamily['function'] = 'Roary'
		OrthologFamily['md5'] = None
		OrthologFamily['protein_translation'] = None

		orthologs = []

		for col in cols:
			row_gff_id = row[col]
			if not pd.isnull(row_gff_id):
				# find if the gene_id is in fact multiple gene_id's tab delimited
				if '\t' in row_gff_id:
					gff_ids = row_gff_id.split('\t')
				else:
					gff_ids = [row_gff_id]

				genome_ref, ID_to_pos, gffid_to_genid = col_to_ref[col]
				for gff_id in gff_ids:
					if gff_id not in gffid_to_genid:
						if '___' in gff_id:
							#chop off extra identifier if it exists
							gff_id = gff_id.split('___')[0]
							if gff_id not in gffid_to_genid:
								raise KeyError("gff ID %s not in file %s (pos 1)"%(gff_id, col))
						else:
							raise KeyError("gff ID %s not in file %s (pos 2)"%(gene_id, col))							
					gene_id = gffid_to_genid[gff_id]
					feature_pos = ID_to_pos[gene_id]
					orthologs.append([gene_id, feature_pos, genome_ref])
		OrthologFamily['orthologs'] = orthologs

		OrthologFamilyList.append(OrthologFamily)

	Pangenome['orthologs'] = OrthologFamilyList

	return Pangenome


def upload_pangenome(cb_url, scratch, Pangenome, workspace_name, pangenome_name):
	"""
	params:
		cb_url         : callback url
		scratch        : folder path to Pangenome object 
		pangenome      : KBaseGenomes.Pangenome like object
		workspace_name : workspace name
		pangenome_name : Pangenome display name
	Returns:
		pangenome_ref: Pangenome workspace reference
		pangenome_info: info on pangenome object
	"""
	dfu = DataFileUtil(cb_url)
	meta = {}
	hidden = 0

	# dump pangenome to scratch for upload
	# data_path = os.path.join(scratch, pangenome_name + '.json')
	# json.dump(pangenome, open(data_path, 'w'))

	if isinstance(workspace_name, int) or workspace_name.isdigit():
		workspace_id = workspace_name
	else:
		workspace_id = dfu.ws_name_to_id(workspace_name)

	save_params = {
		'id': workspace_id,
		'objects': [{
			'type': 'KBaseGenomes.Pangenome',
			'data': Pangenome,
			'name': pangenome_name,
			'meta': meta,
			'hidden': hidden
		}]
	}

	info = dfu.save_objects(save_params)[0]

	ref = "{}/{}/{}".format(info[6], info[0], info[4])
	print("Pangenome saved to {}".format(ref))

	return {
		'pangenome_ref': ref,
		'pangenome_info': info
	}


def roary_report(cb_url, scratch, workspace_name, sum_stats, gene_pres_abs, pangenome_ref, conserved_vs_total_graph, unique_vs_new_graph):
	"""
	params:
		cb_url         : callback url
		workspace_name : name of the workspace
		sum_stats      : summary_statistics.txt file output from Roary
		pangenome_ref  : reference to the pangenome object, or None
	Returns:
		report_name : name of report object  
		report_ref  : reference to report object in workspace
	"""
	report_name = 'Roary_report_'+str(uuid.uuid4())	
	dfu = DataFileUtil(cb_url)

	# Convert output files to HTML
	html_output = format_output_html(sum_stats, gene_pres_abs)

	file_dir = os.path.join(scratch, report_name)
	os.mkdir(file_dir)
	html_path = os.path.join(file_dir, 'output.html')
	with open(html_path, 'w') as f:
		f.write(html_output) 


	html_link = {
		'path': file_dir,
		'name':'output.html',
		# 'label':'Summary_Statistics',
		'description':'Roary Gene Statistics html report'
	}

	csv_link = {
		'path':gene_pres_abs,
		'name':'gene_presence_absence.csv',
		'description':"Data Table of Gene Presence, Gene Absence and other information"
	}

	photo_link_1 = {
		'path': conserved_vs_total_graph,
		'name':'conserved_vs_total_genes.png',
		# 'label':'Conserved_vs_total_genes_graph',
		'description':'Graph of conserved genes vs. total genes'
	}
	photo_link_2 = {
		'path': unique_vs_new_graph,
		'name':'unique_vs_new_genes.png',
		# 'label':'unique_vs_new_genes_graph',
		'description':'Graph of unique genes vs new genes'
	}

	report_client = KBaseReport(cb_url)
	if pangenome_ref is not None:
		report = report_client.create_extended_report({
			'direct_html_link_index':0,
			'html_links':[html_link],
			'file_links':[csv_link, photo_link_1, photo_link_2],
			'workspace_name': workspace_name,
			'report_object_name': report_name,
			'objects_created': [{
				'ref':pangenome_ref,
				'description':"Pangenome Object"
			}]
		})
	else:
		report = report_client.create_extended_report({
			'direct_html_link_index':0,
			'html_links':[html_link],
			'file_links':[csv_link, photo_link_1, photo_link_2],
			'workspace_name': workspace_name,
			'report_object_name': report_name
		})
	return {
		'report_name':report['name'],
		'report_ref': report['ref']
	}