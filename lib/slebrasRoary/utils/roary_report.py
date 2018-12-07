import os
import uuid
import pandas as pd
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WSLargeDataIOClient import WsLargeDataIO

# utils
from .roary_output import format_summary_statistics

def generate_pangenome(gene_pres_abs, path_to_ref_and_ID_pos_dict, pangenome_id, pangenome_name):
	'''
		params:
			gene_pres_abs               : file path to gene_presence_absence.csv output from Roary
			path_to_ref_and_ID_pos_dict : dictionary mapping gff file path to a tuple of (workspace ref, {gene ID :-> file position})
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

	for i, row in df.iterrows():
		OrthologFamily = {}
		OrthologFamily['id'] = row['Gene']
		OrthologFamily['type'] = None
		OrthologFamily['function'] = 'Roary'
		OrthologFamily['md5'] = None
		# figure out how to get this
		OrthologFamily['protein_translation'] = None

		orthologs = []

		for col in cols:
			for path in path_to_ref_and_ID_pos_dict:
				if col in path.split('/')[-1]:
					genome_ref, ID_to_pos = path_to_ref_and_ID_pos_dict[path]
					break
			gene_id = row[col]
			if not pd.isnull(gene_id):
				feature_pos = ID_to_pos[gene_id]
				orthologs.append((gene_id, feature_pos, genome_ref))

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


def roary_report(cb_url, scratch, workspace_name, sum_stats, pangenome_ref, conserved_vs_total_graph, unique_vs_new_graph):
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

	html_sum_stats = format_summary_statistics(sum_stats)
	file_dir = os.path.join(scratch, report_name)
	os.mkdir(file_dir)
	sum_stats_html_path = os.path.join(file_dir, 'sum_stats.html')
	with open(sum_stats_html_path, 'w') as f:
		f.write(html_sum_stats) 

	# get the .txt and .png files into shock
	shock_1 = dfu.file_to_shock({
		'file_path':sum_stats_html_path,
		'make_handle':0,
	})
	shock_2 = dfu.file_to_shock({
		'file_path':conserved_vs_total_graph,
		'make_handle':0,
	})
	shock_3 = dfu.file_to_shock({
		'file_path':unique_vs_new_graph,
		'make_handle':0,
	})

	html_link = {
		'shock_id': shock_1['shock_id'],
		'name':'sum_stats.html',
		'label':'Summary_Statistics',
		'description':'Roary Gene summary statistics'
	}
	photo_link_1 = {
		'shock_id': shock_2['shock_id'],
		'name':'conserved_vs_total_genes.png',
		'label':'html_files',
		'description':'Graph of conserved genes vs. total genes'
	}
	photo_link_2 = {
		'shock_id': shock_3['shock_id'],
		'name':'unique_vs_new_genes.png',
		'label':'html_files',
		'description':'Graph of unique genes vs new genes'
	}

	report_client = KBaseReport(cb_url)
	report = report_client.create_extended_report({
		'direct_html_link_index':0,
		'html_links':[html_link],
		'file_links':[photo_link_1, photo_link_2],
		'workspace_name': workspace_name,
		'report_object_name': report_name,
		'objects_created': [{
			'ref':pangenome_ref,
			'description':"Pangenome Object"
		}]
	})
	# now we send the outputs as a report back to shock

	return {
		'report_name':report['name'],
		'report_ref': report['ref']
	}