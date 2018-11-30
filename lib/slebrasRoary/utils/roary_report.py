import os
import uuid
import pandas as pd
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WSLargeDataIOClient import WsLargeDataIO

def generate_pangenome(gene_pres_abs, path_to_ref, pangenome_id, pangenome_name):
	'''
		params:
			gene_pres_abs:
			path_to_ref: 
			pangenome_id:
			pangenome_name:
		Returns:
			Pangenome : KBaseGenomes.Pangenome like object (see type spec)
	'''
	Pangenome = {}

	Pangenome['genome_refs'] = path_to_ref.values()
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
			for path in path_to_ref:
				if col in path.split('/')[-1]:
					genome_ref = path_to_ref[path]
					break
			gene_id = row[col]
			feature_pos = float(row[col].split('_')[-1])
			orthologs.append((gene_id, feature_pos, genome_ref))


		OrthologFamily['orthologs'] = orthologs

		OrthologFamilyList.append(OrthologFamily)

	Pangenome['orthologs'] = OrthologFamilyList

	return Pangenome


def upload_pangenome(cb_url, scratch, pangenome, workspace_name, pangenome_name):
	"""
	params:
		cb_url:
		scartch: 
		pangenome:
		workspace_name:
		pangenome_name:
	Returns:
		pangenome_ref: 
		pangenome_info:
	"""
	ws_large_data = WsLargeDataIO(cb_url)
	dfu = DataFileUtil(cb_url)
	meta = {}
	if 'hidden' in params and str(params['hidden']).lower() in ('yes', 'true', 't', '1'):
		hidden = 1
	else:
		hidden = 0

	# dump pangenome to scratch for upload
	data_path = os.path.join(scratch, pangenome_name + '.json')
	json.dump(pangenome, open(data_path, 'w'))

	if isinstance(workspace_name, int) or workspace.isdigit():
		workspace_id = workspace_name
	else:
		workspace_id = dfu.ws_name_to_id(workspace_name)

	save_params = {
		'id': workspace_id,
		'objects': [{
			'type': 'KBaseGenomes.Pangenome',
			'data_json_file': data_path,
			'name': pangenome_name,
			'meta': meta,
			'hidden': hidden
		}]
	}

	info = ws_large_data.save_objects(save_params)[0]

	ref = "{}/{}/{}".format(info[6], info[0], info[4])
	print("Pangenome saved to {}".format(ref))

	return {
		'pangenome_ref': ref,
		'pangenome_info': info
	}


def roary_report(cb_url, workspace_name, sum_stats, pangenome_ref):
	"""
	params:
		cb_url         : callback url
		workspace_name : name of the workspace
		sum_stats      :  
		pangenome_ref  : reference to the pangenome object, or None
	Returns:
		report_name : name of report object  
		report_ref  : reference to report object in workspace
	"""
	report_name = 'Roary_report_'+str(uuid.uuid4())

	dfu = DataFileUtil(cb_url)

	# get the .txt file into shock
	shock = dfu.file_to_shock({
		'file_path':sum_stats,
		'make_handle':0,
	})

	results_link = {
		'shock_id': shock['shock_id'],
		'name':'summary_statistics.txt',
		'label':'html_files',
		'description':'Roary TXT file summary statistics'
	}

	report_client = KBaseReport(cb_url)
	report = report_client.create_extended_report({
		'direct_html_link_index':0,
		'html_links':[results_link],
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