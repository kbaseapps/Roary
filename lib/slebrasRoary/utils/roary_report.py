import os
import uuid
import pandas as pd
import json

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport



def generate_pangenome(gene_pres_abs, path_to_ref, pangenome_id, pangenome_name):
	'''
		params:
			gene_pres_abs:
			path_to_ref: 
			pangenome_id:
			pangenome_name:
		Returns:
			Pangenome : json.dump
	'''

	Pangenome = {}

	Pangenome['genome_refs'] = path_to_ref.values()
	Pangenome['id'] = pangenome_id
	Pangenome['name'] =  pangenome_name
	Pangenome['type'] = None

	OrthologFamilyList = []

	consistent_cols  = ["Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences",\
						"Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment",\
						"Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"]
	# load gene_pres_abs data and format it into Pangenome dictionary
	df = pd.read_csv(gene_pres_abs)
	cols = df.drop(consistent_cols).columns.values
	for i, row in df.iterrows():
		OrthologFamily = {}
		OrthologFamily['id'] = row['Gene']
		OrthologFamily['type'] = None
		OrthologFamily['function'] = 'Roary'
		OrthologFamily['md5'] = None
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

	return json.dumps(Pangenome)


def roary_report(cb_url, workspace_name, sum_stats, pangenome):

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
		'html_links',[results_link]
		'workspace_name': workspace_name,
		'report_object_name': report_name
	})
	# now we send the outputs as a report back to shock

	return {
		'report_name':report['name'],
		'report_ref': report['ref']
	}