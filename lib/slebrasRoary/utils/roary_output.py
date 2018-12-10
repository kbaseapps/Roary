import os
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
		loader=PackageLoader('slebrasRoary','utils/templates'),
		autoescape=select_autoescape(['html'])) 

def format_output_html(sum_stats, gene_pres_abs):
	"""
	format summary statistics and gene_presence absence files to html
	"""
	# format txt file
	f = open(sum_stats)
	names = ['core_genes','soft_core_genes','shell_genes','cloud_genes','total_genes']
	results = {}
	i = 0
	for line in f:
		results[names[i]] = line.split()[-1]
		i+=1
	f.close()

	# format csv file
	df = pd.read_csv(gene_pres_abs)
	df.columns = ['_'.join(col.replace('.','').split()) for col in df.columns.values]
	headers = df.columns.values
	d = df.to_dict('index')

	return create_html_tables('output_stats.html', [results], list(d.values()), headers)



def create_html_tables(html_file, formatted_results, formatted_rows, formatted_headers):
	"""
	Render template in 'templates' folder with given inputs
	"""
	template = env.get_template(html_file)
	return template.render(results=formatted_results, headers=formatted_headers, rows=formatted_rows)
