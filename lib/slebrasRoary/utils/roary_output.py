import os
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
		loader=PackageLoader('slebrasRoary','utils/templates'),
		autoescape=select_autoescape(['html'])) 

def format_summary_statistics(sum_stats):
	"""
	format summary statistics to html
	"""
	f = open(sum_stats)
	names = ['core_genes','soft_core_genes','shell_genes','cloud_genes','total_genes']
	results = {}
	i = 0
	for line in f:
		results[names[i]] = line.split()[-1]
		i+=1
	return create_html_tables('sum_stats.html', results)

def create_html_tables(html_file, formatted_results, formatted_headers=None):
	"""
	Render template in 'templates' folder with given inputs
	"""
	template = env.get_template(html_file)
	if formatted_headers:
		return template.render(results=formatted_results, headers=formatted_headers)
	else:
		return template.render(results=formatted_results)

def format_gene_presence_absence(gene_pres_abs):
	df = pd.read_csv(gene_pres_abs)
	df.columns = ['_'.join(col.replace('.','').split()) for col in df.columns.values]
	headers = df.columns.values
	d = df.to_dict('index')
	return create_html_tables('gene_pres_abs.html', list(d.values()), headers)