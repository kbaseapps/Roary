import os
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
		loader=PackageLoader('slebrasRoary','utils/templates'),
		autoescape=select_autoescape(['html'])) 

def format_summary_statistics(sum_stats):
	"""
	implement this to show file as html on results page
	"""
	f = open(sum_stats)
	names = ['core_genes','soft_core_genes','shell_genes','cloud_genes','total_genes']
	results = {}
	i = 0
	for line in f:
		results[names[i]] = line.split()[-1]
		i+=1
	return create_html_tables(results)

def create_html_tables(formatted_results):
	# headers?
	template = env.get_template('sum_stats.html')
	return template.render(results=formatted_results)