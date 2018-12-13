import os
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

# TODO replace all tab characters in this file with 4 spaces (https://www.python.org/dev/peps/pep-0008/#tabs-or-spaces)

env = Environment(
		loader=PackageLoader('slebrasRoary','utils/templates'),
		autoescape=select_autoescape(['html'])) 

def format_output_html(sum_stats, gene_pres_abs):
	"""
	format summary statistics and gene_presence absence files to html
	"""
	# format txt file
        # TODO use context manager
	f = open(sum_stats)
	names = ['core_genes','soft_core_genes','shell_genes','cloud_genes','total_genes']
	results = {}
	i = 0
        # TODO for (idx, line) in enumerate(f):
	for line in f:
            # TODO more comments -- why are you doing this? what does roary's ouput look like?
		results[names[i]] = line.split()[-1]
		i+=1
	f.close()

	# limit number of results to n
	n = 20

	# format csv file
	df = pd.read_csv(gene_pres_abs)
        # TODO comment on why you are doing the below. is it for jinja2?
	df.columns = ['_'.join(col.replace('.','').split()) for col in df.columns.values]
	headers = df.columns.values
	d = df.to_dict('index')

	return create_html_tables('output_template.html', [results], list(d.values())[:n], headers)


def create_html_tables(html_file, formatted_results, formatted_rows, formatted_headers):
	"""
	Render template in 'templates' folder with given inputs
	"""
	template = env.get_template(html_file)
	return template.render(results=formatted_results)#, headers=formatted_headers, rows=formatted_rows)
