import os
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
        loader=PackageLoader('Roary','utils/templates'),
        autoescape=select_autoescape(['html'])) 

def format_output_html(sum_stats, gene_pres_abs, conserved_image_path, unique_image_path):
    """
    format summary statistics and gene_presence absence files to html
    """
    # format txt file
    with open(sum_stats) as f:
        names = ['core_genes','soft_core_genes','shell_genes','cloud_genes','total_genes']
        results = {}
        for i, line in enumerate(f):
            # Format results dictionary to use names for Jinja2 inputs.
            results[names[i]] = line.split()[-1]

    return create_html_tables('output_template.html', [results], conserved_image_path, unique_image_path)


def create_html_tables(html_file, formatted_results, conserved_image_path, unique_image_path):
    """
    Render template in 'templates' folder with given inputs
    """
    template = env.get_template(html_file)
    return template.render(
        results=formatted_results,
        conserved_image=conserved_image_path,
        unique_image=unique_image_path
    )
