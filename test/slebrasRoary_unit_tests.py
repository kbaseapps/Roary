import os
import sys
import pandas as pd
import time
import tempfile
import unittest

from installed_clients.DataFileUtilClient import DataFileUtil

from configparser import ConfigParser

currdir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(currdir, '../lib/slebrasRoary/')))

# import functions to test
from utils.roary_report import generate_pangenome
from utils.roary_proc import run_roary

class TestRoary(unittest.TestCase):

	def test_roary_proc_output(self):
		gff_folder = os.path.join(currdir, 'data/gffs')
		config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
		config = ConfigParser()
		config.read(config_file)
		cfg = {}
		for nameval in config.items('slebrasRoary'):
			cfg[nameval[0]] = nameval[1]
		scratch = cfg['scratch']

		# make sure we have the directory
		self.assertTrue(os.path.isdir(gff_folder))
		files = ['000008285.gff','000021185.gff','000026705.gff',\
				 '000168635.gff','000168815.gff','000196035.gff']

		path_to_ref = {}
		# make sure that all files exist
		j = 0
		for f_name in files:
			f_path = os.path.join(gff_folder,f_name)
			print(f_path)
			self.assertTrue(os.path.isfile(f_path))

			path_to_ref[f_path] = '0000/0/' + str(j)
			j+=1

		output_path = run_roary(scratch, gff_folder, {'blast_p_percentage':98,'max_num_clusters':50000})
		sum_stats = output_path + '/summary_statistics.txt'
		gene_pres_abs = output_path + '/gene_presence_absence.csv'

		pangenome_id = 'kb|test_pangenome'
		pangenome_name = 'Test Pangenome'

		pangenome = generate_pangenome(gene_pres_abs, path_to_ref, pangenome_id, pangenome_name)

		def check_sum_stats(sum_stats):
			ss_file = open(sum_stats)
			compare_nums = [2031, 0, 2505, 0, 4536]
			lines = [l for l in ss_file]
			for i in range(len(lines)):
				line_arr = lines[i].split('\t')
				num = int(line_arr[-1].rstrip())
				self.assertEqual(num, compare_nums[i])
		
		# (row index, column, contents)
		check = [
			('group_10',
				['000008285','000021185','000026705','000168635','000168635','000196035'],
				["GCA_000008285_00337","GCA_000021185_02349","GCA_000026705_00325","GCA_000168635_00315","GCA_000168815_01298","GCA_000196035_00317"]
			),
			('yutF',
				['000008285','000021185','000026705','000168635','000168635','000196035'],
				["GCA_000008285_02418","GCA_000021185_00193","GCA_000026705_02417","GCA_000168635_02482","GCA_000168815_00643","GCA_000196035_02488"]
			)
		]			

		def check_gene_pres_abs(gene_pres_abs):
			# check a few lines
			df_gene = pd.read_csv(gene_pres_abs)
			for gene, cols, vals in check:
				df_row = df_gene[df_gene.Gene == gene][cols]
				for i in range(len(cols)):
					self.assertEqual(df_row[cols[i]].iloc(0), vals[i])

		def check_pangenome(pangenome):
			# check the pangenome output object
			self.assertEqual(len(set(pangenome['genome_refs']).symmetric_difference(set(path_to_ref.values()))), 0)
			self.assertEqual(pangenome['id'], pangenome_id)
			self.assertEqual(pangenome['name'], pangenome_name)
			self.assertEqual(pangenome['type'], None)

			for tup in check:
				for og_fam in pangenome['orthologs']:
					if tup[0] == og_fam['id']:
						gene_ids = set([og[0] for og in og_fam['orthologs']])
						self.assertEqual(len(gene_ids.symmetric_difference(set(tup[2]))),0)
						break

			for og_fam in pangenome['orthologs']:
				# filter 
				og_fam['id']

		# call functions to check the outputs
		check_sum_stats(sum_stats)
		check_gene_pres_abs(gene_pres_abs)
		check_pangenome(pangenome)

if __name__ == '__main__':
	unittest.main()