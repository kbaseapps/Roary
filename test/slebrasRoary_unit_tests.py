import os
import sys
import pandas as pd
import time
import tempfile
import unittest

currdir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(currdir, '../lib/slebrasRoary/')))

# import functions to test
from utils.roary_report import generate_pangenome
from utils.roary_proc import run_roary

class TestRoary(unittest.TestCase):

	def test_roary_proc_output(self):
		gff_folder = os.path.join(currdir, 'data/gffs')

		# make sure we have the directory
		self.asserTrue(os.path.isdir(gff_folder))
		files = ['000008285.gff','000021185.gff','000026705.gff',\
				 '000168635.gff','000168815.gff','000196035.gff']

		path_to_ref = {}
		# make sure that all files exist
		for f_name in files:
			f_path = os.path.join(currdir,f_name)
			self.asserTrue(os.path.isfile(f_path))

			path_to_ref[f_path] = '0000/0/0'

		temp_dir = tempfile.mkdtemp()
		output_path = run_roary(temp_dir, gff_folder)
		sum_stats, gene_pres_abs = roary_outputs(output_path)
		pangenome = generate_pangenome(gene_pres_abs)

		def check_sum_stats(sum_stats):
			ss_file = open(sum_stats)
			compare_nums = [2031, 0, 2505, 0, 4536]
			lines = [l for l in ss_file]
			for i in range(len(lines)):
				line_arr = lines[i].split('\t')
				num = int(line_arr[-1].rstrip())
				self.assertEqual(num, compare_nums[i])

		def check_gene_pres_abs(gene_pres_abs):
			# check a few lines
			# (row index, column, contents)
			check   = [ 
						('group_10',
							['000008285','000021185','000026705','000168635','000168635','000196035'],
							["GCA_000008285_00337","GCA_000021185_02349","GCA_000026705_00325","GCA_000168635_00315","GCA_000168815_01298","GCA_000196035_00317"]
						),
						('yutF',
							['000008285','000021185','000026705','000168635','000168635','000196035'],
							["GCA_000008285_02418","GCA_000021185_00193","GCA_000026705_02417","GCA_000168635_02482","GCA_000168815_00643","GCA_000196035_02488"]
						)
			          ]
			df_gene = pd.read_csv(gene_pres_abs)
			for gene, cols, vals in check:
				df_row = df_gene[df_gene.Gene == gene][cols]
				for i in range(len(cols)):
					self.assertEqual(df_row[cols[i]].iloc(0), vals[i])

		def check_pangenome(pangenome):
			# check the pangenome output object
			pass


		# call functions to check the outputs
		check_sum_stats(sum_stats)
		check_gene_pres_abs(gene_pres_abs)
		check_pangenome(pangenome)

if __name__ == '__main__':
	unittest.main()