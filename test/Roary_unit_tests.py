import os
import sys
import pandas as pd
import numpy as np
import unittest

from configparser import ConfigParser

currdir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(currdir, '../lib/Roary/')))

# import functions to test
from utils.roary_report import generate_pangenome
from utils.roary_proc import run_roary
from utils.roary_inputs import filter_gff


class TestRoary(unittest.TestCase):

    def check_sum_stats(self, sum_stats, compare_nums):
        ''' Function to check the output of the summary_statistics.txt file that Roary produces.
        Reads file line by line and compares to exact values or ranges of values.

        params:
            sum_stats - file path to summary_statistics.txt output file
            comare_nums - list of values to compare summary_statistics values to.
        '''
        ss_file = open(sum_stats)
        lines = [l for l in ss_file]
        for i in range(len(lines)):
            line_arr = lines[i].split('\t')
            num = int(line_arr[-1].rstrip())
            # roary seems to be non-determinisitic so we check for a range of values
            if isinstance(compare_nums[i], tuple):
                self.assertTrue(compare_nums[i][0] <= num <= compare_nums[i][1])
            else:
                self.assertEqual(num, compare_nums[i])
    
    def check_gene_pres_abs(self, gene_pres_abs, check):
        ''' Function to validate the gene_presenece_absence.csv output from Roary.
        Finds the values for the associated genes in 'check' variable and insures that the
        values in the output columns match.

        params: 
            gene_pres_abs - file path to gene_presence_absence.csv
            check - list of tuples in order (gene, [columns], [values])
        '''
        df_gene = pd.read_csv(gene_pres_abs)
        for gene, cols, vals in check:
            df_row = df_gene[df_gene.Gene == gene][cols]
            for i in range(len(cols)):
                if pd.isna(df_row[cols[i]].iloc[0]):
                    self.assertTrue(pd.isna(vals[i]))
                else:
                    self.assertEqual(df_row[cols[i]].iloc[0], vals[i])

    def check_pangenome(self, pangenome, check, path_to_ref, pangenome_id, pangenome_name):
        '''Function validate pangenome output generated from generate_pangenome function.
        Checks for same values that are checked for in check_gene_pres_abs().

        params: 
            pangenome - KBaseGenomes.Pangenome type object to validate
            check -
            path_to_ref - dictionary of {input file path -> (genome reference, {gene ID -> file position})}
        '''
        # check the pangenome output object
        genome_ref_diff = len(set(pangenome['genome_refs']).symmetric_difference(set([tup[0] for tup in path_to_ref.values()])))
        self.assertEqual(genome_ref_diff, 0)
        self.assertEqual(pangenome['id'], pangenome_id)
        self.assertEqual(pangenome['name'], pangenome_name)
        self.assertEqual(pangenome['type'], None)

        for tup in check:
            for og_fam in pangenome['orthologs']:
                if tup[0] == og_fam['id']:
                    gene_ids = set([og[0] for og in og_fam['orthologs']])
                    diff = gene_ids.symmetric_difference(set(tup[2]))
                    # nans are interpreted differently, in the case that the difference is 'nan'
                    # we filter it out in order to not produce an unnecessary error.
                    if len(diff) == 1 and pd.isna(list(diff)[0]):
                        diff = set([])
                    self.assertEqual(len(diff),0)
                    break

    def _dummy_gen_obj(self, gff_file):
        features = []
        with open(gff_file) as f:
            for line in f:
                if '##FASTA' in line:
                    break
                if line[0:2] == "##":
                    continue
    
                ID = line.split('ID=')[-1].split(';')[0]
                features.append({"id":ID})
        return {'features':features, 'id': "dummy_genome"}

    def roary_proc_check(self, files, gff_folder, scratch, sum_stats_check, pres_abs_check):
        '''Run Roary, generate outputs and call necessary tests to validate them.

        params:
            files - list of gff files to run roary on.
            gff_folder - folder where gff files are stored
            scratch - scratch folder
            sum_stats_check - list of values to check summary_statistics.txt against
            pres_abs_check - list of values to check gene_presence_absence.csv against
        '''
        self.assertTrue(os.path.isdir(gff_folder))
        path_to_ref = {}
        # make sure that all files exist
        j = 0
        for f_name in files:
            f_path = os.path.join(gff_folder,f_name)
            self.assertTrue(os.path.isfile(f_path))
            f_path, ID_to_pos, gffid_to_genid, _, _ = filter_gff(f_path, self._dummy_gen_obj(f_path))

            path_to_ref[f_path] = ('0000/0/' + str(j), ID_to_pos, gffid_to_genid)
            j+=1

        params = {'blast_p_percentage':95,'max_num_clusters':50000,\
                  'percent_genes_for_core':99}

        output_path = run_roary(scratch, gff_folder, params)
        sum_stats = os.path.join(output_path, 'summary_statistics.txt')
        gene_pres_abs = os.path.join(output_path, 'gene_presence_absence.csv')

        pangenome_id = 'kb|test_pangenome'
        pangenome_name = 'Test Pangenome'

        pangenome = generate_pangenome(gene_pres_abs, path_to_ref, pangenome_id, pangenome_name)

        self.check_sum_stats(sum_stats, sum_stats_check)
        self.check_gene_pres_abs(gene_pres_abs, pres_abs_check)
        self.check_pangenome(pangenome, pres_abs_check, path_to_ref, pangenome_id, pangenome_name)

    def roary_proc_output_1(self, scratch):
        '''Define input files and validation values for first test

        params:
            scratch - scratch folder path
        '''
        gff_folder = os.path.join(currdir, 'data/gffs')
        # test the following files
        files = ['000008285.gff','000021185.gff','000026705.gff',\
                 '000168635.gff','000168815.gff','000196035.gff']

        # there seems to be a range in values that the output can take.
        # This is reflected as a tuple as (min, max) of range.
        sum_stats_check = [(2011,2013), 0, (2449,2451), 0, 4462]
        pres_abs_check = [
            ('group_10',
                ['000008285','000021185','000026705','000168635','000168815','000196035'],
                ["GCA_000008285_00337","GCA_000021185_02349","GCA_000026705_00325",\
                 "GCA_000168635_00315","GCA_000168815_01298","GCA_000196035_00317"]
            ),
            ('yutF',
                ['000008285','000021185','000026705','000168635','000168815','000196035'],
                ["GCA_000008285_02418","GCA_000021185_00193","GCA_000026705_02417",\
                 "GCA_000168635_02482","GCA_000168815_00643","GCA_000196035_02488"]
            )
        ]
        self.roary_proc_check(files, gff_folder, os.path.join(scratch,'test1'),\
                              sum_stats_check, pres_abs_check)

    def roary_proc_output_2(self, scratch):
        '''Define input files and validation values for second test

        params:
            scratch - scratch folder path
        '''
        gff_folder = os.path.join(currdir,'data/gffs_2')
        # test the following files
        files = ['abietaniphila.gff','abyssi.gff','acidophila.gff','aeruginosa.gff']
        sum_stats_check = [0, 0, 23195, 0, 23195]
        pres_abs_check = [
            ('ihfA',
                ['abietaniphila','abyssi','acidophila','aeruginosa'],
                ['BLS87_RS14525_CDS_1','CNQ84_RS14035_CDS_1',np.nan,'NS85_RS32165_CDS_1']
            ),
            ('rplQ',
                ['abietaniphila','abyssi','acidophila','aeruginosa'],
                [np.nan,'CNQ84_RS18245_CDS_1',np.nan,'NS85_RS06110_CDS_1']
            )
        ]
        self.roary_proc_check(files, gff_folder, os.path.join(scratch, 'test2'),\
                              sum_stats_check, pres_abs_check)

    def test_roary(self):
        '''
        Entrance function for tests. Defines config and scratch folder.
        '''
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        config = ConfigParser()
        config.read(config_file)
        cfg = {}
        for nameval in config.items('Roary'):
            cfg[nameval[0]] = nameval[1]

        scratch = cfg['scratch']
        self.roary_proc_output_1(scratch)
        self.roary_proc_output_2(scratch)

if __name__ == '__main__':
    unittest.main
