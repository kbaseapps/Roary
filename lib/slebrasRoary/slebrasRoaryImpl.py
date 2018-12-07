# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from installed_clients.KBaseReportClient import KBaseReport
from .utils.roary_report import roary_report, generate_pangenome, upload_pangenome
from .utils.roary_inputs import download_gffs
from .utils.roary_proc import run_roary
#END_HEADER


class slebrasRoary:
    '''
    Module Name:
    slebrasRoary

    Module Description:
    A KBase module: slebrasRoary
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/slebras/Roary.git"
    GIT_COMMIT_HASH = "a7d3574191d1e591ddc53ab8485eaac34c0ddd80"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def run_slebrasRoary(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "RoaryParams" (roary input) ->
           structure: parameter "workspace_name" of String, parameter "ref"
           of String, parameter "pangenome_name" of String, parameter
           "blast_p_percentage" of Long, parameter "max_num_clusters" of Long
        :returns: instance of type "RoaryResults" (roary output) ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_slebrasRoary

        # input validation
        # ~conspicuously blank~

        # get input parameters
        workspace_name = params.get('workspace_name')
        pangenome_name = params.get('pangenome_name')
        genome_set_ref = params.get('ref')

        if pangenome_name is not None:
            pangenome_id = "kb|"+pangenome_name
        else:
            pangnome_id = None

        # run meat of operations
        gff_folder_path, path_to_ref_and_ID_pos_dict = download_gffs(self.callback_url, self.shared_folder, genome_set_ref)
        output_path = run_roary(self.shared_folder, gff_folder_path, params)
        sum_stats = os.path.join(output_path,'summary_statistics.txt')
        gene_pres_abs = os.path.join(output_path,'gene_presence_absence.csv')
        conserved_vs_total_graph = os.path.join(output_path,'conserved_vs_total_genes.png')
        unique_vs_new_graph = os.path.join(output_path, 'unique_vs_new_genes.png')

        if pangenome_name:
            pangenome = generate_pangenome(gene_pres_abs, path_to_ref_and_ID_pos_dict, pangenome_id, pangenome_name)
            pangenome_obj = upload_pangenome(self.callback_url, self.shared_folder, pangenome, workspace_name, pangenome_name)
            output = roary_report(self.callback_url, self.shared_folder, workspace_name, sum_stats, pangenome_obj['pangenome_ref'],\
                                    conserved_vs_total_graph, unique_vs_new_graph)
        else:
            output = roary_report(self.callback_url, workspace_name, sum_stats, None, \
                                conserved_vs_total_graph, unique_vs_new_graph)

        #END run_slebrasRoary

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_slebrasRoary return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
