# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
from installed_clients.KBaseReportClient import KBaseReport
from utils.roary_report import roary_report, generate_pangenome
from utils.roary_inputs import download_gffs
from utils.roary_proc import run_roary
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

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
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_slebrasRoary

        # Error Handling

        # ~conspicuously blank~

        # get input parameters
        ws_name = params.get('workspace_name')
        pangenome_name = params.get('pangenome_name')
        pangenome_id = "kb|"+pangenome_name

        # run meat of operations
        gff_folder_path, path_to_ref = download_gffs(params['ref'], self.callback_url, self.shared_folder)
        output_path = run_roary(self.shared_folder, gff_folder_path)
        sum_stats = output_path + '/summary_statistics.txt'
        gene_pres_abs = output_path + '/gene_presence_absence.csv'

        pangenome = generate_pangenome(gene_pres_abs, path_to_ref, pangenome_id, pangenome_name)

        output = roary_report(self.callback_url, ws_name, sum_stats, pangenome)

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
