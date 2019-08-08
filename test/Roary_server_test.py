# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from configparser import ConfigParser

from Roary.RoaryImpl import Roary
from Roary.RoaryServer import MethodContext
from Roary.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

class RoaryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('Roary'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'Roary',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = Roary(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_Roary_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def getGenomeSetRef(self):
        gsr = '22385/60/1'
        return gsr

    def getGenomeRefs(self):
        return ["22385/76/1", "22385/80/1"]

    def getMixedRefs(self):
        return ["22385/82/1", "22385/58/1", "22385/80/1"]

    def initialize(self):
        shutil.rmtree(self.scratch)
        if not os.path.isdir(self.scratch):
            os.mkdir(self.scratch)

    def test_mixed_refs(self):
        self.initialize()
        input_params = {
            'workspace_name': self.getWsName(),
            'ref': self.getMixedRefs(),
            'pangenome_name': "mixed_refs_server_test_pangenome",
            'blast_p_percentage':95,
            'max_num_clusters':50000,
            'percent_genes_for_core':99
        }
        ret = self.getImpl().run_Roary(self.getContext(), input_params)
        print('returned val:', ret)

    def test_genomes(self):
        self.initialize()
        input_params = {
            'workspace_name': self.getWsName(),
            'ref': self.getGenomeRefs(),
            'pangenome_name': "genomes_server_test_pangenome",
            'blast_p_percentage':95,
            'max_num_clusters':50000,
            'percent_genes_for_core':99
        }
        ret = self.getImpl().run_Roary(self.getContext(), input_params)
        print('returned value',ret)

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        self.initialize()
        input_params = {
            'workspace_name': self.getWsName(),
            'ref': [self.getGenomeSetRef()],
            'pangenome_name': "genomeset_server_test_pangenome",
            'blast_p_percentage':95,
            'max_num_clusters':50000,
            'percent_genes_for_core':99    
        }

        ret = self.getImpl().run_Roary(self.getContext(),input_params)
        print('returned value',ret)
        # self.assertTrue('report_name' in ret)
        # self.assertTrue('report_ref' in ret)  
