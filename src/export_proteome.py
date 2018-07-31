#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__date__ = '19 Jun 2018'


import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
import re
import time
from operator import itemgetter
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
from Bio.Seq import Seq
from Bio import SeqIO
from _version import __CHOCOPhlAn_version__

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

class export_proteome:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')

        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.uniref100 = {}
        self.uniref90 = {}
        self.uniref50 = {}

        self.config = config
        if config['verbose']:
            utils.info('Finished.\n')

    def export_core_proteins(self, tax_id, cluster):
        tax_level = self.taxontree.taxid_n[tax_id].rank
        try:
            with open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], tax_level, cluster, tax_id),'rb') as p_handle:
                panproteome = pickle.load(p_handle)
        except FileNotFoundError as ex:
            utils.error('Panproteome {} not found'.format(tax_id))
            raise

        uniref = self.__getattribute__('uniref{}'.format(cluster))
        if not len(uniref):
            all_uref_chunks = filter(re.compile('uniref{}_[0-9]{{1,}}.pkl'.format(cluster)).match, os.listdir('{}/pickled'.format(config['download_base_dir'])))
            for i in all_uref_chunks:
                uniref_chunk = pickle.load(open('{}/pickled/{}'.format(self.config['download_base_dir'], i),'rb'))
                uniref.update({k : v[6] for k,v in uniref_chunk.items()})

        cores = Panproteome.find_core_genes(panproteome)
        cores_seqs = [SeqRecord(id = core, 
                                name = core,
                                description = '',
                                seq = Seq(uniref['UniRef{}_{}'.format(cluster, core)].replace('\n',''), ProteinAlphabet())
                            ) for core in cores]

        os.makedirs('{}{}'.format(self.config['export_dir'], self.config['exportpath_core_proteins']), exist_ok=True)
        SeqIO.write(cores_seqs, '{}{}{}.faa'.format(self.config['export_dir'], self.config['exportpath_core_proteins'], tax_id), 'fasta')