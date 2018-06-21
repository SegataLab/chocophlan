#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

from _version import __version__
__date__ = '11 Apr 2018'


import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
import glob
import time
from operator import itemgetter

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

class export_genomes:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')

        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.config = config
        if config['verbose']:
            utils.info('Finished.\n')

    def elaborate(self):
        with open('export/genomes_list.txt',"wt") as wout:
            [wout.write('{}\t{}\t{}\n'.format(self.proteomes[p]['tax_id'], self.taxontree.print_full_taxonomy(self.proteomes[p]['tax_id']), dict(self.proteomes[p]['ncbi_ids']).get('GCSetAcc',''))) for p in self.proteomes]
