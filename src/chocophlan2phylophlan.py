#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__version__ = '0.01'
__date__ = '11 Apr 2018'


import os
import sys
import argparse as ap
import configparser as cp
import pickle
import resource
import multiprocessing.dummy as dummy
from collections import Counter
from functools import partial
import copy
import glob
import src.process_proteomes as process_proteomes
from src.panproteomes import Panproteome
from operator import itemgetter

if __name__ == '__main__':
    import utils
else:
    import src.utils as utils

CLUSTER = 90

def chocophlan2phylophlan(config):
    taxontree = pickle.load(open(os.path.join(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))
    proteomes = pickle.load(open(os.path.join(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
    d_taxids = taxontree.lookup_by_taxid()

    reference_proteomes = [proteome for proteome in proteomes if proteomes[proteome]['isReference']]

    
    with open(os.path.join(config['download_base_dir'],
                           config['exportpath_phylophlan'])) as phylophlanout:
        for rfid in reference_proteomes:
            tax_id = proteomes[rfid]['tax_id']
            species_name = d_taxids[tax_id].name

            isSpecies = True if d_taxids[tax_id].rank == 'species' else False
            isStrain = True if d_taxids[d_taxids[tax_id].parent_tax_id]
            if 
            panproteome = pickle.load(open(os.path.join(config['download_base_dir'], 
                                                     config['relpath_panproteomes_dir'], 
                                                        'species',
                                                        CLUSTER,
                                                        tax_id,
                                                        '.pkl')))

            core_genes = Panproteome.find_core_genes(panproteome)
            line = '{}\t{}\t{}\t{}\n'.format(tax_id, species_name, rfid, ';'.join(core_genes))