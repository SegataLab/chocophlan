#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__version__ = '0.01'
__date__ = '11 Apr 2018'


import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
import glob
from src.panproteomes import Panproteome
from operator import itemgetter

if __name__ == '__main__':
    import utils
else:
    import src.utils as utils

ranks2code = {'superkingdom': 'k', 'phylum': 'p', 'class': 'c',
                      'order': 'o', 'family': 'f', 'genus': 'g', 'species': 's', 'taxon': 't'}
order = ('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
CLUSTER = 90

def go_up_to_species(taxid):
    father = d_taxids[taxid].parent_tax_id
    if not d_taxids[father].rank == 'species':
        return go_up_to_species(father)
    return father
    

def chocophlan2phylophlan(config):
    taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
    proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
    d_taxids = taxontree.lookup_by_taxid()

    reference_proteomes = [proteome for proteome in proteomes if proteomes[proteome]['isReference']]
    d_out_core = {}
    d_out_refp = {}

    with open(os.path.join(config['download_base_dir'],
                           config['exportpath_phylophlan'])) as phylophlanout:
        for rfid in reference_proteomes:
            tax_id = proteomes[rfid]['tax_id']
            fp = '{}/{}/{}/{}/{}.pkl'.format(config['download_base_dir'], 
                                                     config['relpath_panproteomes_dir'], 
                                                        'species',
                                                        CLUSTER,
                                                        tax_id)
            if os.path.exists(fp):
                panproteome = pickle.load(open(fp,'rb'))
                path = [p for p in d_taxids[1].get_path(d_taxids[tax_id]) if p.rank in ranks2code or (p.rank=='norank' and p.initially_terminal)]
                taxa_str = ''
                hasSpecies = any([True if p.rank == 'species' else False for p in path])
                for i in range(0, len(path)):
                    if path[i].rank in ranks2code:
                        taxa_str += '{}__{}'.format(ranks2code[path[i].rank], path[i].name)
                    elif path[i].rank == 'norank':
                        if not path[i].initially_terminal:
                            if path[i+1].rank in ranks2code and order.index(ranks2code[path[i+1].rank])==order.index(ranks2code[path[i-1].rank])+1:
                                continue
                            elif path[i+1].rank == 'norank':
                                if path[i+1].tax_id != tax_id:
                                    continue 
                            else:
                                taxa_str += '{}__unclassified_{}'.format(order[order.index(ranks2code[path[i-1].rank])+1], path[i-1].name)
                        else:
                            if hasSpecies:
                                taxa_str += '{}__{}'.format(ranks2code['taxon'], path[i].name)   
                    if i < len(path)-1:
                        taxa_str += '|'
                    
                if not d_taxids[tax_id].rank == 'species':

                core_genes = Panproteome.find_core_genes(panproteome)
                d_out_core[tax_id] = (taxa_str, core_genes)
                if tax_id not in d_out_refp:
                    d_out_refp[tax_id] = set()
                d_out_refp[tax_id].add(rfid)

            line = '{}\t{}\t{}\t{}\n'.format(tax_id, species_name, rfid, ';'.join(core_genes))