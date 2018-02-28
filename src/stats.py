#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '03 Jan 2018'
 

import os
import sys
import argparse as ap
import configparser as cp
import pickle
import resource
import glob
import time
if __name__ == '__main__':
    import utils
else:
    import src.utils

def stats(config):
    proteomes = pickle.load(open('{}/pickled/proteomes.pkl'.format(config['download_base_dir']), 'rb'))
    uniprotkb = pickle.load(open('{}/pickled/uniprotkb_idmap.pkl'.format(config['download_base_dir']),'rb'))
    uniref100 = pickle.load(open('{}/pickled/uniref100_idmap.pkl'.format(config['download_base_dir']), 'rb'))
    uniref90 = pickle.load(open('{}/pickled/uniref90_idmap.pkl'.format(config['download_base_dir']), 'rb'))
    uniref50 = pickle.load(open('{}/pickled/uniref50_idmap.pkl'.format(config['download_base_dir']), 'rb'))
    taxontree = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
    d_ranks = taxontree.lookup_by_rank()

    s = '-----Taxonomy-----\n'
          'Total number of species: {}\n'
          '\tBacteria: {}\n'
          '\tArchaea: {}\n'
          'Total number of strains: {}\n'
          '-----UniProtKB-----\n'
          'Total number of proteomes: {}\n'
          'Total number of proteins in proteomes: {}\n'
          'Total number of proteins NOT PRESENT in a proteome: {}\n'
          'Total number of proteins in UniProtKB database: {}\n'
          'Total number of UniRef100 clusters: {}\n'
          '\t\tUniref90 clusters: {}\n'
          '\t\tUniref50 clusters: {}\n'
          'Total RAM usage by ChocoPhlAn indices: {} Gb\n'
          '\n'
          '-----Panproteomes-----\n'
          'Total created panproteomes: {}\n'
          .format(len(d_ranks['species']),
                  sum([len(x) for x in d_ranks['species'] if not x.initially_terminal]),
                  len([x for x in d_ranks['species'] if d_ranks['superkingdom'][0].is_parent_of(x)]),
                  len([x for x in d_ranks['species'] if d_ranks['superkingdom'][1].is_parent_of(x)]),
                  len(proteomes),
                  sum([len(v['members']) for _,v in proteomes.items()]),
                  len(uniprotkb)-sum([len(v['members']) for _,v in proteomes.items()]),
                  len(uniprotkb),
                  len(uniref100),
                  len(uniref90),
                  len(uniref50),
                  resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / pow(2,20),
                  len(glob.glob(config['download_base_dir']+config['relpath_panproteomes_dir']+'/*/*/*')))
    print(s)
    for r, d, f in os.walk(config['download_base_dir']+config['relpath_panproteomes_dir']):
        print(r.replace(config['download_base_dir']+config['relpath_panproteomes_dir']+'/',''), len(f)) if len(f)>0 else None
if __name__ == '__main__':
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']            #UPDATE TO STATS

    stats(config)

    sys.exit(0)
