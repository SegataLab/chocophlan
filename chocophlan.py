#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '28 Sep 2017'


import src.utils as utils
import src.build_taxontree as build_taxontree
import src.download as download
import src.parse_uniprot as parse_uniprot
import src.panproteomes as panproteomes
import src.chocophlan2phylophlan as chocophlan2phylophlan
import src.stats as stats
import src.export_to_metaphlan2 as export_to_metaphlan2
import time
from src.build_taxontree import Nodes as Nodes
import sys
import multiprocessing as mp

def chocophlan():
    utils.info('CHOCOPhlAn: Cluster of HOmologous Cdses fOr PHyLogenetic ANalysis\n')
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    # download.download(config['download'], verbose=config['download']['verbose'])
    # build_taxontree.do_extraction(config['build_taxontree'], verbose=config['build_taxontree']['verbose'])
    # parse_uniprot.parse_uniprot(config['parse_uniprot'])
    # download.download_ncbi_from_proteome_pickle(config['parse_uniprot'])
    # panproteomes.generate_panproteomes(config['panproteomes'])
    #stats.generate_stats(config['stats'])
    #export_to_metaphlan2.run_all(config['export'])
    chocophlan2phylophlan.export_to_phylophlan(config['export'])

if __name__ == '__main__':
    t0 = time.time()
    chocophlan()
    t1 = time.time()
    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
