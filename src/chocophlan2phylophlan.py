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


def chocophlan2phylophlan(config):
    taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
    proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
    d_taxids = taxontree.taxid_n

    reference_proteomes = [proteome for proteome in proteomes if proteomes[proteome]['isReference']]
    d_out_core = {}
    d_out_refp = {}

    
    for rfid in reference_proteomes:
        tax_id = proteomes[rfid]['tax_id']

        if not d_taxids[tax_id].rank == 'species':
            tax_id = taxontree.go_up_to_species(tax_id)

        fp = '{}/{}/{}/{}/{}.pkl'.format(config['download_base_dir'], 
                                                 config['relpath_panproteomes_dir']+"_old", 
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
                

            core_genes = Panproteome.find_core_genes(panproteome)
            
            if tax_id not in d_out_core:
                d_out_core[tax_id] = (taxa_str, set())
            #TODO
            #URL OF CORE UNIREF
            #core_genes = [for ] join con url http://www.uniprot.org/uniref/{}.fasta
            d_out_core[tax_id][1].update(core_genes)

            if tax_id not in d_out_refp:
                d_out_refp[tax_id] = (taxa_str, set())
            #add first reference, non rendundant, redundant
            d_out_refp[tax_id][1].add('http://www.uniprot.org/uniprot/?query=proteome:{}&compress=yes&force=true&format=fasta'.format(rfid))
        else:
            print('Panproteome {} not available'.format(fp))

##CREATE DIRS
##TAXONOMY, URL, VERSION
    with open('{}/{}/taxa2proteomes.txt'.format(config['export_dir'], config['exportpath_phylophlan']), 'w') as t2p_out:
        with open('{}/{}/taxa2core.txt'.format(config['export_dir'], config['exportpath_phylophlan']), 'w') as t2c_out:
            lines_t2p = ['{}\t{}\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_refp.items()]
            lines_t2c = ['{}\t{}\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_core.items()]

            t2p_out.writelines(lines_t2p)
            t2c_out.writelines(lines_t2c)


if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['chocophlan2phylophlan']
    chocophlan2phylophlan(config)
    t1 = time.time()
    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
