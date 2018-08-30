#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__date__ = '11 Apr 2018'


from _version import __CHOCOPhlAn_version__
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

CLUSTER = 90

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class chocophlan2phylophlan:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')

        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.config = config
        self.exportpath = '{}/{}'.format(self.config['export_dir'], self.config['exportpath_phylophlan'])
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)
        if config['verbose']:
            utils.info('Finished.\n')

    def process(self, tax_id):
        if not terminating.is_set():
            fp = '{}/{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], 
                                                     self.config['relpath_panproteomes_dir'], 
                                                        'species',
                                                        CLUSTER,
                                                        tax_id)
            if os.path.exists(fp):
                panproteome = pickle.load(open(fp,'rb'))
                taxa_str = self.taxontree.print_full_taxonomy(tax_id)
                core_genes = Panproteome.find_core_genes(panproteome)
                
                d_out_core = (tax_id, taxa_str, core_genes)
                d_out_refp = (tax_id, taxa_str, [])
                d_out_refg = (tax_id, taxa_str, [])

                #add first reference, non rendundant, redundant
                d_all_prot = {'re':[], 'nr': [], 'rr':[]}
                for p in self.taxontree.get_child_proteomes(self.taxontree.taxid_n[tax_id]):
                    if self.proteomes[p]['isReference']:
                        d_all_prot['re'].append(p)
                    elif self.proteomes[p]['upi']:
                        d_all_prot['rr'].append(p)
                    else:
                        d_all_prot['nr'].append(p)

                d_out_refp[2].extend(d_all_prot['re'] + d_all_prot['rr'] + d_all_prot['nr'])

                try:
                    d_out_refg[2].extend([dict(self.proteomes[p]['ncbi_ids']).get('GCSetAcc','') for p in d_out_refp[2] if p in self.proteomes and 'ncbi_ids' in self.proteomes[p]])
                except Exception as e:
                    utils.info(tax_id+'\n')
                    raise e
                    terminating.set()

                return (d_out_core, d_out_refp, d_out_refg)
            else:
                print('Panproteome {} not available'.format(fp))
        else:
            terminating.set()

    def chocophlan2phylophlan(self):
        d_out_core = {}
        d_out_refp = {}
        d_out_refg = {}

        reference_species = set(
            self.proteomes[proteome]['tax_id'] 
            if self.taxontree.taxid_n[self.proteomes[proteome]['tax_id']].rank == 'species' 
            else self.taxontree.go_up_to_species(self.proteomes[proteome]['tax_id']) 
            for proteome in self.proteomes if self.proteomes[proteome]['isReference']
            )

        if self.config['verbose']:
            utils.info('Started exporting CHOCOPhlAn data for PhyloPhlAn2...')

        terminating = dummy.Event()
        with dummy.Pool(initializer = init_parse, initargs = (terminating, ), processes = self.config['nproc']) as pool:
            d = [x for x in pool.imap_unordered(self.process, reference_species, chunksize=10)]

        if self.config['verbose']:
            utils.info('Done.\n')

        for item in d:
            cor, ref, gen = item
            d_out_core[cor[0]] = (cor[1], cor[2])
            d_out_refp[ref[0]] = (ref[1], ref[2])
            d_out_refg[gen[0]] = (gen[1], gen[2])

        if self.config['verbose']:
            utils.info('Exporting core proteins, proteomes and genomes...\n')

        # REMOVE SPECIES WITHOUT ANY GENOME
        # NO PROTEOMES
        with open('{}/taxa2proteomes.txt'.format(self.exportpath), 'w') as t2p_out:
            with open('{}/taxa2core.txt'.format(self.exportpath), 'w') as t2c_out:
                with open('{}/taxa2genomes.txt'.format(self.exportpath), 'w') as t2g_out:
                    lines_t2p = ['#CHOCOPhlAn version {}\n'.format(__CHOCOPhlAn_version__), '#'+open('data/relnotes.txt').readline(), '#NCBI Taxonomy id\tFull Taxonomy\tList of proteomes\n']
                    lines_t2c = ['#CHOCOPhlAn version {}\n'.format(__CHOCOPhlAn_version__), '#'+open('data/relnotes.txt').readline(), '#NCBI Taxonomy id\tFull Taxonomy\tList of core proteins\n']
                    lines_t2g = ['#CHOCOPhlAn version {}\n'.format(__CHOCOPhlAn_version__), '#'+open('data/relnotes.txt').readline(), '#NCBI Taxonomy id\tFull Taxonomy\tList of genomes\n']

                    lines_t2p.extend(['{}\t{}\thttp://www.uniprot.org/uniprot/?query=proteome:{{}}&compress=yes&force=true&format=fasta\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_refp.items()])
                    lines_t2c.extend(['{}\t{}\thttp://www.uniprot.org/uniref/UniRef90_{{}}.fasta\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_core.items()])
                    lines_t2g.extend(['{}\t{}\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_refg.items()])

                    t2p_out.writelines(lines_t2p)
                    t2c_out.writelines(lines_t2c)
                    t2g_out.writelines(lines_t2g)
        if self.config['verbose']:
            utils.info('Done.\n')

def export_to_phylophlan(config):
    c = chocophlan2phylophlan(config)
    c.chocophlan2phylophlan()

if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['export']
    c = chocophlan2phylophlan(config)
    c.chocophlan2phylophlan()
    t1 = time.time()
    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
