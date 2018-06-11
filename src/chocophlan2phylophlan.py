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

CLUSTER = 90

class chocophlan2phylophlan:
    def __init__(self, config):
        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.config = config

    def process(self, tax_id):
        if not self.taxontree.taxid_n[tax_id].rank == 'species':
            tax_id = self.taxontree.go_up_to_species(tax_id)

        fp = '{}/{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], 
                                                 self.config['relpath_panproteomes_dir'], 
                                                    'species',
                                                    CLUSTER,
                                                    tax_id)
        if os.path.exists(fp):
            panproteome = pickle.load(open(fp,'rb'))
            taxa_str = self.taxontree.print_full_taxonomy(tax_id)
            core_genes = Panproteome.find_core_genes(panproteome)
            
            d_out_core = (tax_id, taxa_str, set())
            d_out_refp = (tax_id, taxa_str, set())
            d_out_refg = (tax_id, taxa_str, set())

            d_out_core[2].update(core_genes)

            #add first reference, non rendundant, redundant
            d_all_prot = {'re':[], 'nr': [], 'rr':[]}
            for p in self.taxontree.get_child_proteomes(self.taxontree.taxid_n[tax_id]):
                if self.proteomes[p]['isReference']:
                    d_all_prot['re'].append(p)
                elif self.proteomes[p]['upi']:
                    d_all_prot['rr'].append(p)
                else:
                    d_all_prot['nr'].append(p)

            d_out_refp[2].add(d_all_prot['re'] + d_all_prot['rr'] + d_all_prot['nr'])

            d_out_refg[2].update([dict(self.proteomes[p]['ncbi_ids']).get('GCSetAcc','') for p in d_out_refp[2]])

            return (d_out_core, d_out_refp, d_out_refg)
        else:
            print('Panproteome {} not available'.format(fp))

    def chocophlan2phylophlan(self):
        d_out_core = {}
        d_out_refp = {}
        d_out_refg = {}

        reference_species = set(
            self.proteomes[proteome]['tax_id'] 
            if self.taxontree.taxid_n[tax_id].rank == 'species' 
            else self.taxontree.go_up_to_species(tax_id) 
            for proteome in self.proteomes if self.proteomes[proteome]['isReference']
            )

        with dummy.Pool(self.config['nproc']) as pool:
            d = [x for x in pool.imap_unordered(self.process, reference_species, chunksize=10)]

        for item in d:
            core, ref, gen = item
            if core[0] not in d_out_core:
                d_out_core[core[0]] = (core[1], set())
            d_out_core[core[0]][1].update(core[2])

            if ref[0] not in d_out_refp:
                d_out_refp[ref[0]] = (ref[1], set())
            d_out_refp[ref[0]][1].update(ref[2])

            if gen[0] not in d_out_refg:
                d_out_refg[gen[0]] = (gen[1], set())
            d_out_refg[gen[0]][1].update(gen[2])

    ##CREATE DIRS
    ##TAXONOMY, URL, VERSION
        with open('{}/{}/taxa2proteomes.txt'.format(self.config['export_dir'], self.config['exportpath_phylophlan']), 'w') as t2p_out:
            with open('{}/{}/taxa2core.txt'.format(self.config['export_dir'], self.config['exportpath_phylophlan']), 'w') as t2c_out:
                with open('{}/{}/taxa2genomes.txt'.format(self.config['export_dir'], self.config['exportpath_phylophlan']), 'w') as t2g_out:
                    lines_t2p = ['#'+open('data/relnotes.txt').readline().strip()]
                    lines_t2c = ['#'+open('data/relnotes.txt').readline().strip()]
                    lines_t2g = ['#'+open('data/relnotes.txt').readline().strip()]

                    lines_t2p.update(['{}\t{}\thttp://www.uniprot.org/uniprot/?query=proteome:\{\}&compress=yes&force=true&format=fasta {}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_refp.items()])
                    lines_t2c.update(['{}\t{}\thttp://www.uniprot.org/uniref/UniRef90_\{\}.fasta {}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_core.items()])
                    lines_t2g.update(['{}\t{}\t{}\n'.format(tax_id, entry[0], ';'.join(entry[1])) for tax_id, entry in d_out_refg.items()])

                    t2p_out.writelines(lines_t2p)
                    t2c_out.writelines(lines_t2c)
                    t2g_out.writelines(lines_t2g)


if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['chocophlan2phylophlan']
    c = chocophlan2phylophlan(config)
    c.chocophlan2phylophlan()
    t1 = time.time()
    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
