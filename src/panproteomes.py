#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__version__ = '0.01'
__date__ = '04 Jan 2018'


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
import src.process_proteomes as process_proteomes
from operator import itemgetter

if __name__ == '__main__':
    import utils
else:
    import src.utils as utils

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class Panproteome:
    ranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    kingdom_to_process = ('Archaea', 'Bacteria')
    coreness = 0.99
    
    def __init__(self, config):
        self.uniparc = {}
        try:
            if config['verbose']:
                utils.info('Loading pickled databases...')
            import gc
            gc.disable()
            self.uniref100 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref100_idmap']),'rb'))
            self.uniref90 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref90_idmap']),'rb'))
            self.uniref50 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref50_idmap']),'rb'))

            self.idmapping = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_uniref_idmap']), 'rb'))
            self.taxontree = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))
            self.proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
            
            if not self.taxontree.are_leaves_trimmed:
                self.taxontree.remove_subtree_taxonomy_by_level('superkingdom', self.kingdom_to_process)
                self.taxontree.tree.root = self.taxontree.tree.root.clades[4]
                self.taxontree.remove_leaves_without_proteomes()
                self.taxontree.are_leaves_trimmed=True
            self.config = config
            gc.enable()
            if config['verbose']:
                utils.info('Finished.\n')
        except Exception as e:
            utils.error(str(e), exit=True)
            raise

    def get_child_proteomes(self, clade):
        if clade.initially_terminal and hasattr(clade,'proteomes'):
            return clade.proteomes
        pp = copy.deepcopy(clade.proteomes) if hasattr(clade,'proteomes') else set()
        for c in clade.clades:
            pp.update(self.get_child_proteomes(c))
        return pp
    
    # Create a dictionary of all proteins and the proteomes in which are found 
    # Find the UniRef cluster at which the protein belongs to and create a dictionary
    # with all the UniRef clusters and the proteomes in which are present
    def process_panproteome(self, item_rank):
        item, rank, cluster = item_rank
        cluster_index = 0 if cluster == 100 else (1 if cluster == 90 else 2)
        proteomes_to_process = self.get_child_proteomes(item)
        panproteome = {}
        uniref_panproteome = {}
        uniref_panproteome['cluster'] = cluster
        uniref_panproteome['tax_id'] = item.tax_id
        uniref_panproteome['rank'] = rank
        uniref_panproteome['number_proteomes'] = len(proteomes_to_process)
        uniref_panproteome['coreness_threshold'] = round(uniref_panproteome['number_proteomes'] * self.coreness)
        uniref_panproteome['members'] = {}

        if len(proteomes_to_process):
            for protein, proteome_id in ((entry,proteome_id) for proteome_id in proteomes_to_process if proteome_id in self.proteomes for entry in self.proteomes[proteome_id]['members']):
                uniref_cluster = self.idmapping.get(protein, None)
                if uniref_cluster is None:
                    urefs = list(set(uref for uref in (self.idmapping.get(upkb) for upkb in self.uniparc.get(protein,['']*7)[6]) if uref is not None))
                    if len(urefs):
                        uniref_cluster = urefs[0][cluster_index]
                else:
                    uniref_cluster = uniref_cluster[cluster_index]
                if uniref_cluster is not None:
                    if uniref_cluster not in uniref_panproteome['members']:
                        uniref_panproteome['members'][uniref_cluster] = { 'coreness': 0,
                                                               'uniqueness': {},
                                                               'copy_number': [],
                                                               'proteomes_present': set(),
                                                               'external_hits' : []
                                                             }
                    uniref_panproteome['members'][uniref_cluster]['proteomes_present'].add(proteome_id)
                    uniref_panproteome['members'][uniref_cluster]['coreness'] = len(uniref_panproteome['members'][uniref_cluster]['proteomes_present'])
                    uniref_panproteome['members'][uniref_cluster]['copy_number'].append(protein)
            

            if len(uniref_panproteome):
                pickle.dump(uniref_panproteome, open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], rank, cluster, item.tax_id),'wb'))

    def create_panproteomes(self, cluster):
        if cluster == 100:
            ranks_to_process = self.ranks[-1::]
        elif cluster == 90:
            ranks_to_process = self.ranks[:3:-1]
        else:
            ranks_to_process = self.ranks
        
        if self.config['verbose']:
            utils.info('Starting creating panproteomes for {}...\n'.format(', '.join(ranks_to_process)))

        for k in ranks_to_process:
            os.makedirs('{}/{}/{}/{}'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], k, cluster), exist_ok=True)
        
        d_ranks = self.taxontree.lookup_by_rank()
        
        elems = [(item, rank, cluster) for rank in ranks_to_process for item in d_ranks[rank]]
        try:
            with dummy.Pool(20) as pool:
                d = [_ for _ in pool.imap_unordered(self.process_panproteome, elems)]

        except Exception as e:
            utils.error(str(e))
            raise


    def calculate_uniqueness(self, panproteome_fp):
        panproteome = pickle.load(open(panproteome_fp,'rb'))
        current_cluster = panproteome['cluster']
        external_clusters = {}
        t_external_clusters = {}

        uniref = getattr(self, 'uniref{}'.format(current_cluster))
        files_to_load = [(prot,uniref['UniRef{}_{}'.format(current_cluster, prot)]) for prot in panproteome['members'] if len(prot)]
        files_to_load.sort(key = lambda x:x[1])

        d_taxids = self.taxontree.lookup_by_taxid()

        item_descendant = [x.tax_id for k in (d_taxids[panproteome['tax_id']].get_terminals(), d_taxids[panproteome['tax_id']].get_nonterminals()) for x in k]

        # Intra cluster uniqueness
        handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], current_cluster, files_to_load[0][1]), 'rb')
        chunk = pickle.load(handle)
        for uniref_id, chunk_id in files_to_load:
            if '{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], current_cluster, chunk_id) != handle.name:
                handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], current_cluster, chunk_id), 'rb')
                chunk = pickle.load(handle)
            
            uniref_entry = process_proteomes.uniref_tuple_to_dict(chunk['UniRef{}_{}'.format(current_cluster, uniref_id)])
            external_clusters[uniref_id] = dict(zip(['UniRef100','UniRef90','UniRef50'],itemgetter(*['UniRef100','UniRef90','UniRef50'])(uniref_entry)))
            
            taxa_is_present = set(x[1] for x in uniref_entry['members'])
            internal_hits = [x for x in taxa_is_present if x in item_descendant]
            external_hits = [x for x in taxa_is_present if x not in item_descendant]

            panproteome['members'][uniref_id]['uniqueness']['{}_{}'.format(current_cluster,current_cluster)] = len(external_hits)

        for cluster, ext_unirefs in external_clusters.items():
            for uniref_c, ids in ext_unirefs.items():
                if uniref_c not in t_external_clusters:
                    t_external_clusters[uniref_c] = {}
                t_external_clusters[uniref_c][cluster] = ids

        for uniref_c in t_external_clusters:
            cluster = int(uniref_c.replace('UniRef',''))
            uniref_map =  getattr(self, 'uniref{}'.format(cluster))

            files_to_load = [(prot,ortholog, uniref_map['{}'.format(ortholog)]) for prot, ortholog in t_external_clusters[uniref_c].items() if len(prot)]
            files_to_load.sort(key = lambda x:x[2])
            
            handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, files_to_load[0][2]), 'rb')            
            chunk = pickle.load(handle)

            for prot, ortholog, chunk_id in files_to_load:
                if '{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, chunk_id) != handle.name:
                    handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, chunk_id), 'rb')
                    chunk = pickle.load(handle)

                uniref_entry = process_proteomes.uniref_tuple_to_dict(chunk['{}'.format(ortholog)])
                taxa_is_present = [x[1] for x in uniref_entry['members']]
                internal_hits = [x for x in taxa_is_present if x in item_descendant]
                external_hits = [x for x in taxa_is_present if x not in item_descendant]

                panproteome['members'][uniref_id]['uniqueness']['{}_{}'.format(current_cluster, cluster)] = len(external_hits)

        pickle.dump(open(panproteome_fp,'wb'))

    @staticmethod
    def find_core_genes(self, panproteome):
        return [gene for gene, _ in filter(lambda gene:gene[1]['coreness'] > panproteome['coreness_threshold'], panproteome['members'].items())]

    def rank_genes(self, panproteome):
        pass

def generate_panproteomes(config):
    p = Panproteome(config)
    # p.create_panproteomes(100)
    d_tax = p.taxontree.lookup_by_taxid()
    p.process_panproteome((d_tax[562],'species',90))

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['panproteomes']
    generate_panproteomes(config)
