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
import glob
from operator import itemgetter

if __name__ == '__main__':
    import utils
    import process_proteomes
else:
    import src.utils as utils
    import src.process_proteomes as process_proteomes

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class Panproteome:
    ranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    kingdom_to_process = ('Archaea', 'Bacteria')
    
    def __init__(self, config):
        try:
            if config['verbose']:
                utils.info('Loading pickled databases...')
            import gc
            gc.disable()
            self.uniparc = {}
            # self.uniref100_tax_id_map = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref100_taxid_idmap']),'rb'))
            self.uniref90_tax_id_map = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref90_taxid_idmap']),'rb'))
            # self.uniref50_tax_id_map = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref50_taxid_idmap']),'rb'))
            # self.uniprotkb = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_idmap']), 'rb'))
            # for i in range(1,196):
                # self.uniparc.update(pickle.load(open('data/pickled/uniparc_{}.pkl'.format(i),'rb')))
            # self.idmapping = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_uniref_idmap']), 'rb'))
            self.taxontree = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))
            self.proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
            self.d_ranks = self.taxontree.lookup_by_rank()
            
            self.config = config
            gc.enable()
            if config['verbose']:
                utils.info('Finished.\n')
        except Exception as e:
            utils.error(str(e), exit=True)
            raise
   
    # Create a dictionary of all proteins and the proteomes in which are found 
    # Find the UniRef cluster at which the protein belongs to and create a dictionary
    # with all the UniRef clusters and the proteomes in which are present
    def process_panproteome(self, item_rank):
        try:
            item, rank, cluster = item_rank
            cluster_index = 0 if cluster == 100 else (1 if cluster == 90 else 2)

            proteomes_to_process = self.taxontree.get_child_proteomes(item)
            panproteome = {}
            uniref_panproteome = {}
            uniref_panproteome['cluster'] = cluster
            uniref_panproteome['tax_id'] = item.tax_id
            uniref_panproteome['rank'] = rank
            uniref_panproteome['number_proteomes'] = len(proteomes_to_process)
            uniref_panproteome['coreness_value'] = 0.5 if uniref_panproteome['number_proteomes'] == 2 else (.66 if uniref_panproteome['number_proteomes'] == 3 else 0.75)
            uniref_panproteome['coreness_threshold'] = round(uniref_panproteome['coreness_value'] * uniref_panproteome['number_proteomes'])
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
                                                                   'copy_number': Counter(),
                                                                   'proteomes_present': set(),
                                                                   'external_hits' : {}
                                                                 }
                        uniref_panproteome['members'][uniref_cluster]['proteomes_present'].add(proteome_id)
                        uniref_panproteome['members'][uniref_cluster]['coreness'] = len(uniref_panproteome['members'][uniref_cluster]['proteomes_present'])
                        uniref_panproteome['members'][uniref_cluster]['copy_number'][(protein, proteome_id)] +=1 

                if len(uniref_panproteome):
                    pickle.dump(uniref_panproteome, open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], rank, cluster, item.tax_id),'wb'))
        except Exception as e:
            utils.error(str(e), exit=True)
            raise

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
        
        try:
            with dummy.Pool(20) as pool:
                d = [_ for _ in pool.imap_unordered(self.process_panproteome, ((item, rank, cluster) for rank in ranks_to_process for item in self.d_ranks[rank]))]

        except Exception as e:
            utils.error(str(e))
            raise

    def calculate_uniqueness(self, panproteome_fp):
        if not terminating.is_set():
            panproteome = pickle.load(open(panproteome_fp, 'rb'))
            panproteome_cluster = panproteome['cluster']
            utils.info(str(panproteome['tax_id'])+'\n')
            item_descendant = [x.tax_id for k in (self.taxontree.taxid_n[panproteome['tax_id']].get_terminals(), self.taxontree.taxid_n[panproteome['tax_id']].get_nonterminals()) for x in k]

            def get_upper_clusters(cluster):
                cs = [100, 90, 50]
                return cs[cs.index(cluster):]

            for cluster in get_upper_clusters(panproteome_cluster):
                if panproteome_cluster == cluster:
                    files_to_load = [(pangene, 'UniRef{}_{}'.format(cluster, pangene)) for pangene in panproteome['members'] if len(pangene)]
                else:
                    self.__setattr('uniref{}_tax_id_map'.format(cluster), pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref{}_taxid_idmap'.format(cluster)]),'rb')))
                    cluster_index = 0 if cluster == 100 else (1 if cluster == 90 else 2)
                    starting_clusters = self.__getattribute__('uniref{}_tax_id_map'.format(panproteome_cluster))
                    files_to_load = [(pangene, starting_clusters['UniRef{}_{}'format(panproteome_cluster, pangene)][1][cluster_index]) 
                                        for pangene in panproteome['members'] if len(pangene)]

                for pangene, uniref_id in files_to_load:
                    taxa_is_present = set(self.uniref90_tax_id_map[uniref_id])
                    external_hits = [x for x in taxa_is_present if x not in item_descendant]

                    panproteome['members'][pangene]['uniqueness']['{}_{}'.format(panproteome_cluster,cluster)] = len(external_hits)
                    panproteome['members'][pangene]['external_hits']['{}_{}'.format(panproteome_cluster,cluster)] = external_hits

                # Intra cluster uniqueness
                # handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, files_to_load[0][2]), 'rb')
                # chunk = pickle.load(handle)
                # for pangene, uniref_id, chunk_id in files_to_load:
                #     if '{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, chunk_id) != handle.name:
                #         handle = open('{}/pickled/uniref{}_{}.pkl'.format(self.config['download_base_dir'], cluster, chunk_id), 'rb')
                #         chunk = pickle.load(handle)
                    
                #     uniref_entry = process_proteomes.uniref_tuple_to_dict(chunk[uniref_id])
                #     external_clusters[pangene] = dict(zip(['UniRef100','UniRef90','UniRef50'],itemgetter(*['UniRef100','UniRef90','UniRef50'])(uniref_entry)))
                    
                #     taxa_is_present = set(x[1] for x in uniref_entry['members'])
                #     # internal_hits = [x for x in taxa_is_present if x in item_descendant]
                #     external_hits = [x for x in taxa_is_present if x not in item_descendant]

                #     panproteome['members'][pangene]['uniqueness']['{}_{}'.format(panproteome_cluster,cluster)] = len(external_hits)
                #     panproteome['members'][pangene]['external_hits']['{}_{}'.format(panproteome_cluster,cluster)] = external_hits
                [panproteome['members'][k].update(v) for file in result for k, v in file.items()]

            pickle.dump(panproteome, open(panproteome_fp, 'wb'))
        else:
            terminating.set()

    @staticmethod
    def find_core_genes(panproteome):
        return [gene for gene, _ in filter(lambda gene:gene[1]['coreness'] >= panproteome['coreness_threshold'], panproteome['members'].items())]


    def rank_genes(self, panproteome):
        pass


    def extract_protein_sequence(self, items):
        if not terminating.is_set():
            chunk, ids = items
            cluster = ids[0].split('_')[0]
            with open('{}/pickled/{}_{}'.format(self.config['download_base_dir'], cluster, chunk), 'rb') as pickled_chunk:
                uniprot_chunk = pickle.load(pickled_chunk)
            entries = itemgetter(*ids)(uniprot_chunk)
            seqs = [itemgetter(0,2)(i) for i in entries]
            with open(outputdir,'w') as w_out:
                wout.writelines(['>{}\n{}\n'.format(s[0],s[1]) for s in seqs])
        else:
            terminating.set()


    @staticmethod
    def export_panproteome_fasta(self, panproteome):
        # with open('export/{}_uniref{}_panproteome.txt'.format(panproteome['tax_id'], panproteome['cluster'])) as export:
        files_to_load = itemgetter(*filter(None,panproteome['members'].keys()))(self.uniref90)
        cluster = panproteome['cluster']
        d_files_to_load = {}
        for up_id, chunk in zip(panproteome['members'].keys(), files_to_load):
            if chunk not in d_files_to_load:
                d_files_to_load[chunk] = []
            d_files_to_load[chunk].append('UniRef{}_{}'.format(cluster, up_id))

        terminating = dummy.Event()
        with dummy.Pool(initializer = init_parse, initargs = (terminating, ), processes = len(d_files_to_load.keys())) as pool:
            result = [pangene for pangene in pool.imap_unordered(extract_protein_sequence, d_files_to_load.items())]

        pass

    def export_pangenome_fasta(self, panproteome):


        pass


def generate_panproteomes(config):
    p = Panproteome(config)

    # with dummy.Pool(config['nproc']) as pool:
    #    d = [_ for _ in pool.imap_unordered(p.create_panproteomes, [100,90,50])]
    
    # terminating = dummy.Event()
    with dummy.Pool(round(config['nproc']/2)) as pool:
        [_ for _ in pool.imap_unordered(p.calculate_uniqueness, (file for file in glob.iglob('{}{}/*/*/*.pkl'.format(p.config['download_base_dir'], p.config['relpath_panproteomes_dir']))))]
      

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['panproteomes']
    generate_panproteomes(config)
    
