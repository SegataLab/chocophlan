#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

from _version import __CHOCOPhlAn_version__
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
import itertools
import re
from operator import itemgetter
import importlib
if __name__ == "__main__":
    import utils
    import stats
    import process_proteomes
else:
    utils = importlib.import_module('src.utils')
    stats = importlib.import_module('src.stats')
    process_proteomes = importlib.import_module('src.process_proteomes')

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class Panproteome:
    ranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    kingdom_to_process = ('Archaea', 'Bacteria')

    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')
        self.uniparc = {}
        for i in filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/{}/uniparc'.format(config['download_base_dir'], config['pickled_dir']))):
            if config['verbose']:
                utils.info('{}\n'.format(i))
            self.uniparc.update({k:v[6] for k,v in pickle.load(open('{}/{}/uniparc/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb')).items() })
        self.idmapping = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_uniref_idmap']), 'rb'))
        self.taxontree = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
        self.d_ranks = self.taxontree.lookup_by_rank()
        
        self.config = config
        if config['verbose']:
            utils.info('Finished.\n')

    def get_upper_clusters(self, cluster):
        cs = [100, 90, 50]
        return cs[cs.index(cluster):]

    # Create a dictionary of all proteins and the proteomes in which are found 
    # Find the UniRef cluster at which the protein belongs to and create a dictionary
    # with all the UniRef clusters and the proteomes in which are present
    def process_panproteome(self, item_rank):
        if not terminating.is_set():
            item, rank, panproteome_cluster = item_rank
            cluster_index = 0 if panproteome_cluster == 100 else (1 if panproteome_cluster == 90 else 2)

            item_descendant = [x.tax_id for x in itertools.chain.from_iterable((item.get_terminals(), item.get_nonterminals()))]
            proteomes_to_process = self.taxontree.get_child_proteomes(item)
            if self.config['discard_low_quality_genomes'] and rank == 'genus':
                proteomes_to_process = [p for p in proteomes_to_process if not self.taxontree.taxid_n[self.proteomes[p]['tax_id']].is_low_quality]

            uniref_panproteome = {}
            uniref_panproteome['cluster'] = panproteome_cluster
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
                        urefs = list(set(uref for uref in 
                            (self.idmapping.get(upkb) for upkb in self.uniparc.get(protein,['']))
                             if uref is not None))
                        if len(urefs):
                            uniref_cluster = urefs[0][cluster_index]
                    else:
                        uniref_cluster = uniref_cluster[cluster_index]
                    if uniref_cluster is not None:
                        if uniref_cluster not in uniref_panproteome['members']:
                            uniref_panproteome['members'][uniref_cluster] = { 'coreness': 0,
                                                                   'uniqueness': {},
                                                                   'uniqueness_nosp' : {},
                                                                   'copy_number': Counter(),
                                                                   'proteomes_present': set(),
                                                                   'external_hits' : {},
                                                                   'external_genomes' : {},
                                                                   'external_species' : {},
                                                                   'external_species_nosp' : {}
                                                                 }
                        uniref_panproteome['members'][uniref_cluster]['proteomes_present'].add(proteome_id)
                        uniref_panproteome['members'][uniref_cluster]['coreness'] = len(uniref_panproteome['members'][uniref_cluster]['proteomes_present'])
                        uniref_panproteome['members'][uniref_cluster]['copy_number'][(protein, proteome_id)] +=1 


                for cluster in self.get_upper_clusters(panproteome_cluster):
                    if panproteome_cluster == cluster:
                        files_to_load = [(pangene, 'UniRef{}_{}'.format(cluster, pangene)) for pangene in uniref_panproteome['members'] if len(pangene)]
                    else:
                        cluster_index = 0 if cluster == 100 else (1 if cluster == 90 else 2)
                        starting_clusters = self.__getattribute__('uniref{}_tax_id_map'.format(panproteome_cluster))
                        files_to_load = [(pangene, starting_clusters.get( 'UniRef{}_{}'.format(panproteome_cluster, pangene),
                                                                          ['',['','','']])[1][cluster_index]
                                                                        )
                                         for pangene in uniref_panproteome['members'] if len(pangene)]

                    destination_clusters = self.__getattribute__('uniref{}_tax_id_map'.format(cluster))
                    for pangene, uniref_id in files_to_load:
                        uniref_id = 'UniRef{}_{}'.format(cluster, uniref_id) if 'UniRef' not in uniref_id else uniref_id
                        taxa_is_present = destination_clusters.get('{}'.format(uniref_id),[''])[0]
                        external_hits = [x for x in taxa_is_present if x[1] not in item_descendant]
                        external_genomes = Counter(self.taxontree.go_up_to_species(taxid) for _, taxid, _ in external_hits)
                        if None in external_genomes: external_genomes.pop(None)
                        external_species = list(external_genomes.keys())
                        external_species_nosp = set(taxid for taxid in external_species if not self.taxontree.taxid_n[taxid].is_low_quality)

                        uniref_panproteome['members'][pangene]['uniqueness']['{}_{}'.format(panproteome_cluster,cluster)] = len(external_species)
                        uniref_panproteome['members'][pangene]['uniqueness_nosp']['{}_{}'.format(panproteome_cluster,cluster)] = len(external_species_nosp)
                        uniref_panproteome['members'][pangene]['external_species']['{}_{}'.format(panproteome_cluster,cluster)] = external_species
                        uniref_panproteome['members'][pangene]['external_species_nosp']['{}_{}'.format(panproteome_cluster,cluster)] = external_species_nosp
                        uniref_panproteome['members'][pangene]['external_genomes']['{}_{}'.format(panproteome_cluster,cluster)] = external_genomes
                        uniref_panproteome['members'][pangene]['external_hits']['{}_{}'.format(panproteome_cluster,cluster)] = external_hits


                if len(uniref_panproteome):
                    stats.Stats.pangenes_stats(uniref_panproteome, self.config)
                    pickle.dump(uniref_panproteome, open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], rank, panproteome_cluster, item.tax_id),'wb'))
        else:
            terminating.set()
            
    def create_panproteomes(self, cluster):
        if cluster == 100:
            ranks_to_process = self.ranks[-1::]
        elif cluster == 90:
            ranks_to_process = self.ranks[-1::]
        else:
            ranks_to_process = self.ranks
        
        if self.config['verbose']:
            utils.info('Starting creating panproteomes for {}...\n'.format(', '.join(ranks_to_process)))

        for k in ranks_to_process:
            os.makedirs('{}/{}/{}/{}'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], k, cluster), exist_ok=True)
        
        try:
            terminating_p = dummy.Event()
            with dummy.Pool(initializer=init_parse, initargs=(terminating_p, ), processes=100) as pool:
                d = [ _ for _ in pool.imap_unordered(self.process_panproteome, 
                    ((item, rank, cluster) for rank in ranks_to_process 
                                           for item in self.d_ranks[rank] 
                                           if self.taxontree.get_child_proteomes(item)),
                    chunksize=5000000)]
        except Exception as e:
            utils.error(str(e))
            raise

    @staticmethod
    def find_core_genes(panproteome):
        return Counter({gene:panproteome['members'][gene]['coreness'] for gene, _ in filter(lambda gene:gene[1]['coreness'] >= panproteome['coreness_threshold'], panproteome['members'].items())})

    def extract_protein_sequence(self, items):
        if not terminating.is_set():
            chunk, ids = items
            cluster = ids[0].split('_')[0]
            with open('{}/{}/{}_{}'.format(self.config['download_base_dir'], config['pickled_dir'], cluster, chunk), 'rb') as pickled_chunk:
                uniprot_chunk = pickle.load(pickled_chunk)
            entries = itemgetter(*ids)(uniprot_chunk)
            seqs = [itemgetter(0,2)(i) for i in entries]
            with open(outputdir,'w') as w_out:
                wout.writelines(['>{}\n{}\n'.format(s[0],s[1]) for s in seqs])
        else:
            terminating.set()


    @staticmethod
    def export_panproteome_fasta(panproteome):
        # with open('export/{}_uniref{}_panproteome.txt'.format(panproteome['tax_id'], panproteome['cluster'])) as export:
        pass
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


def generate_panproteomes(config):
    p = Panproteome(config)
    clusters = [int(x) for x in p.config['uniref_cluster_panproteomes'].split(',')]
    terminating = dummy.Event()
    nproc = config['nproc']
    
    #LOAD ONLY taxid_cluster map of wanted clustering    
    for c in p.get_upper_clusters(max(clusters)):
        p.__setattr__('uniref{}_tax_id_map'.format(c), pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref{}_taxid_idmap'.format(c)]),'rb')))
    
    os.makedirs('{}/{}/'.format(self.config['export_dir'], self.config['panproteomes_stats']), exist_ok=True)
    for cluster in clusters:
        p.create_panproteomes(cluster)

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['panproteomes']
    generate_panproteomes(config)