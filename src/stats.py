#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')

from _version import __CHOCOPhlAn_version__
__date__ = '03 Jan 2018'

import os
import sys
import argparse as ap
import configparser as cp
import pickle
import resource
import glob
import csv
import statistics
import multiprocessing.dummy as dummy
import pandas as pd
import importlib
if __name__ == "__main__":
    import utils
    import panproteomes
else:
    utils = importlib.import_module('src.utils')
    panproteomes = importlib.import_module('src.panproteomes')

def init_parse(terminating_):
    global terminating
    terminating = terminating_
class Stats:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')
        self.proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
        self.uniprotkb = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_idmap']), 'rb'))
        self.uniref100 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref100_idmap']), 'rb'))
        self.uniref90 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref90_idmap']), 'rb'))
        self.uniref50 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref50_idmap']), 'rb'))
        self.taxontree = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.d_ranks = self.taxontree.lookup_by_rank()
        self.config = config

    #species_proteomes = set(x['tax_id'] if taxontree.taxid_n[x['tax_id']].rank =='species' else taxontree.go_up_to_species(x['tax_id']) for x in proteomes.values() )
    def species_stats(self, tax_id):
        if not terminating.is_set():
            species_proteomes = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[tax_id])
            number_proteomes = len(species_proteomes)
            reference_proteomes = len([x for x in species_proteomes if self.proteomes[x]['isReference']])
            redundant_proteomes = len([x for x in species_proteomes if self.proteomes[x]['upi']])
            non_redundant_proteomes = number_proteomes - redundant_proteomes - reference_proteomes
            members_panproteome = [len(self.proteomes[proteome]['members']) for proteome in species_proteomes]
            avg_number_proteins = statistics.mean(members_panproteome)
            median_number_proteins = statistics.median(members_panproteome)
            stdev_number_proteins = statistics.stdev(members_panproteome) if len(members_panproteome) >2 else None
            try:
                panproteome_100 = self.panproteome_stats(tax_id, 100)
                panproteome_90 = self.panproteome_stats(tax_id, 90)
                panproteome_50 = self.panproteome_stats(tax_id, 50)
            except Exception as e:
                utils.error(e)
                terminating.set()
                raise

            return { tax_id: (self.taxontree.print_full_taxonomy(tax_id), number_proteomes, reference_proteomes, redundant_proteomes, non_redundant_proteomes, avg_number_proteins, panproteome_100, panproteome_90, panproteome_50) }
        else:
            terminating.set()
    def panproteome_stats(self, tax_id, cluster):
        try:
            panproteome = pickle.load(open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], 'species', cluster, tax_id),'rb'))
        except FileNotFoundError as ex:
            # utils.error('TAXID {}: Panproteome calculated using UniRef{} does not exist!'.format(tax_id, cluster))
            return {}

        number_members = len(panproteome['members'])
        number_core_proteins = len(panproteomes.Panproteome.find_core_genes(panproteome))

        #coreness
        iter_coreness = [p['coreness'] for p in panproteome['members'].values()]
        core_coreness = [panproteome['members'][p]['coreness'] for p in panproteomes.Panproteome.find_core_genes(panproteome)]
        mean_core_coreness = statistics.mean(core_coreness) if len(core_coreness) else 0
        mean_coreness = statistics.mean(iter_coreness) if len(iter_coreness) else None
        median_coreness = statistics.median(iter_coreness) if len(iter_coreness) else None
        stdev_coreness = statistics.stdev(iter_coreness) if len(iter_coreness) else None

        #uniqueness
        d = {}
        [d.setdefault(k,[]).append(v) for k,v in ((k,v) for p in panproteome['members'].values() for k, v in p['uniqueness'].items())]
        mean_uniqueness = {k:statistics.mean(v) if len(v) else None for k,v in d.items()}
        median_uniqueness = {k:statistics.median(v) if len(v) else None  for k,v in d.items()}
        stdev_uniqueness = {k:statistics.stdev(v) if len(v) else None  for k,v in d.items()}

        #copy_number
        iter_copy_number = [p['copy_number'].values() for p in panproteome['members'].values()]
        mean_copy_number = statistics.mean((statistics.mean(x) for x in iter_copy_number)) if len(iter_copy_number) else None
        median_copy_number = statistics.median((statistics.median(x) for x in iter_copy_number)) if len(iter_copy_number) else None
        ssd = [statistics.stdev(x) for x in iter_copy_number if len(x)>1]
        stdev_copy_number = statistics.stdev(ssd) if len(ssd) > 1 else 0

        #proteomes_present
        iter_proteomes_present = [len(p['proteomes_present']) for p in panproteome['members'].values()]
        mean_proteomes_present = statistics.mean(iter_proteomes_present) if len(iter_proteomes_present) else None
        median_proteomes_present = statistics.median(iter_proteomes_present) if len(iter_proteomes_present) else None
        stdev_proteomes_present = statistics.stdev(iter_proteomes_present) if len(iter_proteomes_present) else None

        ret_d = { 'panproteome_{}_number_members'.format(cluster) : number_members,
               'panproteome_{}_number_core_proteins'.format(cluster) : number_core_proteins,
               'panproteome_{}_mean_cores_coreness'.format(cluster): mean_core_coreness, 
               'panproteome_{}_mean_coreness'.format(cluster): mean_coreness, 
               'panproteome_{}_median_coreness'.format(cluster): median_coreness,
               'panproteome_{}_stdev_coreness'.format(cluster): stdev_coreness,
               'panproteome_{}_mean_copy_number'.format(cluster) : mean_copy_number,
               'panproteome_{}_median_copy_number'.format(cluster) : median_copy_number,
               'panproteome_{}_stdev_copy_number'.format(cluster) : stdev_copy_number,
               'panproteome_{}_mean_proteomes_present'.format(cluster) : mean_proteomes_present, 
               'panproteome_{}_median_proteomes_present'.format(cluster) : median_proteomes_present, 
               'panproteome_{}_stdev_proteomes_present'.format(cluster) : stdev_proteomes_present
              }

        for cc in d:
            ret_d['panproteome_{}_mean_{}_uniqueness'.format(cluster, cc)] = mean_uniqueness[cc]
            ret_d['panproteome_{}_median_{}_uniqueness'.format(cluster, cc)] = median_uniqueness[cc]
            ret_d['panproteome_{}_stdev_{}_uniqueness'.format(cluster, cc)] = stdev_uniqueness[cc]

        return ret_d

    ##TODO: Create panproteme stats folder
    @staticmethod
    def pangenes_stats(panproteome, config):
        species = panproteome['tax_id']
        res = {'{}_{}'.format(species,pangene) : {'tax_id':species,
                   'coreness' : panproteome['members'][pangene]['coreness'],
                   'coreness_perc' : panproteome['members'][pangene]['coreness']/panproteome['number_proteomes'],
                   'external_hits_50' : ';'.join(str(x) for x in panproteome['members'][pangene]['external_hits']['90_50']),
                   'external_hits_90' : ';'.join(str(x) for x in panproteome['members'][pangene]['external_hits']['90_90']),
                   'uniqueness_90' : panproteome['members'][pangene]['uniqueness']['90_90'],
                   'uniqueness_50' : panproteome['members'][pangene]['uniqueness']['90_50'],
                   'uniqueness_nosp_90' : panproteome['members'][pangene]['uniqueness_nosp']['90_90'],
                   'uniqueness_nosp_50' : panproteome['members'][pangene]['uniqueness_nosp']['90_50']
                  } for pangene, v in panproteome['members'].items() if len(pangene)}
        res_df = pd.DataFrame(res).T.sort_values('coreness', ascending=False)
        res_df.to_csv('{}/{}/{}.txt'.format(config['export_dir'], config['panproteomes_stats'], species), sep='\t')


    def stats(self):
        all_stat = { 'Total number of species' : sum(1 for b in self.d_ranks['species'] if self.taxontree.get_child_proteomes(b)),
              'Total number of bacterial species' : sum(1 for b in self.taxontree.taxid_n[2].find_clades() if b.rank == 'species' and self.taxontree.get_child_proteomes(b)),
              'Total number of archaeal species' : sum(1 for b in self.taxontree.taxid_n[2157].find_clades() if b.rank == 'species' and self.taxontree.get_child_proteomes(b)),
              'Total number of eukaryotic species' : sum(1 for b in self.taxontree.taxid_n[2759].find_clades() if b.rank == 'species' and self.taxontree.get_child_proteomes(b)),
              'Total number of strains' : sum(len(b.get_terminals())-1 for b in self.d_ranks['species'] if self.taxontree.get_child_proteomes(b)),
              'Total number of proteomes' : len(self.proteomes),
              'Total number of proteins in proteomes' : len(set(x for _,v in self.proteomes.items() for x in v['members'])),
              'Total number of proteins NOT PRESENT in a proteome' : len(self.uniprotkb)-len(set(x for _,v in self.proteomes.items() for x in v['members'])),
              'Total number of proteins in UniProtKB database' : len(self.uniprotkb),
              'Total number of UniRef100 clusters' : len(self.uniref100),
              'Total number of Uniref90 clusters' : len(self.uniref90),
              'Total number of Uniref50 clusters' : len(self.uniref50),
              'Total created panproteomes' : len(glob.glob(self.config['download_base_dir']+self.config['relpath_panproteomes_dir']+'/*/*/*')),
              'Total RAM usage by ChocoPhlAn indices (Gb)' : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / pow(2,20)
        }

        with open('{}/stats.csv'.format(self.config['export_dir']), 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for k,v in all_stat.items():
                writer.writerow([k,v])

        species_proteomes = set(b.tax_id for b in self.d_ranks['species'] if self.taxontree.get_child_proteomes(b))
        terminating = dummy.Event()
        with dummy.Pool(initializer=init_parse, initargs=(terminating, ), processes=self.config['nproc']) as pool:
            d = {k:v for x in pool.imap_unordered(self.species_stats, species_proteomes, chunksize=20) for k,v in x.items()}

        d_new = {}
        for k,v in d.items():
            v_=dict(zip(['tax_id', 'name', 'number_proteomes', 'number_reference_proteomes', 'number_redundant_proteomes', 'number_non_redundant_proteomes', 'avg_number_proteins'], [k]+list(v[0:6])))
            [v_.update(y) for y in v[6:9]]
            d_new[k] = v_

        with open('{}/panproteomes_stats.csv'.format(self.config['export_dir']), 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=list(d_new[562].keys()), delimiter='\t')
            writer.writeheader()
            for k,v in d_new.items():
                writer.writerow(v)

def generate_stats(config):
    s = Stats(config)
    s.stats()
    
if __name__ == '__main__':
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['stats']            #UPDATE TO STATS

    s = Stats(config)
    s.stats()