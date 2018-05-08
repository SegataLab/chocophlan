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
import statistics
if __name__ == '__main__':
    import utils
else:
    import src.utils

class stats:
    def __init__(config):
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
        species_proteomes = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[tax_id])
        number_proteomes = len(species_proteomes)
        reference_proteomes = len([x for x in species_proteomes if self.proteomes[x]['isReference']])
        redundant_proteomes = len([x for x in species_proteomes if self.proteomes[x]['upi']])
        non_redundant_proteomes = number_proteomes - redundant_proteomes - reference_proteomes
        avg_number_proteins = statistics.mean([len(self.proteomes[proteome]['members']) for proteome in species_proteomes])
        panproteome_100 = self.panproteome_stats(tax_id, 100)
        panproteome_90 = self.panproteome_stats(tax_id, 90)
        panproteome_50 = self.panproteome_stats(tax_id, 50)

        return { tax_id: (reference_proteomes, redundant_proteomes, non_redundant_proteomes, panproteome_100, panproteome_90, panproteome_50) }

    def panproteome_stats(self, tax_id, cluster):
      try:
        panproteome = pickle.load(open('{}{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], 'species', cluster, tax_id),'rb'))
      except FileNotFoundError as ex:
        utils.error('TAXID {}: Panproteome calculated using UniRef{} does not exist!'.format(tax_id, cluster))
        return None
      number_proteomes = panproteome['number_proteomes']
      number_members = len(panproteome['members'])

      #coreness
      mean_coreness = statistics.mean((p['coreness'] for p in panproteome['members'].values()))
      median_coreness = statistics.median((p['coreness'] for p in panproteome['members'].values()))
      stdev_coreness = statistics.stdev((p['coreness'] for p in panproteome['members'].values()))

      #uniqueness
      d = {}
      [d.setdefault(k,[]).append(v) for k,v in ((k,v) for p in panproteome['members'].values() for k, v in p['uniqueness'].items())]
      mean_uniqueness = {k:statistics.mean(v) for k,v in d.items()}
      median_uniqueness = {k:statistics.median(v) for k,v in d.items()}
      stdev_uniqueness = {k:statistics.stdev(v) for k,v in d.items()}

      #copy_number
      mean_copy_number = statistics.mean((len(p['copy_number']) for p in panproteome['members'].values()))
      median_copy_number = statistics.median((len(p['copy_number']) for p in panproteome['members'].values()))
      stdev_copy_number = statistics.stdev((len(p['copy_number']) for p in panproteome['members'].values()))

      #proteomes_present
      mean_proteomes_present = statistics.mean((len(p['proteomes_present']) for p in panproteome['members'].values()))
      median_proteomes_present = statistics.median((len(p['proteomes_present']) for p in panproteome['members'].values()))
      stdev_proteomes_present = statistics.stdev((len(p['proteomes_present']) for p in panproteome['members'].values()))

      return { 'number_proteomes' : number_proteomes,
               'number_members' : number_members,
               'coreness': (mean_coreness, median_coreness, stdev_coreness), 
               'uniqueness': (mean_uniqueness, median_uniqueness, stdev_uniqueness),
               'copy_number': (mean_copy_number, median_copy_number, stdev_copy_number),
               'proteomes_present':(mean_proteomes_present, median_proteomes_present, stdev_proteomes_present) 
              }

    def stats(config):
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
              '-----UniRef-----\n'
              'Total number of UniRef100 clusters: {}\n'
              '\t\tUniref90 clusters: {}\n'
              '\t\tUniref50 clusters: {}\n'
              '-----Panproteomes-----\n'
              'Total created panproteomes: {}\n'
              '-----Memory-----------\n'
              'Total RAM usage by ChocoPhlAn indices: {} Gb\n'
              '\n'
              .format(#Taxonomy
                      len(self.d_ranks['species']),
                      len([x for x in d_ranks['species'] if d_ranks['superkingdom'][0].is_parent_of(x)]),
                      len([x for x in d_ranks['species'] if d_ranks['superkingdom'][1].is_parent_of(x)]),
                      sum([len(x) for x in d_ranks['species'] if not x.initially_terminal]),
                      #UniProtKB
                      len(proteomes),
                      sum([len(v['members']) for _,v in proteomes.items()]),
                      len(uniprotkb)-sum([len(v['members']) for _,v in proteomes.items()]),
                      len(uniprotkb),
                      #UniRef
                      len(uniref100),
                      len(uniref90),
                      len(uniref50),
                      #Panproteomes
                      len(glob.glob(config['download_base_dir']+config['relpath_panproteomes_dir']+'/*/*/*')),
                      #Memory
                      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / pow(2,20))

        s = { 'Total number of species'
              'Total number of bacterial species'
              'Total number of archaeal species'
              'Total number of eukaryotic species'
              'Total number of strains'
              'Total number of proteomes'
              'Total number of proteins in proteomes'
              'Total number of proteins NOT PRESENT in a proteome'
              'Total number of proteins in UniProtKB database'
              'Total number of UniRef100 clusters'
              'Total number of Uniref90 clusters'
              'Total number of Uniref50 clusters'
              'Total created panproteomes'
              'Total RAM usage by ChocoPhlAn indices (Gb)'

        }

        species_proteomes = set(x['tax_id'] if taxontree.taxid_n[x['tax_id']].rank =='species' else taxontree.go_up_to_species(x['tax_id']) for x in proteomes.values() )
        for species in species_proteomes:
            species_stats(self, species)


if __name__ == '__main__':
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']            #UPDATE TO STATS

    stats(config)

    sys.exit(0)
