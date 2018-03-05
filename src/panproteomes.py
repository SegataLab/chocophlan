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
import multiprocessing as mp
from collections import Counter
from functools import partial
import copy

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
    coreness = 0.97
    
    def __init__(self, config):
        self.uniparc = {}
        try:
            if config['verbose']:
                utils.info('Loading pickled databases...')
            import gc
            gc.disable()
            # self.uniref100 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref100_idmap']),'rb'))
            # self.uniref90 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref90_idmap']),'rb'))
            # self.uniref50 = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniref50_idmap']),'rb'))

            # self.uniprotkb_uniref_idmap = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_uniprotkb_idmap']), 'rb'))
            for i in range(1,196):
                self.uniparc.update(pickle.load(open('data/pickled/uniparc_{}.pkl'.format(i),'rb')))
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
        pp = copy.deepcopy(clade.proteomes) if hasattr(clade,'proteomes') else []
        for c in clade.clades:
            pp.extend(self.get_child_proteomes(c))
        return pp
    
    def get_uniref_proteome(self, proteome_id ,cluster):
        cluster = 0 if cluster == 100 else (1 if cluster == 90 else 2)
        uniref_proteome = set()
        proteome = self.proteomes[proteome_id]['members']
        for protein in proteome:
            uniref_cluster = self.idmapping.get(protein, None)
            if uniref_cluster is None:
                urefs = list(set(uref for uref in (self.idmapping.get(upkb) for upkb in self.uniparc.get(protein,['']*7)[6]) if uref is not None))
                if len(urefs):
                    uniref_cluster = urefs[0][cluster]
            else:
                uniref_cluster = self.idmapping[protein][cluster]
            if uniref_cluster is not None:
                uniref_proteome.add(uniref_cluster)
        return list(uniref_proteome)

    def process_panproteome(self, item_rank):
        item, rank, cluster = item_rank
        cluster = 0 if cluster == 100 else (1 if cluster == 90 else 2)
        proteomes_to_process = self.get_child_proteomes(item)
        panproteome = Counter(entry for proteome_id in proteomes_to_process for entry in self.proteomes[proteome_id]['members'])
        uniref_panproteome = Counter()

        for protein, count in panproteome.items():
            uniref_cluster = self.idmapping.get(protein, None)
            if uniref_cluster is None:
                urefs = list(set(uref for uref in (self.idmapping.get(upkb) for upkb in self.uniparc.get(protein,['']*7)[6]) if uref is not None))
                if len(urefs):
                    uniref_cluster = urefs[0][cluster]
            else:
                uniref_cluster = uniref_cluster[cluster]
            if uniref_cluster is not None:
                if uniref_cluster not in uniref_panproteome:
                    uniref_panproteome[(uniref_cluster)] = [0,set()]
                uniref_panproteome[uniref_cluster][0] = uniref_panproteome[uniref_cluster][0] + count
                uniref_panproteome[uniref_cluster][1].add(protein)
                # uniref_panproteome.update([uniref_cluster])
        if len(panproteomes):
            ##ospathjoin
            pickle.dump(panproteomes, open('{}/{}/{}/{}/{}.pkl'.format(self.config['download_base_dir'], self.config['relpath_panproteomes_dir'], rank, cluster, item.tax_id),'wb'))

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
            for elem in elems:
                self.process_panproteome(elem)
        except Exception as e:
            utils.error(str(e))
            raise

        for cluster in uniref_proteomes.values():
            panproteome.update(cluster)

        coreness_threshold = round(len(uniref_proteomes) * self.coreness)
        core_genes = [gene for gene, _ in filter(lambda gene:gene[1] > coreness_threshold, panproteome)]
        accessory_genes = [gene for gene, _ in filter(lambda gene: gene[1] < coreness_threshold, panproteome)]

    # def coreness(self,panproteome):
        

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
