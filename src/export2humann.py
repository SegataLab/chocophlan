#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Fabio Cumbo (fabio.cumbo@unitn.it')

__date__ = '14 Dec 2018'


from _version import *
import os
import argparse as ap
import configparser as cp
import pickle
import glob
import time
import pandas as pd
from operator import itemgetter

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

class chocophlan2humann2:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')

        self.config = config
        self.exportpath = '{}/{}'.format(self.config['export_dir'], self.config['exportpath_phylophlan'])
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)
        
        all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniprotkb/'.format(config['download_base_dir'], config['pickled_dir'])))
        for i in all_uprot_chunks:
            uniprot_chunk = pickle.load(open('{}/{}/uniprotkb/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniprot.update({k : [v[3],v[4],v[1]] for k,v in uniprot_chunk.items()})
        self.uniref90_taxid_map = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_uniref90_taxid_idmap']), 'rb'))
        if config['verbose']:
            utils.info('Finished.\n')

    def get_all_species_panproteomes(self):
        return [x.split('/')[-1][:-4] for x in glob.glob('{}/{}/species/90/*'.format(config['download_base_dir'], config['relpath_panproteomes_dir']))]

    def export_panproteome(self, species_id):
        panproteome = pd.DataFrame.from_csv('{}/{}/{}.txt.bz2'.format(self.config['export_dir'], self.config['panproteomes_stats'], species_id), sep = '\t')
        for pangene in panproteome.index:
            pangene = pangene.split('_')[1]
            if pangene in self.uniref90_taxid_map:
                for upkb, ncbi_tax_id, isRef in self.uniref90_taxid_map[pangene]
                    if str(ncbi_tax_id) == species_id:
                        uniprotkb_entries.append(upkb)
                for upkb in filter(re.compile('^(?!UPI)').match, uniprotkb_entries):
                    
                else:
                    



if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['export']

    c = chocophlan2humann2(config)
    c.chocophlan2humann2()
    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
