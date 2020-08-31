#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '31 Aug 2020'

from chocophlan import *
import chocophlan.utils as utils
import chocophlan.build_taxontree as build_taxontree
import chocophlan.download as download
import chocophlan.parse_uniprot as parse_uniprot
import chocophlan.panproteomes as panproteomes
import chocophlan.export_to_phylophlan as export_to_phylophlan
import chocophlan.stats as stats
import chocophlan.export_to_metaphlan as export_to_metaphlan

def chocophlan():
    logger = utils.setup_logger('.','CHOCOPhlAn_{}'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')))
    
    logger.info('CHOCOPhlAn: Cluster of HOmologous Cdses fOr PHyLogenetic ANalysis\n')
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)

    if not os.path.exists('{}/DONE'.format(config['download']['download_base_dir'])):
        download.download(config['download'], logger, verbose=config['download']['verbose'])
    else:
        logger.info('UniRef and NCBI taxonomy data already downloaded. Moving on.')
    if os.path.exists('{}/DONE'.format(config['download']['download_base_dir'])):
        build_taxontree.do_extraction(config['build_taxontree'], logger, verbose=config['build_taxontree']['verbose'])
    
    if os.path.exists('{}/{}'.format(config['parse_uniprot']['download_base_dir'], config['parse_uniprot']['relpath_pickle_taxontree'])):
        parse_uniprot.parse_uniprot(config['parse_uniprot'])
    
    if os.path.exists('{}/{}'.format(config['parse_uniprot']['download_base_dir'], config['parse_uniprot']['relpath_pickle_taxontree'])) and 
        os.path.exists('{}/{}'.format(config['parse_uniprot']['download_base_dir'], config['parse_uniprot']['relpath_pickle_proteomes'])):
        download.download_ncbi_from_proteome_pickle(config['parse_uniprot'])
    # panproteomes.generate_panproteomes(config['panproteomes'])
    #stats.generate_stats(config['stats'])
    #export_to_metaphlan.run_all(config['export'])
    # chocophlan2phylophlan.export_to_phylophlan(config['export'])

if __name__ == '__main__':
    t0 = time.time()
    chocophlan()
    t1 = time.time()
    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
