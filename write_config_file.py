#!/usr/bin/env python3
__author__ = ('Nicola Segata (nicola.segata@unitn.it),'
              'Francesco Beghini (francesco.beghini@unitn.it),'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '28 Sep 2017'


import src.utils as utils
import os
import argparse as ap
import configparser as cp
import sys


def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('--output', default='settings.cfg',
                   help='Path to config file')
    p.add_argument('--verbose', action='store_true', default=False,
                   help="Prints more stuff")
    p.add_argument('--overwrite', action='store_true', default=False,
                   help=("This flag needs to be set in order to overwrite an "
                         "existing config file"))

    group = p.add_argument_group(title="Download",
                                 description=("Parameters for setting download "
                                              "options"))
    group.add_argument('--refseq_ftp_base', default='ftp.ncbi.nlm.nih.gov')
    group.add_argument('--refseq_bacterial_genomes',
                       default='/refseq/release/bacteria')
    group.add_argument('--refseq_taxonomic_catalogue',
                       default=('/refseq/release/release-catalog/RefSeq-releas'
                                'e84.catalog.gz'))
    group.add_argument('--refseq_taxdump',
                       default=('/pub/taxonomy/taxdump.tar.gz'))

    group.add_argument('--uniprot_ftp_base', default='ftp.uniprot.org')
    group.add_argument('--uniprot_uniref100',
                       default=('/pub/databases/uniprot/current_release/uniref'
                                '/uniref100/uniref100.xml.gz'))
    group.add_argument('--uniprot_reference_proteomes',
                       default=('/pub/databases/uniprot/current_release'
                                '/knowledgebase/reference_proteomes'))
    group.add_argument('--uniprot_pan_proteomes',
                       default=('/pub/databases/uniprot/current_release'
                           '/knowledgebase/pan_proteomes/'))
    group.add_argument('--download_base_dir', default='data/',
                       help='Base directory for raw files to be downloaded to')

    group.add_argument('--relpath_bacterial_genomes',
                       default='/refseq/genomes',
                       help='Directory for genome files')
    group.add_argument('--relpath_taxonomic_catalogue',
                       default='/refseq/catalogue/refseq_catalogue.gz',
                       help='Directory for refseq catalogue file')
    group.add_argument('--relpath_taxdump',
                       default='/refseq/taxdump/refseq_taxdump.tar.gz',
                       help='')
    group.add_argument('--relpath_uniref100',
                       default='/uniprot/uniref/uniref100.xml.gz',
                       help='Directory for uniref100 file')
    group.add_argument('--relpath_reference_proteomes',
                       default='/uniprot/reference_proteomes',
                       help='Directory for the reference proteomes file')
    group.add_argument('--relpath_pan_proteomes',
                       default='/uniprot/pan_proteomes')
    group.add_argument('--relpath_pickle_taxid_contigid',
                       default='/pickled/taxid_contig.pkl',
                       help='')
    group.add_argument('--relpath_pickle_taxid_taxonomy',
                       default='/pickled/taxid_taxonomy.pkl',
                       help='')
    group.add_argument('--relpath_pickle_contigid_filename',
                       default='/pickled/contigid_filename.pkl',
                       help='')
    group.add_argument('--relpath_pickle_proteomes',
                       default='/pickled/proteomes.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniprotdb',
                       default='/pickled/uniprotdb.pkl',
                       help='')
    group.add_argument('--nproc', default=20,
                       help='Number of parallel processes')

    return p.parse_args()


def set_download_options(configparser_object, args, verbose=False):
    configparser_object.add_section('download')
    configparser_object.set('download', 'refseq_ftp_base',
                            args.refseq_ftp_base)
    configparser_object.set('download', 'refseq_bacterial_genomes',
                            args.refseq_bacterial_genomes)
    configparser_object.set('download', 'refseq_taxonomic_catalogue',
                            args.refseq_taxonomic_catalogue)
    configparser_object.set('download', 'refseq_taxdump',
                            args.refseq_taxdump)
    configparser_object.set('download', 'uniprot_ftp_base',
                            args.uniprot_ftp_base)
    configparser_object.set('download', 'uniprot_uniref100',
                            args.uniprot_uniref100)
    configparser_object.set('download', 'uniprot_reference_proteomes',
                            args.uniprot_reference_proteomes)
    configparser_object.set('download', 'uniprot_pan_proteomes',
                            args.uniprot_pan_proteomes)
    configparser_object.set('download', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('download', 'relpath_bacterial_genomes',
                            args.relpath_bacterial_genomes)
    configparser_object.set('download', 'relpath_taxonomic_catalogue',
                            args.relpath_taxonomic_catalogue)
    configparser_object.set('download', 'relpath_taxdump',
                            args.relpath_taxdump)
    configparser_object.set('download', 'relpath_uniref100',
                            args.relpath_uniref100)
    configparser_object.set('download', 'relpath_reference_proteomes',
                            args.relpath_reference_proteomes)
    configparser_object.set('download', 'relpath_pan_proteomes',
                            args.relpath_pan_proteomes)
    configparser_object.set('download', 'verbose', str(verbose))
    configparser_object.set('download', 'nproc', str(args.nproc))

    configparser_object.add_section('extract')
    configparser_object.set('extract', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('extract', 'relpath_taxdump',
                            args.relpath_taxdump)
    configparser_object.set('extract', 'relpath_taxonomic_catalogue',
                            args.relpath_taxonomic_catalogue)
    configparser_object.set('extract', 'relpath_bacterial_genomes',
                            args.relpath_bacterial_genomes)

    configparser_object.set('extract', 'relpath_pickle_taxid_contigid',
                            args.relpath_pickle_taxid_contigid)
    configparser_object.set('extract', 'relpath_pickle_taxid_taxonomy',
                            args.relpath_pickle_taxid_taxonomy)
    configparser_object.set('extract', 'relpath_pickle_contigid_filename',
                            args.relpath_pickle_contigid_filename)

    configparser_object.set('extract', 'verbose', str(verbose))
    configparser_object.set('extract', 'nproc', str(args.nproc))
    
    configparser_object.add_section('process_proteomes')
    configparser_object.set('process_proteomes', 'uniprot_ftp_base',
                            args.uniprot_ftp_base)
    configparser_object.set('process_proteomes', 'uniprot_reference_proteomes', 
                            args.uniprot_reference_proteomes)
    configparser_object.set('process_proteomes', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('process_proteomes', 'relpath_reference_proteomes',
                            args.relpath_reference_proteomes)
    configparser_object.set('process_proteomes', 'relpath_pickle_proteomes',
                            args.relpath_pickle_proteomes)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniprotdb',
                            args.relpath_pickle_uniprotdb)
    configparser_object.set('process_proteomes', 'verbose', str(verbose))
    configparser_object.set('process_proteomes', 'nproc', str(args.nproc))
    return configparser_object


if __name__ == '__main__':
    args = read_params()
    utils.check_config_params(args, verbose=args.verbose)

    config = cp.ConfigParser()
    config = set_download_options(config, args, verbose=args.verbose)

    if os.path.isfile(args.output) and (not args.overwrite):
        sys.exit('Output file "{}" will NOT be overwritten. Exiting...'
                 .format(args.output))
        
    with open(args.output, 'w') as f:
        config.write(f)

    sys.exit(0)
