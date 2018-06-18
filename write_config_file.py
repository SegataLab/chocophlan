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
    group.add_argument('--refseq_taxonomic_catalogue',
                       default=('/refseq/release/release-catalog/'))
    group.add_argument('--refseq_taxdump',
                       default=('/pub/taxonomy/taxdump.tar.gz'))

    group.add_argument('--uniprot_ftp_base', default='ftp.uniprot.org')
    group.add_argument('--ebi_ftp_base', default='ftp.ebi.ac.uk')
    group.add_argument('--uniprot_protomes',
                       default=('/pub/contrib/UniProtKB/'
                               'Proteome_DS/proteome_ds.xml.gz'))
    group.add_argument('--uniprot_uniref100',
                       default=('/pub/databases/uniprot/current_release/uniref'
                                '/uniref100/uniref100.xml.gz'))
    group.add_argument('--uniprot_uniref90',
                       default=('/pub/databases/uniprot/current_release/uniref'
                                '/uniref90/uniref90.xml.gz'))
    group.add_argument('--uniprot_uniref50',
                       default=('/pub/databases/uniprot/current_release/uniref'
                                '/uniref50/uniref50.xml.gz'))
    group.add_argument('--uniprot_reference_proteomes',
                       default=('/pub/databases/uniprot/current_release'
                                '/knowledgebase/reference_proteomes'))
    group.add_argument('--uniprot_sprot',
                       default=('/pub/databases/uniprot/current_release'
                           '/knowledgebase/complete/uniprot_sprot.xml.gz'))
    group.add_argument('--uniprot_trembl',
                       default=('/pub/databases/uniprot/current_release'
                           '/knowledgebase/complete/uniprot_trembl.xml.gz'))
    group.add_argument('--uniparc',
                       default=('/pub/databases/uniprot/current_release'
                           '/uniparc/uniparc_all.xml.gz'))
    group.add_argument('--uniprot_idmapping',
                       default=('/pub/databases/uniprot/current_release'
                           '/knowledgebase/idmapping/idmapping_selected.tab.gz'))
    group.add_argument('--download_base_dir', default='data/',
                       help='Base directory for raw files to be downloaded to')
    
    group.add_argument('--relpath_proteomes_xml',
                       default='/uniprot/proteome_ds.xml.gz',
                       help='')
    group.add_argument('--relpath_genomes',
                       default='/ncbi',
                       help='Directory for genome files')
    group.add_argument('--relpath_taxonomic_catalogue',
                       default='/refseq/catalogue/',
                       help='Directory for refseq catalogue file')
    group.add_argument('--relpath_taxdump',
                       default='/refseq/taxdump/refseq_taxdump.tar.gz',
                       help='')
    group.add_argument('--relpath_uniref100',
                       default='/uniprot/uniref/uniref100.xml.gz',
                       help='Directory for uniref100 file')
    group.add_argument('--relpath_uniref90',
                       default='/uniprot/uniref/uniref90.xml.gz',
                       help='Directory for uniref90 file')
    group.add_argument('--relpath_uniref50',
                       default='/uniprot/uniref/uniref50.xml.gz',
                       help='Directory for uniref50 file')
    group.add_argument('--relpath_uniprot_sprot',
                       default='/uniprot/complete/uniprot_sprot.xml.gz',
                       help='')
    group.add_argument('--relpath_uniprot_trembl',
                       default='/uniprot/complete/uniprot_trembl.xml.gz',
                       help='')
    group.add_argument('--relpath_uniparc',
                       default='/uniprot/uniparc/uniparc_all.xml.gz',
                       help='')
    group.add_argument('--relpath_reference_proteomes',
                       default='/uniprot/reference_proteomes',
                       help='Directory for the reference proteomes file')
    group.add_argument('--relpath_taxon_to_process',
                       default='taxon_to_process',
                       help='')
    
    ### PICKLE FILES
    group.add_argument('--relpath_panproteomes_dir',
                       default='/pickled/panproteomes',
                       help='')
    group.add_argument('--relpath_idmapping',
                       default='/uniprot/idmapping_selected.tab.gz',
                       help='')
    group.add_argument('--relpath_pickle_taxid_contigid',
                       default='/pickled/taxid_contig.pkl',
                       help='')
    group.add_argument('--relpath_pickle_taxid_taxonomy',
                       default='/pickled/taxid_taxonomy.pkl',
                       help='')
    group.add_argument('--relpath_pickle_taxontree',
                       default='/pickled/taxontree.pkl',
                       help='')
    group.add_argument('--relpath_pickle_proteomes',
                       default='/pickled/proteomes.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniprotkb_idmap',
                       default='/pickled/uniprotkb_idmap.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniprotkb_uniref_idmap',
                       default='/pickled/uniprotkb_uniref_idmap.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref100_idmap',
                       default='/pickled/uniref100_idmap.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref90_idmap',
                       default='/pickled/uniref90_idmap.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref50_idmap',
                       default='/pickled/uniref50_idmap.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref100_taxid_idmap',
                       default='/pickled/uniref100_taxid_map.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref90_taxid_idmap',
                       default='/pickled/uniref90_taxid_map.pkl',
                       help='')
    group.add_argument('--relpath_pickle_uniref50_taxid_idmap',
                       default='/pickled/uniref50_taxid_map.pkl',
                       help='')
   
    group.add_argument('--relpath_pickle_contigid_filename',
                       default='/pickled/contigid_filename.pkl',
                       help='')

    group.add_argument('--uniref_cluster_panproteomes',
                      default='90')

    group.add_argument('--discard_low_quality_genomes',
                      action='store_true',
                      default=True)

    group.add_argument('--export_dir',
                       default='/export',
                       help='')
    group.add_argument('--exportpath_phylophlan',
                       default='phylophlan',
                       help='')

    group.add_argument('--nproc', default=7,
                       help='Number of parallel processes')

    return p.parse_args()


def set_download_options(configparser_object, args, verbose=False):
    configparser_object.add_section('download')
    configparser_object.set('download', 'refseq_ftp_base',
                            args.refseq_ftp_base)
    configparser_object.set('download', 'refseq_taxonomic_catalogue',
                            args.refseq_taxonomic_catalogue)
    configparser_object.set('download', 'refseq_taxdump',
                            args.refseq_taxdump)
    configparser_object.set('download', 'ebi_ftp_base',
                            args.ebi_ftp_base)
    configparser_object.set('download', 'uniprot_protomes',
                            args.uniprot_protomes)
    configparser_object.set('download', 'uniprot_ftp_base',
                            args.uniprot_ftp_base)
    configparser_object.set('download', 'uniprot_uniref100',
                            args.uniprot_uniref100)
    configparser_object.set('download', 'uniprot_uniref90',
                            args.uniprot_uniref90)
    configparser_object.set('download', 'uniprot_uniref50',
                            args.uniprot_uniref50)
    configparser_object.set('download', 'uniprot_reference_proteomes',
                            args.uniprot_reference_proteomes)
    configparser_object.set('download', 'uniprot_idmapping',
                            args.uniprot_idmapping)
    configparser_object.set('download', 'uniprot_sprot',
                            args.uniprot_sprot)
    configparser_object.set('download', 'uniprot_trembl',
                            args.uniprot_trembl)
    configparser_object.set('download', 'uniparc',
                            args.uniparc)
    configparser_object.set('download', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('download', 'relpath_genomes',
                            args.relpath_genomes)
    configparser_object.set('download', 'relpath_taxonomic_catalogue',
                            args.relpath_taxonomic_catalogue)
    configparser_object.set('download', 'relpath_taxdump',
                            args.relpath_taxdump)
    configparser_object.set('download', 'relpath_uniref100',
                            args.relpath_uniref100)
    configparser_object.set('download', 'relpath_uniref90',
                            args.relpath_uniref90)
    configparser_object.set('download', 'relpath_uniref50',
                            args.relpath_uniref50)
    configparser_object.set('download', 'relpath_reference_proteomes',
                            args.relpath_reference_proteomes)
    configparser_object.set('download','relpath_proteomes_xml',
                            args.relpath_proteomes_xml)
    configparser_object.set('download', 'relpath_idmapping',
                            args.relpath_idmapping)
    configparser_object.set('download', 'relpath_uniprot_sprot',
                            args.relpath_uniprot_sprot)
    configparser_object.set('download', 'relpath_uniprot_trembl',
                            args.relpath_uniprot_trembl)
    configparser_object.set('download', 'relpath_uniparc',
                            args.relpath_uniparc)
    configparser_object.set('download', 'verbose', str(verbose))
    configparser_object.set('download', 'nproc', str(args.nproc))

    configparser_object.add_section('extract')
    configparser_object.set('extract', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('extract', 'relpath_taxdump',
                            args.relpath_taxdump)
    configparser_object.set('extract', 'relpath_taxonomic_catalogue',
                            args.relpath_taxonomic_catalogue)
    configparser_object.set('extract', 'relpath_pickle_taxid_contigid',
                            args.relpath_pickle_taxid_contigid)
    configparser_object.set('extract', 'relpath_pickle_taxid_taxonomy',
                            args.relpath_pickle_taxid_taxonomy)
    configparser_object.set('extract', 'relpath_pickle_contigid_filename',
                            args.relpath_pickle_contigid_filename)
    configparser_object.set('extract', 'relpath_pickle_taxontree',
                            args.relpath_pickle_taxontree)
    configparser_object.set('extract', 'verbose', str(verbose))
    configparser_object.set('extract', 'nproc', str(args.nproc))
    
    configparser_object.add_section('process_proteomes')
    configparser_object.set('process_proteomes', 'uniprot_ftp_base',
                            args.uniprot_ftp_base)
    configparser_object.set('process_proteomes', 'uniprot_reference_proteomes', 
                            args.uniprot_reference_proteomes)
    configparser_object.set('process_proteomes', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('process_proteomes', 'relpath_genomes',
                            args.relpath_genomes)
    configparser_object.set('process_proteomes', 'relpath_reference_proteomes',
                            args.relpath_reference_proteomes)
    configparser_object.set('process_proteomes', 'relpath_idmapping',
                            args.relpath_idmapping)
    configparser_object.set('process_proteomes', 'relpath_uniprot_sprot',
                            args.relpath_uniprot_sprot)
    configparser_object.set('process_proteomes', 'relpath_uniprot_trembl',
                            args.relpath_uniprot_trembl)
    configparser_object.set('process_proteomes', 'relpath_uniref100',
                            args.relpath_uniref100)
    configparser_object.set('process_proteomes', 'relpath_uniref90',
                            args.relpath_uniref90)
    configparser_object.set('process_proteomes', 'relpath_uniref50',
                            args.relpath_uniref50)
    configparser_object.set('process_proteomes', 'relpath_uniparc',
                            args.relpath_uniparc)
    configparser_object.set('process_proteomes','relpath_proteomes_xml',
                            args.relpath_proteomes_xml)
    configparser_object.set('process_proteomes', 'relpath_pickle_proteomes',
                            args.relpath_pickle_proteomes)
    configparser_object.set('process_proteomes', 'relpath_pickle_taxontree',
                            args.relpath_pickle_taxontree)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniprotkb_uniref_idmap',
                            args.relpath_pickle_uniprotkb_uniref_idmap)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniref100_idmap',
                            args.relpath_pickle_uniref100_idmap)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniref90_idmap',
                            args.relpath_pickle_uniref90_idmap)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniref50_idmap',
                            args.relpath_pickle_uniref50_idmap)

    configparser_object.set('process_proteomes', 'relpath_pickle_uniref100_taxid_idmap',
                            args.relpath_pickle_uniref100_taxid_idmap)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniref90_taxid_idmap',
                            args.relpath_pickle_uniref90_taxid_idmap)
    configparser_object.set('process_proteomes', 'relpath_pickle_uniref50_taxid_idmap',
                            args.relpath_pickle_uniref50_taxid_idmap)

    configparser_object.set('process_proteomes', 'relpath_pickle_uniprotkb_idmap',
                            args.relpath_pickle_uniprotkb_idmap)
    configparser_object.set('process_proteomes', 'relpath_taxon_to_process',
                            args.relpath_taxon_to_process)
    
    configparser_object.set('process_proteomes', 'verbose', str(verbose))
    configparser_object.set('process_proteomes', 'nproc', str(args.nproc))

    configparser_object.add_section('stats')
    configparser_object.set('stats', 'relpath_pickle_proteomes',
                            args.relpath_pickle_proteomes)
    configparser_object.set('stats', 'relpath_pickle_taxontree',
                            args.relpath_pickle_taxontree)

    configparser_object.add_section('panproteomes')
    configparser_object.set('panproteomes', 'relpath_panproteomes_dir',
                            args.relpath_panproteomes_dir)
    configparser_object.set('panproteomes', 'relpath_pickle_taxontree',
                            args.relpath_pickle_taxontree)
    configparser_object.set('panproteomes', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('panproteomes', 'relpath_pickle_uniprotkb_uniref_idmap',
                            args.relpath_pickle_uniprotkb_uniref_idmap)
    configparser_object.set('panproteomes', 'relpath_pickle_uniref100_taxid_idmap',
                            args.relpath_pickle_uniref100_taxid_idmap)
    configparser_object.set('panproteomes', 'relpath_pickle_uniref90_taxid_idmap',
                            args.relpath_pickle_uniref90_taxid_idmap)
    configparser_object.set('panproteomes', 'relpath_pickle_uniref50_taxid_idmap',
                            args.relpath_pickle_uniref50_taxid_idmap)
    configparser_object.set('panproteomes', 'relpath_pickle_uniprotkb_idmap',
                            args.relpath_pickle_uniprotkb_idmap)
    configparser_object.set('panproteomes', 'relpath_pickle_proteomes',
                            args.relpath_pickle_proteomes)
    configparser_object.set('panproteomes', 'uniref_cluster_panproteomes',
                            args.uniref_cluster_panproteomes)
    configparser_object.set('panproteomes', 'discard_low_quality_genomes',
                            args.discard_low_quality_genomes)

    configparser_object.set('panproteomes', 'verbose', str(verbose))
    configparser_object.set('panproteomes', 'nproc', str(args.nproc))


    configparser_object.add_section('export')
    configparser_object.set('export', 'download_base_dir',
                            args.download_base_dir)
    configparser_object.set('export', 'export_dir',
                            args.export_dir)
    configparser_object.set('export', 'exportpath_phylophlan',
                            args.exportpath_phylophlan)
    configparser_object.set('export', 'relpath_panproteomes_dir',
                            args.relpath_panproteomes_dir)
    configparser_object.set('export', 'relpath_pickle_proteomes',
                            args.relpath_pickle_proteomes)
    configparser_object.set('export', 'relpath_pickle_taxontree',
                            args.relpath_pickle_taxontree)

    
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
