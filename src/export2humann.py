#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Fabio Cumbo (fabio.cumbo@unitn.it')

__date__ = '2 Jun 2019'


import argparse as ap
import configparser as cp
import datetime
import glob
import gzip
import itertools
import logging
import math
import multiprocessing.dummy as dummy
import multiprocessing as mp
from functools import partial
import os
import pickle
import re
import resource
import shutil
import time
import bz2
from collections import Counter
from operator import itemgetter
from tempfile import NamedTemporaryFile
import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pyfaidx import Fasta
import errno
import _version as version

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

OUTFILE_PREFIX = 'mpa_v{}_CHOCOPhlAn_{}'.format(version.__MetaPhlAn2_db_version__, version.__CHOCOPhlAn_version__)
shared_variables = type('shared_variables', (object,), {})
log = logging.getLogger(__name__)
ch = logging.FileHandler('CHOCOPhlAn_export2humann_{}.log'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')),'w')
ch.setLevel(logging.INFO)
log.addHandler(ch)

def initialize (config):
    if config['verbose']:
        utils.info('Loading pickled databases...')
    resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))

    shared_variables.config = config
    shared_variables.func_annot = {}
    shared_variables.uniprot = {}
    shared_variables.uniparc = {}
    shared_variables.taxontree = pickle.load(open("{}/{}".format(shared_variables.config['download_base_dir'],shared_variables.config['relpath_pickle_taxontree']), 'rb'))
    shared_variables.proteomes = pickle.load(open("{}{}".format(shared_variables.config['download_base_dir'],shared_variables.config['relpath_pickle_proteomes']), 'rb'))

    os.makedirs('{}/{}'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan']), exist_ok=True)
    os.makedirs('{}/{}/functional_annot/'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2']), exist_ok=True)
    os.makedirs('{}/{}/panproteomes_fna/'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2']), exist_ok=True)
    os.makedirs('{}/{}/gca_upkb_to_nr/'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2']), exist_ok=True)
    all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniprotkb/'.format(shared_variables.config['download_base_dir'], shared_variables.config['pickled_dir'])))
    all_uparc_chunks = filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/{}/uniparc/'.format(shared_variables.config['download_base_dir'], shared_variables.config['pickled_dir'])))
    uniref90_taxid = "{}{}".format(shared_variables.config['download_base_dir'],shared_variables.config['relpath_pickle_uniref90_taxid_idmap'])

    for i in all_uprot_chunks:
        uniprot_chunk = '{}/{}/uniprotkb/{}'.format(shared_variables.config['download_base_dir'], shared_variables.config['pickled_dir'], i)
        shared_variables.uniprot.update({entry[0] : [ entry[3],entry[4],entry[1] ] + list(entry[8:17]) for entry in utils.load_pickle(uniprot_chunk)})        #geneid, proteomeid, tax_id

    for i in all_uparc_chunks:
        uniparc_chunk = '{}/{}/uniparc/{}'.format(shared_variables.config['download_base_dir'], shared_variables.config['pickled_dir'], i)
        shared_variables.uniparc.update({elem[0] : [ elem[1] ]+ list(elem[3:6]) for elem in utils.load_pickle(uniparc_chunk)})

    shared_variables.uniref90_taxid_map = {}
    [shared_variables.uniref90_taxid_map.update(elem) for elem in utils.load_pickle(uniref90_taxid)]
    if shared_variables.config['verbose']:
        utils.info('Finished.\n')

def get_all_species_panproteomes():
    return [x.split('/')[-1][:-4] for x in glob.glob('{}/{}/species/90/*'.format(shared_variables.config['download_base_dir'], shared_variables.config['relpath_panproteomes_dir']))]


def get_uniref90_nucl_seq(item):
    (gca, taxid, taxonomy), t = item
    seqs_dir = '{}/{}/{}/{}/'.format(shared_variables.config['download_base_dir'], shared_variables.config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
    if not os.path.exists(seqs_dir):
        log.error('[{}]\tNo genome and GFF annotations downloaded for {}'.format(taxid, gca))
        return None
    try:
        in_seq_file = glob.glob(seqs_dir+'*.fna.gz')[0]
        in_gff_file = glob.glob(seqs_dir+'*.gff.gz')[0]
    except Exception as e:
        log.error('[{}]Missing genome or GFF annotation for {}'.format(taxid, gca))
        return None

    with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/tmp/chocophlan') as fna_decompressed:
        gff_db = gffutils.create_db(in_gff_file, ':memory:', id_spec='locus_tag', merge_strategy="merge")

        with gzip.open(in_seq_file, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
            shutil.copyfileobj(fna_gz, fna_out)

        with Fasta(fna_decompressed.name) as fna:
            nucls = []
            failed = []
            for NR90, NR50, gene_names in t:
                found = False
                for gene_name in gene_names:
                    c = gff_db.conn.cursor()
                    try:
                        query='{} WHERE featuretype == "gene" AND attributes LIKE ?'.format(gffutils.constants._SELECT), ('%"{}"%'.format(gene_name),)
                        _ = c.execute('{} WHERE featuretype == "gene" AND attributes LIKE ?'.format(gffutils.constants._SELECT), ('%"{}"%'.format(gene_name),))
                    except Exception as ex:
                        log.info('cursor for {}'.format(taxid))
                        log.info("{} WHERE featuretype == 'gene' AND attributes LIKE '%\"{}\"%'".format(gffutils.constants._SELECT, gene_name))
                    results = c.fetchone()
                    if results is not None:
                        feature = gffutils.Feature(**results)
                        found = True
                        break

                if found:  
                    try:
                        nuc_seq = feature.sequence(fna)
                    except Exception as e:
                        print(taxid)
                        raise e
                    mpa_marker = '{}__{}__{}|{}|{}|UniRef50_{}|{}|{}'.format(taxid, NR90[9:], feature['Name'][0], taxonomy, NR90, NR50, len(nuc_seq), gca)
                    record = SeqIO.SeqRecord(seq=Seq(nuc_seq, alphabet=DNAAlphabet()), id=mpa_marker, description='')
                    nucls.append(record)
                else:
                    failed.append(str(NR90))
    if not len(failed): failed = None
    [ gff_db.delete(x) for x in gff_db.all_features() ]
    return (nucls, failed, )

def export_panproteome_centroid(species_id):
    panproteome = pd.read_csv('{}/{}/{}.txt.bz2'.format(shared_variables.config['export_dir'], shared_variables.config['panproteomes_stats'], species_id), sep = '\t', index_col = 0)
    taxonomy = shared_variables.taxontree.print_full_taxonomy(species_id)[0].replace('|','.')
    low_lvls = [c.tax_id for c in shared_variables.taxontree.taxid_n[species_id].find_clades()]
    gc = {}
    func_annot = []
    
    for row in panproteome.itertuples():
        pangene = row.Index
        pangene = 'UniRef90_'+pangene.split('_')[1]
        pangene_nr50 = 'UniRef50_'+shared_variables.uniref90_taxid_map[pangene][1][2]

        for (upkb_id, upid, count) in (x.split(':') for x in row.copy_number.split(';')):
            is_uparc = upkb_id.startswith('UPI')
            upkb_entry = shared_variables.uniprot[upkb_id] if not is_uparc else shared_variables.uniparc[upkb_id]
            if is_uparc:
                gene_names = tuple(x[2] for x in upkb_entry[0] if x[1] in low_lvls)
            else:
                gene_names = tuple(x[1] for x in upkb_entry[0] )
            gca = dict(shared_variables.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
            if not gca:
                continue
            gca = gca.split('.')[0]
            if not is_uparc:
                upkb_fa = shared_variables.uniprot[upkb_id][3:]
            else:
                tmp = [ list(itertools.chain(set(x))) for x in zip(*[shared_variables.uniprot[x[0]][3:] for x in shared_variables.uniref90_taxid_map[pangene][0] 
                if not x[0].startswith('UPI') and x[1] in low_lvls and x[0] in shared_variables.uniprot])]
                upkb_fa = [list(set(itertools.chain.from_iterable(x))) for x in tmp[:6]]
            if upkb_fa: 
                upkb_fa = [ pangene_nr50,
                            ','.join('GO:'+str(x).zfill(7) for x in upkb_fa[0]),
                            ','.join('K'+str(x).zfill(5) for x in upkb_fa[1]),
                            ','.join(str(x) for x in upkb_fa[2]),
                            ','.join('PF'+str(x).zfill(5) for x in upkb_fa[3]),
                            ','.join(str(x) for x in upkb_fa[4]),
                            ','.join(str(x) for x in upkb_fa[5]) ]
                func_annot.append(tuple([pangene, upkb_fa]))

            if len(gca) and len(gene_names):
                if (gca, species_id) not in gc:
                    gc[(gca, species_id, taxonomy)] = set()
                gc[(gca, species_id, taxonomy)].add((pangene, pangene_nr50, gene_names))
                continue
    nproc = max(min(10,len(gc)),1)
    chunksize = 1 if nproc == len(gc) else 1+math.ceil(len(gc)/nproc)

    with dummy.Pool(processes=nproc) as pool_t:
        res = [ _ for _ in pool_t.imap_unordered(get_uniref90_nucl_seq, gc.items(), chunksize = chunksize) if _ is not None]

    failed, pangenes = [], []
    if res:
        pangenes, failed = zip(*res)
        pangenes = list(itertools.chain.from_iterable(filter(None,pangenes)))
    if failed is not None and not all(x is None for x in failed):
        failed = list(filter(None, failed))[0]
        log.warning("[{}]\tFailed to extract nucleotide sequences for the following pangene {}".format(species_id, (','.join(failed))))
    if len(pangenes):
        #save tsv of annots
        with bz2.open('{}/{}/panproteomes_fna/{}.fna.bz2'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2'], species_id), 'wt') as fasta_pangenes:
            [fasta_pangenes.write(marker.format('fasta')) for marker in pangenes]

        with bz2.BZ2File('{}/{}/functional_annot/{}.txt.bz2'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2'], species_id), 'w') as functional_annot:
            pickle.dump(func_annot, functional_annot)

def get_genes_coordinates(item):
    gca, t = item
    coords = []
    seqs_dir = '{}/{}/{}/{}/'.format(shared_variables.config['download_base_dir'], shared_variables.config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
    if not os.path.exists(seqs_dir):
        log.error('[{}]\tNo genome and GFF annotations downloaded for {}'.format(taxid, gca))
        return None
    try:
        in_seq_file = glob.glob(seqs_dir+'*.fna.gz')[0]
        in_gff_file = glob.glob(seqs_dir+'*.gff.gz')[0]
    except Exception as e:
        log.error('[{}]Missing genome or GFF annotation for {}'.format(taxid, gca))
        return None

        try:
            os.symlink(in_seq_file, '{}/{}/{}/panphlan_{}_genomes/{}.fna.gz'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id, species_id, gca))
        except OSError as e:
            if e.errno == errno.EEXIST:
                os.remove('{}/{}/{}/panphlan_{}_genomes/{}.fna.gz'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id, species_id, gca))
                os.symlink(in_seq_file, '{}/{}/{}/panphlan_{}_genomes/{}.fna.gz'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id, species_id, gca))

    with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/tmp/chocophlan') as fna_decompressed:
        gff_db = gffutils.create_db(in_gff_file, ':memory:', id_spec='locus_tag', merge_strategy="merge")

        with gzip.open(in_seq_file, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
            shutil.copyfileobj(fna_gz, fna_out)

        with Fasta(fna_decompressed.name) as fna:
            failed = []
            for NR90, gene_names in t:
                found = False
                for gene_name in gene_names:
                    c = gff_db.conn.cursor()
                    try:
                        _ = c.execute('{} WHERE featuretype == "gene" AND attributes LIKE ?'.format(gffutils.constants._SELECT), ('%"{}"%'.format(gene_name),))
                    except Exception as ex:
                        log.info('cursor for {}'.format(taxid))
                        log.info("{} WHERE featuretype == 'gene' AND attributes LIKE '%\"{}\"%'".format(gffutils.constants._SELECT, gene_name))
                    results = c.fetchone()
                    if results is not None:
                        feature = gffutils.Feature(**results)
                        found = True
                        break

                if found:  
                    coords.append(tuple(NR90, gene_name, gca, feature.chrom, feature.start, feature.end))
                else:
                    failed.append(str(gene_name))

    if not len(failed): failed = None
    [ gff_db.delete(x) for x in gff_db.all_features() ]
    return (coords, failed, )

def export_panphlan_panproteome(species_id):
    panproteome = pd.read_csv('{}/{}/{}.txt.bz2'.format(shared_variables.config['export_dir'], shared_variables.config['panproteomes_stats'], species_id), sep = '\t', index_col = 0)
    os.makedirs('{}/{}/{}'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id), exist_ok=True)
    os.makedirs('{}/{}/{}/panphlan_{}_genomes'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id, species_id), exist_ok=True)
  
    low_lvls = [c.tax_id for c in shared_variables.taxontree.taxid_n[species_id].find_clades()]
    gc = {}
    
    for row in panproteome.itertuples():
        pangene = row.Index
        pangene = 'UniRef90_'+pangene.split('_')[1]

        for (upkb_id, upid, count) in (x.split(':') for x in row.copy_number.split(';')):
            is_uparc = upkb_id.startswith('UPI')
            upkb_entry = shared_variables.uniprot[upkb_id] if not is_uparc else shared_variables.uniparc[upkb_id]
            if is_uparc:
                gene_names = tuple(x[2] for x in upkb_entry[0] if x[1] in low_lvls)
            else:
                gene_names = tuple(x[1] for x in upkb_entry[0] )
            gca = dict(shared_variables.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
            if not gca:
                continue
            gca = gca.split('.')[0]
            if len(gca) and len(gene_names):
                if gca not in gc:
                    gc[gca] = set()
                gc[gca].add((pangene, gene_names))

    nproc = min(10,len(gc))
    chunksize = 1 if nproc == len(gc) else round(nproc/len(gc))
    with dummy.Pool(processes=nproc) as pool_t:
        res = [ _ for _ in pool_t.imap_unordered(get_genes_coordinates, gc.items(), chunksize = chunksize) if _ is not None]

    failed, pangenes_annot = [], []
    if res:
        pangenes_annot, failed = zip(*res)
        pangenes_annot = list(itertools.chain.from_iterable(filter(None,pangenes_annot)))

    if failed is not None and not all(x is None for x in failed):
        failed = list(filter(None, failed))[0]
        log.warning("[{}]\tFailed to extract gene coordinates for the following pangene {}".format(species_id, (','.join(failed))))

    if len(pangenes_annot):
        #save tsv of annots
        with open('{}/{}/{}/panphlan_{}_pangenome.csv'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_panphlan'], species_id, species_id), 'wt') as panphlan_pangenome:
            [panphlan_pangenome.write('\n'.join('\t'.join(gene))) for gene in pangenes_annot]


def export_genome_annotation(panproteome_id):
    gca_id = dict(shared_variables.proteomes[panproteome_id].get('ncbi_ids',{})).get('GCSetAcc',None)
    species_id = shared_variables.proteomes[panproteome_id]['tax_id']
    missing = []
    if gca_id and len(shared_variables.proteomes[panproteome_id]['members']):
        gca_id = gca_id.split('.')[0]
        with bz2.open('{}/{}/gca_upkb_to_nr/{}.txt.bz2'.format(shared_variables.config['export_dir'], shared_variables.config['exportpath_humann2'], gca_id), 'wt') as fasta_pangenes:
            fasta_pangenes.write('#GCA\tUPKB\tNR90\tNR50\n')
            for member in shared_variables.proteomes[panproteome_id]['members']:
                if member in shared_variables.uniprot:
                    NR90, NR50 = shared_variables.uniprot[member][9:11]
                elif member in shared_variables.uniparc:
                    NR90, NR50 = shared_variables.uniparc[member][2:4]
                else:
                    missing.append(member)
                fasta_pangenes.write('{}\t{}\t{}\t{}\n'.format( gca_id, 
                                                                member,
                                                                NR90,
                                                                NR50
                                                            )
                                    )

def run_all(config):
    config['exportpath_metaphlan2'] = os.path.join( config['exportpath_metaphlan2'],OUTFILE_PREFIX)
    config['relpath_gca2taxa'] = config['relpath_gca2taxa'].replace('DATE',version.__UniRef_version__)
    initialize(config)

    mpa_pkl = pickle.load(bz2.BZ2File('{}/{}/{}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'r'))
    
    with mp.Pool(processes = 10, maxtasksperchild = 100) as pll:
        res = [_ for _ in pll.imap_unordered(export_panproteome_centroid, [int(s[0].split('|')[-1]) for s in mpa_pkl['taxonomy'].values()], chunksize=70)]

    all_annot = {}
    OUTFILE_PREFIX = 'HUMAnN2_CHOCOPhlAn_{}'.format(version.__CHOCOPhlAn_version__)
    with bz2.open('{}/{}/{}_functional_annotation_mapping.tsv.bz2'.format(config['export_dir'], config['exportpath_humann2'],OUTFILE_PREFIX), 'wt') as functional_annot:
        for s in glob.glob('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export_201901/humann2/functional_annot/*.txt.bz2'):
            with bz2.BZ2File(s) as p_fa:
                functional_annot.write('NR90\tNR50\tGO\tKO\tKEGG\tPfam\tEC\teggNOG\n')
                functional_annot.write('\n'.join( '{}\t{}'.format( k, '\t'.join(v)) for k, v in pickle.load(p_fa) ) )
    
    with mp.Pool(processes = 10, maxtasksperchild = 100) as pll:
        res = [ _ for _ in pll.imap_unordered(export_genome_annotation, shared_variables.proteomes, chunksize=50)]

    with bz2.open('{}/{}/gca_upkb_to_NR90_NR50.tsv.bz2'.format(config['export_dir'], config['exportpath_humann2']), 'wt') as genome_functional_annot:
        genome_functional_annot.write('NR90\tNR50\tGO\tKO\tKEGG\tPfam\tEC\teggNOG\n')
        for s in glob.iglob('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export_201901/humann2/gca_upkb_to_nr/*.txt.bz2'):
            with bz2.open(s, 'rt') as p_fa:
                genome_functional_annot.writelines( line for line in p_fa if '#GCA' not in line )

if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['export']

    run_all(config)
    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
