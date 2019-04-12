#!/usr/bin/env python3
__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Nicola Segata (nicola.segata@unitn.it), '
            'Nicolai Karcher (karchern@gmail.com),')

__date__ = '21 Mar 2019'

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter
from operator import itemgetter
from pyfaidx import  Fasta
from tempfile import NamedTemporaryFile
from _version import __UniRef_version__
import _version as version
import bz2
import datetime
import gffutils
import glob
import gzip
import hashlib
import importlib
import itertools
import logging
import multiprocessing as mp
import multiprocessing.dummy as dummy
import numpy as np
import os
import pandas as pd
import pickle
import random
import re
import resource
import shutil
import subprocess as sb
import sys
import tarfile
import time

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    utils = importlib.import_module('src.utils')
    from src.panproteomes import Panproteome

OUTFILE_PREFIX = 'mpa_v{}_CHOCOPhlAn_{}'.format(version.__MetaPhlAn2_db_version__, version.__CHOCOPhlAn_version__)

def init_parse(terminating_):
    global terminating
    terminating = terminating_

def grouper(iterable, n, fillvalue=None):
    #"Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

def load_file(f_path):
    with open(f_path, 'rt') as ifn:
        for line in ifn:
            yield line.strip().split()

def isconsecutive(l):
    m = max(l)
    if sum(l) == (m * (m+1))/ 2:
        return True
    return False

def range1(x): return range(1,x)

class CodeTimer:
    def __init__(self, logger, name=None, debug=False):
        self.name = " '"  + name + "'" if name else ''
        self.debug = debug
        self.logger = logger

    def __enter__(self):
        self.start = time.clock()

    def __exit__(self, exc_type, exc_value, traceback):
        self.took = (time.clock() - self.start) * 1000.0
        self.logger.error('Code block' + self.name + ': ' + str(self.took) + ' ms')

class export_to_metaphlan2:
    def __init__(self, config):
        self.debug = config['verbose']
        resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))
        self.log = logging.getLogger(__name__)
        ch = logging.FileHandler('CHOCOPhlAn_export_to_metaphlan_{}.log'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')),'w')
        ch.setLevel(logging.INFO)
        self.log.addHandler(ch)

        if config['verbose']:
            utils.info('Loading pickled databases...')
        self.db = {'taxonomy' : {}, 'markers': {}}
        self.config = config
        self.uniprot = {}
        self.uniparc = {}
        self.uniref90_info_map = {}
        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.exportpath = '{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'])
        all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniprotkb/'.format(config['download_base_dir'], config['pickled_dir'])))
        all_uparc_chunks = filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/{}/uniparc/'.format(config['download_base_dir'], config['pickled_dir'])))
        all_uniref_chunks = filter(re.compile('uniref90_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniref90/'.format(config['download_base_dir'], config['pickled_dir'])))

        for i in all_uprot_chunks:
            uniprot_chunk = '{}/{}/uniprotkb/{}'.format(config['download_base_dir'], config['pickled_dir'], i)
            self.uniprot.update({ entry[0] : [entry[3],entry[4],entry[1], entry[17]] for entry in utils.load_pickle(uniprot_chunk)})

        for i in all_uparc_chunks:
            uniparc_chunk = '{}/{}/uniparc/{}'.format(config['download_base_dir'], config['pickled_dir'], i)
            self.uniparc.update({ entry[0] : entry[1] for entry in utils.load_pickle(uniparc_chunk)})

        for i in all_uniref_chunks:
            uniref90_chunk = '{}/{}/uniref90/{}'.format(config['download_base_dir'], config['pickled_dir'], i)
            self.uniref90_info_map.update({ entry[0] :(entry[2],entry[3:6],len(entry[6])) for entry in utils.load_pickle(uniref90_chunk)})

        ## Check if export folder exists, if not, it creates it
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)

        if config['verbose']:
            utils.info('Finished.\n')

    def extract_markers(self, panproteome, core_coreness_threshold,
                        core_uniqueness_90_threshold = float('Inf'),
                        core_uniqueness_50_threshold = float('Inf'),
                        external_genomes_90_threshold = float('Inf'),
                        external_genomes_50_threshold = float('Inf')
                    ):
        pc_coreness_threshold = (core_coreness_threshold * panproteome['number_proteomes'])/100

        cores_to_use = {k:v['coreness'] for k,v in panproteome['members'].items() if k in self.uniprot and v['coreness'] > pc_coreness_threshold and not self.uniprot[k][3]}     #exclude if it is fragment of a gene

        if panproteome['number_proteomes'] > 1:
            markers = pd.DataFrame.from_dict({'gene' : core,
                        'len' : self.uniref90_info_map['UniRef90_{}'.format(core)][2],
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness_nosp']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness_nosp']['90_50'],
                        'external_genomes_90' : sum((panproteome['members'][core]['external_genomes']['90_90']-self.taxa_to_remove).values()),
                        'external_genomes_50' : sum((panproteome['members'][core]['external_genomes']['90_50']-self.taxa_to_remove).values())}
                        for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness_nosp']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness_nosp']['90_50'] <= core_uniqueness_50_threshold
                            and sum((panproteome['members'][core]['external_genomes']['90_90']-self.taxa_to_remove).values()) <= external_genomes_90_threshold
                            and sum((panproteome['members'][core]['external_genomes']['90_50']-self.taxa_to_remove).values()) <= external_genomes_50_threshold)

        else:
            markers = pd.DataFrame.from_dict({'gene' : core,
                        'len' : self.uniref90_info_map['UniRef90_{}'.format(core)][2], 
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness']['90_50'],
                        'external_genomes_90' : sum((panproteome['members'][core]['external_genomes']['90_90']-self.taxa_to_remove).values()),
                        'external_genomes_50' : sum((panproteome['members'][core]['external_genomes']['90_50']-self.taxa_to_remove).values())}
                        for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness']['90_50'] <= core_uniqueness_50_threshold
                            and sum((panproteome['members'][core]['external_genomes']['90_90']-self.taxa_to_remove).values()) <= external_genomes_90_threshold
                            and sum((panproteome['members'][core]['external_genomes']['90_50']-self.taxa_to_remove).values()) <= external_genomes_50_threshold)
        
        if not markers.empty:
            markers = markers[((markers.len > 150) & (markers.len < 1500))]
            if not markers.empty:
                markers = markers.assign(coreness_score = pd.Series(np.power(markers['coreness_perc'], 1/2)),
                                         uniqueness_50_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_50']))/10000-0.00001))/5),
                                         uniqueness_90_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_90']))/10000-0.00001))/5))
                markers = markers.assign(score = pd.Series((markers['coreness_score'] * markers['uniqueness_50_score'] * markers['uniqueness_90_score'])))
                markers = markers.sort_values('score',ascending=False)

                if panproteome['number_proteomes'] > 1 and not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
                    tiers = []
                    for row in markers.itertuples():
                        if ( (row.coreness_perc >= 0.8) and 
                             (row.uniqueness_90 <= 1) and 
                             (row.uniqueness_50 <= 2) and 
                             ( ((row.external_genomes_90 <= 10 ) and (row.external_genomes_50 <= 5) ) if (external_genomes_90_threshold != float('Inf') and external_genomes_50_threshold != float('Inf')) else True)):
                            tiers.append('A')

                        elif ( (row.coreness_perc >= 0.7) and 
                               (row.uniqueness_90 <= 5) and 
                               (row.uniqueness_50 <= 5) and 
                               ( ((row.external_genomes_90 <= 10 ) and (row.external_genomes_50 <= 10))  if (external_genomes_90_threshold != float('Inf') and external_genomes_50_threshold != float('Inf')) else True)):
                            tiers.append('B')

                        elif ( (row.coreness_perc >= 0.5) and 
                               (row.uniqueness_90 <= 10) and 
                               (row.uniqueness_50 <= 15) and 
                               ( ((row.external_genomes_90 <= 25 ) and (row.external_genomes_50 <= 20))  if (external_genomes_90_threshold != float('Inf') and external_genomes_50_threshold != float('Inf')) else True)):
                            tiers.append('C')
                        else:
                            tiers.append('Z')
                    markers = markers.assign(tier = tiers)
                else:
                    if not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
                        markers = markers.assign(tier = ['A']*len(markers))
                    else:
                        markers = markers.assign(tier = ['U']*len(markers))
        return markers

    def get_p_markers(self, panproteome):
        if panproteome['number_proteomes'] > 1 and not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
            markers = self.extract_markers(panproteome, 80, 1, 2, 10, 5)
            if not markers.empty and markers.tier.value_counts().sum() > 50:
                markers = markers[:150]
            else:
                markers = self.extract_markers(panproteome, 70, 5, 5, 15, 10)
                if not markers.empty and markers.tier.value_counts().sum() > 50 :
                    markers = markers[:150]
                else:
                    markers = self.extract_markers(panproteome, 50, 10, 15, 25, 20)
                    markers = markers[:150]

            if markers.empty:
                markers = pd.DataFrame({'gene':[], 'tier':'Z'})    
        else:
            markers = self.extract_markers(panproteome, 80, 0, 0, 10, 10)
            if markers.empty or (self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality and len(markers) < 100):
                markers = pd.DataFrame({'gene':[], 'tier':'Z'})
            else:
                markers = markers[:100]
        return markers

    '''
        Remove all external hits from Shigella species to improve uniqueness statistic
    '''
    def get_ecoli_markers(self, panproteome):
        shigella_ids = [x.tax_id for x in self.taxontree.taxid_n.values() if 'Shigella' in x.name]
        for pangene in panproteome['members']:
            if len(pangene):
                for cluster in ['90_90','90_50']:
                    external_species_nosp = set(panproteome['members'][pangene]['external_species_nosp'][cluster]).difference(shigella_ids)
                    external_genomes = [taxid for taxid, count in panproteome['members'][pangene]['external_genomes'][cluster].items() if taxid not in shigella_ids for _ in range(count)]

                    panproteome['members'][pangene]['uniqueness_nosp'][cluster] = len(external_species_nosp)
                    panproteome['members'][pangene]['external_genomes'][cluster] = Counter(external_genomes)

        markers = self.extract_markers(panproteome, 50, 5, float('Inf'), 15, 15)
        return markers

    def get_markers_for_species_in_group(self, panproteome):
        group_id = self.taxontree.taxid_n[panproteome['tax_id']].parent_tax_id
        
        group = self.taxontree.taxid_n[group_id]
        same_group_ids = [sp.tax_id for sp in group.find_elements()]

        genus = self.taxontree.taxid_n[group.parent_tax_id]
        same_genus_ids = [sp.tax_id for sp in genus.find_elements()]

        for pangene in panproteome['members']:
            if len(pangene):
                for cluster in ['90_90','90_50']:
                    external_species_nosp = set(panproteome['members'][pangene]['external_species_nosp'][cluster]).difference(same_group_ids)
                    
                    ids_to_remove = set(panproteome['members'][pangene]['external_genomes'][cluster].keys()).intersection(same_group_ids)
                    [panproteome['members'][pangene]['external_genomes'][cluster].pop(taxid) for taxid in ids_to_remove]

                    panproteome['members'][pangene]['uniqueness_nosp'][cluster] = len(external_species_nosp)

        markers = self.extract_markers(panproteome, 50, 20, 20)
        return markers

    '''
        Create marker starting from UniRef90 associated gene DNA sequence
    '''
    def get_uniref90_nucl_seq(self, item):
        (upid, taxid, taxonomy), t = item
        db = {'taxonomy' : {}, 'markers': {}}
        ncbi_ids = dict(self.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
        if ncbi_ids is not None:
            gca = ncbi_ids.split('.')[0]
        else:
            self.log.error('No genome associated to proteome {} available for species {}'.format(upid, taxid))
            return [None, None, None]
        seqs_dir = '{}/{}/{}/{}/'.format(self.config['download_base_dir'], self.config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
        clade = self.taxontree.taxid_n[self.taxontree.go_up_to_species(taxid)]
        if not os.path.exists(seqs_dir):
            self.log.error('For proteome {} / species {} : No genome and GFF annotations downloaded for {}'.format(upid, taxid, gca))
            return [None, None, None]
        try:
            in_seq_file = glob.glob(seqs_dir+'*.fna.gz')[0]
            in_gff_file = glob.glob(seqs_dir+'*.gff.gz')[0]
        except Exception as e:
            self.log.error('For proteome {} / species {} : Missing genome or GFF annotation for {}'.format(upid, taxid, gca))
            return [None, None, None]

        with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/tmp/chocophlan') as fna_decompressed:
            gff_db = gffutils.create_db(in_gff_file, ':memory:', id_spec='locus_tag', merge_strategy="merge")
            with gzip.open(in_seq_file, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
                shutil.copyfileobj(fna_gz, fna_out)

            nucls = []
            failed = []
            taxonomy_db = '{}|t__{}'.format(taxonomy[0], gca)
            for core, gene_names, ext_species in t:
                if type(gene_names) == tuple:
                    gene_names = [x[1] for x in gene_names]
                else:
                    gene_names = (gene_names, )

                found = False
                for gene_name in gene_names:
                    c = gff_db.conn.cursor()
                    try:
                        _ = c.execute('{} WHERE featuretype == "gene" AND attributes LIKE ?'.format(gffutils.constants._SELECT), ('%"{}"%'.format(gene_name),))
                    except Exception as ex:
                        self.log.info('cursor for {}'.format(taxid))
                        self.log.info("{} WHERE featuretype == 'gene' AND attributes LIKE '%\"{}\"%'".format(gffutils.constants._SELECT, gene_name))
                    results = c.fetchone()
                    if results is not None:
                        feature = gffutils.Feature(**results)
                        found = True
                        break

                if found:
                    fna = Fasta(fna_decompressed.name)

                    try:
                        nuc_seq = feature.sequence(fna).upper()
                    except Exception as e:
                        print(taxid)
                        raise e
                    mpa_marker = '{}__{}__{}'.format(taxid, core, feature['Name'][0].replace('|','_'))
                    record = '>{} {}\n{}\n'.format(mpa_marker, 'UniRef90_{};{};{}'.format(core, taxonomy[0], gca), nuc_seq)
                    nucls.append(record)
                    if taxonomy_db not in db['taxonomy']:
                        db['taxonomy'][taxonomy_db] = (taxonomy[1], sum(contig.rlen for contig in fna.faidx.index.values()),)
                    ext_genomes = []
                    for x in ext_species:
                        cp = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[x])
                        if len(cp):
                            ext_genomes.append(dict(self.proteomes[next(iter(cp))].get('ncbi_ids',{})).get('GCSetAcc','').split('.')[0])

                    db['markers'][mpa_marker] = { 'clade': '{}__{}'.format(clade.rank[0], clade.name),
                                                  'ext': list(set(ext_genomes)),
                                                  'len': len(nuc_seq),
                                                  'score': 0,
                                                  'taxon': taxonomy[0]}
                else:
                    failed.append(str(core))
        if not len(failed): failed = None
        return (nucls, db, failed, )

    """
    Use centroid of UniRef90 as plausible marker
    """
    def get_uniref_uniprotkb_from_panproteome(self, panproteome):
        with CodeTimer(name='Panproteome loading', debug=self.debug, logger=self.log):
            panproteome = pickle.load(open(panproteome, 'rb'))
        pp_tax_id = panproteome['tax_id']
        low_lvls = [c.tax_id for c in self.taxontree.taxid_n[pp_tax_id].find_clades()]
        taxonomy = self.taxontree.print_full_taxonomy(pp_tax_id)
        with CodeTimer(name='Marker extraction', debug=self.debug, logger=self.log):
            if pp_tax_id == 562:
                possible_markers = self.get_ecoli_markers(panproteome)
            elif pp_tax_id in [303,     #Pseudomonas putida 
                               76759,   #Pseudomonas monteilii 
                               316,     #Pseudomonas stutzeri
                               1396,    #Bacillus cereus
                               1428,    #Bacillus thuringiensis
                               1405,    #Bacillus mycoides
                               1408,    #Bacillus pumilus
                               1390,    #Bacillus amyloliquefaciens
                               86664,   #Bacillus flexus
                               115862,  #Mycobacterium caprae
                               621      #Shigella boydii
                ]:
                possible_markers = self.extract_markers(panproteome, 50, 14, 10)
                if not possible_markers.empty:
                    possible_markers = possible_markers[:150]
            elif pp_tax_id in [294      #Pseudomonas fluorescens
            ]:
                possible_markers = self.get_markers_for_species_in_group(panproteome)
            else:
                possible_markers = self.get_p_markers(panproteome)
        gc = {}
        upid = ""
        res = []
        with CodeTimer(name='Marker to gene extraction', debug=self.debug, logger=self.log):
            for np90_cluster in possible_markers.gene:
                marker = None
                np90_members = self.uniref90_info_map['UniRef90_{}'.format(np90_cluster)][0]
                # If UniRef90 has member from species or below, use it
                upkb_from_species = list(filter(lambda x:x[1] in low_lvls, np90_members))
                has_repr = any(x[2] for x in upkb_from_species)
                if upkb_from_species:
                    repr_uniprotkb = list(filter(lambda x: (x[2]==has_repr) and 'UPI' not in x[0], upkb_from_species))
                    repr_uniparc = list(filter(lambda x: (x[2]==has_repr) and 'UPI' in x[0], upkb_from_species))
                    if repr_uniprotkb:
                        marker = repr_uniprotkb[0][0]
                    elif repr_uniparc:
                        marker = repr_uniparc[0][0]
                # else:
                #     repr_uniprotkb = list(filter(lambda x: (x[2]==has_repr) and 'UPI' not in x[0], np90_members))
                #     repr_uniparc = list(filter(lambda x: (x[2]==has_repr) and 'UPI' in x[0], np90_members))
                #     if repr_uniprotkb:
                #         marker = repr_uniprotkb[0][0]
                #     elif repr_uniparc:
                #         marker = repr_uniparc[0][0]

                # If entry is from UniProtKB take the gene names from the proteomes in which is present
                if marker in self.uniprot:
                    upids = ["UP{}".format(str(upid).zfill(9)) for upid in self.uniprot[marker][1]]
                    #Some UniProtKB entries (like P00644) do not have any proteome associteda.
                    if not len(upids):
                        tax_id = self.uniprot[marker][2]
                        if hasattr(self.taxontree.taxid_n[tax_id], 'proteomes'):
                            upids = self.taxontree.taxid_n[tax_id].proteomes
                            upids = [up for up in upids if marker in self.proteomes[up]['members']]
                    # Use, if available a reference genome, otherwise, a non-reference and finally a redundant one
                    upid = [x for x in upids if x in self.proteomes and self.proteomes[x]['isReference']]
                    if not len(upid):
                        upid = [x for x in upids if x in self.proteomes and not self.proteomes[x]['upi']]
                        if not len(upid):
                            upid = [x for x in upids if x in self.proteomes and self.proteomes[x]['upi']]
                    if len(upid):
                        upid = upid[0]
                    else:
                        continue
                    gene_names = self.uniprot[marker][0]
                elif marker in self.uniparc:
                        all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['isReference']]
                        if not len(all_genes):
                            all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if not self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                            if not len(all_genes):
                                all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                        if len(all_genes):
                            upid = all_genes[0][0]
                            upid = "UP{}".format(str(upid).zfill(9))
                            gene_names = all_genes[0][1]
                        else:
                            continue
                else:
                    continue

                if len(upid) and len(gene_names):
                    if (upid, panproteome['tax_id'], taxonomy) not in gc:
                        gc[(upid, panproteome['tax_id'], taxonomy)] = set()
                    gc[(upid, panproteome['tax_id'], taxonomy)].add((np90_cluster,
                                                                    gene_names, 
                                                                    tuple(set(x for x in itertools.chain.from_iterable(panproteome['members'][np90_cluster]['external_species'].values())))
                                                                    )
                    )
                else:
                    self.log.warning('No genes or proteomes identified for marker {} in species {}'.format(marker, panproteome['tax_id']))
            
            markers_nucls, db, failed, res = [], {}, [], []
            if len(gc):
                with CodeTimer(name='Sequence extraction', debug=self.debug, logger=self.log):
                    with dummy.Pool(processes=3) as pool:
                        res = [_ for _ in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=10) if _ is not None]
            else:
                self.log.warning('No markers identified for species {}'.format(panproteome['tax_id']))
                return (panproteome['tax_id'], taxonomy)
            if len(res):
                markers_nucls, db, failed = zip(*res)
                for x in db:
                    if x is not None:
                        self.db['taxonomy'].update(x['taxonomy'])
                        self.db['markers'].update(x['markers'])
                markers_nucls = list(itertools.chain.from_iterable(filter(None,markers_nucls)))
            
            #create export/markers_fasta
            if failed is not None and not all(x is None for x in failed):
                failed = list(filter(None, failed))[0]
                self.log.warning("{}: Failed to find features of markers {}".format(pp_tax_id, (','.join(failed))))
            if len(markers_nucls):
                with CodeTimer(name='FNA export', debug=self.debug, logger=self.log):
                    with open('{}/{}/markers_fasta/{}.fna'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wt') as fasta_markers:
                        [fasta_markers.write(marker.format('fasta')) for marker in markers_nucls]
        if len(gc):
            with CodeTimer(name='Pickling markers info', debug=self.debug, logger=self.log):
                pickle.dump(gc, open('{}/{}/markers_info/{}.txt'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wb'))
        with CodeTimer(name='Save marker stats', debug=self.debug, logger=self.log):
            possible_markers.to_csv('{}/{}/markers_stats/{}.txt'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']))

def map_markers_to_genomes(mpa_pkl, taxontree, proteomes, outfile_prefix, config):
    bowtie2 = 'bowtie2'
    BT2OUT_FOLDER = '{}/mpa2_eval/bt2_out'.format(os.path.abspath(config['export_dir']))
    BT2IDX_FOLDER = '{}/mpa2_eval/bt2_idx'.format(os.path.abspath(config['export_dir']))

    os.makedirs('{}/mpa2_eval'.format(config['export_dir']), exist_ok=True)
    os.makedirs(BT2IDX_FOLDER, exist_ok=True)
    os.makedirs(BT2OUT_FOLDER, exist_ok=True)
    
    READLEN = 150
    
    try:
        sb.check_call([bowtie2, "--version"], stdout=sb.DEVNULL)
        sb.check_call(['bowtie2-build', "--version"], stdout=sb.DEVNULL)
    except Exception as e:
        sys.stderr.write("OSError: fatal error running '{}'. Is it in the system path?\n".format(bowtie2))
        sys.exit(1)
    all_genomes = glob.iglob('{}/{}/*/*/*.fna.gz'.format(config['download_base_dir'], config['relpath_genomes']))
    if not os.path.exists('{}/mpa2_eval/bt2_idx/contig2genome.tsv'.format(config['export_dir'])):
        utils.info('\tBuilding a contig->genome map...')
        with dummy.Pool(processes=100) as pool:
            contig2ass = [x for x in pool.imap_unordered(contig2assembly, all_genomes, chunksize=1000)]
        if contig2ass:
            contig2ass = { contig : genome for contigs,genome in contig2ass for contig in contigs }
            with open('{}/mpa2_eval/bt2_idx/contig2genome.tsv'.format(config['export_dir']), 'wt') as fcomp:
                fcomp.write('\n'.join(['{}\t{}'.format(k,v) for k,v in contig2ass.items()]))
        utils.info('\tDone\n')
    else:
        utils.info('\tMapping file (contig->genome) already computed. Loading...')
        contig2ass = { contig : genome for contig,genome in load_file('{}/mpa2_eval/bt2_idx/contig2genome.tsv'.format(config['export_dir'])) }
        utils.info('\tDone.\n')

    contig2ass = { contig : genome[:13] for contig,genome in contig2ass.items() }
    if not os.path.exists('{}/mpa2_eval/bt2_idx/DONE'.format(config['export_dir'])):
        terminating = dummy.Event()
        try:
            utils.info('\tCreating BowTie2 indexes...')
            with mp.Pool(initializer = init_parse, initargs = (terminating, ), processes=30) as pool:
                bt2_indexes = [x for x in pool.imap_unordered(build_bt2_idx, grouper(all_genomes, 1000), chunksize=10)]
                with open('{}/mpa2_eval/bt2_idx/DONE'.format(config['export_dir']), 'wt') as fcomp:
                    fcomp.write('\n'.join(bt2_indexes))
            utils.info('\tDone.\n')
        except Exception as e:
            utils.error('\tFailed to create a BowTie2 indexing...Exiting.')
            exit(1)
    else:
        with open('{}/mpa2_eval/bt2_idx/DONE'.format(config['export_dir']), 'rt') as fcomp:
            bt2_indexes = [x.strip() for x in fcomp]

    mpa_markers_fna = '{}/{}/{}.fna'.format(config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix) 
    mpa_markers_splitted_fna = '{}/{}/{}_splitted.fna'.format(config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix) 
    utils.info('\tSplitting markers in chunks of 150nts...\n')
    split_markers(mpa_markers_fna, mpa_markers_splitted_fna)
    
    if os.path.isfile(mpa_markers_splitted_fna):
        terminating = dummy.Event()
        try:
            utils.info('\tStarted mapping of MetaPhlAn2 markers to reference genomes...')
            bt2_map_args = [ (bt2_idx, mpa_markers_splitted_fna, os.path.join(BT2OUT_FOLDER, bt2_idx.split('/')[-1] + '.sam')) for bt2_idx in bt2_indexes ]
            with mp.Pool(initializer = init_parse, initargs = (terminating, ), processes=20) as pool:
                sams = [x for x in pool.imap_unordered(run_bowtie2, bt2_map_args, chunksize=2)]
        except Exception as e:
            utils.error('\tFailed to map BowTie2...Exiting.')
            exit(1)
        finally:
            utils.remove_file(mpa_markers_splitted_fna)
        utils.info('\tDone.\n')
    else:
        utils.error('\tMetaPhlAn2 markers were not splitted...Exiting.')
        exit(1)

    markers2contig = {}
    
    utils.info('\tLoading SAM files and genomes taxonomy...')
    for spath in sams:
        for line in load_file(spath):
            marker, flag, contig, start, CIGAR = line[0], int(line[1]), line[2], line[3], line[5]
            marker, chunk = marker.split(':::')
            if marker not in markers2contig:
                markers2contig[marker] = {}
            if contig not in markers2contig[marker]:
                markers2contig[marker][contig] = (int(chunk.split('/')[-1]),[])
            markers2contig[marker][contig][1].append((chunk, start, CIGAR, flag))

    gca2taxon = { gca[:-2] : taxid for upid, gca, taxid in ((p, dict(proteomes[p]['ncbi_ids']).get('GCSetAcc',None), proteomes[p]['tax_id']) for p in proteomes if 'ncbi_ids' in proteomes[p]) if gca is not None}
    taxon2gca = {}
    [taxon2gca.setdefault( taxontree.go_up_to_species(taxid) , [ ] ).append(gca) for gca, taxid in gca2taxon.items()]
    utils.info('\tDone.\n')
    marker_present = {}
    utils.info('\tBuilding marker species matrix presence...')
    for marker, contig, total_chunks, chunks in ((marker, contig, total_chunks, chunks) for marker, _ in markers2contig.items() for contig, (total_chunks, chunks) in _.items()):
        flag = set(x[3] for x in chunks)
        reverse_strand = all([hex(f & 0x10) == '0x10' for f in flag])

        chunks.sort(key=lambda x:(int(x[0].split('/')[0]), x[1]))
        idx = 1
        previous_start = 0
        temp_chunks = []
        for c in chunks:
            if previous_start == 0 and c[0].startswith('1/'):
                previous_start = int(c[1]) - (-1 if reverse_strand else +1) * READLEN
            if c[0] == '{}/{}'.format(idx, total_chunks) and previous_start + (-1 if reverse_strand else +1) * READLEN == int(c[1]):
                temp_chunks.append(c)
                idx = idx +1 
                previous_start = int(c[1])

        chunks = temp_chunks
        if [ int(c[1]) for c in chunks[1:] ] == [ int(c[1]) + (-1 if reverse_strand else +1) * READLEN for c in chunks[:-1] ]:
            gca = contig2ass[contig]
            taxon = []
            if gca in gca2taxon:
                taxon = int(gca2taxon[gca])
                if marker not in marker_present:
                     marker_present[marker] = set()
                marker_present[marker].add(taxon)
    utils.info('\tDone.\n')
    markers2externals = {}
    for marker, hits in marker_present.items():
        ext_species = Counter(taxontree.go_up_to_species(h) for h in hits if h !='')
        marker_taxid = int(marker.split('__')[0])
        markers2externals.setdefault(marker, {})
        markers2externals[marker]['coreness'] = ext_species.pop(marker_taxid) if marker_taxid in ext_species else 0
        markers2externals[marker]['external_species'] = [x for x in ext_species.keys()]
        markers2externals[marker]['external_genomes'] = sum(x for x in ext_species.values())
        markers2externals[marker]['ext_species'] = ';'.join('{}:{}'.format(k,v) for k,v in ext_species.items())
    
    utils.info('\tUpdating the MetaPhlAn2 database with new external species...')
    for marker, vals in markers2externals.items():
        ext_genomes = [(taxon2gca[x], taxontree.print_full_taxonomy(x)[0]) for x in vals['external_species']]
        newt = [gca for gcas, taxstr in ext_genomes for gca in gcas if taxstr+'|t__'+gca in mpa_pkl['taxonomy']]
        if newt and marker in mpa_pkl['markers']:
            s_exts = set(mpa_pkl['markers'][marker]['ext'])
            s_exts.update(newt)
            mpa_pkl['markers'][marker]['ext'] = list(s_exts)
    utils.info('\tDone.\n')
    return (mpa_pkl, markers2externals)

def contig2assembly(fn):
    accids = []
    with gzip.open(fn, 'rt', encoding='utf-8') as genome_handle:
        for line in genome_handle:
            if line.startswith('>'):
                accids.append(line[1:].split()[0])

    fn = fn.split('/')[-1]
    return (accids, fn)

def build_bt2_idx(bin_list):
    if not terminating.is_set():
        with NamedTemporaryFile(dir='{}/mpa2_eval/bt2_idx'.format(config['export_dir']), suffix='.fna') as tmp_bin, \
            open(tmp_bin.name, 'wt') as tmp_bin_handle:
            for f in filter(None, bin_list):
                with gzip.open(f, 'rt', encoding='utf-8') as fm:
                    shutil.copyfileobj(fm, tmp_bin_handle, 1024*1024*10)
            idx_name = tmp_bin.name.replace('.fna','_idx')
            bt2_build_comm = ['bowtie2-build', '--threads', '20', tmp_bin.name, idx_name]
            
            try:
                sb.check_call(bt2_build_comm, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
            except sb.CalledProcessError as e:
                terminating.set()
                [ utils.remove_file('{}.{}'.format(idx_name,p)) for p in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"] ]

        return idx_name

def run_bowtie2(args):
    bin_idx, mpa_markers_splitted_fna, out_sam = args
    if not terminating.is_set():
        bt2_mapping_comm = ['bowtie2', '-f', mpa_markers_splitted_fna,
                                       '-x', bin_idx,
                                       '-a', 
                                       '--very-sensitive',
                                       '--no-hd',
                                       '--no-sq',
                                       '--no-unal',
                                       '--threads 20',
                                       '-S', out_sam ]
        try:
            sb.check_call(bt2_mapping_comm, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
        except sb.CalledProcessError as e:
            terminating.set()
            utils.remove_file(out_sam)
            utils.error(str(e))

        return out_sam

def split_markers(ifn, ofn):
    readsize=150
    with open(ofn, "w") as output_handle:
        for record in SeqIO.parse(ifn, "fasta"):
            start = 0
            tlen = len(record.seq)//readsize
            for i in range(1,tlen+1):
                subseq = record.seq[start:readsize*i]
                new_record = SeqIO.SeqRecord(seq=subseq, id='{}:::{}/{}'.format(record.id, i, tlen), description  = '')
                start=readsize*i
                SeqIO.write(new_record, output_handle, "fasta")

def merge_bad_species(sgb_release, gca2taxonomy, config):
    # all_mpa_gca = [x.split('|')[7][3:] for x in mpa_db['taxonomy'].keys()]
    # count_marker_per_clade = Counter(v['taxon'] for v in mpa_db['markers'].values())

    r = re.compile( r"(.*(C|c)andidat(e|us)_.*)|"
                    r"(.*_sp(_.*|$|))|"
                    r"((.*_|^)(b|B)acterium(_.*|))|"
                    r"(.*(eury|)archaeo(n_|te|n$).*)|"
                    r"(.*(endo|)symbiont.*)|"
                    r"(.*genomosp_.*)|"
                    r"(.*unidentified.*)|"
                    r"(.*_bacteria_.*)|"
                    r"(.*_taxon_.*)|"
                    r"(.*_et_al_.*)|"
                    r"(.*_and_.*)|"
                    r"(.*(cyano|proteo|actino)bacterium_.*)")

    bad_genomes = { gca : '|'.join(taxstr.split('|')[:7]) for gca, taxstr in dict(zip(gca2taxonomy.GCA_accession,gca2taxonomy.taxstr)).items() if r.match(taxstr.split('|')[6][3:])}

    gca2sgb = dict((y, (x[1],x[9])) for _, x in sgb_release.iterrows() for y in x[5].split(',') if y in bad_genomes)
    merged = {}
    keep = {}
    for genome, (sgb, taxonomy) in gca2sgb.items():
        nrefgenomes = sgb_release.loc[sgb_release['ID'] == sgb, 'Number of reconstructed genomes'].item()
        #If genome is unique in the SGB, keep it
        if nrefgenomes == 1:
            keep.setdefault(bad_genomes[genome],[]).append(genome)
        #If there are multiple refences, keep only one and the others list in a list of merged genomes
        else:
            merged.setdefault(sgb,[]).append((genome,bad_genomes[genome]))
    
    new_species_merged = {}
    for sgb, gca_tax in merged.items():
        merged_into='|'.join(sgb_release.loc[sgb_release['ID'] == sgb, 'Assigned taxonomy'].item().split('|')[:7])
        species_in_sgb = set(taxstr for (gca, taxstr) in gca_tax)
        # merged_into=max([(s,count_marker_per_clade[s]) for s in species_in_sgb],key = lambda x:x[1])[0]
        new_species_merged.setdefault(merged_into,[]).extend(gca_tax)

    with open('{}/{}/merged_species_spp.tsv'.format(config['export_dir'],config['exportpath_metaphlan2']), 'wt') as fout:
        fout.write('GCA\told_taxonomy\tmerged_into\n')
        fout.write('\n'.join('{}\t{}\t{}'.format(gca, old_t, new_t) for new_t, gca_tax in new_species_merged.items() for gca, old_t in gca_tax))
        fout.write('\n'.join('{}\t{}\t{}'.format(gca, taxa, taxa) for taxa, gca_tax in keep.items() for gca in gca_tax))    

    gca2taxid = dict(zip(gca2taxonomy.GCA_accession,gca2taxonomy.NCBI_taxid))
    taxa_to_remove = [(x[0], gca2taxid[x[0]]) for merged_into, bad_species in new_species_merged.items() if len(bad_species) > 1 for x in bad_species]
    
    # taxa2markers = {}
    # [taxa2markers.setdefault(v['taxon'],[]).append(k) for k,v in mpa_db['markers'].items()]
    
    # for new_t, gca_tax in new_species_merged.items():
    #     for gca, old_t in gca_tax:
    #         if new_t != old_t and old_t in taxa2markers:
    #             [mpa_db['markers'].pop(marker) for marker in taxa2markers[old_t] if marker in mpa_db['markers']]
    #             mpa_db['taxonomy'].pop('{}|t__{}'.format(old_t,gca))

    return taxa_to_remove


def run_all(config):
    #Check if a previous export was performed for the same release
    #If exists but has failed (no DONE file in release folder), create a backup folder and export again
    if os.path.exists('{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)):
       if not os.path.exists('{}/{}/DONE'.format(config['export_dir'],config['exportpath_metaphlan2'])):
           os.rename('{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), 
            '{}/{}/{}_bak_{}'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX, datetime.datetime.today().strftime('%Y%m%d_%H%M')))
    
    for i in ['markers_stats','markers_fasta','markers_info']:
        os.makedirs('{}/{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX, i))

    config['exportpath_metaphlan2'] = os.path.join( config['exportpath_metaphlan2'],OUTFILE_PREFIX)
    config['relpath_gca2taxa'] = config['relpath_gca2taxa'].replace('DATE',__UniRef_version__)

    export = export_to_metaphlan2(config)
    utils.info('Filtering of low quality species...\n')

    gca2taxonomy = pd.read_csv(os.path.join(export.config['export_dir'], export.config['relpath_gca2taxa']), sep='\t')
    # gca2taxonomy = dict(zip(gca2taxonomy.GCA_accession, gca2taxonomy.NCBI_taxid))
    sgb_release = pd.read_csv('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/sgbrepo/releases/SGB.Jan19.txt.bz2', sep='\t')
    sgb_release = sgb_release.loc[sgb_release['# Label'] == 'SGB',]
    self.genomes_to_remove, self.taxa_to_remove = zip(*merge_bad_species(sgb_release, gca2taxonomy))
    
    if not os.path.exists('{}/{}/DONE'.format(export.config['export_dir'],export.config['exportpath_metaphlan2'])):
        species = [ '{}/{}/species/90/{}.pkl'.format(export.config['download_base_dir'], export.config['relpath_panproteomes_dir'], item.tax_id) for item in self.taxontree.lookup_by_rank()['species'] 
                    if self.taxontree.get_child_proteomes(item) and item.tax_id not in self.taxa_to_remove ]
        utils.info('Started extraction of MetaPhlAn2 markers.\n')
        with dummy.Pool(processes=30) as pool:
            failed = [x for x in pool.imap(export.get_uniref_uniprotkb_from_panproteome, species, chunksize=200)]
        utils.info('Done.\n')

        with open('{}/{}/failed_species_markers.txt'.format(export.config['export_dir'], export.config['exportpath_metaphlan2']), 'wt') as ofn_failed:
            ofn_failed.writelines('\n'.join('{}\t{}'.format(taxid, '\t'.join(taxonomy)) for taxid, taxonomy in filter(None, failed)))

        utils.info('Pickling the MetaPhlAn2 database...\n')
        with bz2.BZ2File('{}/{}/{}.orig.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
            pickle.dump(export.db, outfile, pickle.HIGHEST_PROTOCOL)
        utils.info('Done.\n')

        with open('{}/{}/DONE'.format(export.config['export_dir'],export.config['exportpath_metaphlan2']), 'wt') as fcomp:
            fcomp.write(OUTFILE_PREFIX)

        utils.info('Merging all the markers sequences...\n')
        with open('{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wb') as mpa2_markers_fna:
            for f in glob.glob('{}/{}/markers_fasta/*.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'])):
                with open(f, 'rb') as fm:
                    shutil.copyfileobj(fm, mpa2_markers_fna, 1024*1024*10)
        utils.info('Done.\n')
    else:
        utils.info('MetaPhlAn2 markers already extracted. Loading the MetaPhlAn2 database...\n')
        with bz2.open('{}/{}/{}.orig.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX)) as infn:
            export.db = pickle.load(infn)
    utils.info('Done.\n')

    mpa_pkl, markers2ext = map_markers_to_genomes(export.db, export.taxontree, export.proteomes, OUTFILE_PREFIX, export.config)
    
    with open('{}/{}/{}_ext.tsv'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        outfile.write('\n'.join('{}\t{}'.format(marker, ext_species['external_species']) for marker, ext_species in markers2ext.items()))

    with bz2.BZ2File('{}/{}/{}.orig.nomerged.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        pickle.dump(mpa_pkl, outfile, pickle.HIGHEST_PROTOCOL)
  
    try:
        with bz2.open('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wt') as outfile, \
             open('{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'r') as infile:
             outfile.writelines(infile.readlines())
    except:
        utils.remove_file('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX))
        utils.error('Failed to compress the MetaPhlAn2 markers database. Exiting...')
        exit(1)
    finally:
        utils.remove_file('{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX))

    all_gca = [x.split('|')[-1].split('__')[1] for x in mpa_pkl_merged['taxonomy'].keys()]
    for x in mpa_pkl_merged['markers']:
        mpa_pkl_merged['markers'][x]['ext'] = [y for y in set(mpa_pkl_merged['markers'][x]['ext']) if y in all_gca]

    with bz2.BZ2File('{}/{}/{}.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        pickle.dump(mpa_pkl_merged, outfile, pickle.HIGHEST_PROTOCOL)

    with tarfile.TarFile('{}/{}/{}.tar'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as mpa_tar:
        mpa_tar.add('{}/{}/{}.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), recursive=False)
        mpa_tar.add('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), recursive=False)

    md5hash = hashlib.md5()
    with open('{}/{}/{}.tar'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX),'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5hash.update(chunk)

    calc_hash = md5hash.hexdigest()[:32]

    with open('{}/{}/{}.md5'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as tar_md5:
        tar_md5.write(calc_hash)
    
    with open('{}/{}/mpa_latest'.format(export.config['export_dir'], export.config['exportpath_metaphlan2']), 'w') as latest:
        latest.write(OUTFILE_PREFIX)

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    run_all(config['export'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
