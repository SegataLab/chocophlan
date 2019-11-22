#!/usr/bin/env python3
__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Nicola Segata (nicola.segata@unitn.it)')

__date__ = '14 May 2019'

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
        try:
            resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))
        except Exception as ex:
            print(ex)
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

    def get_gca(self, upid):
        return dict(self.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)

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
            markers = markers.set_index('gene')
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
                             (row.uniqueness_90 <= 2) and 
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
            markers = self.extract_markers(panproteome, 80, 2, 2, 10, 5)
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
        (gca, taxid, taxonomy), t = item
        db = {'taxonomy' : {}, 'markers': {}}

        seqs_dir = '{}/{}/{}/{}/'.format(self.config['download_base_dir'], self.config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
        clade = self.taxontree.taxid_n[self.taxontree.go_up_to_species(taxid)]
        if not os.path.exists(seqs_dir):
            self.log.error('[{}]\tNo genome and GFF annotations downloaded for {}'.format(taxid, gca))
            return [None, None, None]
        try:
            in_seq_file = glob.glob(seqs_dir+'*.fna.gz')[0]
            in_gff_file = glob.glob(seqs_dir+'*.gff.gz')[0]
        except Exception as e:
            self.log.error('[{}]\tMissing genome or GFF annotation for {}'.format(taxid, gca))
            return [None, None, None]

        with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/tmp/') as fna_decompressed:
            gff_db = gffutils.create_db(in_gff_file, ':memory:', id_spec='locus_tag', merge_strategy="merge")
            with gzip.open(in_seq_file, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
                shutil.copyfileobj(fna_gz, fna_out)

            nucls = []
            failed = []
            taxonomy_db = '{}|t__{}'.format(taxonomy[0], gca)
            for core, gene_names, ext_species in t:
                # if type(gene_names) == tuple:
                #     gene_names = [x[1] for x in gene_names]
                # else:
                #     gene_names = (gene_names, )

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
                    with Fasta(fna_decompressed.name) as fna:
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
                                                    'taxon': taxonomy[0]
                                                    }
                else:
                    failed.append(str(core))
        if not len(failed): failed = None
        [ gff_db.delete(x) for x in gff_db.all_features() ]
        return (nucls, db, failed, )


    def get_genename_from_gca(self, upkb_from_species, has_repr):
        repr_uniprotkb = list(filter(lambda x: (x[2]==has_repr) and 'UPI' not in x[0] and x[0] in self.uniprot, upkb_from_species))
        repr_uniparc = list(filter(lambda x: (x[2]==has_repr) and 'UPI' in x[0] and x[0] in self.uniparc, upkb_from_species))
        gca = None
        if repr_uniprotkb:
            # If entry is from UniProtKB take the gene names from the proteomes in which is present
            for marker, _, _ in repr_uniprotkb:
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
                    gene_names = self.uniprot[marker][0]
                    for u in upid:
                        ncbi_ids = dict(self.proteomes[u].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
                        if ncbi_ids is not None:
                            gca = ncbi_ids.split('.')[0]
                            return gca, gene_names

        if not gca or repr_uniparc:
            for marker, _, _ in repr_uniparc:
                all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['isReference']]
                if not len(all_genes):
                    all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if not self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                    if not len(all_genes):
                        all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                if len(all_genes):
                    for upid, gene_names in all_genes:
                        upid = "UP{}".format(str(upid).zfill(9))
                        ncbi_ids = dict(self.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
                        if ncbi_ids is not None:
                            gca = ncbi_ids.split('.')[0]
                            return gca, gene_names

        return None, None
    """
    Use centroid of UniRef90 as plausible marker
    """
    def get_uniref_uniprotkb_from_panproteome(self, panproteome):
        try:
            with bz2.open(panproteome, 'r') as p_f:
                panproteome = pickle.load(p_f)
        except:
            self.log.error('[{}]\tCannot load panproteome.'.format(panproteome))
            return 
        pp_tax_id = panproteome['tax_id']
        low_lvls = [c.tax_id for c in self.taxontree.taxid_n[pp_tax_id].find_clades()]
        taxonomy = self.taxontree.print_full_taxonomy(pp_tax_id)
        if pp_tax_id == 562:
            possible_markers = self.get_ecoli_markers(panproteome)
        elif pp_tax_id in self.species_merged_into:
            possible_markers = self.get_p_markers(panproteome)
            if len(possible_markers) < 10:
                possible_markers = self.extract_markers(panproteome, 50, 14, 10)
                if not possible_markers.empty:
                    possible_markers = possible_markers[:150]
        elif pp_tax_id in self.config['taxa_low_markers']:
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
        failed_uniprot_repr = []
        for np90_cluster in possible_markers.index:
            for (upkb_id, upid), count in panproteome['members'][np90_cluster]['copy_number'].items():
                is_uparc = upkb_id.startswith('UPI')
                upkb_entry = self.uniprot[upkb_id] if not is_uparc else self.uniparc[upkb_id]
                if is_uparc:
                    gene_names = tuple(x[2] for x in upkb_entry if x[1] in low_lvls)
                else:
                    gene_names = tuple(x[1] for x in upkb_entry[0] )
                gca = self.get_gca(upid)
                if not gca:
                    continue
                gca = gca.split('.')[0]
                # gca, marker = None, None
                # np90_members = self.uniref90_info_map['UniRef90_{}'.format(np90_cluster)][0]
                # # If UniRef90 has member from species or below, use it
                # upkb_from_species = list(filter(lambda x:int(x[1]) in low_lvls, np90_members))
                # has_repr = any(x[2] for x in upkb_from_species)
                # if upkb_from_species:
                #     for species_list, repr_flag in [(upkb_from_species,has_repr), (upkb_from_species,not has_repr),
                #                                     (np90_members,has_repr), (np90_members,not has_repr)]:
                #         gca, gene_names  = self.get_genename_from_gca(species_list,repr_flag)
                #         if gca is not None and gene_names is not None:
                #             break
                if gca and gene_names:
                    if (gca, panproteome['tax_id'], taxonomy) not in gc:
                        gc[(gca, panproteome['tax_id'], taxonomy)] = set()
                    gc[(gca, panproteome['tax_id'], taxonomy)].add((np90_cluster,
                                                                    gene_names, 
                                                                    tuple(set(x for x in itertools.chain.from_iterable(panproteome['members'][np90_cluster]['external_species'].values())))
                                                                    )
                    )
                    break
                else:
                    failed_uniprot_repr.append(np90_cluster)

        if failed_uniprot_repr:
            self.log.warning('[{}]\tFailed to find an UniProtKB representative for the following UniRef90: {}'.format(panproteome['tax_id'], ','.join(failed_uniprot_repr)))
            possible_markers = possible_markers.loc[~possible_markers.index.isin(failed_uniprot_repr)]

        markers_nucls, db, failed, res = [], {}, [], []
        if len(gc):
                res = list(map(self.get_uniref90_nucl_seq, gc.items()))
        else:
            if len(failed_uniprot_repr) == possible_markers.shape[0]:
                ifnotgenome = not any([True for p in self.taxontree.get_child_proteomes(self.taxontree.taxid_n[pp_tax_id])
                                        if dict(self.proteomes[p].get('ncbi_ids',{})).get('GCSetAcc','')])
                self.log.warning('[{}]\tFailed to find any gene for all possible markers identified.{}'.format(panproteome['tax_id'], "No assemblies are available for the species' proteomes" if ifnotgenome else ''))
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
            possible_markers = possible_markers.loc[~possible_markers.index.isin(failed)]
            self.log.warning("{}: Failed to find features of markers {}".format(pp_tax_id, (','.join(failed))))

        if len(markers_nucls) > 9:
            with open('{}/{}/markers_fasta/{}.fna'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wt') as fasta_markers:
                [fasta_markers.write(marker.format('fasta')) for marker in markers_nucls]
        elif len(markers_nucls) == 0:
            self.log.warning("[{}]\tNo nucleotide sequences extracted".format(pp_tax_id))
        if len(gc):
            pickle.dump(gc, open('{}/{}/markers_info/{}.txt'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wb'))
        possible_markers.to_csv('{}/{}/markers_stats/{}.txt'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']))

def map_markers_to_genomes(mpa_markers_fna, mpa_pkl, taxontree, proteomes, outfile_prefix, config):
    bowtie2 = '/shares/CIBIO-Storage/CM/mir/tools/bowtie2-2.3.4.3/bowtie2'
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
    if not os.path.exists('{}/mpa2_eval/bt2_idx/DONE_{}'.format(config['export_dir'], outfile_prefix)):
        terminating = dummy.Event()
        try:
            utils.info('\tCreating BowTie2 indexes...')
            all_genomes = glob.iglob('{}/{}/*/*/*.fna.gz'.format(config['download_base_dir'], config['relpath_genomes']))
            with mp.Pool(initializer = init_parse, initargs = (terminating, ), processes=30) as pool:
                bt2_indexes = [x for x in pool.imap_unordered(build_bt2_idx, grouper(all_genomes, 1000), chunksize=2)]            
            utils.info('\tDone.\n')
        except Exception as e:
            utils.error('\tFailed to create a BowTie2 indexing...Exiting.')
            return
        with open('{}/mpa2_eval/bt2_idx/DONE_{}'.format(config['export_dir'], outfile_prefix), 'wt') as fcomp:
            fcomp.write('\n'.join(bt2_indexes))
    else:
        with open('{}/mpa2_eval/bt2_idx/DONE_{}'.format(config['export_dir'], outfile_prefix), 'rt') as fcomp:
            bt2_indexes = ['{}/mpa2_eval/bt2_idx/{}'.format(config['export_dir'],x.strip()) for x in fcomp]

    mpa_markers_splitted_fna = mpa_markers_fna.replace('.fna', '_splitted.fna')
    utils.info('\tSplitting markers in chunks of 150nts...')
    split_markers(mpa_markers_fna, mpa_markers_splitted_fna)
    utils.info('\tDone.\n')
    if os.path.isfile(mpa_markers_splitted_fna):
        terminating = dummy.Event()
        try:
            utils.info('\tStarted mapping of MetaPhlAn2 markers to reference genomes...')
            bt2_map_args = [ (bowtie2, bt2_idx, mpa_markers_splitted_fna, os.path.join(BT2OUT_FOLDER, outfile_prefix + "__" + bt2_idx.split('/')[-1] + '.sam')) for bt2_idx in bt2_indexes ]
            with mp.Pool(initializer = init_parse, initargs = (terminating, ), processes=10) as pool:
                sams = [x for x in pool.imap_unordered(run_bowtie2, bt2_map_args, chunksize=2)]
        except Exception as e:
            utils.error(e)
            utils.error('\tFailed to map BowTie2...Exiting.')
            return
        finally:
            utils.remove_file(mpa_markers_splitted_fna)
        utils.info('\tDone.\n')
    else:
        utils.error('\tMetaPhlAn2 markers were not splitted...Exiting.')
        return

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

    gca2taxon = { gca[:-2] : taxontree.go_up_to_species(taxid) for upid, gca, taxid in ((p, dict(proteomes[p]['ncbi_ids']).get('GCSetAcc',None), proteomes[p]['tax_id']) for p in proteomes if 'ncbi_ids' in proteomes[p]) if gca is not None}
    taxon2gca = {}
    [taxon2gca.setdefault( taxid , [ ] ).append(gca) for gca, taxid in gca2taxon.items()]
    utils.info('\tDone.\n')
    marker_present = {}
    utils.info('\tBuilding marker species matrix presence...')
    for marker, contig, total_chunks, chunks in ((marker, contig, total_chunks, chunks) for marker, _ in markers2contig.items() 
                                                                                        for contig, (total_chunks, chunks) in _.items()
                                                ):
        flag = set(x[3] for x in chunks)
        reverse_strand = all([hex(f & 0x10) == '0x10' for f in flag])

        chunks.sort(key=lambda x:(int(x[0].split('/')[0]), x[1]))
        idx = 1
        previous_start = 0
        temp_chunks = []
        for c in chunks:
            if previous_start == 0 and c[0].startswith('1/'):
                previous_start = int(c[1]) - (-1 if reverse_strand else +1) * READLEN
            #Determine if the chunks are consecutive and re-build the marker
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
        with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export_201901/mpa2_eval/bt2_idx', suffix='.fna') as tmp_bin, \
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
    bowtie2, bin_idx, mpa_markers_splitted_fna, out_sam = args
    if not terminating.is_set():
        bt2_mapping_comm = [bowtie2, '-f', mpa_markers_splitted_fna,
                                       '-x', bin_idx,
                                       '-a', 
                                       '--very-sensitive',
                                       '--no-hd',
                                       '--no-sq',
                                       '--no-unal',
                                       '--threads 20',
                                       '-S', out_sam ]
        try:
            # print(' '.join(bt2_mapping_comm))
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
                    r"(.*_CAG_.*)|"
                    r"(.*_sp(_.*|$))|"
                    r"((.*_|^)(b|B)acterium(_.*|))|"
                    r"(.*(eury|)archaeo(n_|te|n$).*)|"
                    r"(.*(endo|)symbiont.*)|"
                    r"(.*genomosp_.*)|"
                    r"(.*unidentified.*)|"
                    r"(.*_bacteria_.*)|"
                    r"(.*_taxon_.*)|"
                    r"(.*_et_al_.*)|"
                    r"(.*_and_.*)|"
                    r"(.*Shigella_.*)|"
                    r"(.*(cyano|proteo|actino)bacterium_.*)")

    bad_genomes = { gca : ('|'.join(taxstr.split('|')[:7]), '|'.join(taxidstr.split('|')[:7]))
                    for _, gca, _, _, taxstr, taxidstr in gca2taxonomy.itertuples()
                    if r.match(taxstr.split('|')[6][3:])
                }

    gca2sgb = dict((y, (x[1], x[9], x[10])) for _, x in sgb_release.iterrows() for y in x['List of reference genomes'].split(','))
    tax2sgb =  dict()
    [tax2sgb.setdefault(taxid.split('|')[-1],[]).append(x['ID']) for _, x in sgb_release.iterrows() for y in x['List of alternative taxonomic IDs'].split(',') for taxid in y.split(':')[::2]]
    merged = {}
    keep = {}
    for genome, (sgb, taxonomy, taxidstr) in gca2sgb.items():
        if genome in bad_genomes:
            taxonomy, taxidstr = bad_genomes[genome]
            nrefgenomes = sgb_release.loc[sgb_release['ID'] == sgb, 'Number of reference genomes'].item()
            #If genome is unique in the SGB, keep it
            if nrefgenomes == 1:
                keep.setdefault((taxonomy, taxidstr,),[]).append(genome)
            #If there are multiple refences, keep only one and the others list in a list of merged genomes
            else:
                merged.setdefault(sgb,[]).append((genome,(taxonomy, taxidstr,)))
    
    #to_remove = set(y.split('|')[-1] for x,(y,z) in mpa_pkl['taxonomy'].items() if 'CAG' in x or 'Shigella' in x)
    # all_species = set((x, y.split('|')[-1]) for x,(y,z) in mpa_pkl['taxonomy'].items())
    report = {}

    # def x(species):
    #     report[species] = {}
    #     report[species]['n_tot'] = len(list(filter(lambda x: x.startswith(species+'_'), mpa_pkl['markers'])))

    #     ext_species_counter = Counter(itertools.chain((int(gca2taxonomy.loc[gca2taxonomy['GCA_accession'] == x].taxidstr.item().split('|')[6]) for marker in filter(lambda x: x.startswith(species+'_'), mpa_pkl['markers'])
    #                                              for x in mpa_pkl['markers'][marker]['ext']
    #                                                 )
    #                                 ))
    #     ext_sgb_counter = Counter(itertools.chain((gca2sgb[x][0] for marker in filter(lambda x: x.startswith(species+'_'), mpa_pkl['markers'])
    #                                              for x in mpa_pkl['markers'][marker]['ext']
    #                                                 )
    #                                 ))
    #     i=0.3
    #     report[species][str(i)+'__all_confounded'] = set(x[0] for x in filter(lambda x: x[1]> i*report[species]['n_tot'], ext_species_counter.items()))
    #     report[species][str(i)+'__all_confounded_len'] = len(report[species][str(i)+'__all_confounded'])
    #     same_sgb_as_species = set(x[0] for x in filter(lambda x: x[1]> i*report[species]['n_tot'], ext_sgb_counter.items()) if species in tax2sgb and x[0] in tax2sgb[species])
    #     report[species][str(i)+'__same_sgb'] = set(x for x in ext_species_counter if str(x) in tax2sgb and set(tax2sgb[str(x)]).intersection(same_sgb_as_species))
        
    #     merged_into = Counter([tax2sgb[str(ext)][0] for ext in report[species][str(i)+'__same_sgb'] for _ in range(ext_sgb_counter[tax2sgb[str(ext)][0]])])
    #     if merged_into:
    #         merged_into = merged_into.most_common(1)[0][0]
    #     report[species][str(i)+'__same_sgb_len'] = len(report[species][str(i)+'__same_sgb'])
    #     report[species]['merged_into'] = merged_into
    #     report[species]['taxonomy'] = taxontree.print_full_taxonomy(int(species))[0]

    # list(map(x, to_remove))
    # report_pd = pd.DataFrame.from_dict(report).T
    # report_pd.to_excel('markers_confounding.xlsx')
    
    #ADD REPORT TO NEW_SPECIES_MERGED

    new_species_merged = {}
    for sgb, gca_tax in merged.items():
        merged_into = ('|'.join(sgb_release.loc[sgb_release['ID'] == sgb, 'Assigned taxonomy'].item().split('|')[:7]),
                       '|'.join(sgb_release.loc[sgb_release['ID'] == sgb, 'Assigned taxonomic ID'].item().split('|')[:7]))
        if r.match(merged_into[0].split('|')[6]) and sgb_release.loc[sgb_release['ID'] == sgb, 'Number of Alternative taxonomies'].item() > 1:
            old = merged_into
            new_merged_into = list(filter(lambda x: not r.match(x[1]), enumerate(sgb_release.loc[sgb_release['ID'] == sgb, 'List of alternative taxonomies'].item().split(','))))
            if new_merged_into:
                merged_into = (new_merged_into[0][1].split(':')[0], sgb_release.loc[sgb_release['ID'] == sgb, 'List of alternative taxonomic IDs'].item().split(',')[new_merged_into[0][0]].split(':')[0], )

        species_in_sgb = set(taxstr for (gca, taxstr) in gca_tax)
        # merged_into=max([(s,count_marker_per_clade[s]) for s in species_in_sgb],key = lambda x:x[1])[0]
        new_species_merged.setdefault(merged_into,set()).update(species_in_sgb)

    # for species, v in report.items():
    #     if v['merged_into']:
    #         merged_into = ('|'.join(sgb_release.loc[sgb_release['ID'] == v['merged_into'], 'Assigned taxonomy'].item().split('|')[:7]),
    #                    '|'.join(sgb_release.loc[sgb_release['ID'] == v['merged_into'], 'Assigned taxonomic ID'].item().split('|')[:7]))
    #         if r.match(merged_into[0].split('|')[6]) and sgb_release.loc[sgb_release['ID'] == sgb, 'Number of Alternative taxonomies'].item() > 1:
    #             old = merged_into
    #             new_merged_into = list(filter(lambda x: not r.match(x[1]), enumerate(sgb_release.loc[sgb_release['ID'] == v['merged_into'], 'List of alternative taxonomies'].item().split(','))))
    #             if new_merged_into:
    #                 merged_into = (new_merged_into[0][1].split(':')[0], sgb_release.loc[sgb_release['ID'] == sgb, 'List of alternative taxonomic IDs'].item().split(',')[new_merged_into[0][0]].split(':')[0], )
    #         new_species_merged.setdefault(merged_into,set()).add( taxontree.print_full_taxonomy(int(species))) 

    with open('{}/{}/merged_species_spp.tsv'.format(config['export_dir'],config['exportpath_metaphlan2']), 'wt') as fout:
        fout.write('old_taxonomy\told_taxonomy_id\tmerged_into\tmerged_into_id\n')
        fout.write('\n'.join('{}\t{}\t{}\t{}'.format(list(old_tax)[0][0], list(old_tax)[0][1], new_t[0], new_t[1]) for new_t, old_tax in new_species_merged.items()))
        fout.write('\n')
        fout.write('\n'.join('{}\t{}\t{}\t{}'.format(taxa[0], taxa[1], taxa[0], taxa[1]) for taxa, gca_tax in keep.items() for gca in gca_tax))

    gca2taxid = dict(zip(gca2taxonomy.GCA_accession,gca2taxonomy.NCBI_taxid))
    taxa_to_remove = [int(x[1].split('|')[6]) for merged_into, bad_species in new_species_merged.items() if len(bad_species) > 1 or (len(bad_species)==1 and list(bad_species)[0] != merged_into) for x in bad_species]
    
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
    config['relpath_gca2taxa'] = config['relpath_gca2taxa'].replace('DATE',version.__UniRef_version__)

    export = export_to_metaphlan2(config)
    utils.info('Filtering of low quality species...\n')

    gca2taxonomy = pd.read_csv(os.path.join(export.config['export_dir'], export.config['relpath_gca2taxa']), sep='\t')
    # gca2taxonomy = dict(zip(gca2taxonomy.GCA_accession, gca2taxonomy.NCBI_taxid))
    sgb_release = pd.read_csv('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/sgbrepo/releases/Jul19/SGB.Jul19.txt.bz2', sep='\t', skiprows=1)
    sgb_release = sgb_release.loc[sgb_release['# Label'] == 'SGB',]
    export.taxa_to_remove = merge_bad_species(sgb_release, gca2taxonomy, export.config)
    export.taxa_to_remove = Counter(export.taxa_to_remove * 50000)
    merged_species = pd.read_csv('{}/{}/merged_species_spp.tsv'.format(config['export_dir'],config['exportpath_metaphlan2']), sep = '\t')
    export.species_merged_into = [int(x.split('|')[-1]) for x in set(merged_species.merged_into_id)]

    if not os.path.exists('{}/{}/DONE'.format(export.config['export_dir'],export.config['exportpath_metaphlan2'])):
        species = [ '{}/{}/species/90/{}.pkl'.format(export.config['download_base_dir'], export.config['relpath_panproteomes_dir'], item.tax_id) 
                    for item in export.taxontree.lookup_by_rank()['species'] 
                    if export.taxontree.get_child_proteomes(item) and item.tax_id not in export.taxa_to_remove ]
        utils.info('Started extraction of MetaPhlAn2 markers.\n')
        with dummy.Pool(processes=30) as pool:
            failed = [x for x in pool.imap(export.get_uniref_uniprotkb_from_panproteome, species, chunksize=200)]
        utils.info('Done.\n')

        with open('{}/{}/failed_species_markers.txt'.format(export.config['export_dir'], export.config['exportpath_metaphlan2']), 'wt') as ofn_failed:
            ofn_failed.writelines('\n'.join('{}\t{}'.format(taxid, '\t'.join(taxonomy)) for taxid, taxonomy in filter(None, failed)))

        utils.info('Pickling the MetaPhlAn2 database...\n')
        with bz2.BZ2File('{}/{}/{}.orig.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
            pickle.dump(export.db, outfile, protocol=2)
        utils.info('Done.\n')

        with open('{}/{}/DONE'.format(export.config['export_dir'],export.config['exportpath_metaphlan2']), 'wt') as fcomp:
            fcomp.write(OUTFILE_PREFIX)

        utils.info('Merging all the markers sequences...\n')
        mpa_fna_db = '{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX)
        with open(mpa_fna_db, 'wb') as mpa2_markers_fna:
            for f in glob.glob('{}/{}/markers_fasta/*.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'])):
                with open(f, 'rb') as fm:
                    mpa2_markers_fna.write(fm.read())
        utils.info('Done.\n')
    else:
        utils.info('MetaPhlAn2 markers already extracted. Loading the MetaPhlAn2 database...\n')
        with bz2.open('{}/{}/{}.orig.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX)) as infn:
            export.db = pickle.load(infn)
    utils.info('Done.\n')

    mpa_pkl, markers2ext = map_markers_to_genomes(mpa_fna_db, export.db, export.taxontree, export.proteomes, OUTFILE_PREFIX, export.config)

    with open('{}/{}/{}_ext.tsv'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        outfile.write('\n'.join('{}\t{}'.format(marker, ext_species['external_species']) for marker, ext_species in markers2ext.items()))

    with bz2.BZ2File('{}/{}/{}.orig.nomerged.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        pickle.dump(mpa_pkl, outfile, protocol=2)

    for tax, (taxid, glen) in mpa_pkl['taxonomy'].items():
        if 'CAG' in tax and not export.taxontree.taxid_n[int(taxid.split('|')[-1])].is_low_quality:
            true_species = '_'.join(tax.split('|')[-2][3:].split('_')[:2]) 
            true_taxid = list(filter(lambda x: x.name==true_species, taxontree.taxid_n.values()))
            if true_taxid:
                true_taxid = true_taxid[0].tax_id
            _taxid = taxid.split('|')[-1]

            if sum(1 for x in ( sum(1 for e in mpa_pkl['markers'][x]['ext'] if taxontree.go_up_to_species(gca2taxon[e]) == true_taxid) / len(mpa_pkl['markers'][x]['ext']) 
                if len(mpa_pkl['markers'][x]['ext']) else 0 
                for x in filter(lambda x:x.startswith(_taxid), mpa_pkl['markers'])
            ) if x > 0.5) /  len(list(filter(lambda x:x.startswith(_taxid), mpa_pkl['markers']))) > 0.5:
                torem = list(filter(lambda x:x.startswith(_taxid), mpa_pkl['markers']))
                [mpa_pkl['markers'].pop(x) for x in torem]
    
        ##remove c__Coriobacteriia and 92706 Brevibacterium flavum
    try:
        with bz2.open('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wt') as outfile, \
             open('{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                if record.id in mpa_pkl['markers']:
                    outfile.write(record.format('fasta'))
        with bz2.open('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'rt') as outfile:
            ids = [x[1:].split()[0].strip() for x in outfile if x[0] == '>']
        for marker in set(mpa_pkl['markers']).difference(ids):
            mpa_pkl['markers'].pop(marker)
                
    except:
        utils.remove_file('{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX))
        utils.error('Failed to compress the MetaPhlAn2 markers database. Exiting...')
        exit(1)
    finally:
        utils.remove_file('{}/{}/{}.fna'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX))

    all_gca = set(x.split('|')[-1].split('__')[1] for x in mpa_pkl['taxonomy'].keys())

    for x in mpa_pkl['markers']:
        mpa_pkl['markers'][x]['ext'] = list(all_gca.intersection(set(mpa_pkl['markers'][x]['ext'])))

    with bz2.BZ2File('{}/{}/{}.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as outfile:
        pickle.dump(mpa_pkl, outfile, protocol=2)

    with bz2.open('{}/{}/{}_marker_info.txt.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wt') as outfile:
        outfile.write('\n'.join('{}\t{}'.format(m, mpa_pkl['markers'][m]) for m in mpa_pkl['markers']))

    with tarfile.TarFile('{}/{}/{}.tar'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as mpa_tar:
        mpa_tar.add(name = '{}/{}/{}.pkl'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX),
                    arcname = '{}.pkl'.format(OUTFILE_PREFIX),
                    recursive = False )
                    
        mpa_tar.add(name = '{}/{}/{}.fna.bz2'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX),
                    arcname = '{}.fna.bz2'.format(OUTFILE_PREFIX),
                    recursive = False)

    md5hash = hashlib.md5()
    with open('{}/{}/{}.tar'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX),'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5hash.update(chunk)

    calc_hash = md5hash.hexdigest()[:32]

    with open('{}/{}/{}.md5'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'w') as tar_md5:
        tar_md5.write(calc_hash)
    
    with open('{}/{}/mpa_latest'.format(export.config['export_dir'], export.config['exportpath_metaphlan2']), 'w') as latest:
        latest.write(OUTFILE_PREFIX)

    all_markers_stats = []

    for s in species:
        tax_id = int(s.split('/')[-1].split('.')[0])
        try:
            marker_stat = pd.read_csv('{}/{}/markers_stats/{}.txt'.format(export.config['export_dir'], export.config['exportpath_metaphlan2'], tax_id))
            cc = marker_stat.tier.value_counts().to_dict()
        except FileNotFoundError:
            cc = {'A': 0, 'B': 0, 'C' : 0, 'U' : 0, 'Z' : 0}
        cc['n_markers'] = sum(cc.values())
        cc['tax_id'] = tax_id
        cc['tax_str'] = export.taxontree.print_full_taxonomy(tax_id)[0] if tax_id in export.taxontree.taxid_n else '' 
        cc['has_genome'] = any([True for p in export.taxontree.get_child_proteomes(export.taxontree.taxid_n[tax_id])
                                        if dict(export.proteomes[p].get('ncbi_ids',{})).get('GCSetAcc','')]) if tax_id in export.taxontree.taxid_n else False
        cc['excluded<10markers'] = True if cc['n_markers'] < 10 else False
        all_markers_stats.append(cc)

    all_markers_stats_df = pd.DataFrame.from_dict(all_markers_stats).fillna(0).astype({'A': int, 'B': int, 'C' : int, 'U' : int, 'Z' : int})
    all_markers_stats_df.to_csv('{}/{}/all_markers_stats.txt'.format(export.config['export_dir'], export.config['exportpath_metaphlan2']), sep = '\t', index = False)


def init_fuzzy(ncbi_):
    global ncbi_db
    ncbi_db = ncbi_

def get_fuzzy_taxonomy(k):
    return (k, ncbi_db.get_fuzzy_name_translation(k.split('|')[6][3:].replace('_',' '), sim=0.7))

def build_virus_db(config, taxontree, proteomes):
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    new_db_fp = '{0}/{1}/{2}/{2}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)
    v20_db_fp = '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/metaphlan2/metaphlan_databases/mpa_v20_m200.pkl'
    new_db_pkl = pickle.load(bz2.open(new_db_fp))
    v20_db_pkl = pickle.load(bz2.open(v20_db_fp))
    
    virus_markers = list(filter(lambda x: v20_db_pkl['markers'][x]['taxon'].startswith('k__V'), v20_db_pkl['markers']))
    tmp_pkl = { 'markers' : {}, 'taxonomy' : {} }
    tmp_pkl['markers'] = {k : v20_db_pkl['markers'][k] for k in virus_markers}
    tmp_pkl['taxonomy'] = {k:v for k,v in v20_db_pkl['taxonomy'].items() if k.startswith('k__V')}

    str2taxa = {}
    to_proc = []
    for x in tmp_pkl['taxonomy']:
        if x.startswith('k__Vir'):
            taxname = x.split('|')[6][3:].replace('_',' ')
            taxid = ncbi.get_name_translator([taxname])
            if not taxid:
                to_proc.append(x)
            else:
                taxid = taxid[taxname][0]
            if taxid:
                str2taxa[x] = taxid


    with mp.Pool(initializer = init_fuzzy, initargs = (ncbi, ), processes=30) as pool:
        d = {x[0] : x[1][0] for x in pool.imap_unordered(get_fuzzy_taxonomy, to_proc, chunksize=10)}

    # partial_fuzzy = partial(ncbi.get_fuzzy_name_translation,sim=0.7 )
    # d = [x[0] for x in map(partial_fuzzy, (k for k in to_proc)) ]
    
    for x, taxid in d:
        str2taxa[x] = taxid

    newtaxa = {}
    for taxstr, taxid in str2taxa.items():
        if 's__Marseillevirus' in taxstr:
            new_taxstr, new_taxid = ('k__Viruses|p__Viruses_unclassified|c__Viruses_unclassified|o__Viruses_unclassified|f__Marseilleviridae|g__Marseillevirus|s__Marseillevirus_marseillevirus',
                                     '10239||||944644|1513458|694581')
        elif 's__Mavirus' in taxstr:
            new_taxstr, new_taxid = ('k__Viruses|p__Viruses_unclassified|c__Viruses_unclassified|o__Viruses_unclassified|f__Lavidaviridae|g__Mavirus|s__Cafeteriavirus_dependent_mavirus',
                                     '10239||||1914302|993034|1932923')
        else:
            new_taxstr, new_taxid = taxontree.print_full_taxonomy(taxid)
        newtaxa[taxstr] = ((new_taxstr, new_taxid), tmp_pkl['taxonomy'].pop(taxstr))

    for old_taxstr, ((new_taxstr, new_taxid), genome_size) in newtaxa.items():
        clade = old_taxstr.split('|')[-1]
        tmp_pkl['taxonomy']['|'.join([new_taxstr, clade])] = (new_taxid, genome_size)

    all_m = list(tmp_pkl['markers'].keys())
    for marker in all_m:
        if len(tmp_pkl['markers'][marker]['ext']):
            tmp_pkl['markers'][marker]['ext'] = set(filter(lambda x: not x.startswith('GC'), tmp_pkl['markers'][marker]['ext']))
        old_taxon = tmp_pkl['markers'][marker]['taxon']
        newtaxa_key = list(filter(lambda x: x.startswith(old_taxon), newtaxa))
        if newtaxa_key:
            old_taxon = newtaxa_key[0]
            if old_taxon not in newtaxa:
                if 't__' in old_taxon:
                    #Keep the old clade and update the taxonomy till the t__
                    new_clade = tmp_pkl['markers'][marker]['clade']
                    old_taxon = '|'.join(old_taxon.split('|')[:-1])
                    new_taxon = '{}|{}'.format(newtaxa[old_taxon][0][0], new_clade)
                    taxid = newtaxa[old_taxon][0][1].split('|')[-1]
                elif 'g__' in old_taxon:
                    #Update the whole taxonomy
                    last_level = old_taxon.split('|')[-1][:3]
                    level_member = list(filter(lambda x: x.startswith(old_taxon), newtaxa.keys()))[0]
                    new_taxon = re.search( '.*{}.*\|'.format(last_level), newtaxa[level_member][0][0]).group()[:-1]
                    new_clade = new_taxon.split('|')[-1]
                    taxid = newtaxa[level_member][0][1].split('|')[new_taxon.count('|')]
            else:
                new_taxon = newtaxa[old_taxon][0][0]
                new_clade = new_taxon.split('|')[-1]
                taxid = newtaxa[old_taxon][0][1] .split('|')[-1]
            
        tmp_pkl['markers'][marker]['clade'] = new_clade
        tmp_pkl['markers'][marker]['taxon'] = new_taxon
        tmp_pkl['markers']['{}__{}'.format(taxid, marker)] = tmp_pkl['markers'].pop(marker)
    
    tmp_pkl['taxonomy'].update(new_db_pkl['taxonomy'])

    mpa_v294_viruses_fp = '{0}/{1}/{2}_viruses/{2}.fna'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)
    mpa_v294_fp = '{0}/{1}/{2}/{2}.fna.bz2'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)
    with bz2.open('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/metaphlan2/metaphlan_databases/mpa_v20_m200.fna.bz2', 'rb') as mpa_v20_fna:
        with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/tmp/') as fna_decompressed,  open(fna_decompressed.name, 'wb') as fna_out:
            fna_out.write(mpa_v20_fna.read())

            with Fasta(fna_decompressed.name) as v20_markers, \
                open(mpa_v294_viruses_fp, 'wt') as fout:
                    for m in tmp_pkl['markers']:
                        fout.write('>{}\n{}\n'.format(m, v20_markers[m.split('__')[1]][:].seq ))

    mpa_pkl, markers2ext = map_markers_to_genomes(mpa_v294_viruses_fp, tmp_pkl, taxontree, proteomes, OUTFILE_PREFIX+"_viruses", config)
    mpa_pkl['markers'].update(new_db_pkl['markers'])
    with open(mpa_v294_viruses_fp, 'at') as mpa_v29_virus, bz2.open(mpa_v294_fp, 'rb') as mpa_v20_fna:
        mpa_v29_virus.write(mpa_v20_fna.read())

    with bz2.open('{0}/{1}/{2}_viruses/{2}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wb') as mpa_viruses_pkl:
        pickle.dump(mpa_pkl, mpa_viruses_pkl, protocol=2)

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
