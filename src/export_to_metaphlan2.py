#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

__date__ = '11 Apr 2018'

import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
import _version as version
import glob
import time
import re
import gzip
import itertools
import shutil
import logging
import importlib
import traceback
import pandas as pd
import numpy as np
import gffutils
import resource
import bz2
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from pyfaidx import  Fasta
from tempfile import NamedTemporaryFile
from operator import itemgetter

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    utils = importlib.import_module('src.utils')
    from src.panproteomes import Panproteome

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class CodeTimer:
    def __init__(self, logger, name=None, debug=False):
        self.name = " '"  + name + "'" if name else ''
        self.debug = debug
        self.logger = logger

    def __enter__(self):
        if self.debug:
            self.start = time.clock()

    def __exit__(self, exc_type, exc_value, traceback):
        if self.debug:
            self.took = (time.clock() - self.start) * 1000.0
            self.logger.info('Code block' + self.name + ': ' + str(self.took) + ' ms')

class export_to_metaphlan2:
    def __init__(self, config):
        self.debug = config['verbose']
        resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))
        ch = logging.FileHandler('CHOCOPhlAn_export_to_metaphlan.log')
        ch.setLevel(logging.INFO)
        self.log = logging.getLogger('export_to_metaphlan2')
        self.log.addHandler(ch)

        if config['verbose']:
            utils.info('Loading pickled databases...')
        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        # self.uniref90_taxid_map = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_uniref90_taxid_idmap']), 'rb'))
        self.db = {'taxonomy' : {}, 'markers': {}}
        self.config = config
        self.uniprot = {}
        self.uniparc = {}
        self.uniref90_taxid_map = {}
        self.exportpath = '{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'])
        all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniprotkb/'.format(config['download_base_dir'], config['pickled_dir'])))
        all_uparc_chunks = filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/{}/uniparc/'.format(config['download_base_dir'], config['pickled_dir'])))
        all_uniref_chunks = filter(re.compile('uniref90_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniref90/'.format(config['download_base_dir'], config['pickled_dir'])))

        for i in all_uprot_chunks:
            uniprot_chunk = pickle.load(open('{}/{}/uniprotkb/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniprot.update({k : [v[3],v[4],v[1]] for k,v in uniprot_chunk.items()})

        for i in all_uparc_chunks:
            uniparc_chunk = pickle.load(open('{}/{}/uniparc/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniparc.update({k : v[1] for k,v in uniparc_chunk.items()})

        for i in all_uniref_chunks:
            uniref90_chunk = pickle.load(open('{}/{}/uniref90/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniref90_taxid_map.update({k:(set(t[:3] for t in v[2]), v[3:6], len(v[6])) for k,v in uniref90_chunk.items()})

        ## Check if export folder exists, if not, it creates it
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)

        if config['verbose']:
            utils.info('Finished.\n')

    def extract_markers(self, panproteome, core_coreness_threshold=80, core_uniqueness_90_threshold=1, core_uniqueness_50_threshold=10):
        pc_coreness_threshold = (core_coreness_threshold * panproteome['number_proteomes'])/100
        cores_to_use = {k:v['coreness'] for k,v in panproteome['members'].items() if v['coreness'] >= pc_coreness_threshold}

        if panproteome['number_proteomes'] > 1:
            markers = [{'gene' : core,
                        'len' : self.uniref90_taxid_map['UniRef90_{}'.format(core)][2], 
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness_nosp']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness_nosp']['90_50']} for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness_nosp']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness_nosp']['90_50'] <= core_uniqueness_50_threshold]

        else:
            markers = [{'gene' : core,
                        'len' : self.uniref90_taxid_map['UniRef90_{}'.format(core)][2], 
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness']['90_50']} for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness']['90_50'] <= core_uniqueness_50_threshold]

        markers = pd.DataFrame.from_dict(markers)
        if len(markers):
            markers = markers[((markers.len > 100) & (markers.len < 1000))]
        return markers

    def rank_markers(self, markers):
        if len(markers):
            markers = markers.assign(coreness_score = pd.Series(np.power(markers['coreness_perc'], 1/2)))
            markers = markers.assign(uniqueness_50_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_50']))/10000-0.00001))/5))
            markers = markers.assign(uniqueness_90_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_90']))/10000-0.00001))/5))
            markers = markers.assign(score = pd.Series((markers['coreness_score'] * markers['uniqueness_50_score'] * markers['uniqueness_90_score'])))
            markers = markers.sort_values('score',ascending=False)
            return markers

    def get_p_markers(self, panproteome):
        if panproteome['number_proteomes'] > 1 and not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
            markers = self.extract_markers(panproteome)
            if len(markers) >= 200:
                markers = self.rank_markers(markers)
                markers = markers[:150]
            else:
                markers = self.extract_markers(panproteome, 70, 5, 30)
                if len(markers) >= 150:
                    markers = self.rank_markers(markers)
                    markers = markers[:150]
                else:
                    markers = self.extract_markers(panproteome, 50, float('Inf'), float('Inf') )
                    if len(markers) > 0:
                        markers = self.rank_markers(markers)
                        markers = markers[:150]

            if not len(markers):
                markers = pd.DataFrame({'gene':[], 'tier':'Z'})
            else:
                tiers = []
                for row in markers.itertuples():
                    if ((row.coreness_perc >= 0.8) & (row.uniqueness_90 <= 1) & (row.uniqueness_50 <= 5)):
                        tiers.append('A')
                    elif ((row.coreness_perc >= 0.7) & (row.uniqueness_90 <= 5) & (row.uniqueness_50 <= 10)):
                        tiers.append('B')
                    elif (row.coreness_perc >= 0.5):
                        tiers.append('C')
                    else:
                        print('aaaa')
                        break
                markers = markers.assign(tier = tiers)
                
        else:
            markers = self.extract_markers(panproteome, 80, 0, 0)
            if not len(markers) or (self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality and len(markers) < 100):
                markers = pd.DataFrame({'gene':[], 'tier':'Z'})
            else:
                markers = self.rank_markers(markers)
                markers = markers[:100]
                markers = markers.assign(tier = ['U']*len(markers))

        return markers

    '''
        Remove all external hits from Shigella species to improve uniqueness statistic
    '''
    def get_ecoli_markers(self, panproteome):
        shigella_ids = [x.tax_id for x in self.taxontree.taxid_n.values() if 'Shigella' in x.name]
        for pangene in panproteome['members']:
            if len(pangene):
                for cluster in panproteome['members'][pangene]['external_hits']:
                    external_species_nosp = set(taxid for taxid in panproteome['members'][pangene]['external_species_nosp'][cluster] if taxid not in shigella_ids)
                    panproteome['members'][pangene]['uniqueness_nosp'][cluster] = len(external_species_nosp)
        return self.get_p_markers(panproteome)

    # def exctract_gene_from_genome(self, item):
    #     if not terminating.is_set():
    #         upid, t  = item
    #         gca = dict(self.proteomes[upid]['ncbi_ids']).get('GCSetAcc','')
    #         seqs_dir = 'data/ncbi/{}/{}/'.format(gca[4:][:6], gca[4:][6:-2])
    #         in_seq_file = glob.glob(seqs_dir+'*fna*')[0]
    #         in_gff_file = glob.glob(seqs_dir+'*gff*')[0]

    #         file_handle, gff_decompressed = tempfile.mkstemp(prefix="gff")
    #         seq_parser = SeqIO.parse(gzip.open(in_seq_file,'rt'), "fasta")
    #         seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(in_seq_file,'rt'), "fasta"))

    #         with gzip.open(in_gff_file, 'rb') as gff_gz:
    #             shutil.copyfileobj(gff_gz, open(gff_decompressed, 'wb'))

    #         gff_parser = list(GFF.parse(open(gff_decompressed), base_dict = seq_dict))
    #         os.remove(gff_decompressed)

    #         nucls = []
    #         for core_pangene, upkb in t:
    #             if 'UPI' not in upkb:
    #                 gene_name = self.uniprot[upkb]
    #             else:
    #                 all_genes = dict([(x[0],x[2]) for x in self.uniparc[upkb]])
    #                 upid_ = int(upid[2:])
    #                 gene_name = all_genes[upid_]

    #             if type(gene_name) == tuple:
    #                 gene_names = [x[1] for x in gene_name]
    #             else:
    #                 gene_names = (gene_name, )

    #             for gene_name in gene_names:
    #                 for chromosome in gff_parser:
    #                     for feature in chromosome.features:
    #                         if feature.type == 'gene' and gene_name in feature.qualifiers.get('Name'):
    #                             nuc_seq = feature.location.extract(chromosome)
    #                             nuc_seq.id = nuc_seq.features[0].qualifiers.get('Name',[''])[0]
    #                             nuc_seq.description = 'UniRef90_{}'.format(core_pangene)
    #                             nucls.append(nuc_seq)
    #         return nucls
    #     else:
    #         terminating.set()

    '''
        Create marker starting from UniRef90 associated gene DNA sequence
    '''
    def get_uniref90_nucl_seq(self, item):
        (upid, taxid, taxonomy), t = item
        db = {'taxonomy' : {}, 'markers': {}}
        ncbi_ids = dict(self.proteomes[upid]['ncbi_ids']).get('GCSetAcc',None)
        if ncbi_ids is not None:
            gca = ncbi_ids.split('.')[0]
        else:
            return []
        seqs_dir = 'data/ncbi/{}/{}/'.format(gca[4:][:6], gca[4:][6:])
        clade = self.taxontree.taxid_n[self.taxontree.go_up_to_species(taxid)]
        if not os.path.exists(seqs_dir):
            self.log.error('Missing genome and GFF annotations for {}'.format(gca))
            return []
        try:
            in_seq_file = glob.glob(seqs_dir+'*fna*')[0]
            in_gff_file = glob.glob(seqs_dir+'*gff*')[0]
        except Exception as e:
            self.log.error('Missing genome or GFF annotation for {} / {}. {}'.format(gca, upid, e))
            return []

        with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/tmp/chocophlan') as fna_decompressed:
            gff_db = gffutils.create_db(in_gff_file, ':memory:', id_spec='locus_tag', merge_strategy="merge")
            with gzip.open(in_seq_file, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
                shutil.copyfileobj(fna_gz, fna_out)

            nucls = []
            failed = []
            taxonomy_db = '{}|t__{}'.format(taxonomy, gca)
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
                        nuc_seq = feature.sequence(fna)
                    except Exception as e:
                        print(taxid)
                        raise e
                    mpa_marker = '{}__{}__{}'.format(taxid, core, feature['Name'][0])
                    record = SeqIO.SeqRecord(seq=Seq(nuc_seq, alphabet=DNAAlphabet()), id=mpa_marker, description='UniRef90_{};{};{}'.format(core, taxonomy, gca))
                    nucls.append(record)
                    if taxonomy_db not in db['taxonomy']:
                        db['taxonomy'][taxonomy_db] = sum(record.rlen for record in fna.faidx.index.values())
                    ext_genomes = []
                    for x in ext_species:
                        cp = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[x])
                        if len(cp):
                            ext_genomes.append(dict(self.proteomes[next(iter(cp))].get('ncbi_ids',{})).get('GCSetAcc','').split('.')[0])

                    db['markers'][mpa_marker] = { 'clade': '{}__{}'.format(clade.rank[0], clade.name),
                                              'ext': list(set(ext_genomes)),
                                              'len': len(record),
                                              'score': 0,
                                              'taxon': taxonomy }
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

        with CodeTimer(name='Panproteome taxonomy extraction', debug=self.debug, logger=self.log):
            try:
                taxonomy = self.taxontree.print_full_taxonomy(pp_tax_id)
                # print(taxonomy)
            except Exception as ex:
                try:
                    taxonomy = self.taxontree.print_full_taxonomy(pp_tax_id)
                except Exception as e:
                    self.log.error('CRITICAL ERROR: {}'.format(pp_tax_id))
                    self.log.error(traceback.format_exc())
                    raise ex
        with CodeTimer(name='Marker extraction', debug=self.debug, logger=self.log):
            if pp_tax_id == 562:
                possible_markers = self.get_ecoli_markers(panproteome)
            else:
                possible_markers = self.get_p_markers(panproteome)
        gc = {}
        upid = ""
        res = []
        with CodeTimer(name='Marker to gene extraction', debug=self.debug, logger=self.log):
            for np90_marker in possible_markers.gene:
                np90_members = self.uniref90_taxid_map['UniRef90_{}'.format(np90_marker)][0]
                for m in np90_members:
                    # If UniRef90 has member from species or below, use it
                    if m[2] == True and (m[1] == pp_tax_id or m[1] in self.taxontree.taxid_n[pp_tax_id].find_clades(tax_id=pp_tax_id)):
                        marker = m[0]
                        if marker in panproteome['members']:
                            break
                    # Otherwise, take the repr one
                    elif m[2] == True:
                        marker = m[0]
                        if marker in panproteome['members']:
                            break
                try:
                    if 'UPI' not in marker:
                        '''
                        If entry is from UniProtKB take the gene names from the proteomes in which is present
                        '''
                        if marker in self.uniprot:
                            upids = ["UP{}{}".format("0"*(9-len(str(upid))),upid) for upid in self.uniprot[marker][1]]
                            '''
                                Some UniProtKB entries (like P00644) do not have any proteome associteda.

                            '''
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
                        else:
                            continue
                    else:
                        if marker in self.uniparc:
                            all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['isReference']]
                            if not len(all_genes):
                                all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if not self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                            if len(all_genes):
                                upid = all_genes[0][0]
                                upid = "UP{}{}".format("0"*(9-len(str(upid))),upid)   
                                gene_names = all_genes[0][1]
                            else:
                                continue
                        else:
                            continue
                except:
                    self.log.info('proteome: {}\nmarker: {}\n'.format(panproteome['tax_id'], marker))

                if len(upid) and len(gene_names):
                    if (upid, panproteome['tax_id'], taxonomy) not in gc:
                        gc[(upid, panproteome['tax_id'], taxonomy)] = set()
                    gc[(upid, panproteome['tax_id'], taxonomy)].add((marker, gene_names, tuple(panproteome['members'][marker]['external_species_nosp']['90_90'])))
            try:
                markers_nucls, db, failed, res = [], {}, [], []
                if len(gc):
                    with CodeTimer(name='Sequence extraction', debug=self.debug, logger=self.log):
                        #res = [self.get_uniref90_nucl_seq(x) for x in gc.items()]
                        with dummy.Pool(processes=3) as pool:
                            res = [_ for _ in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=10) if _ is not None]
                        if len(res):
                            markers_nucls, db, failed = zip(*res)
                            for x in db:
                                self.db['taxonomy'].update(x['taxonomy'])
                                self.db['markers'].update(x['markers'])
                            markers_nucls = list(itertools.chain.from_iterable(markers_nucls))
                else:
                    self.log.error('No markers identified for species {}'.format(panproteome['tax_id']))
                    return panproteome['tax_id']
            except Exception as ex:
                self.log.error('Cannot parse panproteome {}'.format(panproteome['tax_id']))
                return panproteome['tax_id']
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

def run_all(): 
    config = utils.read_configs('settings.cfg')
    config = utils.check_configs(config)['export']
    outfile_prefix = 'mpa_v{}_CHOCOPhlAn_{}'.format(version.__MetaPhlAn2_db_version__, version.__CHOCOPhlAn_version__)

    if os.path.exists('{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix)):
        os.rename('{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix), 
            '{}/{}/{}_bak_{}'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix, datetime.datetime.today().strftime('%Y%m%d_%H%M')))
    for i in ['markers_stats','markers_fasta','markers_info']:
        os.makedirs('{}/{}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix, i))

    config['exportpath_metaphlan2'] = config['exportpath_metaphlan2'] + '/' + outfile_prefix

    species = glob.glob('{}/{}/species/90/*'.format(config['download_base_dir'], config['relpath_panproteomes_dir']))
    export = export_to_metaphlan2(config)


    with dummy.Pool(processes=30) as pool:
        failed = [x for x in pool.imap(export.get_uniref_uniprotkb_from_panproteome, species, chunksize=200000)]
    
    with open('{}/{}/failed_markers.txt'.format(config['export_dir'], config['exportpath_metaphlan2']), 'wt') as ofn_failed:
        ofn_failed.writelines('\n'.join(str(x) for x in filter(None, failed)))

    all_gca = [x.split('|')[-1].split('__')[1] for x in export.db['taxonomy'].keys()] 
    for x in export.db['markers']:
        export.db['markers'][x]['ext'] = [y for y in set(export.db['markers'][x]['ext']) if y in all_gca]
                

    with bz2.BZ2File('{}/{}/{}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix), 'w') as outfile:
        pickle.dump(export.db, outfile, pickle.HIGHEST_PROTOCOL)

    os.system('cat {}/{}/markers_fasta/*.fna > {}/{}/{}.fna'.format(config['export_dir'], config['exportpath_metaphlan2'],config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix))
    os.system('bowtie2-build {}/{}/{}.fna {}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix, config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix))