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
import datetime
import bz2
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from pyfaidx import  Fasta
from tempfile import NamedTemporaryFile
from operator import itemgetter
import subprocess as sb

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    utils = importlib.import_module('src.utils')
    from src.panproteomes import Panproteome

def init_parse(terminating_):
    global terminating
    terminating = terminating_

def log_result(retval):
    # This is called whenever foo_pool(i) returns a result.
    failed.append(retval)

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
            self.uniref90_taxid_map.update({k:(v[2],v[3:6],len(v[6])) for k,v in uniref90_chunk.items()})

        ## Check if export folder exists, if not, it creates it
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)

        if config['verbose']:
            utils.info('Finished.\n')

    def extract_markers(self, panproteome, core_coreness_threshold, core_uniqueness_90_threshold, core_uniqueness_50_threshold, external_genomes_90_threshold, external_genomes_50_threshold):
        pc_coreness_threshold = (core_coreness_threshold * panproteome['number_proteomes'])/100
        cores_to_use = {k:v['coreness'] for k,v in panproteome['members'].items() if v['coreness'] >= pc_coreness_threshold}

        if panproteome['number_proteomes'] > 1:
            markers = pd.DataFrame.from_dict({'gene' : core,
                        'len' : self.uniref90_taxid_map['UniRef90_{}'.format(core)][2],
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness_nosp']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness_nosp']['90_50'],
                        'external_genomes_90' : sum(panproteome['members'][core]['external_genomes']['90_90'].values()),
                        'external_genomes_50' : sum(panproteome['members'][core]['external_genomes']['90_50'].values())}
                        for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness_nosp']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness_nosp']['90_50'] <= core_uniqueness_50_threshold
                            and sum(panproteome['members'][core]['external_genomes']['90_90'].values()) <= external_genomes_90_threshold
                            and sum(panproteome['members'][core]['external_genomes']['90_50'].values()) <= external_genomes_50_threshold)

        else:
            markers = pd.DataFrame.from_dict({'gene' : core,
                        'len' : self.uniref90_taxid_map['UniRef90_{}'.format(core)][2], 
                        'coreness_perc': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                        'coreness': panproteome['members'][core]['coreness'],
                        'uniqueness_90' : panproteome['members'][core]['uniqueness']['90_90'],
                        'uniqueness_50' : panproteome['members'][core]['uniqueness']['90_50'],
                        'external_genomes_90' : sum(panproteome['members'][core]['external_genomes']['90_90'].values()),
                        'external_genomes_50' : sum(panproteome['members'][core]['external_genomes']['90_50'].values())}
                        for core in cores_to_use if len(core)
                            if panproteome['members'][core]['uniqueness']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness']['90_50'] <= core_uniqueness_50_threshold
                            and sum(panproteome['members'][core]['external_genomes']['90_90'].values()) <= external_genomes_90_threshold
                            and sum(panproteome['members'][core]['external_genomes']['90_50'].values()) <= external_genomes_50_threshold)
        
        if not markers.empty:
            markers = markers[((markers.len > 100) & (markers.len < 1500))]
            if not markers.empty:
                markers = markers.assign(coreness_score = pd.Series(np.power(markers['coreness_perc'], 1/2)),
                                         uniqueness_50_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_50']))/10000-0.00001))/5),
                                         uniqueness_90_score = pd.Series(-np.log(1-((10000-np.minimum(10000,markers['uniqueness_90']))/10000-0.00001))/5))
                markers = markers.assign(score = pd.Series((markers['coreness_score'] * markers['uniqueness_50_score'] * markers['uniqueness_90_score'])))
                markers = markers.sort_values('score',ascending=False)

                if panproteome['number_proteomes'] > 1 and not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
                    tiers = []
                    for row in markers.itertuples():
                        if ((row.coreness_perc >= 0.8) & (row.uniqueness_90 <= 1) & (row.uniqueness_50 <= 5) & (row.external_genomes_90 <= 10 ) & (row.external_genomes_50) <= 10 ):
                            tiers.append('A')
                        elif ((row.coreness_perc >= 0.7) & (row.uniqueness_90 <= 5) & (row.uniqueness_50 <= 10) & (row.external_genomes_90 <= 10 ) & (row.external_genomes_50) <= 10 ):
                            tiers.append('B')
                        elif (row.coreness_perc >= 0.5):    
                            tiers.append('C')
                    markers = markers.assign(tier = tiers)
                else:
                    if not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
                        markers = markers.assign(tier = ['A']*len(markers))
                    else:
                        markers = markers.assign(tier = ['U']*len(markers))
        return markers

    def get_p_markers(self, panproteome):
        if panproteome['number_proteomes'] > 1 and not self.taxontree.taxid_n[panproteome['tax_id']].is_low_quality:
            markers = self.extract_markers(panproteome, 80, 1, 5, 10, 10)
            if not markers.empty and markers.tier.value_counts().sum() > 50:
                markers = markers[:150]
            else:
                markers = self.extract_markers(panproteome, 70, 5, 10, 10, 10)
                if not markers.empty and markers.tier.value_counts().sum() > 50 :
                    markers = markers[:150]
                else:
                    markers = self.extract_markers(panproteome, 50, float('Inf'), float('Inf'), 10, 10)
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
                for cluster in panproteome['members'][pangene]['external_hits']:
                    external_species_nosp = set(taxid for taxid in panproteome['members'][pangene]['external_species_nosp'][cluster] if taxid not in shigella_ids)
                    panproteome['members'][pangene]['uniqueness_nosp'][cluster] = len(external_species_nosp)
        return self.get_p_markers(panproteome)

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
            self.log.error('No genome associated to proteome {} available for species {}'.format(upid, taxid))
            return [None, None, None]
        seqs_dir = 'data/ncbi/{}/{}/'.format(gca[4:][:6], gca[4:][6:])
        clade = self.taxontree.taxid_n[self.taxontree.go_up_to_species(taxid)]
        if not os.path.exists(seqs_dir):
            self.log.error('For proteome {} / species {} : No genome and GFF annotations downloaded for {}'.format(upid, taxid, gca))
            return [None, None, None]
        try:
            in_seq_file = glob.glob(seqs_dir+'*fna*')[0]
            in_gff_file = glob.glob(seqs_dir+'*gff*')[0]
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
                        nuc_seq = feature.sequence(fna)
                    except Exception as e:
                        print(taxid)
                        raise e
                    mpa_marker = '{}__{}__{}'.format(taxid, core, feature['Name'][0])
                    record = SeqIO.SeqRecord(seq=Seq(nuc_seq, alphabet=DNAAlphabet()), id=mpa_marker, description='UniRef90_{};{};{}'.format(core, taxonomy[0], gca))
                    nucls.append(record)
                    if taxonomy_db not in db['taxonomy']:
                        db['taxonomy'][taxonomy_db] = (taxonomy[1], sum(record.rlen for record in fna.faidx.index.values()),)
                    ext_genomes = []
                    for x in ext_species:
                        cp = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[x])
                        if len(cp):
                            ext_genomes.append(dict(self.proteomes[next(iter(cp))].get('ncbi_ids',{})).get('GCSetAcc','').split('.')[0])

                    db['markers'][mpa_marker] = { 'clade': '{}__{}'.format(clade.rank[0], clade.name),
                                              'ext': list(set(ext_genomes)),
                                              'len': len(record),
                                              'score': 0,
                                              'taxon': taxonomy[0] }
                else:
                    failed.append(str(core))
        if not len(failed): failed = None
        return (nucls, db, failed, )

    """
    Use centroid of UniRef90 as plausible marker
    """
    def get_uniref_uniprotkb_from_panproteome(self, panproteome):
        panproteome = pickle.load(open(panproteome, 'rb')) 
        pp_tax_id = panproteome['tax_id']

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
        if pp_tax_id == 562:
            possible_markers = self.get_ecoli_markers(panproteome)
        else:
            possible_markers = self.get_p_markers(panproteome)
        gc = {}
        upid = ""
        res = []
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
                        if not len(all_genes):
                            all_genes = [(x[0],x[2]) for x in self.uniparc[marker] if self.proteomes["UP{}{}".format("0"*(9-len(str(x[0]))),x[0])]['upi']]
                    if len(all_genes):
                        upid = all_genes[0][0]
                        upid = "UP{}{}".format("0"*(9-len(str(upid))),upid)
                        gene_names = all_genes[0][1]
                    else:
                        continue
                else:
                    continue

            if len(upid) and len(gene_names):
                if (upid, panproteome['tax_id'], taxonomy) not in gc:
                    gc[(upid, panproteome['tax_id'], taxonomy)] = set()
                gc[(upid, panproteome['tax_id'], taxonomy)].add((marker, gene_names, tuple(panproteome['members'][marker]['external_species_nosp']['90_90'])))

            markers_nucls, db, failed, res = [], {}, [], []
            if len(gc):
                # for item in gc.items():
                #     res = self.get_uniref90_nucl_seq(item)
                #     if res is not None:
                #         res_markers_nucls, res_db, res_failed = res
                #         failed.extend(res_failed)
                #         for x in res_db:
                #             if x is not None:
                #                 self.db['taxonomy'].update(x['taxonomy'])
                #                 self.db['markers'].update(x['markers'])
                #         markers_nucls.extend([marker.format('fasta') for marker in filter(None,res_markers_nucls)])
                    # res.append(self.get_uniref90_nucl_seq(item))
                # with dummy.Pool(processes=3) as pool:
                #     for res in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items()):
                #         if res is not None:
                #             res_markers_nucls, res_db, res_failed = res
                #             failed.extend(res_failed)
                #             for x in res_db:
                #                 if x is not None:
                #                     self.db['taxonomy'].update(x['taxonomy'])
                #                     self.db['markers'].update(x['markers'])
                #             markers_nucls.extend([marker.format('fasta') for marker in filter(None,res_markers_nucls)])
                with dummy.Pool(processes=3) as pool:
                    res = [_ for _ in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=10) if _ is not None]
                if len(res):
                    markers_nucls, db, failed = zip(*res)
                    for x in db:
                        self.db['taxonomy'].update(x['taxonomy'])
                        self.db['markers'].update(x['markers'])
                    markers_nucls = [marker.format('fasta') for marker in itertools.chain.from_iterable(markers_nucls)]

            else:
                self.log.warning('No markers identified for species {}'.format(panproteome['tax_id']))
                return panproteome['tax_id']
            
            #create export/markers_fasta
            if failed is not None and not all(x is None for x in failed):
                failed = list(filter(None, failed))[0]
                self.log.warning("{}: Failed to find features of markers {}".format(pp_tax_id, (','.join(failed))))
            if len(markers_nucls):
                with open('{}/{}/markers_fasta/{}.fna'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wt') as fasta_markers:
                    fasta_markers.write(''.join(markers_nucls))
        if len(gc):
            pickle.dump(gc, open('{}/{}/markers_info/{}.txt'.format(self.config['export_dir'], self.config['exportpath_metaphlan2'], panproteome['tax_id']), 'wb'))
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

    failed = []
    # for s in species:
    #     failed.append(export.get_uniref_uniprotkb_from_panproteome(s))
    with dummy.Pool(processes=30) as pool:
        failed = [x for x in pool.imap(export.get_uniref_uniprotkb_from_panproteome, species, chunksize=200000)]

    # pool = dummy.Pool(processes=150)
    # res = [ pool.apply_async(export.get_uniref_uniprotkb_from_panproteome, args = (s,), callback=log_result) for s in species ]
    # pool.close()
    # pool.join()

    with open('{}/{}/failed_markers.txt'.format(config['export_dir'], config['exportpath_metaphlan2']), 'wt') as ofn_failed:
        ofn_failed.writelines('\n'.join(str(x) for x in filter(None, failed)))

    all_gca = [x.split('|')[-1].split('__')[1] for x in export.db['taxonomy'].keys()] 
    for x in export.db['markers']:
        export.db['markers'][x]['ext'] = [y for y in set(export.db['markers'][x]['ext']) if y in all_gca]
                

    with bz2.BZ2File('{}/{}/{}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], outfile_prefix), 'w') as outfile:
        pickle.dump(export.db, outfile, pickle.HIGHEST_PROTOCOL)


    #######
    merge_comm = 'cat {}/{}/markers_fasta/*.fna > {}/{}/{}.fna'.format(config['export_dir'], config['exportpath_metaphlan2'],config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix)
    bwt2b_comm = 'bowtie2-build {}/{}/{}.fna {}/{}/{}'.format(config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix, config['export_dir'], config['exportpath_metaphlan2'],outfile_prefix)
    sb.run(merge_comm.split())
    sb.run(bwt2b_comm.split(), stderr = sb.DEVNULL, stdout = sb.DEVNULL, check = True)