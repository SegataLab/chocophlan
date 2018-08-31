#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

from _version import __CHOCOPhlAn_version__
__date__ = '11 Apr 2018'


import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
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
# from BCBio import GFF
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

class export_to_metaphlan2:
    def __init__(self, config):
        resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))
        logging.basicConfig(level=logging.INFO, filename='CHOCOPhlAn_export_to_metaphlan.log', filemode='w')
        ch = logging.StreamHandler()
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
            self.uniref90_taxid_map.update({k:(set(t[:3] for t in v[2]), v[3:6]) for k,v in uniref90_chunk.items()})

        ## Check if export folder exists, if not, it creates it
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)

        if config['verbose']:
            utils.info('Finished.\n')



    @staticmethod
    def get_p_markers(panproteome):

        def rank_markers(markers):
            markers = pd.DataFrame(markers)
            markers['coreness_rank'] = np.power(markers.coreness, 1/6)
            markers['uniqueness_50_rank'] = -np.log(1-((10000-np.minimum(10000,markers.uniqueness_50))/10000-0.00001))/5
            markers['uniqueness_90_rank'] = -np.log(1-((10000-np.minimum(10000,markers.uniqueness_90))/10000-0.00001))/5
            markers['ranking'] = markers.coreness_rank * markers.uniqueness_50_rank * markers.uniqueness_90_rank
            markers = markers.sort_values('ranking',ascending=False)
            return markers

        def extract_markers(panproteome, core_coreness_threshold=80, core_uniqueness_90_threshold=1, core_uniqueness_50_threshold=10):
            pc_coreness_threshold = (core_coreness_threshold * panproteome['number_proteomes'])//100
            cores_to_use = {k:v['coreness'] for k,v in panproteome['members'].items() if v['coreness'] >= pc_coreness_threshold}

            try:
                return [{'gene' : core,
                         'coreness': (panproteome['members'][core]['coreness'] / panproteome['number_proteomes']),
                         'uniqueness_90' : panproteome['members'][core]['uniqueness_nosp']['90_90'],
                         'uniqueness_50' : panproteome['members'][core]['uniqueness_nosp']['90_50']} for core in cores_to_use if len(core)
                                if panproteome['members'][core]['uniqueness_nosp']['90_90'] <= core_uniqueness_90_threshold 
                            and panproteome['members'][core]['uniqueness_nosp']['90_50'] <= core_uniqueness_50_threshold]
            except:
                utils.info('Something went wrong in the extractiion of markers for panproteome {}\n'.format(panproteome['tax_id']))

        markers = extract_markers(panproteome)
        if not len(markers):
            tier = 'Z'
        else:
            if len(markers) >= 200:
                markers = rank_markers(markers)
                markers = markers[:200]
                tier = 'A'
            else:
                markers = extract_markers(panproteome, 70, 5, 30)
                if len(markers) >= 150:
                    markers = rank_markers(markers)
                    markers = markers[:150]
                    tier = 'B'
                else:
                    markers = extract_markers(panproteome, 50, float('Inf'), float('Inf') )
                    markers = rank_markers(markers)
                    markers = markers[:50]
                    tier = 'C'
        return (markers, tier, )

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
        return export_to_metaphlan2.get_p_markers(panproteome)

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
                    nuc_seq = feature.sequence(fna)
                    mpa_marker = '{}__{}__{}'.format(taxid, core, feature['Name'][0])
                    record = SeqIO.SeqRecord(seq=Seq(nuc_seq, alphabet=DNAAlphabet()), id=mpa_marker, description='UniRef90_{};{};{}'.format(core, taxonomy, gca))
                    nucls.append(record)
                    if taxonomy_db not in self.db['taxonomy']:
                        self.db['taxonomy'][taxonomy_db] = sum(record.rlen for record in fna.faidx.index.values())
                    ext_genomes = []
                    for x in ext_species:
                        cp = self.taxontree.get_child_proteomes(self.taxontree.taxid_n[x])
                        if len(cp):
                            ext_genomes.append(dict(self.proteomes[next(iter(cp))].get('ncbi_ids',{})).get('GCSetAcc','').split('.')[0])

                    self.db['markers'][mpa_marker] = { 'clade': '{}__{}'.format(clade.rank[0], clade.name),
                                              'ext': ext _genomes,
                                              'len': len(record),
                                              'score': 0,
                                              'taxon': taxonomy }
                else:
                    failed.append(core)
        return (nucls, failed, )

    """
    Use centroid of UniRef90 as plausible marker
    """
    def get_uniref_uniprotkb_from_panproteome(self, panproteome):
        try:
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
                possible_markers, tier = self.get_ecoli_markers(panproteome)
            else:
                possible_markers, tier = export_to_metaphlan2.get_p_markers(panproteome)
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
                markers_nucls, failed = zip(*map(self.get_uniref90_nucl_seq, gc.items()))
            except TypeError:
                self.log.error('Cannot parse panproteome {}'.format(panproteome['tax_id']))
                return panproteome['tax_id']
            else:
                #create export/metaphlan2
                if len(failed):
                    self.log.warning("Failed to find features of markers {}".format(','.join(failed)))
                if len(markers_nucls):
                    with open('export/metaphlan2/{}.fasta'.format(panproteome['tax_id']), 'wt') as fasta_markers:
                        [fasta_markers.write(marker.format('fasta')) for marker in markers_nucls]
        except Exception as ex:
            self.log.error("FAILED TO PRECESS PANPROTOEME: {}".format(pp_tax_id))
            raise ex
     

def run_all():
    config = utils.read_configs('settings.cfg')
    config = utils.check_configs(config)['export']
    export = export_to_metaphlan2(config)
    species = glob.glob('data/pickled/panproteomes/species/90/*')
    with dummy.Pool(processes=80) as pool:
        failed = [x for x in pool.imap(export.get_uniref_uniprotkb_from_panproteome, species, chunksize=200)]
    outfile_prefix = 'mpa_v21_CHOCOPhlAn_{}'.format(__CHOCOPhlAn_version__)
    with bz2.BZ2File('export/metaphlan2/{}'.format(outfile_prefix), 'w') as outfile:
        pickle.dump(self.db, outfile, pickle.HIGHEST_PROTOCOL)

    os.system('cat export/metaphlan2/*.fasta > export/metaphlan2/{}.markers.fasta && bowtie2-build export/metaphlan2/{}.markers.fasta {}'.format(outfile, outfile))
