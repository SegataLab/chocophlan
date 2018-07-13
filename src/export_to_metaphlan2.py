#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)')

from _version import __version__
__date__ = '11 Apr 2018'


import os
import argparse as ap
import configparser as cp
import pickle
import multiprocessing.dummy as dummy
import glob
import time
import tempfile
import re
import gzip
import itertools
import shutil
import gc
from BCBio import GFF
from Bio import SeqIO
from operator import itemgetter

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

def init_parse(terminating_):
    global terminating
    terminating = terminating_

class export_to_metaphlan2:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')
        gc.disable()
        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.config = config
        self.uniprot = {}
        self.uniparc = {}

        all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/pickled'.format(config['download_base_dir'])))
        all_uparc_chunks = filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/pickled'.format(config['download_base_dir'])))

        for i in all_uprot_chunks:
            uniprot_chunk = pickle.load(open('{}/pickled/{}'.format(config['download_base_dir'], i),'rb'))
            self.uniprot.update({k : v[3:5] for k,v in uniprot_chunk.items()})

        for i in all_uparc_chunks:
            uniparc_chunk = pickle.load(open('{}/pickled/{}'.format(config['download_base_dir'], i),'rb'))
            self.uniparc.update({k : v[1] for k,v in uniparc_chunk.items()})
        gc.enable()
        if config['verbose']:
            utils.info('Finished.\n')

    CORE_CORENESS_THRESOLD = 0.9
    CORE_UNIQUENESS_90_THRESHOLD = 1
    CORE_UNIQUENESS_50_THRESHOLD = 10
    ## Coreness greater or equal than  90%
    ## Uniqueness 90_50 less or equal than 10
    ## Uniqueness 90_90 less or equal than 1, no sps
    ## Copy number less than 1.5

    @staticmethod
    def rank_genes(panproteome, config):
        pc_coreness_thresold = 0.9 * panproteome['number_proteomes']
        cores_to_use = {k:v for k,v in Panproteome.find_core_genes(panproteome).items() if v > pc_coreness_thresold}

        markers = [core for core in cores_to_use 
                        if panproteome['members'][core]['uniqueness']['90_90'] <= 1 
                    and panproteome['members'][core]['uniqueness']['90_50'] <= 10]
        return markers


    def elaborate(self):
        with open('export/genomes_list.txt',"wt") as wout:
            [wout.write('{}\t{}\t{}\n'.format(self.proteomes[p]['tax_id'], self.taxontree.print_full_taxonomy(self.proteomes[p]['tax_id']), dict(self.proteomes[p]['ncbi_ids']).get('GCSetAcc',''))) for p in self.proteomes]

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

    def get_uniref90_nucl_seq(self, item):
        upid, t = item
        gca = dict(self.proteomes[upid]['ncbi_ids']).get('GCSetAcc','')
        seqs_dir = 'data/ncbi/{}/{}/'.format(gca[4:][:6], gca[4:][6:-2])
        in_seq_file = glob.glob(seqs_dir+'*fna*')[0]
        in_gff_file = glob.glob(seqs_dir+'*gff*')[0]

        file_handle, gff_decompressed = tempfile.mkstemp(prefix="gff")
        try:
            seq_parser = SeqIO.parse(gzip.open(in_seq_file,'rt'), "fasta")

            seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(in_seq_file,'rt'), "fasta"))                

            with gzip.open(in_gff_file, 'rb') as gff_gz:
                shutil.copyfileobj(gff_gz, open(gff_decompressed, 'wb'))
        except Exception as e:
            utils.error('Missing genome or GFF annotation for {} / {} '.format(gca, upid))
            return None

        try:
            gff_parser = list(GFF.parse(open(gff_decompressed), base_dict = seq_dict))
        except:
            return None
        os.remove(gff_decompressed)

        nucls = []
        for core, gene_names in t:
            if type(gene_names) == tuple:
                gene_names = [x[1] for x in gene_names]
            else:
                gene_names = (gene_name, )

            try:
                for gene_name in gene_names:
                    for chromosome in gff_parser:
                        for feature in chromosome.features:
                            if feature.type == 'gene' and gene_name in feature.qualifiers.get('Name'):
                                nuc_seq = feature.location.extract(chromosome)
                                nuc_seq.id = nuc_seq.features[0].qualifiers.get('Name',[''])[0]
                                nuc_seq.description = 'UniRef90_{}'.format(core)
                                nucls.append(nuc_seq)
            except:
                utils.info(upid+'\n')
        return nucls

    """
    Use centroid of UniRef90 as plausible marker
    """
    def get_uniref_uniprotkb_from_panproteome(self, panproteome):
        if not terminating.is_set():
            panproteome = pickle.load(open(panproteome, 'rb'))
            possible_markers = Panproteome.rank_genes(panproteome)
            gc = {}
            upid = ""
            nucls = []
            for marker in possible_markers:
                try:
                    if 'UPI' not in marker:
                        if marker in self.uniprot:
                            upids = ["UP{}{}".format("0"*(9-len(str(upid))),upid) for upid in self.uniprot[marker][1]]
                            upid = [x for x in upids if x in self.proteomes and self.proteomes[x]['isReference']]
                            if not len(upid):
                                upid = [x for x in upids if x in self.proteomes and not self.proteomes[x]['upi']]
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
                    utils.info('proteome: {}\nmarker: {}\n'.format(panproteome['tax_id'], marker))

                if len(upid):
                    if upid not in gc:
                        gc[upid] = []
                    gc[upid].append((marker, gene_names))

            try:
                for item in gc.items():
                    res = self.get_uniref90_nucl_seq(item)
                    if res:
                        nucls.extend(res)
                # terminating_mp = dummy.Event()
                # with dummy.Pool(initializer=init_parse, initargs=(terminating_mp, ), processes=len(gc)) as pool:
                #     nucls = [x for x in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=50)]
            except TypeError as e:
                utils.error('Cannot parse panproteome {}'.format(panproteome['tax_id']))
                return panproteome['tax_id']
            else:
                #create export/metaphlan2
                if len(nucls):
                    with open('export/metaphlan2/{}.fasta'.format(panproteome['tax_id']), 'wt') as fasta_markers:
                        [fasta_markers.write(marker.format('fasta')) for marker in nucls]
        else:
            terminating.set():
    # def get_uniprotkb_from_panproteome(self, uniprotkb_id):
    #     pass
    #     cores = Panproteome.find_core_genes(panproteome)
    #     gc = {}
    #     [gc.setdefault(upid, []).append((core_pangene, upkb)) for core_pangene in cores for upkb, upid in panproteome['members'][core_pangene]['copy_number']]

    #     terminating = dummy.Event()
    #     with dummy.Pool(initializer=init_parse, initargs=(terminating, ), processes=100) as pool:
    #         nucls = [x for x in pool.imap_unordered(self.exctract_gene_from_genome, (x for x in gc.items()), chunksize=50)]


    def run_all():
        species = glob.glob('data/pickled/panproteomes/species/90/*')
        terminating = dummy.Event()
        with dummy.Pool(initializer=init_parse, initargs=(terminating, ), processes=150) as pool:
            failed = [x for x in pool.imap_unordered(self.get_uniref_uniprotkb_from_panproteome, species, chunksize=50) if x is not None]