#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Fabio Cumbo (fabio.cumbo@unitn.it')

__date__ = '14 Dec 2018'


import argparse as ap
import configparser as cp
import datetime
import glob
import gzip
import itertools
import logging
import multiprocessing.dummy as dummy
import os
import pickle
import re
import resource
import shutil
import time
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

from _version import *

if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    import src.utils as utils
    from src.panproteomes import Panproteome

class chocophlan2humann2:
    def __init__(self, config):
        if config['verbose']:
            utils.info('Loading pickled databases...')
        resource.setrlimit(resource.RLIMIT_NOFILE, (131072, 131072))
        self.log = logging.getLogger(__name__)
        ch = logging.FileHandler('CHOCOPhlAn_export2humann_{}.log'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')),'w')
        ch.setLevel(logging.INFO)
        self.log.addHandler(ch)

        self.config = config
        self.uniprot = {}
        self.uniparc = {}
        self.taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']), 'rb'))
        self.proteomes = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']), 'rb'))
        self.exportpath = '{}/{}'.format(self.config['export_dir'], self.config['exportpath_humann2'])
        if not os.path.exists(self.exportpath):
            os.makedirs(self.exportpath)
        
        all_uprot_chunks = filter(re.compile('(trembl|sprot)_[0-9]{1,}.pkl').match, os.listdir('{}/{}/uniprotkb/'.format(config['download_base_dir'], config['pickled_dir'])))
        for i in all_uprot_chunks:
            uniprot_chunk = pickle.load(open('{}/{}/uniprotkb/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniprot.update({k : [v[3],v[4],v[1]] + list(v[8:14]) for k,v in uniprot_chunk.items()})        #geneid, proteomeid, tax_id

        all_uparc_chunks = filter(re.compile('uniparc_[0-9]{4}.pkl').match, os.listdir('{}/{}/uniparc/'.format(config['download_base_dir'], config['pickled_dir'])))
        for i in all_uparc_chunks:
            uniparc_chunk = pickle.load(open('{}/{}/uniparc/{}'.format(config['download_base_dir'], config['pickled_dir'], i),'rb'))
            self.uniparc.update({k : v[1] for k,v in uniparc_chunk.items()})

        self.uniref90_taxid_map = pickle.load(open("{}{}".format(config['download_base_dir'],config['relpath_pickle_uniref90_taxid_idmap']), 'rb'))
        if config['verbose']:
            utils.info('Finished.\n')

    def get_all_species_panproteomes(self):
        return [x.split('/')[-1][:-4] for x in glob.glob('{}/{}/species/90/*'.format(config['download_base_dir'], config['relpath_panproteomes_dir']))]


    def get_uniref90_nucl_seq(self, item):
        (upid, taxid, taxonomy), t = item
        ncbi_ids = dict(self.proteomes[upid].get('ncbi_ids',[(None,None)])).get('GCSetAcc',None)
        if ncbi_ids is not None:
            gca = ncbi_ids.split('.')[0]
        else:
            self.log.error('No genome associated to proteome {} available for species {}'.format(upid, taxid))
            return [None, None, None]
        seqs_dir = '{}/{}/{}/{}/'.format(config['download_base_dir'], config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
        if not os.path.exists(seqs_dir):
            self.log.error('For proteome {} / species {} : No genome and GFF annotations downloaded for {}'.format(upid, taxid, gca))
            return [None, None, None]
        try:
            in_seq_file = glob.glob(seqs_dir+'*.fna*')[0]
            in_gff_file = glob.glob(seqs_dir+'*.gff*')[0]
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
                else:
                    failed.append(str(core))
        if not len(failed): failed = None
        return (nucls, failed, )

    def export_panproteome(self, species_id):
        panproteome = pd.DataFrame.from_csv('{}/{}/{}.txt.bz2'.format(self.config['export_dir'], self.config['panproteomes_stats'], species_id), sep = '\t')
        taxonomy = self.taxontree.print_full_taxonomy(species_id)
        gc = {}
        func_annot = {}
        for pangene in panproteome.index:
            uniprotkb_entries = []
            pangene = 'UniRef90_'+pangene.split('_')[1]
            if pangene in self.uniref90_taxid_map:
                for upkb, ncbi_tax_id, isRef in self.uniref90_taxid_map[pangene][0]:
                    if self.taxontree.go_up_to_species(ncbi_tax_id) == species_id and isRef:
                        uniprotkb_entries.append(upkb)

                for upkb in filter(re.compile('^(?!UPI)').match, uniprotkb_entries):
                    gene_id, upid, tax_id = self.uniprot[upkb]
                    upid = "UP{}".format(str(upid).zfill(upid))
                    if (upid, panproteome['tax_id'], taxonomy) not in gc:
                        gc[(upid, panproteome['tax_id'], taxonomy)] = set()
                    gc[(upid, panproteome['tax_id'], taxonomy)].add((upkb, gene_id))
                    func_annot.setdefault(upkb,[]).append(self.uniprot[upkb][3:])
                else:
                    upids = [("UP{}".format(upid.zfill(upid)), gene_id) for upid, tax_id, gene_id in self.uniparc[upkb] if tax_id == species_id ]
                    for upid, gene_id in upids:
                        if (upid, panproteome['tax_id'], taxonomy) not in gc:
                            gc[(upid, panproteome['tax_id'], taxonomy)] = set()
                        gc[(upid, panproteome['tax_id'], taxonomy)].add((upkb, gene_id))

                with dummy.Pool(processes=3) as pool:
                    res = [_ for _ in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=10) if _ is not None]

                pangenes, failed = zip(*res)
                pangenes = list(itertools.chain.from_iterable(filter(None,pangenes)))
            
            if failed is not None and not all(x is None for x in failed):
                failed = list(filter(None, failed))[0]
                self.log.warning("{}: Failed to find features for pangene {}".format(pp_tax_id, (','.join(failed))))
            if len(pangenes):
                #save tsv of annots
                func_annot
                with open('{}/{}/markers_fasta/{}.fna'.format(self.config['export_dir'], self.config['exportpath_humann2'], panproteome['tax_id']), 'wt') as fasta_pangenes:
                    [fasta_pangenes.write(marker.format('fasta')) for marker in pangenes]


if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['export']

    c = chocophlan2humann2(config)
    c.chocophlan2humann2()
    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
