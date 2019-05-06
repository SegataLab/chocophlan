#!/usr/bin/env python3
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Fabio Cumbo (fabio.cumbo@unitn.it')

__date__ = '3 Apr 2019'


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
import _version as version

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
        os.makedirs('{}/{}/functional_annot/'.format(self.config['export_dir'], self.config['exportpath_humann2']), exist_ok=True)
        os.makedirs('{}/{}/panproteomes_fna/'.format(self.config['export_dir'], self.config['exportpath_humann2']), exist_ok=True)
        
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
        seqs_dir = '{}/{}/{}/{}/'.format(self.config['download_base_dir'], self.config['relpath_genomes'], gca[4:][:6], gca[4:][6:])
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
            for NR90, NR50, gene_names in t:
                # if isinstance(gene_names,tuple):
                #     gene_names = [x[1] for x in itertools.chain.from_iterable(gene_names)]
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
                    fna = Fasta(fna_decompressed.name)
                    try:
                        nuc_seq = feature.sequence(fna)
                    except Exception as e:
                        print(taxid)
                        raise e
                    mpa_marker = '{}__{}__{}|{}|{}|UniRef50_{}|{}'.format(taxid, NR90, feature['Name'][0], taxonomy[0], NR90, NR50, len(nuc_seq))
                    record = SeqIO.SeqRecord(seq=Seq(nuc_seq, alphabet=DNAAlphabet()), id=mpa_marker, description='')
                    nucls.append(record)
                else:
                    failed.append(str(NR90))
        if not len(failed): failed = None
        return (nucls, failed, )

    def export_panproteome(self, species_id):
        panproteome = pd.DataFrame.from_csv('{}/{}/{}.txt.bz2'.format(self.config['export_dir'], self.config['panproteomes_stats'], species_id), sep = '\t')
        taxonomy = self.taxontree.print_full_taxonomy(species_id)
        low_lvls = [c.tax_id for c in self.taxontree.taxid_n[species_id].find_clades()]
        gc = {}
        func_annot = {}
        for pangene in panproteome.index:
            pangene = 'UniRef90_'+pangene.split('_')[1]
            pangene_nr50 = self.uniref90_taxid_map[pangene][1][2]
            if pangene in self.uniref90_taxid_map:
                nr90_members = list(filter(lambda x:x[1] in low_lvls, self.uniref90_taxid_map[pangene][0]))
                if not nr90_members:
                    nr90_members = list(filter(lambda x:x[2] and x[0] in self.uniprot, self.uniref90_taxid_map[pangene][0]))
                has_repr = any(x[2] for x in nr90_members)
                nr90_upkb_entry = list(filter(lambda x:x[0] in self.uniprot and x[2] == has_repr, nr90_members))
                if nr90_upkb_entry:
                    upkb_id = nr90_upkb_entry[0][0]
                    upkb_entry = self.uniprot[upkb_id]
                    gene_names = tuple(x[1] for x in upkb_entry[0])
                    upkb_fa = self.uniprot[upkb_id][3:]
                    upkb_fa = [ ','.join('GO:'+str(x).zfill(7) for x in upkb_fa[0]),
                                ','.join('KO:'+str(x).zfill(5) for x in upkb_fa[1]),
                                ','.join(str(x) for x in upkb_fa[2]),
                                ','.join('PF'+str(x).zfill(5) for x in upkb_fa[3]),
                                ','.join(str(x) for x in upkb_fa[4]),
                                ','.join(str(x) for x in upkb_fa[5]) ]
                    func_annot[pangene] = upkb_fa

                    gene_id, upid, tax_id = upkb_entry[:3]
                    upids = ["UP{}".format(str(upid).zfill(9)) for upid in upkb_entry[1]]
                    #Some UniProtKB entries (like P00644) do not have any proteome associteda.
                    if not len(upids):
                        tax_id = upkb_entry[2]
                        if hasattr(self.taxontree.taxid_n[tax_id], 'proteomes'):
                            upids = [up for up in self.taxontree.taxid_n[tax_id].proteomes if upkb_id in self.proteomes[up]['members']]
                            if not len(upids):
                                possible_entry = [entry for entry in nr90_members if entry[1] == tax_id and entry[0] != upkb_id]
                                if (possible_entry):
                                    upids = self.taxontree.taxid_n[tax_id].proteomes
                                    gene_names = tuple( (x[1] for x in self.uniprot[upkb_id][0]) if upkb in self.uniprot 
                                                  else 
                                                  ( tuple( x[2] for x in self.uniparc[upkb]) if upkb in self.uniparc 
                                                   else None ) 
                                                  for upkb, taxid, isRef in possible_entry )
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
                    if len(upid) and len(gene_names):
                        if (upid, species_id, taxonomy) not in gc:
                            gc[(upid, species_id, taxonomy)] = set()
                        gc[(upid, species_id, taxonomy)].add((pangene, pangene_nr50, gene_names))

        with dummy.Pool(processes=100) as pool:
            res = [_ for _ in pool.imap_unordered(self.get_uniref90_nucl_seq, gc.items(), chunksize=10) if _ is not None]

        pangenes, failed = zip(*res)
        pangenes = list(itertools.chain.from_iterable(filter(None,pangenes)))
        
        if failed is not None and not all(x is None for x in failed):
            failed = list(filter(None, failed))
            self.log.warning("{}: Failed to find features for pangene {}".format(species_id, (','.join(failed))))
        if len(pangenes):
            #save tsv of annots
            with bz2.open('{}/{}/functional_annot/{}.tsv.bz2'.format(self.config['export_dir'], self.config['exportpath_humann2'], species_id), 'wt') as functional_annot:
                functional_annot.write('pangene\tGO\tKO\tKEGG\tPfam\tEC\teggNOG\n')
                functional_annot.write('\n'.join( '{}\t{}'.format( k, '\t'.join(v)) for k, v in func_annot.items() ) )

            with bz2.open('{}/{}/panproteomes_fna/{}.fna.bz2'.format(self.config['export_dir'], self.config['exportpath_humann2'], species_id), 'wt') as fasta_pangenes:
                [fasta_pangenes.write(marker.format('fasta')) for marker in pangenes]


if __name__ == '__main__':
    t0 = time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config, verbose=args.verbose)
    config = config['export']

    export_humann2 = chocophlan2humann2(config)
    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)