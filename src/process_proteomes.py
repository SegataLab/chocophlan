#!/usr/bin/env python3

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '7 Nov 2017'

import utils
import os
import gzip
import pickle
import multiprocessing as mp
import glob
import time
import math
import sys
import numpy as np
import h5py
from lxml import etree
import copyreg
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from functools import partial

kingdom_to_process = ['Bacteria','Archaea']

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def init_parse(terminating_, tree_, table_):
    global terminating
    terminating = terminating_
    global tree
    tree = tree_
    global table
    table = table_
def element_unpickler(data):
    return etree.fromstring(data)

def element_pickler(element):
    data = etree.tostring(element)
    return element_unpickler, (data,)


def parse_xml_elem(chunk, config):
    if not terminating.is_set():
        try:
            tag_to_parse = ["accession","name","dbReference","sequence","organism","sequence"]
            dbReference = ["EMBL","Ensembl","GeneID","GO","KEGG","KO","Pfam","Refseq","Proteomes"]
            elem = etree.fromstring(chunk)
            d_prot = {}
            del chunk
            sequence = ""
            for children in tag_to_parse:
                tag_children = "{http://uniprot.org/uniprot}"+children
                if children == "sequence":
                    sequence = [x for x in elem.iterchildren(tag=tag_children)][0].text
                if children == "name" or children == "accession":
                    d_prot[children] = [x for x in elem.iterchildren(tag=tag_children)][0].text
                if children == "organism":
                    org = [x for x in elem.iterchildren(tag_children)][0]
                    lineage = [x for x in org.iterchildren("{http://uniprot.org/uniprot}lineage")][0]
                    kingdom = lineage.getchildren()[0].text
                    d_prot['tax_id'] = [x.get('id') for x in org.iterchildren(tag='{http://uniprot.org/uniprot}dbReference')][0]
                    d_prot['kingdom'] = kingdom
                if children == "dbReference":
                    if children not in d_prot:
                        d_prot[children] = {}
                    for ref in elem.iterchildren(tag=tag_children):
                        if ref.get('type') in dbReference:
                            if ref.get('type') not in d_prot[children]:
                                d_prot[children][ref.get('type')] = []
                            d_prot[children][ref.get('type')].append(ref.get('id').encode('utf-8'))
            elem.clear()
            
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
            
            uniprotkb = table['/uniprotkb'] 
            if not uniprotkb.__contains__(d_prot['accession']):
                protein = uniprotkb.create_group(d_prot['accession'])
            
            protein = uniprotkb[d_prot['accession']]

            protein.attrs['accession'] = d_prot['accession']
            protein.attrs['name'] = d_prot['name']
            protein.attrs['sequence'] = sequence
            protein.attrs['tax_id'] = d_prot['tax_id']

            for ref in d_prot['dbReference'].keys():
                str_size = max(list(map(len,d_prot['dbReference'][ref])))
                if not protein.__contains__(ref):
                    reference = protein.create_dataset(ref, (len(d_prot['dbReference'][ref]),), dtype=np.dtype("S{}".format(str_size)))
                
                protein[ref][:] = d_prot['dbReference'][ref]
            table.flush()
        except Exception as e:
            terminating.set()
            utils.error(str(e))
            raise
    else:
        terminating.set()


def filt(e):
    return ('end' in e[0]) and ('entry' in e[1].tag)

def yield_filtered_xml(tree):
    for _, elem in filter(filt,tree):
        yield etree.tostring(elem)
        elem.clear()

def parse_xml(table, l_input, config):
    for input in l_input:
        if config['verbose']:
            utils.info('Starting processing {} file...\n'.format(input))
        tree = etree.iterparse(input)
        terminating = mp.Event()
        chunksize = int(config['nproc'])
        
        parse_xml_partial = partial(parse_xml_elem, config=config)
        with mp.Pool(initializer=init_parse, initargs=(terminating,tree,table,), processes=int(config['nproc'])) as pool:
            try:
                [ _ for _ in pool.imap(parse_xml_partial, yield_filtered_xml(tree), chunksize=chunksize)]
                del tree
            except Exception as e:
                utils.error(str(e))
                utils.error('Processing failed',  exit=True)
                raise
        if config['verbose']:
            utils.info('Done processing {} file!\n'.format(input))
    if config['verbose']:
        utils.info("Done\n")

def process(input, is_reference=True):
    if not terminating.is_set():
        prot = {}
        with gzip.open(input, 'rt') as idmap:
            for line in idmap:
                uniprotkbid, k, v = line.split()
                if uniprotkbid not in prot:
                    prot[uniprotkbid] = {k:[v]}
                else:
                    if k not in prot[uniprotkbid]:
                        prot[uniprotkbid][k] = [v]
                    else:
                        prot[uniprotkbid][k].append(v)
        protid, taxid = os.path.split(input)[1].split('.')[0].split('_')
        prot = {'proteome_id': protid
               ,'taxid': taxid
               ,'isreference' : is_reference
               ,'proteins': prot }
        return prot
    else:
        terminating.set()
            
def parse_reference_proteomes(config, verbose=False):
    terminating = mp.Event()
    chunksize = math.floor(int(config['nproc']) / 2)

    if not os.path.exists("{}/{}".format(config['download_base_dir'],config['relpath_reference_proteomes'])):
        utils.error('Required files does not exist! Exiting...', exit=True)
    
    if verbose:
        utils.info("Starting parallel proteomes processing...\n",)

    for k in kingdom_to_process:
        basepath = "{}/{}/{}/*.idmapping.gz".format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)

        chunksize = config['nproc']
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=int(config['nproc'])) as pool:
            try:
                proteomes = { item['proteome_id']:item for item in [x for x in pool.imap(process, ls, chunksize = chunksize)]}
            except Exception as ex:
                utils.error(str(ex))
                utils.error('Processing failed',  exit=True)
                raise
    pickle.dump(proteomes, open(config['download_base_dir'] + config['relpath_pickle_proteomes'], "wb" ))
    if verbose:
        utils.info("Done\n")

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']

    choco = utils.createDataset(config['download_base_dir']+config['relpath_chocophlan_database'])
    if not choco.__contains__("/uniprotkb"):
        choco.create_group('/uniprotkb')
    parse_xml(choco, ["data/uniprot/complete/uniprot_sprot.xml", "data/uniprot/complete/uniprot_trembl.xml"], config)
    #parse_xml(choco,["data/uniprot/complete/uniprot_sprot_example.xml"], config)
    choco.close()
    
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
