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
from lxml import etree
import copyreg

kingdom_to_process = ['Bacteria','Archaea']

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def init_parse(terminating_, tree_):
    global terminating
    terminating = terminating_
    global queue
    queue = tree_

def element_unpickler(data):
        return etree.fromstring(data)

def element_pickler(element):
    data = etree.tostring(element)
    return element_unpickler, (data,)


def parse_xml_elem(chunk):
    if not terminating.is_set():
        try:
            tag_to_parse = ["accession","name","dbReference","sequence","organism"]
            _ , elem = chunk
            d_prot = {}
            for children in tag_to_parse:
                tag_children = "{http://uniprot.org/uniprot}"+children
                if children == "sequence" or children == "name" or children == "accession":
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
                        if ref.get('type') not in d_prot[children]:
                            d_prot[children][ref.get('type')] = []
                        d_prot[children][ref.get('type')].append(ref.get('id'))
            return(d_prot)
        except Exception as e:
            terminating.set()
            utils.error(str(e))
            raise
    else:
        terminating.set()

def parse_xml(l_input, config):
    kbdb = []
    filt = lambda e : 'end' in e[0] and 'entry' in e[1].tag
    for input in l_input:
        tree = etree.iterparse(input)
        terminating = mp.Event()
        chunksize = int(config['nproc'])
        queue = mp.Queue()
    
        copyreg.pickle(etree._Element, element_pickler, element_unpickler)

        with mp.Pool(initializer=init_parse, initargs=(terminating,queue,), processes=int(config['nproc'])) as pool:
            try:
                kbdb.extend([x for x in pool.imap(parse_xml_elem, filter(filt,tree), chunksize = chunksize)])
            except Exception as e:
                utils.error(str(e), init_new_line=True)
                utils.error('Processing failed', init_new_line=True, exit=True)
                raise
    kbdb = {item['accession']: item for item in kbdb}
    pickle.dump(kbdb, open(config['download_base_dir'] + config['relpath_pickle_uniprotdb'], "wb" ))
    if verbose:
        utils.info("Done\n", init_new_line=True)

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
        utils.info("Starting parallel proteomes processing...\n", init_new_line=True)

    for k in kingdom_to_process:
        basepath = "{}/{}/{}/*.idmapping.gz".format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)

        chunksize = config['nproc']
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=int(config['nproc'])) as pool:
            try:
                proteomes = { item['proteome_id']:item for item in [x for x in pool.imap(process, ls, chunksize = chunksize)]}
            except Exception as ex:
                utils.error(str(ex), init_new_line=True)
                utils.error('Processing failed', init_new_line=True, exit=True)
                raise
    pickle.dump(proteomes, open(config['download_base_dir'] + config['relpath_pickle_proteomes'], "wb" ))
    if verbose:
        utils.info("Done\n", init_new_line=True)

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    
    #parse_reference_proteomes(config['process_proteomes'],config['process_proteomes']['verbose'])
    parse_xml(["../data/uniprot/complete/uniprot_sprot.xml", "../data/uniprot/complete/uniprot_trembl.xml"], config['process_proteomes'])

    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
