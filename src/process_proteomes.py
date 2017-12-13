#!/usr/bin/env python3

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '7 Nov 2017'

if __name__ == '__main__':
    import utils
else:
    import src.utils as utils
import os
import gzip
import pickle
import re
import multiprocessing as mp
import glob
import time
import math
import sys
import numpy as np
import h5py
import shutil
from itertools import zip_longest
from lxml import etree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Clade as BClade
from Bio.Phylo.BaseTree import Tree as BTree
from functools import partial

kingdom_to_process = ['Bacteria','Archaea']
dbReference = ['GeneID','Proteomes']
#dbReference = ['EMBL','Ensembl','GeneID','GO','KEGG','KO','Pfam','Refseq','Proteomes']

def filt(e):
    return ('end' in e[0]) and ('entry' in e[1].tag)

def yield_filtered_xml_string(tree):
    for _, elem in filter(filt,tree):
        yield etree.tostring(elem)
        elem.clear()

def yield_filtered_xml(tree):
    for _, elem in filter(filt,tree):
        yield elem
        elem.clear()

def grouper(iterable, n, fillvalue=None):
    #"Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def init_parse(terminating_, tree_,counter_):
    global terminating
    terminating = terminating_
    global tree
    tree = tree_
    global counter
    counter = counter_

def createDataset(filepath):
    try:
        f = h5py.File(filepath, 'a', driver='sec2')#, backing_store = True, block_size=1024)
        return f
    except Exception as e:
        utils.error(str(e), exit=True)

def parse_uniprotkb_xml_elem(elem, config, db):
    if not terminating.is_set() and elem is not None:
        tag_to_parse = ['accession','name','dbReference','sequence','organism','sequence']
        try:
            d_prot = {}
            sequence = ''
            elem = etree.fromstring(elem)
            org = [x for x in elem.iterchildren('{http://uniprot.org/uniprot}organism')][0]
            kingdom = [x for x in org.iterchildren('{http://uniprot.org/uniprot}lineage')][0].getchildren()[0].text
            proteome = [x.get("value") for x in elem.iterchildren('{http://uniprot.org/uniprot}dbReference') if x.get('id') == "Proteomes" ]
            
            if kingdom in kingdom_to_process:
                global counter
                with counter.get_lock():
                    counter.value += 1
                for children in tag_to_parse:
                    tag_children = '{http://uniprot.org/uniprot}'+children
                    if children == 'sequence':
                        sequence = [x for x in elem.iterchildren(tag=tag_children)][0].text
                    if children == 'name' or children == 'accession':
                        d_prot[children] = [x for x in elem.iterchildren(tag=tag_children)][0].text
                    if children == 'organism':
                        org = [x for x in elem.iterchildren(tag_children)][0]
                        d_prot['tax_id'] = [x.get('id') for x in org.iterchildren(tag='{http://uniprot.org/uniprot}dbReference')][0]
                    if children == 'dbReference':
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
                return(d_prot)
        except Exception as e:
            utils.error(str(e))
            raise
    else:
        terminating.set()
    
def parse_uniprotkb_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = config['nproc']
    from multiprocessing import Value
    counter = Value('i',0)

    db=os.path.splitext(os.path.basename(xml_input))[0].split('.')[0]

    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))
    
    tree = etree.iterparse(gzip.GzipFile(xml_input))
    if(os.path.exists("{}/pickled/idmap.pkl".format(config['download_base_dir']))):
        idmap=pickle.load(open("{}/pickled/idmap.pkl".format(config['download_base_dir']),'rb'))
    else:
        idmap={}
    
    parse_uniprotkb_xml_elem_partial = partial(parse_uniprotkb_xml_elem, config=config, db=db)
    with mp.Pool(initializer=init_parse, initargs=(terminating, tree,counter,), processes=chunksize) as pool:
        try:
            file_chunk = 1
            for group in grouper(yield_filtered_xml_string(tree), 1000000):
                d={x['accession']:x for x in pool.imap(parse_uniprotkb_xml_elem_partial, group, chunksize=chunksize) if x is not None}
                idmap.update(dict.fromkeys(d.keys(), file_chunk))
                pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],db, file_chunk),'wb'), -1)
                file_chunk+=1
            pickle.dump(idmap, open("{}/pickled/idmap.pkl".format(config['download_base_dir']),'wb'), -1)
            print(db+str(counter.value))
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        if config['verbose']:
            utils.info('Done processing {} file!\n'.format(xml_input))
    if config['verbose']:
        utils.info('Done processing UniProtKB\n')

def create_proteomes(database, config, verbose = False):
    if verbose:
        utils.info('Starting proteomes processing...\n')
    g_uniprotkb = database['/uniprotkb']
    d_proteomes = {}
    for g_protein in g_uniprotkb.values():
        if 'Proteomes' in g_protein:
            proteome = g_protein['Proteomes'][0]
            accession = g_protein.attrs['accession']
            if proteome not in d_proteomes:
                d_proteomes[proteome] = []
            d_proteomes[proteome].append(accession.encode('utf-8'))
    
    if 'proteomes' not in database:
        database.create_group('/proteomes')
    
    g_proteomes = database['/proteomes']
    
    for proteome in d_proteomes:
        str_size = max((map(len,d_proteomes[proteome])))
        if proteome not in g_proteomes:
            g_proteomes.create_dataset(proteome, (len(d_proteomes[proteome]),), dtype=np.dtype('S{}'.format(str_size)))
        g_proteomes[proteome][:] = d_proteomes[proteome]
        
        g_proteomes[proteome].attrs['tax_id'] = database['/uniprotkb/{}'.format(d_proteomes[proteome][0].decode())].attrs['tax_id']
        g_proteomes[proteome].attrs['isReference'] = False
        
    if not os.path.exists('{}/{}'.format(config['download_base_dir'],config['relpath_reference_proteomes'])):
        utils.error('Required files does not exist! Exiting...', exit=True)
    
    if verbose:
        utils.info('Starting reference proteomes processing...\n',)

    for k in kingdom_to_process:
        basepath = '{}/{}/{}/*.idmapping.gz'.format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)
        for protid, taxid in [os.path.split(file)[1].split('.')[0].split('_') for file in ls]:
            protid = protid.encode('utf-8')
            if protid in g_proteomes:
                g_proteomes[protid].attrs['tax_id'] = taxid
                g_proteomes[protid].attrs['isReference'] = True
    if verbose:
        utils.info('Done\n')

def parse_uniref_xml_elem(elem, config):
    if not terminating.is_set() and elem is not None:
        tag_to_parse = ['member', 'representativeMember', 'property'] 
        global counter
        with counter.get_lock():
            counter.value += 1
        try:
            elem = etree.fromstring(elem)
            d_uniref = {}
            d_uniref['id'] = elem.get('id')
            d_uniref['members'] =[]
            for tag in tag_to_parse:
                tag_children = '{http://uniprot.org/uniref}'+tag
                if tag == 'property':
                    for pro in elem.iterchildren(tag_children):
                        if pro.get('type') == 'common taxon ID':
                            d_uniref['common_taxid'] = pro.get('value')
                if tag == 'representativeMember' or tag == 'member':
                    member = [x for x in elem.iterchildren(tag_children)]
                    properties = [x.iterchildren('{http://uniprot.org/uniref}property') for y in member for x in y.iterchildren('{http://uniprot.org/uniref}dbReference')]
                    for p in properties:
                        for pro in p:
                            if pro.get('type') == 'UniProtKB accession':
                                accession = pro.get('value')
                                isRepr = b"True" if tag == 'representativeMember' else b"False"
                                d_uniref['members'].append((accession.encode('utf-8'),isRepr))
                            if pro.get('type') == 'UniRef100 ID':
                                d_uniref['UniRef100'] = pro.get('value')
                            if pro.get('type') == 'UniRef90 ID':
                                d_uniref['UniRef90'] = pro.get('value')
                            if pro.get('type') == 'UniRef50 ID':
                                d_uniref['UniRef50'] = pro.get('value')
                
            return d_uniref
            
        except Exception as e:
            utils.error(str(e))
            with open("{}/uniref/FAILED".format(config['temp_folder']),'w+') as out:
                out.write(elem.get('id')+"\n")
            raise
    else:
        terminating.set()

def create_uniref_dataset(xml, config):
    uniref_xml = etree.iterparse(gzip.GzipFile(xml))
    
    terminating = mp.Event()
    chunksize = config['nproc']
    from multiprocessing import Value
    counter = Value('i',0)
    cluster = os.path.basename(xml).split('.')[0]

    if(os.path.exists("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'], cluster))):
        idmap=pickle.load(open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster), 'rb'))
    else:
        idmap={}

    file_chunk = 1
    with mp.Pool(initializer=init_parse, initargs=(terminating, uniref_xml, counter,), processes=chunksize) as pool:
        try:
            if config['verbose']:
                utils.info("Starting processing UniRef {} database\n".format(cluster))
            parse_uniref_xml_elem_partial = partial(parse_uniref_xml_elem, config=config)
            for group in grouper(yield_filtered_xml_string(uniref_xml), 1000000):
                d={x['id']:x for x in pool.imap(parse_uniref_xml_elem_partial, group, chunksize=chunksize) if x is not None}
                pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],cluster, file_chunk),'wb'), -1)
                idmap.update(dict.fromkeys(d.keys(), file_chunk))

                file_chunk+=1

        except Exception as e:
            utils.error(str(e))
            raise
    pickle.dump(idmap, open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster),'wb'), -1)
    if config['verbose']:
        utils.info('UniRef {} database processed successfully.\n'.format(cluster))

         

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']

    #parse_uniprotkb_xml(config['download_base_dir']+config['relpath_uniprot_sprot'], config)
#    parse_uniprotkb_xml(config['download_base_dir']+config['relpath_uniprot_trembl'], config)
    processes = [#mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_sprot'], config,)),
    #             mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_trembl'], config,)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref100'],config,)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref90'],config,)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref50'],config,)),

                ]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    #    #create_proteomes(choco,config)
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(float(t1 - t0)))

    sys.exit(0)
