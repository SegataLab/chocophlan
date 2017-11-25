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
dbReference = ["EMBL","Ensembl","GeneID","GO","KEGG","KO","Pfam","Refseq","Proteomes"]

def filt(e):
    return ('end' in e[0]) and ('entry' in e[1].tag)

def yield_filtered_xml(tree):
    for _, elem in filter(filt,tree):
        yield elem
        elem.clear()

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def init_parse(terminating_, tree_):
    global terminating
    terminating = terminating_
    global tree
    tree = tree_

def createDataset(filepath):
    try:
        #f = h5py.File(filepath, 'a', driver='stdio')
        f = h5py.File(filepath, 'a', driver='core', backing_store = True, block_size=1024)

        return f
    except Exception as e:
        utils.error(str(e), exit=True)

def parse_uniprotkb_xml_elem(elem, config, database):
    tag_to_parse = ["accession","name","dbReference","sequence","organism","sequence"]
    try:
        d_prot = {}
        sequence = ""
        org = [x for x in elem.iterchildren('{http://uniprot.org/uniprot}organism')][0]
        kingdom = [x for x in org.iterchildren("{http://uniprot.org/uniprot}lineage")][0].getchildren()[0].text
        if kingdom in kingdom_to_process:
            for children in tag_to_parse:
                tag_children = "{http://uniprot.org/uniprot}"+children
                if children == "sequence":
                    sequence = [x for x in elem.iterchildren(tag=tag_children)][0].text
                if children == "name" or children == "accession":
                    d_prot[children] = [x for x in elem.iterchildren(tag=tag_children)][0].text
                if children == "organism":
                    org = [x for x in elem.iterchildren(tag_children)][0]
                    d_prot['tax_id'] = [x.get('id') for x in org.iterchildren(tag='{http://uniprot.org/uniprot}dbReference')][0]
                if children == "dbReference":
                    if children not in d_prot:
                        d_prot[children] = {}
                    for ref in elem.iterchildren(tag=tag_children):
                        if ref.get('type') in dbReference:
                            if ref.get('type') not in d_prot[children]:
                                d_prot[children][ref.get('type')] = []
                            d_prot[children][ref.get('type')].append(ref.get('id').encode('utf-8'))
        
            uniprotkb = database['/uniprotkb'] 
            if not uniprotkb.__contains__(d_prot['accession']):
                protein = uniprotkb.create_group(d_prot['accession'])
            
            protein = uniprotkb[d_prot['accession']]
            protein.attrs['accession'] = d_prot['accession']
            protein.attrs['name'] = d_prot['name']
            protein.attrs['sequence'] = sequence
            protein.attrs['tax_id'] = d_prot['tax_id']
    
            for ref in d_prot['dbReference'].keys():
                str_size = max(map(len,d_prot['dbReference'][ref]))
                if not protein.__contains__(ref):
                    reference = protein.create_dataset(ref, (len(d_prot['dbReference'][ref]),), dtype=np.dtype("S{}".format(str_size)), compression="gzip")
                
                protein[ref][:] = d_prot['dbReference'][ref]

        elem.clear()
        
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    except Exception as e:
        utils.error(str(e))
        raise

def parse_uniprotkb_xml(database, l_input, config):
    for input in l_input:
        if config['verbose']:
            utils.info('Starting processing {} file...\n'.format(input))
        tree = etree.iterparse(input)
        
        try:
            [parse_uniprotkb_xml_elem(b,config,database) for b in yield_filtered_xml(tree)]
            del tree
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        if config['verbose']:
            utils.info('Done processing {} file!\n'.format(input))
    if config['verbose']:
        utils.info("Done\n")

def create_proteomes(database, config, verbose = False):
    if verbose:
        utils.info("Starting proteomes processing...\n")
    g_uniprotkb = database['/uniprotkb']
    d_proteomes = {}
    for g_protein in g_uniprotkb.values():
        if g_protein.__contains__('Proteomes'):
            proteome = g_protein['Proteomes'][0]
            accession = g_protein.attrs['accession']
            if proteome not in d_proteomes:
                d_proteomes[proteome] = []
            d_proteomes[proteome].append(accession.encode('utf-8'))
    
    if not database.__contains__('proteomes'):
        database.create_group('/proteomes')
    
    g_proteomes = database['/proteomes']
    
    for proteome in d_proteomes:
        str_size = max((map(len,d_proteomes[proteome])))
        if not g_proteomes.__contains__(proteome):
            g_proteomes.create_dataset(proteome, (len(d_proteomes[proteome]),), dtype=np.dtype("S{}".format(str_size)))
        g_proteomes[proteome][:] = d_proteomes[proteome]
        

        g_proteomes[proteome].attrs['tax_id'] = database['/uniprotkb/{}'.format(d_proteomes[proteome][0])].attrs['tax_id']
        g_proteomes[proteome].attrs['isReference'] = False
        
    if not os.path.exists("{}/{}".format(config['download_base_dir'],config['relpath_reference_proteomes'])):
        utils.error('Required files does not exist! Exiting...', exit=True)
    
    if verbose:
        utils.info("Starting reference proteomes processing...\n",)

    for k in kingdom_to_process:
        basepath = "{}/{}/{}/*.idmapping.gz".format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)
        for protid, taxid in [os.path.split(file)[1].split('.')[0].split('_') for file in ls]:
            protid = protid.encode('utf-8')
            if g_proteomes.__contains__(protid):
                g_proteomes[protid].attrs['tax_id'] = taxid
                g_proteomes[protid].attrs['isReference'] = True
    if verbose:
        utils.info("Done\n")

def parse_uniref_xml_elem(elem, config, database):
    tag_to_parse = ['member', 'representativeMember', 'property'] 
    try:
        d_uniref = {}
        d_uniref = {}
        d_uniref["id"] = elem.get('id')
        d_uniref['members'] = {}
        for tag in tag_to_parse:
            tag_children = "{http://uniprot.org/uniref}"+tag
            if tag == "property":
                for pro in elem.iterchildren(tag_children):
                    if pro.get('type') == 'common taxon ID':
	    	        d_uniref["common_taxid"] = pro.get('value')
	    if tag == "representativeMember" or tag == "member":
	        member = [x for x in elem.iterchildren(tag_children)]
	        properties = [x.iterchildren('{http://uniprot.org/uniref}property') for y in member for x in y.iterchildren("{http://uniprot.org/uniref}dbReference")]
	        for p in properties:
	    	    for pro in p:
	    	        if pro.get('type') == "UniProtKB accession":
	    		    accession = pro.get('value')
	    		    d_uniref['members'][accession] = {}
	    		    d_uniref['members'][accession]["isRepr"] = True if tag == "representativeMember" else False
                        if pro.get('type') == "UniRef100 ID":
	    		    d_uniref['members'][accession]["UniRef100"] = pro.get('value')
	    		if pro.get('type') == "UniRef90 ID":
	    		    d_uniref['members'][accession]["UniRef90"] = pro.get('value')
	    		if pro.get('type') == "UniRef50 ID":
	    		    d_uniref['members'][accession]["UniRef50"] = pro.get('value')

    except Exception as e:
        utils.error(str(e))
        raise

def create_uniprot_dataset(database, config, verbose=False):
    uniprot100_xml = etree.iterparse(config['download_base_dir']+config['relpath_uniref100'])
    uniprot90_xml = etree.iterparse(config['download_base_dir']+config['relpath_uniref90'])
    uniprot50_xml = etree.iterparse(config['download_base_dir']+config['relpath_uniref50'])
    
    uniprot = BTree()
    uniprot_entries = {}
    
    

    entry = BClade(clades=[], name=name)

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']

    with createDataset(config['download_base_dir']+config['relpath_chocophlan_database']) as choco:
        if not choco.__contains__("/uniprotkb"):
            choco.create_group('/uniprotkb')
        parse_uniprotkb_xml(choco, ["data/uniprot/complete/uniprot_sprot.xml", "data/uniprot/complete/uniprot_trembl.xml"], config)
        #parse_xml(choco,["data/uniprot/complete/uniprot_sprot_example.xml"], config)
        #parse_uniprotkb_xml(choco,["data/uniprot/complete/uniprot_sprot.xml"], config)
        create_proteomes(choco,config)
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(float(t1 - t0)))

    sys.exit(0)
