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
import multiprocessing as mp
import glob
import time
import math
import sys
import numpy as np
import h5py
from lxml import etree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Clade as BClade
from Bio.Phylo.BaseTree import Tree as BTree
from functools import partial

kingdom_to_process = ['Bacteria','Archaea']
dbReference = ['EMBL','Ensembl','GeneID','GO','KEGG','KO','Pfam','Refseq','Proteomes']

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
        f = h5py.File(filepath, 'a', driver='core', backing_store = True, block_size=1024)
        return f
    except Exception as e:
        utils.error(str(e), exit=True)

def parse_uniprotkb_xml_elem(elem, database):
    tag_to_parse = ['accession','name','dbReference','sequence','organism','sequence']
    try:
        d_prot = {}
        sequence = ''
        org = [x for x in elem.iterchildren('{http://uniprot.org/uniprot}organism')][0]
        kingdom = [x for x in org.iterchildren('{http://uniprot.org/uniprot}lineage')][0].getchildren()[0].text
        if kingdom in kingdom_to_process:
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
        
            uniprotkb = database['/uniprotkb'] 
            if d_prot['accession'] not in uniprotkb:
                protein = uniprotkb.create_group(d_prot['accession'])
            
            protein = uniprotkb[d_prot['accession']]
            protein.attrs['accession'] = d_prot['accession']
            protein.attrs['name'] = d_prot['name']
            protein.attrs['sequence'] = sequence
            protein.attrs['tax_id'] = d_prot['tax_id']
    
            for ref in d_prot['dbReference'].keys():
                str_size = max(map(len,d_prot['dbReference'][ref]))
                if ref not in protein:
                    reference = protein.create_dataset(ref, (len(d_prot['dbReference'][ref]),), dtype=np.dtype('S{}'.format(str_size)), compression='gzip')
                
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
            [parse_uniprotkb_xml_elem(b,database) for b in yield_filtered_xml(tree)]
            del tree
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        if config['verbose']:
            utils.info('Done processing {} file!\n'.format(input))
    if config['verbose']:
        utils.info('Done\n')

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
    if not terminating.is_set():
        tag_to_parse = ['member', 'representativeMember', 'property'] 
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
            if not os.path.exists("{}/uniprot".format(config['temp_folder'])):
                os.makedirs("{}/uniprot/".format(config['temp_folder']),exist_ok = True)
                
            accession = d_uniref['id']
            with h5py.File("{}/uniprot/{}".format(config['temp_folder'],accession)) as database:
                #if 'uniref' not in database:
                #    database.create_group('/uniref')
                #for c in ['100','90','50']:
                #    if '/uniref/{}'.format(c) not in database:
                #        database.create_group('/uniref/{}'.format(c))

                #g_uniref = database['/uniref']
            
                accession = d_uniref['id']
                cluster, pid = accession.split('_')
                cluster = cluster.split('UniRef')[1]

                if cluster not in database:
                    database.create_group('/{}'.format(cluster))
                if accession not in database['/{}'.format(cluster)]:
                    database.create_group('/{}/{}'.format(cluster, accession))

                entry = database['/{}/{}'.format(cluster,accession)]
                if 'UniRef100' in d_uniref:
                    entry['UniRef100'] = h5py.SoftLink('/uniref/100/{}'.format(d_uniref['UniRef100']))
                if 'UniRef90' in d_uniref:
                    entry['UniRef90'] = h5py.SoftLink('/uniref/90/{}'.format(d_uniref['UniRef90']))
                if 'UniRef50' in d_uniref:
                    entry['UniRef50'] = h5py.SoftLink('/uniref/50/{}'.format(d_uniref['UniRef50']))

                str_len = max([max(len(k),len(v)) for k,v in d_uniref['members']])
                members = entry.create_dataset('members', (len(d_uniref['members']),2), dtype = np.dtype('S{}'.format(str_len)))
                members[:] = d_uniref['members']

        except Exception as e:
            utils.error(str(e))
            with open("{}/uniprot/FAILED".format(config['temp_folder']),'w+') as out:
                out.write(elem.get('id')+"\n")
            raise
    else:
        terminating.set()

def create_uniref_dataset(database, config, verbose=False):
    uniprot100_xml = etree.iterparse(gzip.GzipFile(config['download_base_dir']+config['relpath_uniref100']))
    uniprot90_xml = etree.iterparse(gzip.GzipFile(config['download_base_dir']+config['relpath_uniref90']))
    uniprot50_xml = etree.iterparse(gzip.GzipFile(config['download_base_dir']+config['relpath_uniref50']))
    
    terminating = mp.Event()
    chunksize = config['nproc']
    uniprot_entries=[]
    
    for tree in [uniprot100_xml, uniprot50_xml, uniprot90_xml]:
        with mp.Pool(initializer=init_parse, initargs=(terminating, tree,), processes=chunksize) as pool:
            try:
                if verbose:
                    util.info("Starting processing UniProt database\n")
                parse_uniref_xml_elem_partial = partial(parse_uniref_xml_elem, config=config)
                [_ for _ in pool.imap(parse_uniref_xml_elem_partial, yield_filtered_xml_string(tree), chunksize= chunksize)]
                #uniprot_entries.append([parse_uniref_xml_elem(x) for x in yield_filtered_xml(tree)])
                if verbose:
                    utils.info("Done\n")
            except Exception as e:
                utils.error(str(e))
                raise
        if verbose:
            utils.info('Merging UniRef files...\n')
        for _, _, f in os.walk("{}/uniprot/".format(config['temp_folder'])):
            with h5py.File(f) as temp:
                cluster = list(temp)[0]
                entry = list(p[cluster])[0]
                temp.copy('{}/{}'.format(cluster,entry), database['/uniref/{}'.format(cluster)])
                os.unlink(f)
        if verbose:
            utils.info('Done\n')
    if verbose:
        utils.info('UniProt database processed successfully.\n')

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']

    with createDataset(config['download_base_dir']+config['relpath_chocophlan_database']) as choco:
        if '/uniprotkb' not in choco:
            choco.create_group('/uniprotkb')
        #parse_uniprotkb_xml(choco, ['data/uniprot/complete/uniprot_sprot.xml', 'data/uniprot/complete/uniprot_trembl.xml'], config)
        #parse_uniprotkb_xml(choco,['data/uniprot/complete/uniprot_sprot_example.xml'], config)
        #parse_uniprotkb_xml(choco,['data/uniprot/complete/uniprot_sprot.xml'], config)
        #create_proteomes(choco,config)
        create_uniref_dataset(None, config)
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(float(t1 - t0)))

    sys.exit(0)
