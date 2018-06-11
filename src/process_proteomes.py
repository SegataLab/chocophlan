#!/usr/bin/env python3

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '7 Nov 2017'

if __name__ == '__main__':
    import utils
    import extract
else:
    import src.utils as utils
    import src.extract as extract
import os
import gzip
import pickle
import re
import multiprocessing as mp
import multiprocessing.dummy as mpdummy
import glob
import time
import sys
import requests
import asyncio
import concurrent.futures
from itertools import zip_longest
from lxml import etree
from functools import partial
# from blist import sorteddict
import traceback

kingdom_to_process = ['Bacteria','Archaea', 'Eukaryota']
genus_to_process = ['Candida', 'Blastocystis', 'Saccharomyces']
dbReference = ['EMBL','EnsemblBacteria','GeneID','GO','KEGG','KO','Pfam','Proteomes', 'EC', 'eggNOG']
GROUP_CHUNK = 1000000
API_SIZE = 10

taxid_to_process = set()
uniprotkb_uniref_idmap = {}


def grouper(iterable, n, fillvalue=None):
    #"Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def clear_node(elem):
    elem.clear()
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]

def yield_filtered_xml_string(tree):
    for _, elem in tree:
        yield etree.tostring(elem)
        elem.clear()

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def taxon_to_process(config, taxontree):
    global taxid_to_process
    d_ranks = taxontree.lookup_by_rank()
    taxons = []
    with open(config['relpath_taxon_to_process'],'rt') as taxa:
        for line in taxa:
            name, rank = line.strip().split('\t')
            taxons.extend([x.tax_id for x in d_ranks[rank] if name == x.name])

    [taxid_to_process.update([x.tax_id for x in taxontree.taxid_n[t].get_nonterminals()]) for t in taxons]
    [taxid_to_process.update([x.tax_id for x in taxontree.taxid_n[t].get_terminals()]) for t in taxons]

def is_taxon_processable(tax_id):
    global taxid_to_process
    return tax_id in taxid_to_process

def uniprot_tuple_to_dict(v):
    return {'accession' : v[0],
            'tax_id' : v[1],
            'sequence' : v[2],
             'RefSeq' : v[3],
             'Proteomes' : tuple("UP{}{}".format("0"*(9-len(str(x))),x) for x in v[4]),
             'EMBL' : v[5],
             'EnsemblBacteria' : v[6],
             'GeneID' : v[7],
             'GO' : tuple("GO:{}{}".format('0'*(7-len(str(x))),x) for x in v[8]),
             'KO' : tuple("K{}{}".format("0"*(5-len(str(x))),x) for x in v[9]),
             'KEGG' : v[10],
             'Pfam' : tuple("PF{}{}".format('0'*(5-len(str(x))),x) for x in v[11]),
             'EC' : v[12],
             'EggNOG': v[13],
             'UniRef100' : v[14],
             'UniRef90' : v[15],
             'UniRef50' : v[16]
            }

def parse_uniref_xml_elem(elem, config):
    if not terminating.is_set() and elem is not None:
        elem = etree.fromstring(elem)
        tag_to_parse = ['member', 'representativeMember', 'property'] 
        try:
            d_uniref = {}
            d_uniref['id'] = elem.get('id')
            d_uniref['members'] =[]
            d_uniref['sequence'] = elem.get('sequence')
            for tag in tag_to_parse:
                tag_children = '{http://uniprot.org/uniref}'+tag
                if tag == 'property':
                    for pro in elem.iterchildren(tag_children):
                        if pro.get('type') == 'common taxon ID':
                            d_uniref['common_taxid'] = int(pro.get('value'))
                if tag == 'representativeMember' or tag == 'member':
                    member = [x for x in elem.iterchildren(tag_children)]
                    ref = [x for y in member for x in y.iterchildren('{http://uniprot.org/uniref}dbReference')]
                    for r in ref:
                        tax_id = 0
                        if r.get('type') == 'UniProtKB ID':
                            accession = [y.get('value') for y in r if y.get('type') == 'UniProtKB accession'][0]
                            uniparc_cluster = [y.get('value') for y in r if y.get('type') == 'UniParc ID'][0]
                        elif r.get('type') == 'UniParc ID':
                            accession = r.get('id')
                            uniparc_cluster = None
                        isRepr = True if tag == 'representativeMember' else False

                        for pro in r:
                            if pro.get('type') == 'UniRef100 ID':
                                d_uniref['UniRef100'] = pro.get('value')[10:]
                            elif pro.get('type') == 'UniRef90 ID':
                                d_uniref['UniRef90'] = pro.get('value')[9:]
                            elif pro.get('type') == 'UniRef50 ID':
                                d_uniref['UniRef50'] = pro.get('value')[9:]
                            elif pro.get('type') == 'NCBI taxonomy':
                                tax_id = int(pro.get('value'))
                        d_uniref['members'].append((accession,tax_id, isRepr))
                        if uniparc_cluster is not None:
                            d_uniref['members'].append((uniparc_cluster,tax_id, isRepr))
            t_uniref = (d_uniref['id'],
                        d_uniref.get('common_taxid',None),
                        tuple(d_uniref['members']),
                        d_uniref.get('UniRef100',''),
                        d_uniref.get('UniRef90', ''),
                        d_uniref.get('UniRef50',''),
                        d_uniref.get('sequence','')
                        )
                    
            elem.clear()
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
            return t_uniref
            
        except Exception as e:
            utils.error('Failed to elaborate item: '+ elem.get('id'))
            raise
    else:
        terminating.set()

def create_uniref_dataset(xml, config):  
    uniref_xml = etree.iterparse(gzip.GzipFile(xml), events = ('end',), tag = '{http://uniprot.org/uniref}entry', huge_tree = True)

    terminating = mp.Event()
    chunksize = config['nproc']
    cluster = os.path.basename(xml).split('.')[0]
    if config['verbose']:
        utils.info("Starting processing UniRef {} database\n".format(cluster))

    idmap={}
    taxon_map = {}
    upids = []

    parse_uniref_xml_elem_partial = partial(parse_uniref_xml_elem, config=config)
    with mp.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool:
        try:
            for file_chunk, group in enumerate(grouper(yield_filtered_xml_string(uniref_xml), GROUP_CHUNK),1):
                d={x[0]:x for x in pool.imap_unordered(parse_uniref_xml_elem_partial, group, chunksize=chunksize) if x is not None}
                upids.extend([(m[0],c[0]) for c in d.values() for m in c[2] if 'UPI' in m[0]])
                pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],cluster, file_chunk),'wb'), -1)
                idmap.update(dict.fromkeys(d.keys(), file_chunk))
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing of {} failed.'.format(xml_input))
            raise
    with open("{}/pickled/{}_uniparc_idmap.pkl".format(config['download_base_dir'],cluster),'wb') as pickle_uniref_uniparc_idmap:
        pickle.dump(upids, pickle_uniref_uniparc_idmap, -1)

    with open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster),'wb') as pickle_uniref_idmap:
        pickle.dump(idmap, pickle_uniref_idmap, -1)

    with open("{}/pickled/{}_taxid_idmap.pkl".format(config['download_base_dir'],cluster),'wb') as pickle_taxid_map:
        pickle.dump(taxon_map, pickle_taxid_map, -1)
    if config['verbose']:
        utils.info('UniRef {} database processed successfully.\n'.format(cluster))

# UniProtKB elemtes are mapped in a int-based id tuple
#  0 accession
#  1 tax_id
#  2 sequence
#  3 RefSeq
#  4 Proteomes
#  5 EMBL
#  6 Ensembl
#  7 GeneID
#  8 GO
#  9 KO
# 10 KEGG
# 11 Pfam
def parse_uniprotkb_xml_elem(elem, config):
    if not terminating.is_set() and elem is not None:
        global uniprotkb_uniref_idmap, taxid_to_process
        elem = etree.fromstring(elem)
        nsprefix = '{http://uniprot.org/uniprot}'
        tag_to_parse = ['accession','dbReference','sequence','organism','sequence']
        try:
            d_prot = {}
            sequence = ''
            taxid = int(elem.find('.//{}organism/{}dbReference[@type="NCBI Taxonomy"]'.format(nsprefix,nsprefix)).get('id'))
            if is_taxon_processable(taxid):
                d_prot['accession'] =  elem.find('.//{}accession'.format(nsprefix)).text
                d_prot['sequence'] = elem.find('.//{}sequence'.format(nsprefix)).text
                d_prot['gene'] = [(name.get('type') ,name.text) for gene in elem.findall('.//{}gene'.format(nsprefix)) if gene is not None for name in gene]
                d_prot['tax_id'] = taxid
                for ref in elem.iterchildren(tag='{http://uniprot.org/uniprot}dbReference'):
                    if ref.get('type') in dbReference:
                        if ref.get('type') not in d_prot:
                            d_prot[ref.get('type')] = []
                        d_prot[ref.get('type')].append(ref.get('id'))

                t_prot = (d_prot.get('accession'),  #0
                          d_prot.get('tax_id'),      #1
                          d_prot.get('sequence',''),   #2
                          tuple(d_prot.get('gene','')),     #3
                          tuple(int(x[2:]) for x in d_prot.get('Proteomes',[]) if x is not None),  #4 UP000005640
                          tuple(d_prot.get('EMBL','')),       #5
                          tuple(d_prot.get('EnsemblBacteria','')),    #6
                          tuple(d_prot.get('GeneID','')),     #7
                          tuple(int(x[3:]) for x in d_prot.get('GO',[]) if x is not None),         #8 GO:0016847
                          tuple(int(x[1:]) for x in d_prot.get('KO',[]) if x is not None),         #9 K09972
                          tuple(d_prot.get('KEGG','')),       #10
                          tuple(int(x[2:]) for x in d_prot.get('Pfam',[]) if x is not None),        #11
                          tuple(d_prot.get('EC','')),
                          tuple(d_prot.get('EggNOG','')),
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[0],
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[1],
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[2]
                         )

                elem.clear()
                for ancestor in elem.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]

                return(t_prot)

        except Exception as e:
            utils.error('Failed to elaborate item: '+ elem.find('.//{}accession'.format(nsprefix)).text)
            raise
    else:
        terminating.set()
    
def parse_uniprotkb_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = int(config['nproc'])

    db= -1 if 'uniprot_sprot' in xml_input else +1
    db_name = 'sprot' if db == -1 else 'trembl'

    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))
    
    tree = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = '{http://uniprot.org/uniprot}entry', huge_tree = True)

    idmap = {}
    d = {}
    parse_uniprotkb_xml_elem_partial = partial(parse_uniprotkb_xml_elem, config=config)
    
    with mpdummy.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool:
        for file_chunk, group in enumerate(grouper(yield_filtered_xml_string(tree), GROUP_CHUNK),1):
            d = {x[0]:x for x in pool.imap_unordered(parse_uniprotkb_xml_elem_partial, group, chunksize=int(GROUP_CHUNK/chunksize)) if x is not None}
            try:
                if len(d):
                    idmap.update(dict.fromkeys(d.keys(), db*file_chunk))
                    with open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],db_name, file_chunk),'wb') as pickle_chunk:
                        pickle.dump(d, pickle_chunk, -1)
            except Exception as e:
                utils.error(str(e))
                print(traceback.format_exc())
                utils.error('Processing of {} failed.'.format(xml_input),  exit=True)
                raise
    if len(idmap):
        with open("{}/pickled/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'],db_name),'wb') as pickle_idmap:
            pickle.dump(idmap, pickle_idmap, -1)
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

def parse_uniparc_xml_elem(elem, config):
    global uniprotkb_uniref_idmap, taxid_to_process
    nsprefix = '{http://uniprot.org/uniparc}'
    if elem is not None:
        elem = etree.fromstring(elem)
        try:
            d_prot = {}
            d_prot['accession'] = elem.find('.//{}accession'.format(nsprefix)).text
            d_prot['sequence'] = elem.find('.//{}sequence'.format(nsprefix)).text.replace('\n','')
            active_entries = elem.findall('.//{}dbReference[@active="Y"]'.format(nsprefix))
            d_prot['uniprotkb_ids'] = list(set(x.get('id') for x in active_entries if 'UniProtKB' in x.get('type')))
            entry_with_protid = [x for x in active_entries for y in x.iterchildren('{http://uniprot.org/uniparc}property') if y.get('type') == 'proteome_id']
            d_prot['Proteomes'] = set()
            for entry in entry_with_protid:
                gene_name, proteome_id = '', ''
                tax_id = None
                for p in entry:
                    if p.get('type') == 'proteome_id':
                        proteome_id = int(p.get('value')[2:])
                    elif p.get('type') == 'NCBI_taxonomy_id':
                        tax_id = int(p.get('value'))
                    elif p.get('type') == 'gene_name':
                        gene_name = p.get('value')
                if tax_id is not None and is_taxon_processable(tax_id):
                    d_prot['Proteomes'].add((proteome_id, tax_id, gene_name))

            t_prot = (d_prot.get('accession'),  
                      tuple(d_prot.get('Proteomes',[])),
                      d_prot.get('sequence'),
                      uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[0],
                      uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[1],
                      uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[2],
                      d_prot.get('uniprotkb_ids',[])
                     )
            elem.clear()
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]

            ## If entry does not belong to a proteome, just throw it, wont need.
            if not len(t_prot[1]):
                return None
                
            return(t_prot)
        except Exception as e:
            utils.error('Failed to elaborate item: '+ elem.findall('.//{}accession'.format(nsprefix)).text)
            raise

def parse_uniparc_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = config['nproc']
    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))

    # tree = etree.iterparse(xml_input, events = ('end',), tag = '{http://uniprot.org/uniparc}entry', huge_tree = True)
    tree = etree.iterparse(gzip.open(xml_input), events = ('end',), tag = '{http://uniprot.org/uniparc}entry', huge_tree = True)
    idmap = {}
    parse_uniparc_xml_elem_partial = partial(parse_uniparc_xml_elem, config=config)
    for file_chunk, group in enumerate(grouper(yield_filtered_xml_string(tree), GROUP_CHUNK),1000):
        with mpdummy.Pool(initializer=initt, initargs=(terminating,), processes=int(chunksize)) as pool:
            d={x[0]:x for x in pool.imap_unordered(parse_uniparc_xml_elem_partial, group, chunksize=int(GROUP_CHUNK/chunksize)) if x is not None}
        try:
            if len(d):
                idmap.update(dict.fromkeys(d.keys(), file_chunk))
                with open("{}/pickled/uniparc_{}.pkl".format(config['download_base_dir'], file_chunk),'wb') as pickle_chunk:
                    pickle.dump(d, pickle_chunk, -1)
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing of {} failed.'.format(xml_input),  exit=True)
            raise
    pickle.dump(idmap, open("{}/pickled/uniparc_idmap.pkl".format(config['download_base_dir']),'wb'), -1)
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

def get_members(line):
    if not terminating.is_set():
        line = line.split('\t')
        u100 = line[7].split('_')[1] if len(line[7]) else '' 
        u90 =  line[8].split('_')[1] if len(line[8]) else '' 
        u50 =  line[9].split('_')[1] if len(line[9]) else '' 
        return (line[0], (u100, u90, u50))
    else:
        terminating.set()

### 0    UniProtKB-AC
### 1    UniProtKB-ID
### 2    GeneID (EntrezGene)
### 3    RefSeq
### 4    GI
### 5    PDB
### 6    GO
### 7    UniRef100
### 8    UniRef90
### 9    UniRef50
### 10   UniParc
### 11   PIR
### 12   NCBI-taxon
### 13   MIM
### 14   UniGene
### 15   PubMed
### 16   EMBL
### 17   EMBL-CDS
### 18   Ensembl
### 19   Ensembl_TRS
### 20   Ensembl_PRO
### 21   Additional PubMed
def parse_uniprotkb_uniref_idmapping(config):
    if config['verbose']:
        utils.info('Starting processing idmapping.gz\n')
    input_file = config['download_base_dir'] + config['relpath_idmapping']
    terminating = mp.Event()
    chunksize = config['nproc']
    idmapping = {}

    with gzip.open(input_file, 'rt', encoding='utf-8') as f:
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool:
            idmapping = {x[0]:x[1] for x in pool.imap_unordered(get_members, f)}
    with open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb') as pickle_idmapping:
        pickle.dump(idmapping, pickle_idmapping, -1)
    if config['verbose']:
        utils.info('Done processing idmapping.gz\n')

def get_uniprotkb_entry(idmap, accession,config):
    fn = idmap[accession]
    if fn < 0:
        db = 'sprot'
    elif fn < 1000:
        db = 'trembl'
    else:
        db = 'uniparc'
    chunk = pickle.load(open('{}/pickled/{}_{}.pkl'.format(config['download_base_dir'],db,abs(fn)),'rb'))
    return uniprot_tuple_to_dict(chunk[accession])

def process(f):
    if not terminating.is_set():
        x = pickle.load(open(f,'rb'))
        d_proteome = {}
        try:
            if 'uniparc' not in f:
                for k,v in x.items():
                    if v[4] is not None:
                        proteomes = tuple("UP{}{}".format("0"*(9-len(str(x))),x) for x in v[4])
                        accession = k
                        taxid = v[1]
                        for proteome in proteomes:
                            if proteome not in d_proteome:
                                d_proteome[proteome] = {"isReference" : False, "members" : [], 'tax_id' : taxid, 'upi' : False}
                            d_proteome[proteome]['members'].append(accession)
            else:
                for entry in x.values():
                    for proteome, taxid in (("UP{}{}".format("0"*(9-len(str(upi))),upi),taxid) for upi, taxid, _ in entry[1]):
                        if proteome not in d_proteome:
                            d_proteome[proteome] = {'members' : [], 'isReference' : False, 'tax_id' : taxid, 'upi' : True}
                        d_proteome[proteome]['members'].append(entry[0])
        except Exception as e:
            utils.error(f)
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        return d_proteome
    else:
        terminating.set()

def parse_proteomes_xml_elem(elem, config):
    if not terminating.is_set() and elem is not None:
        elem = etree.fromstring(elem)
        try:
            d_prot = {}
            taxid = int(elem.attrib['taxonomy'])
            upi = True if len(list(elem.iterchildren(tag='redundantTo'))) else False
            if is_taxon_processable(taxid):
                accession = elem.attrib['upid']
                d_prot[accession] = {}
                d_prot[accession]['members'] = []
                d_prot[accession]['upi'] = upi
                d_prot[accession]['tax_id'] = taxid
                d_prot[accession]['isReference'] = True if elem.getchildren()[2].text == 'true' else False
                d_prot[accession]['ncbi_ids'] = [(c.get('type'),c.get('id')) for c in elem.iterchildren(tag='dbReference')]
                if not upi:
                    d_prot[accession]['members'] = [x.get('accession') for k in elem.iterchildren(tag='component') for x in k.iterchildren(tag='protein')]
                elem.clear()
                for ancestor in elem.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]
                return(d_prot)
        except Exception as e:
            utils.error('Failed to elaborate item: '+str(e))
            utils.error(elem.attrib['upid'])
            raise
    else:
        terminating.set()

def get_request(requestURL, session=None):
    r = session.get(requestURL, timeout = None)
    return r

# def get_upp(accession):
#     baseURL = "https://www.ebi.ac.uk/proteins/api/uniparc/proteome/{}?offset={}&size={}&rfActive=true"
#     requestURL = baseURL.format(accession, 0, API_SIZE)
#     entries = []
#     print('Downloading ' + accession)
#     with requests.Session() as s:
#         s.headers.update({ "Accept" : "application/xml"})
#         r = s.get(requestURL, timeout = None)
#         entries.append(r.text)
        
#         total_records = int(r.headers['x-pagination-totalrecords'])
#         step = int(r.links['self']['url'].split('?')[1].split('&')[1].split('=')[1])
#         first_offset = int(r.links['self']['url'].split('?')[1].split('&')[0].split('=')[1])
#         last_offset = int(r.links['last']['url'].split('?')[1].split('&')[0].split('=')[1])

#         if last_offset != first_offset:
#             requestURLs = [baseURL.format(accession,offset,step) for offset in range(step, last_offset+step, step)]
#             get_request_p = partial(get_request, session=s)
#             with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
#                 [entries.append(r.text) for r in executor.map(get_request_p, requestURLs, chunksize=10)]
#         d = []
#         if sum([e.count('<entry') for e in entries]) == total_records:
#             for entry in entries:
#                 entry = etree.fromstring(entry.encode())
#                 d.extend([parse_uniparc_xml_elem(etree.tostring(e), config) for e in entry])
#         else:
#             pickle.dump(entries, open('{}.pkl'.format(accession)))
#             utils.error('Failed to elaborate proteome: '+accession)
#             raise
#     print('Downloaded ' + accession)
#     return (accession, d)
    
# async def download(proteomes):
#     with concurrent.futures.ThreadPoolExecutor(max_workers=15) as executor:
#         loop = asyncio.get_event_loop()
#         futures = [loop.run_in_executor(executor, get_upp, k) for k,v in proteomes.items() if v['upi']]
#         results = [x for x in await asyncio.gather(*futures)]
#         utils.info('Downloaded.')
#         return results

def create_proteomes(xml_input, config):
    if config['verbose']:
        utils.info('Starting proteomes processing...\n')
    d_proteomes = {}
    
    terminating = mp.Event()
    chunksize = config['nproc']
    
    # r = re.compile('.*(trembl|sprot|uniparc).*')
    r = re.compile('.*(uniparc).*.pkl')
    chunks = []

    tree = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = 'proteome', huge_tree = True)
    parse_proteomes_xml_elem_partial = partial(parse_proteomes_xml_elem, config = config)
    try:
        with mpdummy.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool:
            chunks=[x for x in pool.imap_unordered(parse_proteomes_xml_elem_partial, yield_filtered_xml_string(tree), chunksize=chunksize)]
    except Exception as e:
        utils.error(str(e))
        utils.error('Processing failed')
        raise
       
    d_proteomes = { k: {'members' : v.get('members', set()), 
                        'isReference' : v['isReference'], 
                        'tax_id' : v['tax_id'], 
                        'upi' : v['upi'], 
                        'ncbi_ids' : v['ncbi_ids']}
                    for k,v in ( i for chunk in chunks if chunk is not None for i in chunk.items() )}
                    
    try:
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool:
            chunks = [x for x in pool.imap_unordered(process, 
                        (x for x in filter(r.match,glob.iglob("{}/pickled/*".format(config['download_base_dir']))))
                          , chunksize = chunksize)]
    except Exception as e:
        utils.error(str(e))
        utils.error('Processing failed')
        raise

    for k,v in (i for chunk in chunks if chunk is not None for i in chunk.items()):
        if k not in d_proteomes:
            d_proteomes[k] = {'members' : set(), 
                              'isReference' : v['isReference'], 
                              'tax_id' : v['tax_id'], 
                              'upi' : v['upi']}
        if d_proteomes[k]['upi'] and v['upi']:
            d_proteomes[k]['members'] = set(d_proteomes[k]['members'])
            d_proteomes[k]['members'].update(v['members'])
    
    with open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'wb') as pickle_proteomes:
        pickle.dump(d_proteomes, pickle_proteomes, -1)
    if config['verbose']:
        utils.info('Done\n')

def annotate_taxon_tree(config):
    utils.info('Starting annotation of taxonomic tree with proteome ids\n')

    with open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'rb') as pickle_proteomes:
        proteomes = pickle.load(pickle_proteomes)
    with open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb') as pickle_taxontree:
        taxontree = pickle.load(pickle_taxontree)

    for protid, v in proteomes.items():
        taxid = v['tax_id']
        try:
            clade = taxontree.taxid_n[int(taxid)]
            if not hasattr(clade,'proteomes'):
                clade.proteomes = set()
            clade.proteomes.add(protid)
        except:
            print(taxid)
    # taxontree.are_leaves_trimmed = False
    with open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_taxontree']),'wb') as pickle_taxontree:
        pickle.dump(taxontree, pickle_taxontree, -1)
    utils.info('Done.')

def uniref_tuple_to_dict(v):
    return {'id': v[0],
            'common_taxid': v[1],
            'members': v[2],
            'UniRef100': 'UniRef100_{}'.format(v[3]) if len(v[3])>0 else v[0],
            'UniRef90': 'UniRef90_{}'.format(v[4]) if len(v[4])>0 else v[0],
            'UniRef50': 'UniRef50_{}'.format(v[5]) if len(v[5])>0 else v[0]
           }

def merge_uniparc_idmapping(config):
    if config['verbose']:
        utils.info('Started merging UniParc-UniRef idmapping.\n')
    uniparc_idmapping = {}
    clusters = ('uniref100','uniref90','uniref50')
    for c in clusters:
        cluster_idmapping = pickle.load(open("{}/pickled/{}_uniparc_idmap.pkl".format(config['download_base_dir'],c),'rb'))
        
        for upi, cluster in cluster_idmapping:
            if upi not in uniparc_idmapping:
                uniparc_idmapping[upi] = []
            uniparc_idmapping[upi].append(cluster)
        os.unlink("{}/pickled/{}_uniparc_idmap.pkl".format(config['download_base_dir'],c))

    for upi, clusters in uniparc_idmapping.items():
        ur100 = [x.split('_')[1] for x in clusters if 'UniRef100_' in x]
        ur90 = [x.split('_')[1] for x in clusters if 'UniRef90_' in x]
        ur50 = [x.split('_')[1] for x in clusters if 'UniRef50_' in x]
        clusters = [ur100[0] if len(ur100) else '', ur90[0] if len(ur90) else '', ur50[0] if len(ur50) else '']
                   
        uniparc_idmapping[upi] = tuple(clusters)
    
    idmapping = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'rb'))
    idmapping.update(uniparc_idmapping)
    #update to global
    with open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb') as f_out:
        pickle.dump(idmapping, f_out, -1)

    if config['verbose']:
        utils.info('Done.\n')

def get_uniref_entry(config, accession):
    if 'UniRef100' in accession:
        cluster = 'uniref100'
    elif 'UniRef90' in accession:
        cluster = 'uniref90'
    else:
        cluster = 'uniref50'
    idmap = pickle.load(open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster),'rb'))
    fn = idmap[accession]
    chunk = pickle.load(open('{}_{}.pkl'.format(cluster,fn),'rb'))
    return uniref_tuple_to_dict(chunk[accession])

def merge_idmap(config):
    x = ["{}/pickled/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'], 'sprot'), 
         "{}/pickled/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'], 'trembl'),
         "{}/pickled/{}_idmap.pkl".format(config['download_base_dir'], 'uniparc')]
    utils.info('Merging UniProtKB and UniParc ids...\n')

    idmap = {}
    for fin in x:
        with open(fin, 'rb') as pickled:
            idmap.update(pickle.load(pickled))
        os.unlink(fin)

    pickle.dump(idmap,open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_uniprotkb_idmap']), 'wb'), -1)
    
    utils.info('Done merging UniProtKB ids.\n')

def extract_uniref_taxonomy(fp):
    u = pickle.load(open(fp,'rb'))
    utils.info(fp+'\n')
    ret = {}
    for k, v in u.items():
        ret[k] = set(x[1] for x in v[2])
    return ret

def export_uniref_taxonomy():
    # from src.utils import *
    # from src.panproteomes import *
    # from blist import sorteddict

    with mpdummy.Pool(30) as pool:
        # u50 = [x for x in pool.imap_unordered(extract_uniref_taxonomy, ['data/pickled/uniref50_{}.pkl'.format(x) for x in range(1,32)], chunksize = 10)]
        # u90 = [x for x in pool.imap_unordered(extract_uniref_taxonomy, ['data/pickled/uniref90_{}.pkl'.format(x) for x in range(1,75)], chunksize = 10)]
        u100 = [x for x in pool.imap_unordered(extract_uniref_taxonomy, ['data/pickled/uniref100_{}.pkl'.format(x) for x in range(1,144)], chunksize = 10)]
        [u.update(x) for x in u100]
def process_proteomes(config):
    os.makedirs('{}/pickled'.format(config['download_base_dir']), exist_ok=True)
    os.makedirs('{}/pickled/uniparc_entries'.format(config['download_base_dir']), exist_ok=True)

    # create_uniref_dataset(config['download_base_dir']+config['relpath_uniref100'],config)
    # create_uniref_dataset(config['download_base_dir']+config['relpath_uniref90'],config)
    # create_uniref_dataset(config['download_base_dir']+config['relpath_uniref50'],config)

    global uniprotkb_uniref_idmap

    utils.info('Loading NCBI taxonomic tree...')
    taxontree_path = "{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree'])
    if os.path.exists(taxontree_path):
        taxontree = pickle.load(open(taxontree_path,'rb'))
        taxon_to_process(config, taxontree)
        utils.info('Done.\n')
    else:
        utils.error('NCBI taxonomic tree not found. Exiting...', exit = True)
        
    # parse_uniprotkb_uniref_idmapping(config)
    # if all([os.path.exists("{}/pickled/{}_uniparc_idmap.pkl".format(config['download_base_dir'],c)) for c in ('uniref100','uniref90','uniref50')]) and os.path.exists('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap'])):
    #     merge_uniparc_idmapping(config)
    # else:
    #     utils.error('Failed to process UniRef database. Exiting...', exit = True)
    
    utils.info('Loading UniProtKB-UniRef mapping...')
    uniprotkb_uniref_idmap = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'rb'))
    utils.info('Done.\n')
    
    # parse_uniprotkb_xml(config['download_base_dir']+config['relpath_uniprot_sprot'], config)
    parse_uniprotkb_xml(config['download_base_dir']+config['relpath_uniprot_trembl'], config)
    # parse_uniparc_xml(config['download_base_dir']+config['relpath_uniparc'],config)

    # merge_idmap(config)
    # create_proteomes(config['download_base_dir']+config['relpath_proteomes_xml'], config)
    # annotate_taxon_tree(config)

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']
    process_proteomes(config)
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(float(t1 - t0)))
