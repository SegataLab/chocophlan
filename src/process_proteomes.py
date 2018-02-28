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
import glob
import time
import math
import sys
import shutil
from itertools import zip_longest
from lxml import etree
from functools import partial

kingdom_to_process = ['Bacteria','Archaea']
genus_to_process = ['Candida', 'Blastocystis', 'Saccharomyces']
#dbReference = ['GeneID','Proteomes']
dbReference = ['EMBL','EnsemblBacteria','GeneID','GO','KEGG','KO','Pfam','RefSeq','Proteomes']
group_chunk = 1000000

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

def init_parse(terminating_, tree_,uniprotkb_uniref_idmap_):
    global terminating
    terminating = terminating_
    global tree
    tree = tree_
    global uniprotkb_uniref_idmap
    uniprotkb_uniref_idmap = uniprotkb_uniref_idmap_


def taxon_to_process(config, taxontree):
    d_ranks = taxontree.lookup_by_rank()
    global taxid_to_process
    taxid_to_process = []
    with open(config['relpath_taxon_to_process'],'rt') as taxa:
        for line in taxa:
            name, rank = line.strip().split('\t')
            taxid_to_process.extend([x.tax_id for x in d_ranks[rank] if name == x.name])

def is_taxon_processable(taxid, taxon_to_process):
    for x in taxon_to_process:
        if d_taxids[x].is_parent_of(d_taxids[taxid]):
            return True
    return False

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
             'UniRef100' : v[12],
             'UniRef90' : v[13],
             'UniRef50' : v[14]
            }

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
        tag_to_parse = ['accession','dbReference','sequence','organism','sequence']
        try:
            d_prot = {}
            sequence = ''
            elem = etree.fromstring(elem)
            org = [x for x in elem.iterchildren('{http://uniprot.org/uniprot}organism')][0]
            taxid = int([x.get('id') for x in org.iterchildren(tag='{http://uniprot.org/uniprot}dbReference')][0])

            if is_taxon_processable(taxid,taxid_to_process): #or len([t.text for t in taxon if t.text in genus_to_process]):
                for children in tag_to_parse:
                    tag_children = '{http://uniprot.org/uniprot}'+children
                    if children == 'name' or children == 'sequence' or children == 'accession':
                        d_prot[children] = [x for x in elem.iterchildren(tag=tag_children)][0].text
                    if children == 'organism':
                        org = [x for x in elem.iterchildren(tag_children)][0]
                        d_prot['tax_id'] = int([x.get('id') for x in org.iterchildren(tag='{http://uniprot.org/uniprot}dbReference')][0])
                    if children == 'dbReference':
                        for ref in elem.iterchildren(tag=tag_children):
                            if ref.get('type') in dbReference:
                                if ref.get('type') not in d_prot:
                                    d_prot[ref.get('type')] = []
                                d_prot[ref.get('type')].append(ref.get('id'))


                elem.clear()
                t_prot = (d_prot.get('accession'),  #0
                          d_prot.get('tax_id'),      #1
                          d_prot.get('sequence',''),   #2
                          tuple(d_prot.get('RefSeq','')),     #3
                          tuple(int(x[2:]) for x in d_prot.get('Proteomes',[]) if x is not None),  #4 UP000005640
                          tuple(d_prot.get('EMBL','')),       #5
                          tuple(d_prot.get('EnsemblBacteria','')),    #6
                          tuple(d_prot.get('GeneID','')),     #7
                          tuple(int(x[3:]) for x in d_prot.get('GO',[]) if x is not None),         #8 GO:0016847
                          tuple(int(x[1:]) for x in d_prot.get('KO',[]) if x is not None),         #9 K09972
                          tuple(d_prot.get('KEGG','')),       #10
                          tuple(int(x[2:]) for x in d_prot.get('Pfam',[]) if x is not None),        #11
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[0],
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[1],
                          uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[2]
                         )
                for ancestor in elem.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]
                return(t_prot)
        except Exception as e:
            print(d_prot)
            utils.error(str(e))
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
    global uniprotkb_uniref_idmap
    uniprotkb_uniref_idmap = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'rb'))
    taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb'))
    taxon_to_process(config, taxontree)
    global d_taxids
    d_taxids = taxontree.lookup_by_taxid()

    parse_uniprotkb_xml_elem_partial = partial(parse_uniprotkb_xml_elem, config=config)

    for file_chunk, group in enumerate(grouper((etree.tostring(v) for _, v in tree),group_chunk),1):
        try:
            with mp.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool:
                d={x[0]:x for x in pool.imap_unordered(parse_uniprotkb_xml_elem_partial, group, chunksize=chunksize) if x is not None}
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        idmap.update(dict.fromkeys(d.keys(), db*file_chunk))
        pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],db_name, file_chunk),'wb'), -1)
    pickle.dump(idmap, open("{}/pickled/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'],db_name),'wb'), -1)
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

def parse_uniparc_xml_elem(elem, config, uniprotkb_uniref_idmap):
    #if not terminating.is_set() and elem is not None:
    tag_to_parse = ['accession','dbReference', 'sequence']
    try:
        d_prot = {}
        # elem = etree.fromstring(elem)
        for children in tag_to_parse:
            tag_children = '{http://uniprot.org/uniparc}'+children
            if children == 'accession' or children == 'sequence':
                d_prot[children] = [x for x in elem.iterchildren(tag=tag_children)][0].text.replace('\n','')
            if children == 'dbReference':
                active_entries = [x for x in elem.iterchildren(tag=tag_children) if x.get('active') == 'Y']
                d_prot['uniprotkb_ids'] = [x.get('id') for x in active_entries if 'UniProtKB' in x.get('type') ]
                entry_with_protid = [x for x in active_entries for y in x.iterchildren('{http://uniprot.org/uniparc}property') if y.get('type') == 'proteome_id']
                d_prot['Proteomes'] = set()
                for entry in entry_with_protid:
                    for p in entry:
                        if p.get('type') == 'proteome_id':
                            proteome_id = int(p.get('value')[2:])
                        elif p.get('type') == 'NCBI_taxonomy_id':
                            tax_id = int(p.get('value'))
                    if is_taxon_processable(tax_id,taxid_to_process):
                        d_prot['Proteomes'].add((proteome_id, tax_id))
        elem.clear()
        t_prot = (d_prot.get('accession'),  
                  tuple(d_prot.get('Proteomes',[])),
                  d_prot.get('sequence'),
                  uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[0],
                  uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[1],
                  uniprotkb_uniref_idmap.get(d_prot.get('accession'),['','',''])[2],
                  d_prot.get('uniprotkb_ids',[])
                 )

        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]

        ## If entry does not belong to a proteome, just throw it, wont need.
        if not len(t_prot[1]):
            return None
            
        return(t_prot)
    except Exception as e:
        print(d_prot)
        utils.error(str(e))
        raise
    #else:
    #    terminating.set()

def parse_uniparc_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = config['nproc']*6

    if config['verbose']:
        utils.info('Loading UniParc-UniRef cross-reference map...\n')
    t0 = time.time()
    uniprotkb_uniref_idmap = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'rb'))
    t1 = time.time()
    if config['verbose']:
        utils.info('Cross-reference map loaded in {} seconds\n'.format(int(t1-t0)))

    tree = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = '{http://uniprot.org/uniparc}entry', huge_tree = True)
    idmap = {}
    #parse_uniparc_xml_elem_partial = partial(parse_uniparc_xml_elem, config=config)
    d = []
    for file_chunk, group in enumerate(tree):  #Offset for identify UniParc entries in the id:chunk map file
        if file_chunk == 0 or file_chunk % 1000000:
            res = parse_uniparc_xml_elem(group[1], config, uniprotkb_uniref_idmap)
            if res is not None:
                d.append(res)
        else:
            d = {x[0]:x for x in d}
        #try:
        #    with mp.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool:
        #        print('In pool')
        #        d={x[0]:x for x in pool.imap_unordered(parse_uniparc_xml_elem_partial, group, chunksize=chunksize) if x is not None}
        #except Exception as e:
        #    utils.error(str(e))
        #    utils.error('Processing failed',  exit=True)
        #    raise
            print('{} entry processed.'.format(file_chunk, flush = True))
            idmap.update(dict.fromkeys(d.keys(), int(file_chunk/1000000)+1000))
            pickle.dump(d, open("{}/pickled/uniparc_{}.pkl".format(config['download_base_dir'], int(file_chunk/1000000)+1000),'wb'), -1)
            d = []
    pickle.dump(idmap, open("{}/pickled/uniparc_idmap.pkl".format(config['download_base_dir']),'wb'), -1)
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

def get_members(line):
    line = line.split('\t')
    u100 = line[7].split('_')[1] if len(line[7]) else '' 
    u90 =  line[8].split('_')[1] if len(line[8]) else '' 
    u50 =  line[9].split('_')[1] if len(line[9]) else '' 
    return {line[0]: (u100, u90, u50)}

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
            idmapping = {k:v for x in pool.imap_unordered(get_members, f) for k,v in x.items()}
    pickle.dump(idmapping, open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb'), -1)
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
                    for proteome, taxid in (("UP{}{}".format("0"*(9-len(str(upi))),upi),taxid) for upi, taxid in entry[1]):
                        if proteome not in d_proteome:
                            d_proteome[proteome] = {'members' : [], 'isReference' : False, 'tax_id' : taxid, 'upi' : True}
                        d_proteome[proteome]['members'].append(entry[0])
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        return d_proteome
    else:
        terminating.set()

def create_proteomes(config):
    if config['verbose']:
        utils.info('Starting proteomes processing...\n')
    d_proteomes = {}
    
    terminating = mp.Event()
    chunksize = config['nproc']
    
    
    r = re.compile('.*(trembl|sprot|uniparc).*')
    chunks = []

    try:
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool:
            chunks = [x for x in pool.imap_unordered(process, 
                        (x for x in filter(r.match,glob.iglob("{}/pickled/*".format(config['download_base_dir']))))
                          , chunksize = chunksize)]
    except Exception as e:
        utils.error(str(e))
        utils.error('Processing failed',  exit=True)
        raise

    for chunk in chunks:
        for k,v in chunk.items():
            if k not in d_proteomes:
                d_proteomes[k] = {'members' : set(), 'isReference' : False, 'tax_id' : v['tax_id'], 'upi' : v['upi']}
            d_proteomes[k]['members'].add(v['members'])

    if config['verbose']:
        utils.info('Done processing\n')

    if not os.path.exists('{}/{}'.format(config['download_base_dir'],config['relpath_reference_proteomes'])):
        utils.error('Required files does not exist! Exiting...', exit=True)
    
    if config['verbose']:
        utils.info('Starting reference proteomes processing...\n',)

    for k in kingdom_to_process:
        basepath = '{}/{}/{}/*.idmapping.gz'.format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)
        for protid, taxid in [os.path.split(file)[1].split('.')[0].split('_') for file in ls]:
            if protid in d_proteomes:
                d_proteomes[protid]['isReference'] = True
            else:
                utils.info(protid+"\n")
                
    pickle.dump(d_proteomes, open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'wb'), -1)
    if config['verbose']:
        utils.info('Done\n')

def annotate_taxon_tree(config):
    utils.info('Starting annotation of taxonomic tree with proteome ids\n')
    proteomes = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'rb'))
    taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb'))
    taxids = taxontree.lookup_by_taxid()
    for protid, v in proteomes.items():
        taxid = v['tax_id']
        try:
            clade = taxids[int(taxid)]
            if not hasattr(clade,'proteomes'):
                clade.proteomes = set()
            clade.proteomes.add(protid)
        except:
            print(taxid)
    taxontree.are_leaves_trimmed = False
    pickle.dump(taxontree, open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_taxontree']),'wb'), -1)
    utils.info('Done.')

def uniref_tuple_to_dict(v):
    return {'id': v[0],
            'common_taxid': v[1],
            'members': v[2],
            'UniRef100': 'UniRef100_{}'.format(v[3]) if len(v[3])>0 else '',
            'UniRef90': 'UniRef90_{}'.format(v[4]) if len(v[4])>0 else '',
            'UniRef50': 'UniRef50_{}'.format(v[5]) if len(v[5])>0 else ''
           }

def parse_uniref_xml_elem(elem, config):
    if not terminating.is_set() and elem is not None:
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
                            d_uniref['common_taxid'] = int(pro.get('value'))
                if tag == 'representativeMember' or tag == 'member':
                    member = [x for x in elem.iterchildren(tag_children)]
                    ref = [x for y in member for x in y.iterchildren('{http://uniprot.org/uniref}dbReference')]
                    for r in ref:
                        if r.get('type') == 'UniProtKB ID':
                            accession = [y.get('value') for y in r if y.get('type') == 'UniProtKB accession'][0]
                            uniparc_cluster = [y.get('value') for y in r if y.get('type') == 'UniParc ID'][0]
                        elif r.get('type') == 'UniParc ID':
                            accession = r.get('id')
                            uniparc_cluster = None
                        isRepr = True if tag == 'representativeMember' else False
                        d_uniref['members'].append((accession,isRepr))
                        if uniparc_cluster is not None:
                            d_uniref['members'].append((uniparc_cluster,isRepr))
                        for pro in r:
                            if pro.get('type') == 'UniRef100 ID':
                                d_uniref['UniRef100'] = pro.get('value')[10:]
                            if pro.get('type') == 'UniRef90 ID':
                                d_uniref['UniRef90'] = pro.get('value')[9:]
                            if pro.get('type') == 'UniRef50 ID':
                                d_uniref['UniRef50'] = pro.get('value')[9:]
            t_uniref = (d_uniref['id'],
                        d_uniref.get('common_taxid',None),
                        tuple(d_uniref['members']),
                        d_uniref.get('UniRef100',''),
                        d_uniref.get('UniRef90', ''),
                        d_uniref.get('UniRef50','')
                        )
            return t_uniref
            
        except Exception as e:
            utils.error(str(e))
            with open("{}/uniref/FAILED".format(config['temp_folder']),'w+') as out:
                out.write(elem.get('id')+"\n")
            raise
    else:
        terminating.set()

def create_uniref_dataset(xml, config):
    uniref_xml = etree.iterparse(gzip.GzipFile(xml), events = ('end',), tag = '{http://uniprot.org/uniref}entry', huge_tree = True)
    
    terminating = mp.Event()
    chunksize = config['nproc']
    cluster = os.path.basename(xml).split('.')[0]

    idmap={}
    upids = []
    if config['verbose']:
        utils.info("Starting processing UniRef {} database\n".format(cluster))

    parse_uniref_xml_elem_partial = partial(parse_uniref_xml_elem, config=config)

    for file_chunk, group in enumerate(grouper((etree.tostring(v) for _, v in uniref_xml),group_chunk),1):
        try:
            with mp.Pool(initializer=init_parse, initargs=(terminating, uniref_xml, None,), processes=chunksize) as pool:
                d={x[0]:x for x in pool.imap_unordered(parse_uniref_xml_elem_partial, group, chunksize=chunksize) if x is not None}
        except Exception as e:
            utils.error(str(e))
            raise
        upids.extend([(m[0],c[0]) for c in d.values() for m in c[2] if 'UPI' in m[0]])
        pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],cluster, file_chunk),'wb'), -1)
        idmap.update(dict.fromkeys(d.keys(), file_chunk))

    pickle.dump(upids, open("{}/pickled/{}_uniparc_idmap.pkl".format(config['download_base_dir'],cluster),'wb'), -1)
    pickle.dump(idmap, open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster),'wb'), -1)
    if config['verbose']:
        utils.info('UniRef {} database processed successfully.\n'.format(cluster))

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
    pickle.dump(idmapping, open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb'), -1)

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
        
def process_proteomes(config):
    os.makedirs('{}/pickled'.format(config['download_base_dir']), exist_ok=True)
    step1     = [mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref100'],config)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref90'],config)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref50'],config)),
                 mp.Process(target=parse_uniprotkb_uniref_idmapping, args=(config,))]

    step2     = [mp.Process(target=merge_uniparc_idmapping, args=(config,))]

    step3     = [mp.Process(target=parse_uniparc_xml, args=(config['download_base_dir']+config['relpath_uniparc'],config)),
                 mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_sprot'], config)),
                 mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_trembl'], config))
                ]
                

    for p in step1:
       p.start()

    for p in step1:
       p.join()

    for p in step2:
       p.start()

    for p in step2:
       p.join()

    for p in step3:
        p.start()

    for p in step3:
        p.join()

    merge_idmap(config)
    create_proteomes(config)
    annotate_taxon_tree(config)

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

    sys.exit(0)
