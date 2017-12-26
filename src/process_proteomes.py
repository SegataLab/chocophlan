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
#dbReference = ['GeneID','Proteomes']
dbReference = ['EMBL','EnsemblBacteria','GeneID','GO','KEGG','KO','Pfam','RefSeq','Proteomes']

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

def init_parse(terminating_, tree_,idmapping_):
    global terminating
    terminating = terminating_
    global tree
    tree = tree_
    global idmapping
    idmapping = idmapping_

def uniprot_tuple_to_dict(v):
    return {'accession' : v[0],
            'tax_id' : v[1],
            'sequence' : v[2],
             'RefSeq' : v[3],
             'Proteomes' : ["UP{}{}".format("0"*(9-len(x)),x) for x in v[4]],
             'EMBL' : v[5],
             'EnsemblBacteria' : v[6],
             'GeneID' : v[7],
             'GO' : ["GO:{}{}".format('0'*(7-len(x)),x) for x in v[8]],
             'KO' : ["K{}{}".format("0"*(5-len(x)),x) for x in v[9]],
             'KEGG' : v[10],
             'Pfam' : ["PF{}{}".format('0'*(5-len(x)),x) for x in v[11]],
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
            kingdom = [x for x in org.iterchildren('{http://uniprot.org/uniprot}lineage')][0].getchildren()[0].text
            
            if kingdom in kingdom_to_process:
                for children in tag_to_parse:
                    tag_children = '{http://uniprot.org/uniprot}'+children
                    if children == 'name' or children == 'accession' or children == 'sequence':
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
                          d_prot.get('sequence'),   #2
                          tuple(d_prot.get('RefSeq')),     #3
                          (int(x[2:]) for x in d_prot.get('Proteomes',[]) if x is not None),  #4 UP000005640
                          tuple(d_prot.get('EMBL')),       #5
                          tuple(d_prot.get('EnsemblBacteria')),    #6
                          tuple(d_prot.get('GeneID')),     #7
                          (int(x[3:]) for x in d_prot.get('GO',[]) if x is not None),         #8 GO:0016847
                          (int(x[1:]) for x in d_prot.get('KO',[]) if x is not None),         #9 K09972
                          tuple(d_prot.get('KEGG')),       #10
                          (int(x[2:]) for x in d_prot.get('Pfam',[]) if x is not None),        #11
                          idmapping[d_prot.get('accession')][0],
                          idmapping[d_prot.get('accession')][1],
                          idmapping[d_prot.get('accession')][2]
                         )
                for ancestor in elem.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]
                return(t_prot)
        except Exception as e:
            utils.error(str(e))
            raise
    else:
        terminating.set()
    
def parse_uniprotkb_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = int(config['nproc'])

    db= -1 if 'uniprot_sprot' in xml_input else +1
    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))
    
    tree = etree.iterparse(gzip.GzipFile(xml_input))
    if(os.path.exists("{}/pickled/uniprotkb_idmap.pkl".format(config['download_base_dir']))):
        idmap=pickle.load(open("{}/pickled/uniprotkb_idmap.pkl".format(config['download_base_dir']),'rb'))
    else:
        idmap={}
    
    idmapping = pickle.load(open(config['download_base_dir']+'/pickled/idmapping.pkl','rb'))

    parse_uniprotkb_xml_elem_partial = partial(parse_uniprotkb_xml_elem, config=config)
    with mp.Pool(initializer=init_parse, initargs=(terminating, tree,idmapping,), processes=chunksize) as pool:
        try:
            file_chunk = 1
            for group in grouper(yield_filtered_xml_string(tree),1000000):
                d={x[0]:x for x in pool.imap(parse_uniprotkb_xml_elem_partial, group, chunksize=chunksize) if x is not None}
                idmap.update(dict.fromkeys(d.keys(), db*file_chunk))
                pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],db, file_chunk),'wb'), -1)
                file_chunk+=1
            pickle.dump(idmap, open("{}/pickled/uniprotkb_idmap.pkl".format(config['download_base_dir']),'wb'), -1)
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        if config['verbose']:
            utils.info('Done processing {} file!\n'.format(xml_input))
    if config['verbose']:
        utils.info('Done processing UniProtKB\n')

def get_members(line):
    line = line.split('\t')
    u100 = line[7].split('_')[1] if len(line[7]) else None
    u90 =  line[8].split('_')[1] if len(line[8]) else None
    u50 =  line[9].split('_')[1] if len(line[9]) else None
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
def parse_idmapping(config):
    if config['verbose']:
        utils.info('Starting processing idmapping.gz\n')
    input_file = config['download_base_dir'] + config['relpath_idmapping']
    terminating = mp.Event()
    chunksize = config['nproc']
    idmapping = {}

    with gzip.open(input_file, 'rt', encoding='utf-8') as f:
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool:
            idmapping = {k:v for x in pool.imap(get_members, f) for k,v in x.items()}
    pickle.dump(idmapping, open(config['download_base_dir']+'/pickled/idmapping.pkl','wb'))
    if config['verbose']:
        utils.info('Done processing idmapping.gz\n')

def get_uniprotkb_entry(config, accession):
    idmap = pickle.load(open("{}/pickled/uniprotkb_idmap.pkl".format(config['download_base_dir']),'rb'))
    fn = idmap[accession]
    if fn < 0:
        db = 'sprot'
    else:
        db = 'trembl'
    chunk = pickle.load(open('uniprot_{}_{}.pkl'.format(db,abs(fn)),'rb'))
    return uniprot_tuple_to_dict(chunk[accession])


def create_proteomes(config, verbose = False):
    if verbose:
        utils.info('Starting proteomes processing...\n')
    d_proteomes = {}
    
    terminating = mp.Event()
    chunksize = config['nproc']
    
    def process(f):
        if not terminating.is_set():
            x = pickle.load(open(f,'rb'))
            d_proteome = {}
            for k,v in x.items():
                if v[4] is not None:
                    proteomes = v[4]
                    accession = k
                    taxid = v[1]
                    for proteome in proteomes:
                        if proteome not in d_proteome:
                            d_proteome[proteome] = {"isReference" : False, "members" : [], 'tax_id' : taxid}
                        d_proteome[proteome]['members'].append(accession)
            return d_proteome
        else:
            terminating.set()
    
    ls = glob.glob("{}/pickled/uniprot_*".format(config['download_base_dir']))
    with mp.Pool(initializer=init_parse, initargs=(terminating,None,None), processes=chunksize) as pool:
        try:
            chunks = [x for x in pool.imap(process, ls, chunksize = chunksize)]
            for chunk in chunks:
                for k,v in chunk.items():
                    if k not in d_proteomes:
                        d_proteomes[k] = {'members' : [], 'isReference' : False, 'tax_id' : v['tax_id']}
                    d_proteomes[k]['members'].extend(v['members'])

        except Exception as e:
            utils.error(str(e))
            utils.error('Processing failed',  exit=True)
            raise
        if config['verbose']:
            utils.info('Done processing\n')

    if not os.path.exists('{}/{}'.format(config['download_base_dir'],config['relpath_reference_proteomes'])):
        utils.error('Required files does not exist! Exiting...', exit=True)
    
    if verbose:
        utils.info('Starting reference proteomes processing...\n',)

    for k in kingdom_to_process:
        basepath = '{}/{}/{}/*.idmapping.gz'.format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)
        for protid, taxid in [os.path.split(file)[1].split('.')[0].split('_') for file in ls]:
            if protid in d_proteomes:
                d_proteomes[protid]['isReference'] = True
            else:
                utils.info(protid+"\n")
    pickle.dump(d_proteomes, open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'wb'))
    if verbose:
        utils.info('Done\n')

def annotate_taxon_tree(config):
    proteomes = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'rb'))
    taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb'))
    taxids = extract.Nodes.lookup_by_taxid(taxontree) 
    i=0
    l=len(proteomes)
    for protid, v in proteomes.items():
        i+=1
        taxid = v['tax_id']
        try:
            clade = taxids[int(taxid)]
            clade.proteome = protid
        except:
            continue
        print("{}/{}".format(i,l))
    pickle.dump(taxontree, open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_taxontree']),'wb'))

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
                    properties = [x.iterchildren('{http://uniprot.org/uniref}property') for y in member for x in y.iterchildren('{http://uniprot.org/uniref}dbReference')]
                    for p in properties:
                        for pro in p:
                            if pro.get('type') == 'UniProtKB accession' or pro.get('type') == 'UniParc ID':
                                accession = pro.get('value')
                                isRepr = True if tag == 'representativeMember' else False
                                d_uniref['members'].append((accession,isRepr))
                            if pro.get('type') == 'UniRef100 ID':
                                d_uniref['UniRef100'] = pro.get('value')[10:]
                            if pro.get('type') == 'UniRef90 ID':
                                d_uniref['UniRef90'] = pro.get('value')[9:]
                            if pro.get('type') == 'UniRef50 ID':
                                d_uniref['UniRef50'] = pro.get('value')[9:]
            if d_uniref.get('common_taxid') is None:
                print(d_uniref['id'])
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

def uniref_tuple_to_dict(v):
    return {'id': v[0],
            'common_taxid': v[1],
            'members': v[2],
            'UniRef100': 'UniRef100_{}'.format(v[3]) if len(v[3])>0 else '',
            'UniRef90': 'UniRef90_{}'.format(v[4]) if len(v[4])>0 else '',
            'UniRef50': 'UniRef50_{}'.format(v[5]) if len(v[5])>0 else ''
           }

def create_uniref_dataset(xml, config):
    uniref_xml = etree.iterparse(gzip.GzipFile(xml))
    
    terminating = mp.Event()
    chunksize = config['nproc']
    cluster = os.path.basename(xml).split('.')[0]

    if(os.path.exists("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'], cluster))):
        idmap=pickle.load(open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster), 'rb'))
    else:
        idmap={}

    file_chunk = 1
    with mp.Pool(initializer=init_parse, initargs=(terminating, uniref_xml, None,), processes=chunksize) as pool:
        try:
            if config['verbose']:
                utils.info("Starting processing UniRef {} database\n".format(cluster))
            parse_uniref_xml_elem_partial = partial(parse_uniref_xml_elem, config=config)
            for group in grouper(yield_filtered_xml_string(uniref_xml), 1000000):
                d={x[0]:x for x in pool.imap(parse_uniref_xml_elem_partial, group, chunksize=chunksize) if x is not None}
                pickle.dump(d, open("{}/pickled/{}_{}.pkl".format(config['download_base_dir'],cluster, file_chunk),'wb'), -1)
                idmap.update(dict.fromkeys(d.keys(), file_chunk))

                file_chunk+=1

        except Exception as e:
            utils.error(str(e))
            raise
    pickle.dump(idmap, open("{}/pickled/{}_idmap.pkl".format(config['download_base_dir'],cluster),'wb'), -1)
    if config['verbose']:
        utils.info('UniRef {} database processed successfully.\n'.format(cluster))

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

        
def process_proteomes(config):
    os.makedirs('{}/pickled'.format(config['download_base_dir']), exist_ok=True)
    step1     = [#mp.Process(target=parse_idmapping, args=(config,)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref100'],config)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref90'],config)),
                 mp.Process(target=create_uniref_dataset, args=(config['download_base_dir']+config['relpath_uniref50'],config)),
                ]
    setp2     = [mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_sprot'], config)),
                 mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_trembl'], config)),
                ]

    for p in step1:
        p.start()

    for p in step1:
        p.join()

    for p in step2:
        p.start()

    for p in step2:
        p.join()

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
