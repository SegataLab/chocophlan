#!/usr/bin/env python3

__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Fabio Cumbo (fabio.cumbo@unitn.it),'
              'Francesco Asnicar (f.asnicar@unitn.it),'
              'Nicola Segata (nicola.segata@unitn.it)')
from _version import __CHOCOPhlAn_version__
__date__ = '25 Mar 2019'

if __name__ == '__main__':
    import utils
    import extract
else:
    import src.utils as utils
    import src.extract as extract
import os
import gzip
import re
import multiprocessing as mp
import multiprocessing.dummy as mpdummy
import glob
import time
import sys
import concurrent.futures
import pickle
import itertools
from lxml import etree
from functools import partial
import traceback

ns = {'upkb' : 'http://uniprot.org/uniprot', 'nr' : 'http://uniprot.org/uniref', 'up' : 'http://uniprot.org/uniparc'}
GROUP_CHUNK = 1000000

taxid_to_process = set()
uniprotkb_uniref_idmap = {}


def grouper(iterable, n, fillvalue=None):
    #"Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

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
             'UniRef50' : v[16],
             'isFragment' : v[17]
            }

def parse_uniref_xml_elem(elem):
    if not terminating.is_set() and elem is not None:
        elem = etree.fromstring(elem)
        members = []
        t_uniref = None
        try:
            nr_id = ''.join(elem.xpath("/nr:entry/@id", namespaces=ns, smart_strings=False))
            common_taxid = ''.join(elem.xpath("/nr:entry/nr:property[@type='common taxon ID']/@value", namespaces=ns, smart_strings=False))
            common_taxid = int(common_taxid) if common_taxid else common_taxid
            #TODO Fix extraction of upper UniRef level, now returns a collapsed list
            UniRef100 = set(elem.xpath("(/nr:entry/nr:member|nr:representativeMember)/nr:dbReference/nr:property[@type='UniRef100 ID']/@value", namespaces = ns, smart_strings=False))
            UniRef90 = set(elem.xpath("(/nr:entry/nr:member|nr:representativeMember)/nr:dbReference/nr:property[@type='UniRef90 ID']/@value", namespaces = ns, smart_strings=False))
            UniRef50 = set(elem.xpath("(/nr:entry/nr:member|nr:representativeMember)/nr:dbReference/nr:property[@type='UniRef50 ID']/@value", namespaces = ns, smart_strings=False))
            xml_all_members = elem.xpath("/nr:entry/nr:member|nr:representativeMember",namespaces=ns, smart_strings=False)

            for member in xml_all_members:
                accession = member.xpath("./nr:dbReference/nr:property[@type='UniProtKB accession']/@value|./nr:dbReference[@type='UniParc ID']/@id", namespaces = ns, smart_strings=False)[0]
                uniparc_cluster = ''.join(member.xpath("./nr:dbReference/nr:property[@type='UniParc ID']/@value", namespaces = ns, smart_strings=False))
                tax_id = ''.join(member.xpath("./nr:dbReference/nr:property[@type='NCBI taxonomy']/@value", namespaces = ns, smart_strings=False))
                tax_id = int(tax_id) if len(tax_id) else common_taxid
                isRepr = True if member.xpath('name()', namespaces=ns, smart_strings=False) == 'representativeMember' else False
                members.append((accession,int(tax_id), isRepr))
                if uniparc_cluster:
                    members.append((uniparc_cluster,int(tax_id), isRepr))

            t_uniref = (nr_id,
                        common_taxid,
                        tuple(members),
                        UniRef100 if len(UniRef100) > 1 else ''.join(UniRef100)[10:],
                        UniRef90 if len(UniRef90) > 1 else ''.join(UniRef90)[9:],
                        UniRef50 if len(UniRef50) > 1 else ''.join(UniRef50)[9:],
                        ''.join(elem.xpath("/nr:entry/nr:representativeMember/nr:sequence/text()", namespaces=ns, smart_strings=False)).replace('\n','')
                        )
        except Exception as e:
            utils.error('Failed to elaborate item: '+ elem.get('id'))
            terminating.set()
        finally:
            elem.clear()
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
        
        return t_uniref

def parse_uniref_xml(xml_input, config):  
    uniref_xml = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = '{http://uniprot.org/uniref}entry', huge_tree = True)

    terminating = mp.Event()
    chunksize = config['nproc']
    cluster = os.path.basename(xml_input).split('.')[0]
    if config['verbose']:
        utils.info("Starting processing UniRef {} database\n".format(cluster))

    if not os.path.exists("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], cluster)):
        os.makedirs("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], cluster))

    with mp.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool,\
         open("{}/{}/{}_uniparc_idmap.pkl".format(config['download_base_dir'],config['pickled_dir'], cluster),'wb') as pickle_uniref_uniparc_idmap,\
         open("{}/{}/{}_idmap.pkl".format(config['download_base_dir'],config['pickled_dir'], cluster),'wb') as pickle_uniref_idmap,\
         open("{}/pickled/{}_taxid_idmap.pkl".format(config['download_base_dir'],cluster),'wb') as pickle_taxid_map:
        try:
            file_chunk = 1
            pickle_chunk = open("{}/{}/{}/{}_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], cluster, cluster, file_chunk),'wb')
            for index, elem in enumerate(pool.imap_unordered(parse_uniref_xml_elem, yield_filtered_xml_string(uniref_xml), chunksize=chunksize)):
                if not index % GROUP_CHUNK:
                    pickle_chunk.close()
                    file_chunk = 1 + (index//GROUP_CHUNK )
                    pickle_chunk = open("{}/{}/{}/{}_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], cluster, cluster, file_chunk),'wb')
                utils.optimized_dump(pickle_chunk, elem)
                utils.optimized_dump(pickle_uniref_uniparc_idmap, list(set((m[0],elem[0]) for m in elem[2] if 'UPI' in m[0])))
                utils.optimized_dump(pickle_uniref_idmap, (elem[0], file_chunk, ))
                utils.optimized_dump(pickle_taxid_map, { elem[0] : (set(t[:3] for t in elem[2]), elem[3:6]) })
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing of {} failed.'.format(xml_input))
            raise
        finally:
            pickle_chunk.close()
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
# 12 EC
# 13 eggNOG
# 14 NR100
# 15 NR90
# 16 NR50
# isFragment
def parse_uniprotkb_xml_elem(elem):
    if not terminating.is_set() and elem is not None:
        global uniprotkb_uniref_idmap, taxid_to_process
        elem = etree.fromstring(elem)
        t_prot = None
        try:
            tax_id = ''.join(elem.xpath("/upkb:entry/upkb:organism/upkb:dbReference[@type='NCBI Taxonomy']/@id", namespaces=ns, smart_strings=False))
            tax_id = int(tax_id) if tax_id else tax_id
            if is_taxon_processable(tax_id):
                accession = ''.join(elem.xpath("/upkb:entry/upkb:accession[1]/text()", namespaces=ns, smart_strings=False))
                gene = elem.xpath("/upkb:entry/upkb:gene/upkb:name/@type|/upkb:entry/upkb:gene/upkb:name/text()", namespaces=ns, smart_strings=False)
                gene = tuple(zip(gene[0::2], gene[1::2]))
                nr_ids = uniprotkb_uniref_idmap.get(accession,['','',''])
                t_prot = (accession,  #0
                          tax_id,      #1
                          ''.join(elem.xpath("/upkb:entry/upkb:sequence/text()", namespaces=ns, smart_strings=False)),   #2
                          gene,     #3
                          tuple(int(x[2:]) for x in elem.xpath("/upkb:entry/upkb:dbReference[@type='Proteomes']/@id", namespaces=ns, smart_strings=False)),  #4 UP000005640
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='EMBL']/@id", namespaces=ns, smart_strings=False)),       #5
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='EnsemblBacteria']/@id", namespaces=ns, smart_strings=False)),    #6
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='GeneID']/@id", namespaces=ns, smart_strings=False)),     #7
                          tuple(int(x[3:]) for x in elem.xpath("/upkb:entry/upkb:dbReference[@type='GO']/@id", namespaces=ns, smart_strings=False)),         #8 GO:0016847
                          tuple(int(x[1:]) for x in elem.xpath("/upkb:entry/upkb:dbReference[@type='KO']/@id", namespaces=ns, smart_strings=False)),         #9 K09972
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='KEGG']/@id", namespaces=ns, smart_strings=False)),       #10
                          tuple(int(x[2:]) for x in elem.xpath("/upkb:entry/upkb:dbReference[@type='Pfam']/@id", namespaces=ns, smart_strings=False)),        #11
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='EC']/@id", namespaces=ns, smart_strings=False)),
                          tuple(elem.xpath("/upkb:entry/upkb:dbReference[@type='eggNOG']/@id", namespaces=ns, smart_strings=False)),
                          nr_ids[0],
                          nr_ids[1],
                          nr_ids[2],
                          True if elem.xpath("/upkb:entry/upkb:sequence/@fragment", namespaces=ns, smart_strings=False) else False
                        )
        except Exception as e:
            utils.error('Failed to elaborate item: ' + ''.join(elem.xpath("/upkb:entry/upkb:accession[1]/text()", namespaces=ns, smart_strings=False)))
            terminating.set()
        finally:
            elem.clear()
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
        
        return t_prot
    
def parse_uniprotkb_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = int(config['nproc'])

    db = -1 if 'uniprot_sprot' in xml_input else +1
    db_name = 'sprot' if db == -1 else 'trembl'

    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))
    
    tree = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = '{http://uniprot.org/uniprot}entry', huge_tree = True)

    with mpdummy.Pool(initializer=initt, initargs=(terminating, ), processes=chunksize) as pool, \
         open("{}/{}/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'], config['pickled_dir'], db_name),'ab') as pickle_idmap:
        try:
            file_chunk = 1
            pickle_chunk = open("{}/{}/uniprotkb/{}_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], db_name, file_chunk),'ab')
            for index, elem in enumerate(pool.imap_unordered(parse_uniprotkb_xml_elem, yield_filtered_xml_string(tree)),1):
                if not index % GROUP_CHUNK:
                    pickle_chunk.close()
                    file_chunk = 1 + (index//GROUP_CHUNK )
                    pickle_chunk = open("{}/{}/uniprotkb/{}_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], db_name, file_chunk),'ab')
                if elem:
                    utils.optimized_dump(pickle_chunk, elem)
                    utils.optimized_dump(pickle_idmap, (elem[0], file_chunk*db,))
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing of {} failed.'.format(xml_input),  exit=True)
            raise
        finally:
            pickle_chunk.close()
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

def parse_uniparc_xml_elem(elem):
    global uniprotkb_uniref_idmap, taxid_to_process
    if elem is not None:
        elem = etree.fromstring(elem)
        t_prot = None
        try:
            accession = ''.join(elem.xpath("/up:entry/up:accession/text()", namespaces=ns, smart_strings=False))
            uniprotkb_ids = list(set(elem.xpath("/up:entry/up:dbReference[@active='Y' and starts-with(@type, 'UniProtKB')]/@id", namespaces=ns, smart_strings=False)))
            entry_with_protid = [x for x in elem.xpath("/up:entry/up:dbReference[@active='Y']", namespaces=ns, smart_strings=False) for y in x.iterchildren('{http://uniprot.org/uniparc}property') if y.get('type') == 'proteome_id']
            proteomes = set()
            for entry in entry_with_protid:
                proteome_id = int(''.join(entry.xpath("./up:property[@type = 'proteome_id']/@value", namespaces=ns, smart_strings=False))[2:])
                tax_id = ''.join(entry.xpath("./up:property[@type = 'NCBI_taxonomy_id']/@value", namespaces=ns, smart_strings=False))
                gene_name = ''.join(entry.xpath("./up:property[@type = 'gene_name']/@value", namespaces=ns, smart_strings=False))
                if tax_id and is_taxon_processable(int(tax_id)):
                    proteomes.add((proteome_id, int(tax_id), gene_name))

            nr_ids = uniprotkb_uniref_idmap.get(accession,['','',''])
            t_prot = ( accession, 
                       tuple(proteomes),
                       ''.join(elem.xpath("/up:entry/up:sequence/text()", namespaces=ns, smart_strings=False)).replace('\n',''),
                       nr_ids[0],
                       nr_ids[1],
                       nr_ids[2],
                       uniprotkb_ids
                     )
        except Exception as e:
            utils.error('Failed to elaborate item: '+ ''.join(elem.xpath("/up:entry/up:accession/text()", namespaces=ns, smart_strings=False)))
            raise
        finally:
            elem.clear()
            for ancestor in elem.xpath('ancestor-or-self::*'):
                while ancestor.getprevious() is not None:
                    del ancestor.getparent()[0]
        ## If entry does not belong to a proteome, just throw it, wont need.
        if not t_prot[1]:
            return None
            
        return t_prot

def parse_uniparc_xml(xml_input, config):
    terminating = mp.Event()
    chunksize = int(config['nproc'])

    if not os.path.exists("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], 'uniparc')):
        os.makedirs("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], 'uniparc'))
    if config['verbose']:
        utils.info('Starting processing {} file...\n'.format(xml_input))

    tree = etree.iterparse(gzip.open(xml_input), events = ('end',), tag = '{http://uniprot.org/uniparc}entry', huge_tree = True)
    
    with mpdummy.Pool(initializer=initt, initargs=(terminating,), processes=chunksize) as pool,\
         open("{}/{}/uniparc_idmap.pkl".format(config['download_base_dir'], config['pickled_dir']),'ab') as pickle_idmap:
        try:
            file_chunk = 1000
            pickle_chunk = open("{}/{}/uniparc/uniparc_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], file_chunk),'wb')
            for index, elem in enumerate(pool.imap_unordered(parse_uniparc_xml_elem, yield_filtered_xml_string(tree)),1):
                if not index % GROUP_CHUNK:
                    pickle_chunk.close()
                    file_chunk = 1000 + (index//GROUP_CHUNK)
                    pickle_chunk = open("{}/{}/uniparc/uniparc_{}.pkl".format(config['download_base_dir'], config['pickled_dir'], file_chunk),'wb')
                if elem :
                    utils.optimized_dump(pickle_chunk, elem)
                    utils.optimized_dump(pickle_idmap, (elem[0], file_chunk,))
        except Exception as e:
            utils.error(str(e))
            utils.error('Processing of {} failed.'.format(xml_input),  exit=True)
            raise
        finally:
            pickle_chunk.close()
    if config['verbose']:
        utils.info('Done processing {} file!\n'.format(xml_input))

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
def extract_uniref_map(line):
    line = line.split('\t')
    u100 = line[7].split('_')[1] if line[7] else '' 
    u90 =  line[8].split('_')[1] if line[8] else '' 
    u50 =  line[9].split('_')[1] if line[9] else '' 
    return (line[0], (u100, u90, u50))

def parse_uniprotkb_uniref_idmapping(config):
    if config['verbose']:
        utils.info('Starting processing idmapping.gz\n')
    input_file = config['download_base_dir'] + config['relpath_idmapping']
    chunksize = int(config['nproc'])

    with gzip.open(input_file, 'rt', encoding='utf-8') as fin,\
         open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb') as pickle_idmapping:
            for x in fin:
                utils.optimized_dump(pickle_idmapping, extract_uniref_map(x))
    if config['verbose']:
        utils.info('Done processing idmapping.gz\n')

def extract_proteomes_from_pickle(pickle_chunk):
    end = False
    entries = []
    with open(pickle_chunk,'rb') as pickle_file:
        while not end:
            try:
                entry = pickle.load(pickle_file)
                if entry is not None:
                    entries.append(entry)
            except EOFError:
                end=True

    d_proteome = {}
    for entry in entries:
        for upi, taxid, gene_name in entry[1]:
            proteome = "UP{}".format(str(upi).zfill(9))
            if proteome not in d_proteome:
                d_proteome[proteome] = {'members' : [], 'isReference' : False, 'tax_id' : taxid, 'upi' : True}
            d_proteome[proteome]['members'].append(entry[0])

    return d_proteome

def parse_proteomes_xml_elem(elem):
    if elem is not None:
        elem = etree.fromstring(elem)
        t_prot = None
        try:
            taxid = int(''.join(elem.xpath('@taxonomy', smart_strings=False)))
            upi = True if len(elem.xpath('/proteome/redundantTo', smart_strings=False)) else False
            if is_taxon_processable(taxid):
                accession = ''.join(elem.xpath('@upid', smart_strings=False))
                members = set(elem.xpath('/proteome/component/protein/@accession', smart_strings=False))
                isReference = elem.xpath('/proteome/isReferenceProteome/text()',smart_strings=False) == ['true']
                ncbi_ids = elem.xpath('/proteome/dbReference/@type|/proteome/dbReference/@id',smart_strings=False)
                ncbi_ids = tuple(zip(ncbi_ids[0::2],ncbi_ids[1::2]))

                t_prot = ( accession,
                           taxid,
                           upi,
                           isReference,
                           ncbi_ids,
                           members
                        )
        except Exception as e:
            utils.error('Failed to elaborate item: '+str(e))
            utils.error(elem.attrib['upid'])
        finally:
                elem.clear()
                for ancestor in elem.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]
        
        return t_prot

def create_proteomes_pkl(xml_input, config):
    if config['verbose']:
        utils.info('Starting proteomes processing...\n')
    d_proteomes = {}
    
    chunksize = config['nproc']
    
    # r = re.compile('.*(trembl|sprot|uniparc).*')
    r = re.compile('.*(uniparc).*.pkl')
    chunks = []

    tree = etree.iterparse(gzip.GzipFile(xml_input), events = ('end',), tag = 'proteome', huge_tree = True)
    try:
        with mpdummy.Pool(processes=chunksize) as pool:
            for prot in pool.imap_unordered(parse_proteomes_xml_elem, yield_filtered_xml_string(tree), chunksize=chunksize):
                if prot is not None:
                    d_proteomes[prot[0]] = {'tax_id' : prot[1],
                                            'upi' : prot[2],
                                            'isReference' : prot[3],
                                            'ncbi_ids' : prot[4],
                                            'members' : prot[5]
                                        }
    except Exception as e:
        utils.error(str(e))
        utils.error('Processing failed')
        raise

    pickle_to_process = filter(r.match,glob.iglob("{}/{}/uniparc/*".format(config['download_base_dir'], config['pickled_dir'])))
       
    try:
        with mp.Pool(processes=chunksize) as pool:
            for item in pool.imap_unordered(extract_proteomes_from_pickle, pickle_to_process, chunksize = chunksize):
                for k,v in item.items():
                    if k not in d_proteomes:
                        d_proteomes[k] = {'members' : set(), 
                                          'isReference' : v['isReference'], 
                                          'tax_id' : v['tax_id'], 
                                          'upi' : v['upi']}
                    if d_proteomes[k]['upi'] and v['upi']:
                        d_proteomes[k]['members'] = set(d_proteomes[k]['members'])
                        d_proteomes[k]['members'].update(v['members'])
        
    except Exception as e:
        utils.error(str(e))
        utils.error('Processing failed')
        raise
  
    with open("{}{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'wb') as pickle_proteomes:
        pickle.dump(d_proteomes, pickle_proteomes, -1)
    if config['verbose']:
        utils.info('Done\n')

def annotate_taxon_tree(config):
    utils.info('Starting annotation of taxonomic tree with proteome ids\n')

    with open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_proteomes']),'rb') as pickle_proteomes, \
         open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb') as pickle_taxontree:
        proteomes = pickle.load(pickle_proteomes)
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
    uniparc_idmapping = {}
    ur_clusters = ('uniref100','uniref90','uniref50')

    if all([os.path.exists("{}/{}/{}_uniparc_idmap.pkl".format(config['download_base_dir'], config['pickled_dir'], c)) for c in ur_clusters]) and \
    os.path.exists('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap'])):
        if config['verbose']:
            utils.info('Started merging UniParc-UniRef idmapping.\n')
        for c in ur_clusters:
            cluster_idmapping = "{}/{}/{}_uniparc_idmap.pkl".format(config['download_base_dir'],config['pickled_dir'], c)
            
            for upi, cluster in itertools.chain.from_iterable(utils.load_pickle(cluster_idmapping)):
                if upi not in uniparc_idmapping:
                    uniparc_idmapping[upi] = set()
                uniparc_idmapping[upi].add(cluster)
            os.unlink("{}/{}/{}_uniparc_idmap.pkl".format(config['download_base_dir'],config['pickled_dir'], c))

        for upi, clusters in uniparc_idmapping.items():
            ur100 = [x.split('_')[1] for x in clusters if 'UniRef100_' in x]
            ur90 = [x.split('_')[1] for x in clusters if 'UniRef90_' in x]
            ur50 = [x.split('_')[1] for x in clusters if 'UniRef50_' in x]
            clusters = [ur100[0] if len(ur100) else '', ur90[0] if len(ur90) else '', ur50[0] if len(ur50) else '']
            uniparc_idmapping[upi] = tuple(clusters)
        
        idmapping = '{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap'])

        for upkb, cluster in utils.load_pickle(idmapping):
            uniparc_idmapping[upkb] = cluster

        #update to global
        with open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'wb') as f_out:
            pickle.dump(uniparc_idmapping, f_out, -1)

        if config['verbose']:
            utils.info('Done.\n')

    return uniparc_idmapping

def merge_idmap(config):
    x = ["{}/{}/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'],config['pickled_dir'],  'sprot'), 
         "{}/{}/uniprotkb_{}_idmap.pkl".format(config['download_base_dir'], config['pickled_dir'], 'trembl'),
         "{}/{}/{}_idmap.pkl".format(config['download_base_dir'], config['pickled_dir'], 'uniparc')]
    utils.info('Merging UniProtKB and UniParc ids...\n')

    with open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_uniprotkb_idmap']), 'wb') as idmap_pkl:
        [utils.optimized_dump(idmap_pkl, elem) for fin in x for elem in utils.load_pickle(fin)]

    utils.info('Done merging UniProtKB ids.\n')

def annotate_with_sgb(config):
    # taxontree = pickle.load(open("{}/{}".format(config['download_base_dir'],config['relpath_pickle_taxontree']),'rb'))
    # proteomes
    pass

def process_proteomes(config):
    os.makedirs('{}/{}'.format(config['download_base_dir'],config['pickled_dir']), exist_ok=True)
    if not os.path.exists("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], 'uniprotkb')):
        os.makedirs("{}/{}/{}".format(config['download_base_dir'], config['pickled_dir'], 'uniprotkb'))

    step1 = [ mp.Process(target=parse_uniref_xml, args=(config['download_base_dir']+config['relpath_uniref100'],config)),
              mp.Process(target=parse_uniref_xml, args=(config['download_base_dir']+config['relpath_uniref90'],config)),
              mp.Process(target=parse_uniref_xml, args=(config['download_base_dir']+config['relpath_uniref50'],config))
            ]
    # for p in step1:
    #     p.start()

    # for p in step1:
    #     p.join()
        
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
    uniprotkb_uniref_idmap = merge_uniparc_idmapping(config)
    
    if not uniprotkb_uniref_idmap:
        utils.info('Loading UniProtKB-UniRef mapping...')
        uniprotkb_uniref_idmap = pickle.load(open('{}{}'.format(config['download_base_dir'],config['relpath_pickle_uniprotkb_uniref_idmap']),'rb'))
        utils.info('Done.\n')
    
    step2 = [ mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_sprot'], config)),
              mp.Process(target=parse_uniprotkb_xml, args=(config['download_base_dir']+config['relpath_uniprot_trembl'], config)),
              mp.Process(target=parse_uniparc_xml, args=(config['download_base_dir']+config['relpath_uniparc'],config))
            ]
    for p in step2:
        p.start()

    for p in step2:
        p.join()

    merge_idmap(config)
    create_proteomes_pkl(config['download_base_dir']+config['relpath_proteomes_xml'], config)
    annotate_taxon_tree(config)
    annotate_with_sgb(config)

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
