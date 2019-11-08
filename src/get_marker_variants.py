__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicola Segata (nicola.segata@unitn.it)')

__date__ = '18 Oct 2019'

import numpy
import itertools
import pandas as pd
import pickle
import bz2
import glob
import os
import re
import errno
from Bio import SeqIO
from functools import partial
from pyfaidx import Fasta
from pyfaidx import FastaIndexingError
import multiprocessing as mp
import multiprocessing.dummy as dummy
import importlib
from _version import __UniRef_version__
from tempfile import NamedTemporaryFile
import _version as version
import logging
import datetime
import shutil
import subprocess as sb
import pybedtools
pybedtools.set_tempdir('/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/tmp')
VSEARCH_PATH = '/shares/CIBIO-Storage/CM/mir/tools/vsearch-2.13.6/bin/vsearch'
if __name__ == '__main__':
    import utils
    from panproteomes import Panproteome
else:
    utils = importlib.import_module('src.utils')
    from src.panproteomes import Panproteome

def load_file(f_path):
    with open(f_path, 'rt') as ifn:
        for line in ifn:
            yield line.strip().split()

OUTFILE_PREFIX = 'mpa_v{}_CHOCOPhlAn_{}'.format(version.__MetaPhlAn2_db_version__, version.__CHOCOPhlAn_version__)

log = logging.getLogger(__name__)
ch = logging.FileHandler('CHOCOPhlAn_get_marker_variants_{}.log'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')),'w')
ch.setLevel(logging.INFO)
log.addHandler(ch)

def run_vsearch(args):
    vsearch, cluster_fast, uc = args
    vsearch_comm = [vsearch, '--cluster_fast', os.path.abspath(cluster_fast),
                            '--id' , '0.97',
                            '--uc', os.path.abspath(uc)
                    ]
    try:
        # print(' '.join(bt2_mapping_comm))
        xx = sb.check_call(vsearch_comm, stdout=sb.DEVNULL, stderr=sb.PIPE)
    except sb.CalledProcessError as e:
        utils.error(str(e))

    return xx


def load_file(f_path):
    with open(f_path, 'rt') as ifn:
        for line in ifn:
            yield line.strip().split('\t')

def build_pangenome(taxid, config):
    species_genomes = glob.iglob('{0}/{1}/{2}/panphlan_{2}_genomes/*'.format(config['export_dir'], config['exportpath_panphlan'], taxid))
    with mp.Pool(processes=10) as pool:
        panphlan_csv = pd.read_csv('{0}/{1}/{2}/panphlan_{2}_pangenome.csv'.format(config['export_dir'], config['exportpath_panphlan'], taxid), sep='\t', names = ['nr90','gene_name','gca','contig','start','end'])
        partial_get_ffn_from_fna = partial(get_ffn_from_fna, panphlan_csv=panphlan_csv)
        d = [x for x in pool.imap_unordered(partial_get_ffn_from_fna, species_genomes, chunksize=1)]
    species_genomes = glob.iglob('{0}/{1}/{2}/panphlan_{2}_genomes/*'.format(config['export_dir'], config['exportpath_panphlan'], taxid))

    with open('{0}/{1}/{2}/panphlan_{2}_pangenome.ffn'.format(config['export_dir'], config['exportpath_panphlan'], taxid), 'wt') as ffn_out:
        for genome in species_genomes:
            genome = os.readlink(genome)
            try:
                with open(genome.replace('fna','ffn')[:-3], 'rt') as ffn_in:
                    ffn_out.write(ffn_in.read())
                    ffn_out.write('\n')
                os.unlink(genome.replace('fna','ffn')[:-3])
            except FileNotFoundError as ex:
                log.error('[{}]\tffn for {} does not exist.'.format(taxid, genome))
    
    with open('{0}/{1}/{2}/panphlan_{2}_pangenome.tochange'.format(config['export_dir'], config['exportpath_panphlan'], taxid), 'wt') as tochange:
        tochange.write(
            '\n'.join(
                ['{}:{}-{}\t{}__{}__{}'.format(line[3],int(line[4])-1, line[5], taxid, line[0][9:], line[1]) 
                    for line in load_file('{0}/{1}/{2}/panphlan_{2}_pangenome.csv'.format(config['export_dir'], config['exportpath_panphlan'], taxid))
                ]
            )
        )

    uc = run_vsearch([VSEARCH_PATH, 
                '{0}/{1}/{2}/panphlan_{2}_pangenome.ffn'.format(config['export_dir'], config['exportpath_panphlan'], taxid), 
                '{0}/{1}/{2}/panphlan_{2}_pangenome.clusters'.format(config['export_dir'], config['exportpath_panphlan'], taxid)
    ])

    if uc==0:
        log.info('[{}] VSEARCH clustering has beeen performed.\n'.format(taxid))

def get_ffn_from_fna(genome, panphlan_csv):
    genome = os.readlink(genome)
    gff = pybedtools.BedTool(genome.replace('fna','gff'))
    with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/tmp/') as fna_decompressed:
        with gzip.open(genome, 'rb') as fna_gz, open(fna_decompressed.name, 'wb') as fna_out:
            fna_out.write(fna_gz.read())

        a = gff.sequence(fi=fna_decompressed.name, s=True, split=True)
        genome_id = re.search('GCA_[0-9]{9}', genome)
        panphlan_csv = panphlan_csv.loc[panphlan_csv.gca == genome_id.group(),]
        with Fasta(a.seqfn, duplicate_action="first", key_function = lambda x: x.split('(')[0]) as ffn_in, open(genome.replace('fna','ffn')[:-3], 'wt') as ffn_out:
            ids_to_extract = set('{}:{}-{}'.format(x.contig, x.start-1, x.end) for x in panphlan_csv.itertuples())
            ids_to_extract = ids_to_extract.intersection(ffn_in.keys())
            ffn_out.write('\n'.join(
                ['>{}\n{}'.format(x, ffn_in[x]) for x in ids_to_extract]
            ))

def load_clustering(config, taxid, contig2ass):
    uc_fp = '{}/{}/{}/panphlan_{}_pangenome.clusters'.format(config['export_dir'], config['exportpath_panphlan'], taxid, taxid)
    clusters = pd.read_csv( uc_fp,
                            sep = '\t', 
                            names = ['type', 'c_number', 'length', 'perc_id', 'strand', 'none1', 'none2', 'comp_aln', 'qseqid', 'tseqid']
                        )
    
    marker_map = pd.read_csv( uc_fp.replace('clusters', 'tochange'),
                            sep = ' ',
                            names = ['uc_name', 'mpa_name']
                        )

    contigs = marker_map.uc_name.str.split(':', expand=True).iloc[:,0]
    marker_map = marker_map.assign(contigs = contigs)
    marker_map = marker_map.merge(contig2ass, on = ['contigs', 'contigs'])
    clusters = clusters.merge(marker_map, left_on = 'qseqid', right_on = 'uc_name')

    return clusters

def get_hnn_variants_from_uc(taxid, contig2ass, config):
    humann_centroids_fp = '{}/{}/panproteomes_fna/{}.fna.bz2'.format(config['export_dir'], config['exportpath_humann2'], taxid)
    with bz2.open(humann_centroids_fp, 'rt') as hnn_cf:
        hnn_centroids = [ line[1:].strip() for line in hnn_cf if line[0] == '>']
    hnn_centroids = {x.split('|')[0] : x for x in hnn_centroids}

    try:
        clusters = load_clustering(config, taxid, contig2ass)
    except Exception as e:
        log.error("[{}] {}".format(taxid,e))
        return

    def get_centroid_variant(clusters, centroid):
        cname, length = centroid.split('|')[0],centroid.split('|')[-1]
        retval = (cname, )
        marker_cluster = clusters.loc[(clusters.mpa_name == cname) & (clusters.length.between(int(length)-1,int(length)+1) ),'c_number']
        if marker_cluster.empty:
            marker_cluster = clusters.loc[(clusters.mpa_name.str.startswith(cname.rsplit('__',1)[0])) & (clusters.length.between(int(length)-1,int(length)+1) ),'c_number']
        if not marker_cluster.empty:
            marker_cluster = marker_cluster.unique()[0]
            marker_hits = clusters.loc[clusters.c_number == marker_cluster,]
            
            if marker_hits.shape[0] > 2:
                marker_hits = marker_hits.sort_values(by='perc_id').loc[marker_hits.perc_id != '100.0', ]
                variants_id = tuple(marker_hits.groupby('comp_aln').head(n=1).loc[marker_hits.comp_aln != '*', ].qseqid)
                if len(variants_id) > 1:
                    retval = (cname, variants_id,)

        return retval

    variants = list(filter(None, (get_centroid_variant(clusters, c) for c in hnn_centroids.values())))

    with NamedTemporaryFile(dir='/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/tmp/') as fna_decompressed:
        with bz2.BZ2File(humann_centroids_fp, 'rb') as fna_bz2, open(fna_decompressed.name, 'wb') as fna_out:
            fna_out.write(fna_bz2.read())
            try:
                with Fasta('{}/{}/{}/panphlan_{}_pangenome.ffn'.format(config['export_dir'], config['exportpath_panphlan'], taxid, taxid), duplicate_action="first") as ffn_panphlan, \
                    Fasta(fna_decompressed.name, key_function = lambda x: x.split('|')[0]) as ffn_humann:
                    with bz2.open('{}/{}/panproteome_variants/{}.fna.bz2'.format(config['export_dir'], config['exportpath_humann2'], taxid), 'wt') as fout:
                        for v in variants:
                            fout.write('>{}\n{}\n'.format(hnn_centroids[v[0]], ffn_humann[v[0].split(' ')[0]][:].seq ))
                            if len(v) == 2:
                                fout.write('\n'.join('>{}/{} {}\n{}'.format(hnn_centroids[v[0]],i,x, ffn_panphlan[x][:].seq) for i, x in enumerate(v[1],1) ))
                                fout.write('\n')
            except FastaIndexingError as fie:
                log.error(fie)
                return

def get_variant_name(clusters, marker):
    retval = (marker, )
    marker_cluster = clusters.loc[(clusters.mpa_name == marker[0]) & (clusters.genome == marker[1]),'c_number']
    if marker_cluster.empty:
        marker_cluster = clusters.loc[(clusters.mpa_name.str.startswith(marker[0].rsplit('__',1)[0])) & (clusters.genome == marker[1]),'c_number']
    if not marker_cluster.empty:
        marker_cluster = marker_cluster.unique()[0]
        marker_hits = clusters.loc[clusters.c_number == marker_cluster,]
        if marker_hits.shape[0] > 2:
            marker_hits = marker_hits.sort_values(by='perc_id').loc[marker_hits.perc_id != '100.0', ]
            variants_id = tuple(marker_hits.groupby('comp_aln').head(n=1).loc[marker_hits.comp_aln != '*', ].qseqid)
            if len(variants_id) > 1:
                retval = (marker, variants_id,)
    
    return retval


def get_mpa_variants_from_uc(taxid, all_markers, contig2ass, config): 
    markers_descr = [ marker for marker in all_markers if marker[0].startswith(taxid+'__')]
    if not any([len(x[0].split('__'))==2 for x in markers_descr]):
        markers = [(x[0], x[1].split(';')[-1], ) if len(x)==2 else (x[0], x[0].split('__')[0].rsplit('_', 1), ) for x in markers_descr]

        markers_descr = dict(markers_descr) if all(len(x)==2 for x in markers_descr) else dict(zip(itertools.chain.from_iterable(markers_descr), ['']*len(markers_descr)))

        try:
            clusters = load_clustering(config, taxid, contig2ass)
        except Exception as e:
            log.error('[{}]\t{}'.format(taxid, e))
            return

        try:
            variants = list(filter(None, (get_variant_name(clusters, marker) for marker in markers)))
        except Exception as ex:
            print('{}\t{}'.format(taxid, ex))
            return

        variants = list(filter(lambda x: len(x) == 2, variants))
        if variants:
            with Fasta('{0}/{1}/{2}/panphlan_{2}_pangenome.ffn'.format(config['export_dir'], config['exportpath_panphlan'], taxid), duplicate_action="first") as ffn:
                with bz2.open('{0}/{1}/{2}_variants/marker_variants/{3}.fna'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX, taxid), 'wt') as fout:
                    for v in variants:
                        names = clusters.loc[(clusters.mpa_name == v[0][0]) & (clusters.genome == v[0][1]),['uc_name', 'mpa_name']]
                        if names.empty:
                            names = clusters.loc[(clusters.mpa_name.str.startswith(v[0][0].rsplit('__',1)[0])) & (clusters.genome == v[0][1]),['uc_name', 'mpa_name']]
                        if not names.empty:
                            uc_name, mp_name = names.iloc[0,]
                            fout.write('\n'.join('>{}/{} {}\n{}'.format(mp_name,i,x, ffn[x][:].seq) for i, x in enumerate(v[1],1) ))
                            fout.write('\n')
            log.info('Species {} processed'.format(taxid))
        else:
            log.info('No variants found for species {}'.format(taxid))
    else:
        log.info('Species {} with no standard marker name'.format(taxid))
    
    

def run_all(config):
    os.makedirs('{0}/{1}/panproteome_variants'.format(config['export_dir'],config['exportpath_humann2']), exist_ok=True)
    os.makedirs('{0}/{1}/{2}_variants/marker_variants'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), exist_ok=True)
    mpa_pkl = pickle.load(bz2.open('{}/{}/{}/{}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX, OUTFILE_PREFIX)))
    all_markers = [x[1:].strip().split() for x in bz2.open('{}/{}/{}/{}.fna
    .bz2'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX, OUTFILE_PREFIX), 'rt') if x[0] == '>' ]

    all_species_mpa = set(marker[0].split('_')[0] for marker in all_markers)
    all_species_hnn = set(x.split('/')[-1].split('.')[0] for x in glob.iglob(os.path.join(config['export_dir'],config['exportpath_humann2'],'panproteomes_fna','*')))
    
    exported_pangenomes_dir = '{}/{}'.format(config['export_dir'], config['exportpath_panphlan'])
    
    contig2ass = { contig : genome for contig,genome in load_file('{}/mpa2_eval/bt2_idx/contig2genome.tsv'.format(config['export_dir'])) }
    contig2ass = { contig : genome[:13] for contig,genome in contig2ass.items() }
    contig2ass = pd.DataFrame.from_dict(contig2ass, orient='index')
    contig2ass.index.name = 'contigs'
    contig2ass.reset_index(inplace = True)
    contig2ass.columns = ['contigs', 'genome']
    

    partial_build_pangenome = partial(build_pangenome, config=config)

    with dummy.Pool(processes=30) as pool:
        d = [x for x in pool.imap_unordered(partial_build_pangenome, all_species_hnn, chunksize=10)]

    partial_get_mpa_variants_from_uc = partial(get_mpa_variants_from_uc, all_markers=all_markers, contig2ass=contig2ass, config=config)
    partial_get_hnn_variants_from_uc = partial(get_hnn_variants_from_uc, contig2ass=contig2ass, config=config)

    with mp.Pool(processes=30) as pool:
        d = [x for x in pool.imap_unordered(partial_get_mpa_variants_from_uc, all_species_mpa, chunksize=10)]

    with bz2.open('{0}/{1}/{2}_variants/{2}_variants.fna.bz2'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'wt') as mpa_fna_variants:
        for variant in glob.iglob('{0}/{1}/{2}_variants/marker_variants/*.fna.bz2'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)):
            with bz2.open(variant, 'rt') as v_in:
                mpa_fna_variants.write(v_in.read())
        
        with bz2.open('{0}/{1}/{2}/{2}.fna.bz2'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX), 'rt') as mpa_fna:
            mpa_fna_variants.write(mpa_fna.read())

    try:
        os.symlink(os.path.abspath('{0}/{1}/{2}/{2}.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)), 
                   os.path.abspath('{0}/{1}/{2}_variants/{2}_variants.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX))
        )
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(os.path.abspath('{0}/{1}/{2}_variants/{2}_variants.pkl'.format(config['export_dir'], config['exportpath_metaphlan2'], OUTFILE_PREFIX)))


    with mp.Pool(processes=30) as pool:
        d = [x for x in pool.imap_unordered(partial_get_hnn_variants_from_uc, all_species_hnn, chunksize=1)]

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    run_all(config['export'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
