#!/usr/bin/env python

""" Given an gca and an output folder, download the corresponding gff and fna files

To run:
$ ./download_gff_fna.py id_file.txt downloads_folder

Will download the gff and fna files for all ids in the input file to a
folder named "downloads_folder" in the current working directory.
"""
author__ = ('Nicola Segata (nicola.segata@unitn.it), '
            'Francesco Beghini (francesco.beghini@unitn.it)'
            'Nicolai Karcher (karchern@gmail.com),'
            'Francesco Asnicar (f.asnicar@unitn.it)'
            'Lauren McIver (lauren.j.mciver@gmail.com)')

__version__ = '0.01'
__date__ = '04 Jan 2018'

import os
import sys
import tempfile
import time
import datetime
import pickle
from functools import partial
import multiprocessing.dummy as dummy

if __name__ == '__main__':
    import utils
    from process_proteomes import initt
else:
    import src.utils as utils
    from src.process_proteomes import initt

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
    from urllib.request import urlcleanup
except ImportError:
    from urllib import urlretrieve
    from urllib import urlcleanup

NCBI_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
GCA_COLUMN=0
BSA_COLUMN=2
TAX_COLUMN=5
FTP_COLUMN=19
GFF_EXTENSION="_genomic.gff.gz"
FNA_EXTENSION="_genomic.fna.gz"

def download_file(url,file):
    """ Download the url to the file location """

    # try to download the file
    try:
        file, headers = urlretrieve(url,file)
        status = 0
    except:
        status = 1

    # clean the cache, fixes issue with downloading
    # two ftp files one after the next in python2 (some versions)
    urlcleanup()
    return status

def get_ncbi_assembly_info():
    """ Download the assembly data from NCBI """

    # create a tempfile for the download
    file_handle, new_file = tempfile.mkstemp(prefix="ncbi_download")

    # try to download the file
    utils.info("Downloading the assembly data info from NCBI\n")
    download_file(NCBI_URL,new_file)
    
    # read in the file, ignoring headers
    utils.info("Loading assembly data...\n")
    data = []
    for line in open(new_file):
        if not line.startswith("#"):
            data.append(line.rstrip().split("\t"))

    # remove the temp file
    os.remove(new_file)

    return data

def get_gca_ftp(name,data=None):
    """ Get the full ftp path to the basename for all downloads for a given gca """

    if not data:
        data = get_ncbi_assembly_info()

    # allow for different versions of the gca
    name=name.split(".")[0]+"."

    #u tils.info("Searching for ftp: {}\n".format(name))
    ftp=None
    for line in data:
        if line[GCA_COLUMN].startswith(name):
            ftp=line[FTP_COLUMN]+"/"+os.path.basename(line[FTP_COLUMN])
            break
    if not ftp or not ftp.startswith("ftp"):
        ftp=None
        utils.error("ERROR: Unable to find ftp: {}".format(name))

    return ftp

def file_size(file):
    """ Return the size of the file in MB """

    try:
        size = os.path.getsize(file) / (1024.0**2)
    except OSError:
        size = 0

    return size

def download_gff_fna(gca_name, folder, data=None, delete=None, total_size=0):
    """ Download the gff and fna for the given gca to the folder specified """

    # get the url for the gff and fna
    base_ftp=get_gca_ftp(gca_name, data=data)
    if not base_ftp:
        return
    gff_ftp=base_ftp+GFF_EXTENSION
    fna_ftp=base_ftp+FNA_EXTENSION

    # download the files
    gff_file=os.path.join(folder,os.path.basename(gff_ftp))
    fna_file=os.path.join(folder,os.path.basename(fna_ftp))

    status = 0
    for url,file in zip([gff_ftp,fna_ftp],[gff_file,fna_file]):
        status += download_file(url,file)
        if delete:
            os.remove(file)
            utils.info("Deleted: " + file)
    return status

def process(item, data, config):
    if not terminating.is_set():
        k,v = item
        ncbi_ids = dict(v['ncbi_ids'])
        id_ = ''
        if 'GCSetAcc' in ncbi_ids:
            id_ = ncbi_ids['GCSetAcc']
        elif 'Biosample' in ncbi_ids:   
            id_ = [line[GCA_COLUMN] for line in data if line[BSA_COLUMN] == ncbi_ids['Biosample'] and int(line[TAX_COLUMN]) == v['tax_id']]
            if len(id_): id_=id_[0]
        if len(id_):
            ncbi_ids['GCSetAcc'] = id_
            try:
                download_folder = '{}/{}/{}/{}/{}'.format(config['download_base_dir'], config['relpath_genomes'], id_.split('_')[1][0:3],id_.split('_')[1][3:6],id_.split('_')[1][6:9])
            except:
                utils.info('{}\n'.format(id_))
            os.makedirs(download_folder, exist_ok=True)
            status = download_gff_fna(id_,download_folder,data)
            if status:
                utils.info("Download for {} has failed!\n".format(id_))
                return id_
    else:
        terminating.set()

def download_ncbi(config):
    # download the assembly data once
    data = get_ncbi_assembly_info()
    terminating = dummy.Event()
    utils.info("Loading proteomes data...\n")
    proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))

    partial_process = partial(process, data=data, config=config)
    with dummy.Pool(initializer=initt, initargs=(terminating, ), processes=config['nproc']) as pool:
        failed = [f for f in pool.imap_unordered(partial_process, [(k,v) for k, v in proteomes.items() if 'ncbi_ids' in v], chunksize=config['nproc'])]

    with open('failed_GCA.txt','w') as f:
        f.writelines(faield)
            
if __name__ == "__main__":
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['process_proteomes']
    config['relpath_genomes'] = '/ncbi'

    download_ncbi(config)

