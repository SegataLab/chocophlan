#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '01 Oct 2017'


import multiprocessing as mp
import time
import ftplib
import math
import sys
import os
if __name__ == '__main__':
    import utils
else:
    import src.utils as utils
import re
from functools import partial

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def do_download(inputs, verbose):
    if not terminating.is_set():
        try:
            ftp_base, full_link, output_path = inputs
            
            ftp = ftplib.FTP(ftp_base)  # Login to ftp server
            ftp.login()

            dir_path = os.path.dirname(output_path)

            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)

            if not os.path.exists(output_path):
                with open(output_path, "wb") as fileout:
                    if verbose:
                        utils.info("Downloading {}\n".format(full_link))
                    ftp.retrbinary("RETR " + full_link, fileout.write)
                    ftp.quit()
            else:
                if verbose:
                    utils.info("File {} already present\n".format(output_path))
        except Exception as e:
            terminating.set()
#            utils.remove_file(output_path)
            utils.error(str(e))
            #utils.error('Download failed for\n    {}'.format(full_link))
            raise
    else:
        terminating.set()


def download(config, verbose=False):
    # TODO: MD5 checksum testing!
    # TODO: Check if folders exist!
    # TODO: Incoorporate the download of the large uniref file!
    # TODO: The refseq catalogue needs to have the version included in its name!

    ### uniprot XML ###
#    argument_list = [(config['uniprot_ftp_base'],
#                      config['uniprot_uniref100'],
#                      config['download_base_dir'] + config['relpath_uniref100'])]
#
#    ### Bacterial refseq genomes ###
#    ftp = ftplib.FTP(config['refseq_ftp_base'])
#    ftp.login()
#    ftp.cwd(config['refseq_bacterial_genomes'])
#    ls = ftp.nlst()
#    ftp.quit()
#    argument_list += [(config['refseq_ftp_base'],
#                       '/'.join([config['refseq_bacterial_genomes'], entry]),
#                       '/'.join([config['download_base_dir'], config['relpath_bacterial_genomes'],
#                                 entry]))
#                      for entry in ls if "genomic.fna.gz" in entry]
#    ### RefSeq catalogue ###
#    argument_list.append((config['refseq_ftp_base'],
#                          config['refseq_taxonomic_catalogue'],
#                          config['download_base_dir'] + config['relpath_taxonomic_catalogue']))
#
#    ### Refseq taxdump ###
#    argument_list.append((config['refseq_ftp_base'],
#                          config['refseq_taxdump'],
#                          config['download_base_dir'] + config['relpath_taxdump']))
#
#    ### UniProt Reference Proteomes ###
#    ftp = ftplib.FTP(config['uniprot_ftp_base'])
#    ftp.login()
#    ftp.cwd(config['uniprot_reference_proteomes'])
#    ls = ftp.nlst()
#    ftp.quit()
#
#    r = re.compile("Reference_Proteomes_.*\.tar\.gz")
#    ref_prot = [x for x in filter(r.match, ls)][0]
#
#    argument_list.append((config['uniprot_ftp_base'],
#                          "{}/{}".format(config['uniprot_reference_proteomes'],ref_prot),
#                          config['download_base_dir'] + config['relpath_reference_proteomes']))
    ### Pan proteome download ###
    ftp = ftplib.FTP(config['uniprot_ftp_base'])
    ftp.login()
    ftp.cwd(config['uniprot_pan_proteomes'])
    ls = ftp.nlst()
    ftp.quit()
    argument_list = []
    argument_list += [(config['uniprot_ftp_base'],
                       '/'.join([config['uniprot_pan_proteomes'], entry]),
                       '/'.join([config['download_base_dir'], config['relpath_pan_proteomes'],
                                entry]))
                       for entry in ls if "fasta.gz" in entry]
    
    terminating = mp.Event()
    #chunksize = math.floor(len(argument_list) / (int(config['nproc']) * 2))
    with mp.Pool(initializer=initt, initargs=(terminating,),
                 processes=config['nproc']) as pool:
        try:
            if verbose:
                utils.info("Starting parallel download\n")
            
            # Need to define a partial function because function passed via imap require only one argument
            do_download_partial = partial(do_download,verbose=config['verbose'])

            [_ for _ in pool.imap(do_download_partial, argument_list,
                                  chunksize=1)]
        except Exception as e:
            utils.error(str(e))
            utils.error('Download failed', exit=True)
            raise

    if verbose:
        utils.info('Download succesful\n')


if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)
    
    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    print(config)
    download(config['download'], verbose=config['download']['verbose'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
