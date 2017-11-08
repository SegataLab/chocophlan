#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
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


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def do_download(inputs):
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
                    ftp.retrbinary("RETR " + full_link, fileout.write)
            else:
                utils.info("File {} already present\n".format(output_path))
        except Exception as e:
            terminating.set()
            utils.remove_file(output_path)
            utils.error(str(e), init_new_line=True)
            utils.error('Download failed for\n    {}'.format(full_link),
                        init_new_line=True)
            raise
    else:
        terminating.set()


def download(config, verbose=False):
    # TODO: MD5 checksum testing!
    # TODO: Check if folders exist!
    # TODO: Incoorporate the download of the large uniref file!
    # TODO: The refseq catalogue needs to have the version included in its name!

    ### uniprot XML ###
    argument_list = [(config['uniprot_ftp_base'],
                      config['uniprot_uniref100'],
                      config['download_base_dir'] + config['relpath_uniref100'])]

    ### uniprot reference proteomes ###
    argument_list += [(config['uniprot_ftp_base'],
                      config['uniprot_reference_proteomes'],
                      config['download_base_dir'] + config['relpath_reference_proteomes'])]

    ### Bacterial refseq genomes ###
    ftp = ftplib.FTP(config['refseq_ftp_base'])
    ftp.login()
    ftp.cwd(config['refseq_bacterial_genomes'])
    ls = ftp.nlst()
    argument_list += [(config['refseq_ftp_base'],
                       '/'.join([config['refseq_bacterial_genomes'], entry]),
                       '/'.join([config['download_base_dir'], config['relpath_bacterial_genomes'],
                                 entry]))
                      for entry in ls if "genomic.fna.gz" in entry]

    ### RefSeq catalogue ###
    argument_list.append((config['refseq_ftp_base'],
                          config['refseq_taxonomic_catalogue'],
                          config['download_base_dir'] + config['relpath_taxonomic_catalogue']))

    ### Refseq taxdump ###
    argument_list.append((config['refseq_ftp_base'],
                          config['refseq_taxdump'],
                          config['download_base_dir'] + config['relpath_taxdump']))

    ### UniProt Reference Proteomes ###
    ftp = ftplib.FTP(config['uniprot_ftp_base'])
    ftp.login()
    ftp.cwd(config['uniprot_ref_proteomes'])
    ls = ftp.nlst()

    r = re.compile("Reference_Proteomes_.*\.tar\.gz")
    ref_prot = [x for x in filter(r.match, ls)][0]

    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_ref_proteomes'],
                          "/" + ref_prot,
                          config['download_base_dir'] + config['relpath_refprot']))

    terminating = mp.Event()
    chunksize = math.floor(len(argument_list) / (int(config['nproc']) * 2))

    with mp.Pool(initializer=initt, initargs=(terminating,),
                 processes=config['nproc']) as pool:
        try:
            if verbose:
                utils.info("Starting parallel download\n", init_new_line=True)

            [_ for _ in pool.imap(do_download, argument_list,
                                  chunksize=chunksize if chunksize else 1)]
        except Exception as e:
            utils.error(str(e), init_new_line=True)
            utils.error('Download failed', init_new_line=True, exit=True)
            raise

    if verbose:
        utils.info('Download succesful\n', init_new_line=True)


if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    download(config['download'], verbose=config['download']['verbose'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
