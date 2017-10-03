#!/usr/bin/env python3

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '01 Oct 2017'


if __name__ == '__main__':
    import utils
else:
    import src.utils as utils
import multiprocessing as mp
import time
import ftplib
import math
import sys

# The initt function as well as most of the do_download functionality have been taken from Francesco Asnicar's PhyloPhlAn2.
def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def do_download(inputs):
    ftp_base, full_link, output_path = inputs
    # Download
    if not terminating.is_set():
        try:
            # Code to download!
            # Login to ftp server
            ftp = ftplib.FTP(ftp_base)
            ftp.login()
            # Download
            with open(output_path, "wb") as fileout:
                ftp.retrbinary("RETR " + full_link, fileout.write)
        except Exception as e:
            terminating.set()
            utils.error(str(e), init_new_line=True)
            utils.error('Download failed for\n    {}'.format(ftp_link), init_new_line=True)
            raise
    else:
        terminating.set()

def download(config):
    # TODO: MD5 checksum testing!
    # TODO: Check if folders exist! 
    # TODO: Incoorporate the download of the large uniref file!
    # TODO: The refseq catalogue needs to have the version included in its name!

    # Create argument list for parallel download of files
    argument_list = []

    ### uniprot XML ###
    argument_list.append((config['uniprot_ftp_base'], config['uniprot_uniref100'], config['download_dir'] + '/uniprot/uniref/uniref100.xml.gz'))

    ### Bacterial refseq genomes ###
    ftp = ftplib.FTP(config['refseq_ftp_base'])
    ftp.login() 
    ftp.cwd(config['refseq_bacterial_genomes'])
    ls = ftp.nlst()
    refseq_arguments = [(config['refseq_ftp_base'], config['refseq_bacterial_genomes'] + '/' + entry, config['download_dir'] + '/refseq/genomes/' + entry) for entry in ls if "genomic.fna.gz" in entry]
    argument_list.extend(refseq_arguments) # Use extend instead of append, as append would append refseq_arguments as a list.

    ### RefSeq catalogue ###
    argument_list.append((config['refseq_ftp_base'], config['refseq_taxonomic_catalogue'], config['download_dir'] + '/refseq/catalogue/refseq_catalogue.gz'))

    terminating = mp.Event()
    chunksize = math.floor(len(argument_list) / (int(config['nproc']) * 2))
    with mp.Pool(initializer=initt, initargs=(terminating,),
             processes=int(config['nproc'])) as pool:
        try:
            utils.info("Starting parallel download.", init_new_line=True)
            [_ for _ in pool.imap(do_download, argument_list, chunksize=chunksize if chunksize else 1)]
        except Exception as e:
            utils.error(str(e), init_new_line=True)
            utils.error('Download failed.', init_new_line=True,
                  exit=True)
            raise
    utils.info('Download succesful.', init_new_line=True)
if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['download'] 

    download(config)

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
