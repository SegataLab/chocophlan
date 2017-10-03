#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '01 Oct 2017'


from src.utils import *


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

            with open(output_path, "wb") as fileout:
                ftp.retrbinary("RETR " + full_link, fileout.write)
        except Exception as e:
            terminating.set()
            error(str(e), init_new_line=True)
            error('Download failed for\n    {}'.format(full_link),
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
                      config['download_dir'] + '/uniprot/uniref/uniref100.xml.gz')]

    ### Bacterial refseq genomes ###
    ftp = ftplib.FTP(config['refseq_ftp_base'])
    ftp.login()
    ftp.cwd(config['refseq_bacterial_genomes'])
    ls = ftp.nlst()
    argument_list += [(config['refseq_ftp_base'],
                       '/'.join([config['refseq_bacterial_genomes'], entry,
                                 config['download_dir'], 'refseq/genomes',
                                 entry]))
                      for entry in ls if "genomic.fna.gz" in entry]

    ### RefSeq catalogue ###
    argument_list.append((config['refseq_ftp_base'],
                          config['refseq_taxonomic_catalogue'],
                          config['download_dir'] + '/refseq/catalogue/refseq_catalogue.gz'))

    terminating = mp.Event()
    chunksize = math.floor(len(argument_list) / (int(config['nproc']) * 2))

    with mp.Pool(initializer=initt, initargs=(terminating,),
                 processes=config['nproc']) as pool:
        try:
            if verbose:
                info("Starting parallel download.", init_new_line=True)

            [_ for _ in pool.imap(do_download, argument_list,
                                  chunksize=chunksize if chunksize else 1)]
        except Exception as e:
            error(str(e), init_new_line=True)
            error('Download failed', init_new_line=True, exit=True)
            raise


if __name__ == '__main__':
    t0 = time.time()

    args = read_params()
    check_params(args, verbose=args.verbose)

    config = read_configs(args.config_file, verbose=args.verbose)
    config = check_configs(config)

    download(config['download'], verbose=config['verbose'])

    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
