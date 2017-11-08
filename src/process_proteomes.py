#!/usr/bin/env python3

__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '7 Nov 2017'

import utils
import os
import gzip
import pickle
import multiprocessing as mp
import glob


kingdom_to_process = ['Bacteria','Archaea']

def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_

def process(input, is_reference=True):
    if not terminating.is_set():
        prot = {}
        with gzip.open(input, 'rt') as idmap:
            for line in idmap:
                uniprotkbid, k, v = line.split()
                if uniprotkbid not in prot:
                    prot[uniprotkbid] = {k:[v]}
                else:
                    if k not in prot[uniprotkbid]:
                        prot[uniprotkbid][k] = [v]
                    else:
                        prot[uniprotkbid][k].append(v)
        protid, taxid = os.path.split(input)[1].split('.')[0].split('_')
        prot = {'proteome_id': protid
               ,'taxid': taxid
               ,'isreference' : is_reference
               ,'proteins': prot }
        return prot
    else:
        terminating.set()
            
def parse_reference_proteomes(config, verbose=False):
    terminating = mp.Event()
    chunksize = math.floor(len(argument_list) / (int(config['nproc']) * 2))

    for k in kingdom_to_process:
        basepath = "{}/{}/{}/*.idmapping.gz".format(config['download_base_dir'],config['relpath_reference_proteomes'], k)
        ls = glob.glob(basepath)

        chunksize = config['nproc']
        with mp.Pool(initializer=initt, initargs=(terminating,), processes=config['nproc']) as pool:
            try:
                if verbose:
                    utils.info("Starting parallel proteomes processing...\n", init_new_line=True)
                proteomes = { item['proteome_id']:item for item in [x for x in pool.imap(process, ls, chunksize = chunksize)]}
            except Exception as ex:
                utils.error(str(e), init_new_line=True)
                utils.error('Processing failed', init_new_line=True, exit=True)
                raise
        if verbose:
                utils.info("Done\n", init_new_line=True)
    pickle.dump(proteomes, open(config['download_base_dir'] + config['relpath_pickle_proteomes'], "wb" ))

if __name__ == '__main__':
    t0=time.time()
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    parse_reference_proteomes(config)
    
    t1=time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
