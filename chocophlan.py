#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '28 Sep 2017'

from src.utils import *
from src import download

def chocophlan():
    args = read_params()
    check_params(args, verbose=args.verbose)

    config = read_configs(args.config_file, verbose=args.verbose)
    config = check_configs(configs, verbose=args.verbose)

    download.download(config['download'])


if __name__ == '__main__':
    t0 = time.time()
    chocophlan()
    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
