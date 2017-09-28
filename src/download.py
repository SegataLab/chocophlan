#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '28 Sep 2017'


import sys
import argparse as ap
import time


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = ap.ArgumentParser(description="")

    p.add_argument('--verbose', action='store_true', default=False,
                   help="Makes ChocoPhlAn verbose")
    p.add_argument('-v', '--version', action='store_true', default=False,
                   help="Prints the current ChocoPhlAn version")

    return p.parse_args()


def check_params(args, verbose=False):
    if args.version:
        info('download.py version {} ({})\n'.format(__version__, __date__),
             exit=True)


if __name__ == '__main__':
    t0 = time.time()

    args = read_params()
    check_params(args, verbose=args.verbose)

    t1 = time.time()
    info('Total elapsed time {}s\n'.format(int(t1 - t0)))
    sys.exit(0)
