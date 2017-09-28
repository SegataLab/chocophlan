#!/usr/bin/env python3


import os
import sys
import argparse as ap
import configparser as cp


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
                   help="Prints more stuff")

    return p.parse_args()


def check_params(args):
    pass


# AVAILABLE OPTIONS:


if __name__ == '__main__':
    args = read_params()
    check_params(args, verbose=args.verbose)
    config = cp.ConfigParser()

    if (os.path.isfile(args.output)) and args.overwrite and args.verbose:
        info('Output file "{}" will be overwritten\n'.format(args.output))

    with open(args.output, 'w') as f:
        config.write(f)

    sys.exit(0)
