#!/usr/bin/env python3


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '03 Oct 2017'


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


def remove_file(filename, path=None, verbose=False):
    to_remove = path if path else ''
    to_remove += filename if filename else ''

    if to_remove:
        if os.path.exists(to_remove):
            if verbose:
                error('removing "{}"'.format(to_remove))

            os.remove(to_remove)
        elif verbose:
            error('file "{}" not found'.format(to_remove))


def read_params():
    p = ap.ArgumentParser(description="")

    p.add_argument('-f', '--config_file', type=str, default=None,
                   help="The configuration file to load")
    p.add_argument('--verbose', action='store_true', default=False,
                   help="Makes ChocoPhlAn verbose")
    p.add_argument('-v', '--version', action='store_true', default=False,
                   help="Prints the current ChocoPhlAn version")

    return p.parse_args()


def check_params(args, verbose=False):
    if args.version:
        info('ChocoPhlAn version {} ({})\n'.format(__version__, __date__),
             exit=True)
    # checking configuration file
    if not args.config_file:
        error('-f (or --config_file) must be specified', exit=True)
    elif not os.path.isfile(args.config_file):
        error('configuration file "{}" not found'.format(args.config_file),
              exit=True)


def check_config_params(args, verbose=False):
	# I have created a check_config_params function for write_config_file.py
	# The reason for this is that check_params checks the version, which I'm not sure the
	# write_config_file.py script really needs. Alternatively, we can leave remove the version check from
	# the check_params function here and decorate it with a version check if needed.
	pass


def read_configs(config_file, verbose=False):
    configs = {}
    config = cp.ConfigParser()
    config.read(config_file)

    if verbose:
        info('Reading configuration file "{}"\n'.format(config_file))

    for section in config.sections():  # "DEFAULT" section not included!
        configs[section.lower()] = {}

        for option in config[section]:
            configs[section.lower()][option.lower()] = config[section][option]

    return configs


def check_configs(config, verbose=False):
    for section_key, section_dic in config.items():
        for section, value in section_dic.items():
            if 'nproc' in section:
                try:
                    section_dic[section] = int(value)
                except Exception as e:
                    error(str(e), init_new_line=True)
                    error('nproc is not an int!\n    {}'
                          .format('{}: {}\n    '.join(zip(['section_key', 'section', 'value'],
                                                          [section_key, section, value]))),
                          init_new_line=True, exit=True)
            elif 'verbose' in section:
                try:
                    section_dic[section] = bool(value)
                except Exception as e:
                    error(str(e), init_new_line=True)
                    error('verbose is not a bool!\n    {}'
                          .format('{}: {}\n    '.join(zip(['section_key', 'section', 'value'],
                                                          [section_key, section, value]))),
                          init_new_line=True, exit=True)
            else:
                section_dic[section] = trim_trailing_slashes(value)

    return config


def trim_trailing_slashes(input_string):
    try:
        return input_string.strip().rstrip("/")
    except Exception as e:
        error(str(e), init_new_line=True)
        error('Supplied object not a string\n    {}'.format(input_string),
              init_new_line=True, exit=True)
        raise
