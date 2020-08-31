#!/usr/bin/env python3

__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)')
__date__ = '20 Feb 2019'

from _version import __CHOCOPhlAn_version__
import importlib
import pickle
from Bio.Phylo.BaseTree import Clade as BClade
if __name__ == '__main__':
    import utils
else:
    utils = importlib.import_module('src.utils')

def main(config):
    utils.info("Loading proteomes data...\n")
    proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
    utils.info("Loading NCBI taxontree...\n")
    taxontree = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))
    
    nodes_sgb = read_file('export/sgb_info.tsv')
    gca2taxon = { line[0]: line[1] for line in read_file('gca2taxon.tsv') }
    annotate_sgb(taxontree, proteomes, nodes_sgb, gca2taxon)


if __name__ == '__main__':
    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    config = config['export']

    main(config)

def read_file(fp):
    with open(fp) as handle:
        header = handle.readline()
        r = [line.strip().split('\t') for line in handle]
    return r

def annotate_sgb(taxontree, proteomes, nodes_sgb, gca2taxon):
    for line in nodes_sgb:
        ( GCA_id, SGB_id, tax_id, ncbi_taxonomy, n_reconstructed_genomes, 
          n_reference, isuSGB, lvl_extimated_taxonomy, estimated_taxonomy, 
          avg_dist_from_reference_genome, taxonomy_closest_genome, 
          family_taxonomy_16SrRNA, genus_taxonomy_16SrRNA, 
          estimated_species, estimated_species_taxid ) = line

        if estimated_taxonomy != 'NA':
            tax_id = int(tax_id)
            estimated_species_taxid = int(estimated_species_taxid)
            if tax_id in taxontree.taxid_n and estimated_species_taxid in taxontree.taxid_n:
                node = taxontree.taxid_n[tax_id]
                node_parent_species = taxontree.taxid_n[estimated_species_taxid]

                if 'SGB'+SGB_id not in taxontree.taxid_n:
                    sgb_node = BClade(clades=[], name='SGB'+SGB_id)
                    sgb_node.tax_id = 'SGB'+SGB_id
                    sgb_node.parent_tax_id = node_parent_species.tax_id
                    sgb_node.rank = 'sgb'
                    sgb_node.initially_terminal = False
                    taxontree.taxid_n['SGB'+SGB_id] = sgb_node
                    node_parent_species.clades.append(sgb_node)
                
                if node.rank != 'species':
                    if node in taxontree.taxid_n[node.parent_tax_id].clades:
                        taxontree.taxid_n[node.parent_tax_id].clades.remove(node)
                    node.parent_tax_id = 'SGB'+SGB_id
                    taxontree.taxid_n['SGB'+SGB_id].clades.append(node)
                else:
                    proteome = gca2taxon.get(GCA_id,'')
                    if proteome:
                        proteomes[proteome]['tax_id'] = 'SGB'+SGB_id
