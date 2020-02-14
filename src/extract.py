#! /usr/bin/env python3
# 20171006: This script takes as input the downloaded refseq/uniprot files
# and extracts information in such a manner that it can be conveniently
# accessed.

# Current plan: In order to extract sequences belonging to taxonomic units
# (i.e. species), we create a contig_ID to tax_ID mapping from the refseq
# catalogue file. We then read the bacterial genome files (data/refseq/genomes/*),
#  use the contig_ID to tax_ID mapping to write contig sequences to tax_ID files.

__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it),'
              'Fabio Cumbo (fabio.cumbo@unitn.it),'
              'Nicola Segata (nicola.segata@unitn.it)')

from _version import __CHOCOPhlAn_version__
__date__ = '25 Mar 2019'

if __name__ == '__main__':
    import utils   
else:
    import src.utils as utils
import time
import sys
import os
import tarfile
import re
import copy
import gzip
import pickle
import glob
import itertools
from collections import defaultdict
from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree as BTree
from Bio.Phylo.BaseTree import Clade as BClade


class Nodes:
    #
    # Format of nodes.dmp from RefSeq documentation
    #
    # ---------
    #
    # This file represents taxonomy nodes. The description for each node includes
    # the following fields:
    #
    #   tax_id                  -- node id in GenBank taxonomy database
    #   parent tax_id               -- parent node id in GenBank taxonomy database
    #   rank                    -- rank of this node (superkingdom, kingdom, ...)
    #   embl code               -- locus-name prefix; not unique
    #   division id             -- see division.dmp file
    #   inherited div flag  (1 or 0)        -- 1 if node inherits division from parent
    #   genetic code id             -- see gencode.dmp file
    #   inherited GC  flag  (1 or 0)        -- 1 if node inherits genetic code from parent
    #   mitochondrial genetic code id       -- see gencode.dmp file
    #   inherited MGC flag  (1 or 0)        -- 1 if node inherits mitochondrial gencode from parent
    #   GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
    #   hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
    #   comments                -- free-text comments and citations
    #
    reduced_tax_levels = ['superkingdom', 'phylum',
                      'class', 'order', 'family', 'genus', 'species']
    def __init__(self):
        pass

    def __init__(self, nodes_dmp, tax_ids_to_names=None):
        tmp_nodes = {}
        # Go through every line of Nodes file to construct tree. tmp_nodes will
        # be a dictionary pointing from the taxid to its clade
        # with open( nodes_dmp_file ) as inpf:
        for line in nodes_dmp:
            (tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div_flag,
             genetic_code_id, inherited_GC_flag, mitochondrial_genetic_code, inherited_MGC_flag,
             GenBank_hidden_flag, hidden_subtree_root_flag, comments) = line[::2]

    # For every entry in Nodes (every location in the tree) create clade containing the scientific name and pointer to the parent node.
    # Specify the rank of the clade and the taxonomic ID of the root.
            name = (tax_ids_to_names[int(tax_id)]
                    if tax_ids_to_names else None)

            clade = BClade(clades=[], name=name)
            clade.parent_tax_id = int(parent_tax_id)
            clade.rank = re.sub(r'\W+', '', rank).strip("_")
            clade.tax_id = int(tax_id)
            clade.initially_terminal = False
            #clade.accession = accessions[clade.tax_id] if clade.tax_id in accessions else []
            tmp_nodes[clade.tax_id] = clade

            # can add any other info in node.dmp

    # Build the tree using all the clades (iterate through clades using
    # tmp_nodes)
        self.tree = BTree()
        for node in tmp_nodes.values():
            # node = parent is the trick from NCBI to identify the root
            if node.tax_id == node.parent_tax_id:
                self.tree.root = node
                continue
            parent = tmp_nodes[node.parent_tax_id]
            parent.clades.append(node)
           
        self.taxid_n = tmp_nodes
        self.leaves_taxids = []
        self.num_nodes = 0


        # Determine initial leaves of the tree. This function is called once after loading of the tree and should NOT be called at any time later, as
        # many logical concepts of functions here are based on this assumption.
        self.determine_initial_leaves()
        self.get_leave_ids()
        self.get_nr_nodes()

    # Recursively goes through all clades in the tree. Each clade root gets
    # list of all accessions in the clade.
    def add_internal_accessions(self, clade=None):
        if not clade:
            clade = self.tree.root

        clade.all_accessions = [] + \
            ([clade.accession] if clade.accession else [])

        for child in clade.clades:
            clade.all_accessions += self.add_internal_accessions(child)
        return clade.all_accessions

    # Recursively go through tree, remove references to clades that have no
    # accession information in any of their nodes.
    def remove_subtrees_without_accessions(self, clade=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if len(c.all_accessions)]
        for c in clade.clades:
            self.remove_subtrees_without_accessions(c)
    
    def remove_subtree_taxonomy_by_level(self, rank, names, clade=None):
        if not clade:
            clade = self.tree.root
        if clade.rank == rank:
            clade.clades = [c for c in clade.clades if clade.name in names]
        for c in clade.clades:
            self.remove_subtree_taxonomy_by_level(rank,names,c)
    
    def remove_leaves_without_proteomes(self,clade=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if not (c.initially_terminal and not hasattr(c,'proteome'))]
        for c in clade.clades:
            self.remove_leaves_without_proteomes(c)

    # Recursively go through the tree, and at each node remove references to
    # child clades pertaining to plasmid DNA.
    def remove_plasmids(self, clade=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if 'plasmid' not in c.name]
        for c in clade.clades:
            self.remove_plasmids(c)

    def lookup_by_rank(self):
        rankid={}
        for clade in self.tree.find_clades():
            if clade.rank:
                if clade.rank not in rankid:
                    rankid[clade.rank] = []
                rankid[clade.rank].append(clade)
        return rankid

    def determine_initial_leaves(self, clade=None):
        if not clade:
            clade = self.tree.root
        if clade.is_terminal():
            clade.initially_terminal=True
        for c in clade.clades:
            self.determine_initial_leaves(c)

    # Recursively go through the tree and remove each node that is initially terminal and whose taxid is not in the list of taxids to keep.
    # Beware that this function will only "spare" a taxon if it is initially terminal. This might be problematic in later usecases, since I'm
    # not sure whether we can garantuee that all species taxids (corresponding to proteomes) are actually ALL leave nodes in the initial tree.
    def remove_leaves(self, clade=None, taxids_to_keep=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if (not c.initially_terminal) or (c.initially_terminal and c.tax_id in taxids_to_keep)]
        for c in clade.clades:
            self.remove_leaves(c, taxids_to_keep=taxids_to_keep)

    def get_nr_nodes(self, clade=None, num_nodes=None):
        if not clade:
            clade = self.tree.root
        if not num_nodes:
            self.num_nodes = 0
        clade.clades = [c for c in clade.clades]
        for c in clade.clades:
            self.num_nodes += 1
            self.get_nr_nodes(c, num_nodes = self.num_nodes)

    # This function removes"stub" subtrees that (can) remain after pruning unwanted leaves with the remove_leaves function.
    def remove_stub(self, clade=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if (c.is_terminal() and c.initially_terminal) or (not c.is_terminal())]
        for c in clade.clades:
            self.remove_stub(c)

    # This function removes stubs until there are none left, returning only the tree of interest.
    def remove_subtrees(self, clade=None):
        itera = 0
        while True:
            self.get_nr_nodes()
            nds_bf = self.num_nodes
            self.remove_stub()
            self.get_nr_nodes()
            nds_af = self.num_nodes
            # Pruning is complete if no nodes were removed between last and current iteration.
            if (nds_af == nds_bf): 
                print("Subtree pruning complete")
                break
            itera += 1
            print("Removing stub subtrees. Iteration ", str(itera))

    def get_leave_ids(self, clade=None, recur=False):
        if not recur:
            self.leaves_taxids = []
        if not clade:
            clade = self.tree.root
        if clade.is_terminal():
            self.leaves_taxids.append(clade.tax_id)
        for c in clade.clades:
            self.get_leave_ids(c, recur=True)

    def go_up_to_species(self, taxid, up_to_group=False):
        if taxid in self.taxid_n and taxid > 1:
            rank = self.taxid_n[taxid].rank
            rank_index = self.reduced_tax_levels.index(rank) if rank in self.reduced_tax_levels else float('infinity')
            if rank_index > self.reduced_tax_levels.index('species'):
                if up_to_group and rank == 'speciesgroup':
                    return taxid
                else:
                    father = self.taxid_n[taxid].parent_tax_id
                    return self.go_up_to_species(father, up_to_group)
            elif rank_index == self.reduced_tax_levels.index('species'):
                if up_to_group:
                    father = self.taxid_n[taxid].parent_tax_id
                    if self.taxid_n[father].rank == 'speciessubgroup':
                        return self.go_up_to_species(father, up_to_group)
                    if self.taxid_n[father].rank == 'speciesgroup':
                        return father
                return taxid

    def get_child_proteomes(self, clade):
        if clade.initially_terminal and hasattr(clade,'proteomes'):
            return copy.deepcopy(clade.proteomes)
        pp = copy.deepcopy(clade.proteomes) if hasattr(clade,'proteomes') else set()
        for c in clade.clades:
            pp.update(self.get_child_proteomes(c))
        return pp
     
    def print_full_taxonomy(self, tax_id, include_groups=False):
        ranks2code = {'superkingdom': 'k', 'phylum': 'p', 'class': 'c',
                      'order': 'o', 'family': 'f', 'genus': 'g', 'species': 's','speciesgroup': 's', 'taxon': 't', 'sgb' : 't'}
        order = ('k', 'p', 'c', 'o', 'f', 'g', 's', 't')

        # path = [p for p in self.tree.root.get_path(self.taxid_n[tax_id]) if p.rank in ranks2code or (p.rank=='norank')]
        parent_tax_id = tax_id
        path = []
        if self.taxid_n[tax_id].rank == 'speciesgroup':
            path.append(self.taxid_n[tax_id])
            parent_tax_id = self.taxid_n[tax_id].parent_tax_id

        while(parent_tax_id != 1):
            curr_tax = self.taxid_n[parent_tax_id]
            if curr_tax.rank in ranks2code or (curr_tax.rank=='norank'):
                if curr_tax.rank == 'speciesgroup' and not include_groups:
                    pass
                else:
                    path.append(curr_tax)
            parent_tax_id = curr_tax.parent_tax_id

        path.reverse()
        if path[0].name == 'cellular_organisms':
            _ = path.pop(0)
        taxa_str, taxa_ids = [], []

        hasSpecies = any([True if p.rank == 'species' else False for p in path])
        hasTaxon = any([True if (p.rank == 'norank' or p.rank == 'taxon') and p.initially_terminal else False for p in path])

        if hasSpecies and hasTaxon:
            i=-1
            isSpecies=False
            while not isSpecies:
                if not path[i].initially_terminal and path[i].rank=='norank':
                    path.remove(path[i])
                if path[i].rank=='species':
                    isSpecies=True
                    path[-1].rank = 'taxon'
                i-=1
        path = [p for p in path if p.rank != 'norank']
        taxa_str = ['{}__{}'.format(ranks2code[path[x].rank], path[x].name) for x in range(len(path)) if path[x].rank != 'norank']
        taxa_ids = ['{}__{}'.format(ranks2code[path[x].rank], path[x].tax_id) for x in range(len(path)) if path[x].rank != 'norank']

        for x in range(len(order)-1) if not hasTaxon else range(len(order)):
            if x < len(taxa_str):
                if not taxa_str[x].startswith(order[x]):
                    end_lvl = order.index(taxa_str[x][0])
                    missing_levels_str = ['{}__{}_unclassified'.format(order[i], taxa_str[x-1][3:]) for i in range(x, end_lvl)]
                    missing_levels_ids = ['{}__'.format(order[i]) for i in range(x, end_lvl)]
                    for i in range(len(missing_levels_str)):
                        taxa_str.insert(x+i, missing_levels_str[i])
                        taxa_ids.insert(x+i, missing_levels_ids[i])
        
        return ('|'.join(taxa_str), '|'.join([t.split('__')[1] for t in taxa_ids]), )
        
    def print_tree(self, out_file_name, reduced=False):

        tree = self.reduced_tree if reduced else self.tree

        #to_print = tree.find_clades({"sequence_data": True})

        ranks2code = {'superkingdom': 'k', 'phylum': 'p', 'class': 'c',
                      'order': 'o', 'family': 'f', 'genus': 'g', 'species': 's', 'taxon': 't'}

        def red_rank(rank):
            if reduced and rank in ranks2code:
                return ranks2code[rank]
            return rank

        with open(out_file_name, "w") as outf:

            def trac_print_t(clade, names=None):
                if names is None:
                    if clade.name == 'root':
                        names = ""
                    else:
                        names = red_rank(clade.rank) + '__' + clade.name
                else:
                    names += ("|" if names else "") + \
                        red_rank(clade.rank) + '__' + clade.name

                # if clade.is_terminal():
                if clade.tax_id is not None and clade.name != 'root':
                    outf.write("\t".join([clade.name,
                                          # t.accession['status'],
                                          #",".join(t.accession['gen_seqs']),
                                          str(clade.tax_id),
                                          # t.accession['code'],
                                          #",".join(t.accession['accession']),
                                          # str(t.accession['len']),
                                          names
                                          ]) + "\n")

                if not clade.is_terminal():
                    for c in clade.clades:
                        trac_print_t(c, names)
            trac_print_t(tree.root)

            """
            for t in tree.get_terminals():
                tax = "|".join([red_rank(p.rank)+'__'+p.name for p in tree.get_path( t )])
                outf.write("\t".join( [ t.name,
                                        #t.accession['status'],
                                        #",".join(t.accession['gen_seqs']),
                                        str(t.tax_id),
                                        #t.accession['code'],
                                        #",".join(t.accession['accession']),
                                        #str(t.accession['len']),
                                        tax
                                        ]    )+"\n")
            """

    def get_tree_with_reduced_taxonomy(self, superkingdom="Bacteria"):
        reduced_tax_levels = ['superkingdom', 'phylum',
                              'class', 'order', 'family', 'genus', 'species']
        self.reduced_tree = copy.deepcopy(self.tree)

        def add_noranks(clade):
            if not clade.rank in self.reduced_tax_levels:
                clade.rank = 'norank'
            for c in clade.clades:
                add_noranks(c)

        def remove_noranks(clade):

            run = True

            while run:
                run = False
                new_clades = []
                for c in clade.clades:
                    # if len(c.clades) and c.rank not in reduced_tax_levels:
                    # if not hasattr(c,"sequence_data") and c.rank not in
                    # reduced_tax_levels:
                    if len(c.clades) and c.rank not in reduced_tax_levels:
                        run = True
                        c.rank = "norank"
                        new_clades += c.clades
                    else:
                        new_clades.append(c)
                    # if hasattr(c,"sequence_data") and c.rank not in
                    # reduced_tax_levels:
                    if c.rank not in reduced_tax_levels:
                        c.rank = "norank"
                clade.clades = new_clades
            for c in clade.clades:
                if len(c.clades):
                    remove_noranks(c)

        def add_taxa(clade):
            # if clade.rank == "norank" and hasattr(clade,"sequence_data"):
            if clade.rank == "norank" and clade.is_terminal():
                clade.rank = "taxon"

            """
            if not len(clade.clades) and clade.accession:
                if clade.rank == 'species':
                    newclade = copy.deepcopy( clade )
                    clade.accession = []
                    clade.sequence_data = False
                    newclade.rank = "taxon"
                    clade.clades = [newclade]
            """

            for c in clade.clades:
                add_taxa(c)

        def add_internal_missing_levels(clade, lev=1):
            if clade.rank == "taxon":
                return
            cur_lev = self.reduced_tax_levels.index(clade.rank) if lev > 0 else 0

            jumps, to_add = [], ""
            for i, c in enumerate(clade.clades):
                if c.rank == 'taxon':
                    continue
                next_lev = self.reduced_tax_levels.index(c.rank)
                if next_lev == cur_lev + 1 or c.rank == 'superkingdom':
                    add_internal_missing_levels(c)
                    continue

                for i, l in enumerate(self.reduced_tax_levels[:-1]):
                    if clade.rank == l and c.rank != self.reduced_tax_levels[i + 1]:
                        jumps.append(c)
                        to_add = self.reduced_tax_levels[i + 1]
            if jumps:
                children_ok = [c for c in clade.clades if c not in jumps]
                newclade = copy.deepcopy(clade)
                newclade.clades = jumps
                clade.clades = [newclade] + children_ok
                newclade.rank = to_add
                newclade.name = clade.name if "_noname" in clade.name else clade.name + "_noname"
                add_internal_missing_levels(newclade)

        def reduce_double_taxa(clade):
            if clade.rank == 'species' and len(clade.clades):
                torem = []
                for c in clade.clades:
                    if c.rank == 'taxon':
                        if len(c.clades):
                            clade.clades += c.clades
                            torem.append(c)
                clade.clades = [c for c in clade.clades if c not in torem]
                return
            for c in clade.clades:
                reduce_double_taxa(c)

        #add_noranks( self.reduced_tree.root )

        utils.info('Removing noranks from the taxonomy\n')
        remove_noranks(self.tree.root)
        #add_noranks( self.reduced_tree.root)
        utils.info('Adding taxon names to the taxonomy\n')
        add_taxa(self.tree.root)
        # utils.info('Adding internal missing taxonomic levels\n')
        # add_internal_missing_levels(self.reduced_tree.root, lev=-1)
        utils.info('Removing duplicated taxa\n')
        reduce_double_taxa(self.tree.root)

class Names:
    #
    # Format of names.dmp from RefSeq documentation
    #
    # ---------
    #
    # Taxonomy names file has these fields:
    #
    #   tax_id                  -- the id of node associated with this name
    #   name_txt                -- name itself
    #   unique name             -- the unique variant of this name if name not unique
    #   name class              -- (synonym, common name, ...)
    #

    def __init__(self, names_dmp):
        # Read from file names.dmp, get information in every field
        self.tax_ids_to_names = {}
        for line in names_dmp:
            tax_id, name_txt, unique, name_class = line[::2]

            # extracting scientific names only (at least for now) which are unique!
        # tax_ids_to_names relates taxid to the sole scientific name of the
        # organism
            if name_class == "scientific name":
                name = re.sub(r'\W+', '_', name_txt).strip("_")
                self.tax_ids_to_names[int(tax_id)] = name

    def get_tax_ids_to_names(self):
        return self.tax_ids_to_names

def do_extraction(config, verbose=False):
    taxdump_dir = config['download_base_dir'] + config['relpath_taxdump']

    names_buf, nodes_buf = read_taxdump(taxdump_dir)

    utils.info("Starting extraction of taxonomic names... ")
    names = Names(names_buf)
    utils.info("Finished\n")

    utils.info("Creating taxonomy files from nodes.dmp... ")
    tax_tree = Nodes(nodes_buf, names.get_tax_ids_to_names())
    utils.info("Finished\n")

    tax_tree.tree.root.clades = tax_tree.tree.root.clades[0:2] + tax_tree.tree.root.clades[4::]

    utils.info("Refining tree... ")
    # tax_tree.get_tree_with_reduced_taxonomy()
    r = re.compile( r"(.*(C|c)andidat(e|us)_.*)|"
                    r"(.*_sp(_.*|$).*)|"
                    r"((.*_|^)(b|B)acterium(_.*|))|"
                    r"(.*(eury|)archaeo(n_|te|n$).*)|"
                    r"(.*(endo|)symbiont.*)|"
                    r"(.*genomosp_.*)|"
                    r"(.*unidentified.*)|"
                    r"(.*_bacteria_.*)|"
                    r"(.*_taxon_.*)|"
                    r"(.*_et_al_.*)|"
                    r"(.*_and_.*)|"
                    r"(.*(cyano|proteo|actino)bacterium_.*)")
    merged_low_quality = [tax_tree.taxid_n[x].tax_id for x in tax_tree.taxid_n if r.match(tax_tree.taxid_n[x].name)]
    [tax_tree.taxid_n[x].__setattr__('is_low_quality', True) for x in merged_low_quality]
    [tax_tree.taxid_n[x].__setattr__('is_low_quality', False) for x in set(tax_tree.taxid_n.keys()).difference(merged_low_quality)]
    utils.info('Finished postprocessing the taxonomy\n')

    utils.info("Pickling tree... ")
    os.makedirs(config['download_base_dir'] + '/pickled', exist_ok=True)
    pickle.dump(tax_tree, 
                open(config['download_base_dir'] + config['relpath_pickle_taxontree'], "wb" ))
    utils.info("Finished\n")


# Additional idea: we can also save the position of each contig within each file like so:
# Final_dir = {"contig_id_1": ("file_name_abc", 252624)}, where the second entry of the tuple is the number of the contig within the file.
# Can be easily added while reading and creating the contigid to filename mapping. This would later safe us a LOT of time, right?
def read_taxdump(taxdump_dir, verbose=False):
    utils.info('Reading the NCBI taxdump file from {}\n'.format(taxdump_dir))
    try:
        tarf = None

        if os.path.isfile(taxdump_dir):
            tarf = tarfile.open(taxdump_dir, "r:gz")
        else:
            utils.error('{} does not exists. Exiting...'.format(taxdump_dir))
            sys.exit()
        for m in tarf.getmembers():
            if m.name == "names.dmp":
                names_buf = (l.decode("utf-8").strip().split('\t')
                             for l in tarf.extractfile(m).readlines())
            if m.name == "nodes.dmp":
                nodes_buf = (l.decode("utf-8").strip().split('\t')
                             for l in tarf.extractfile(m).readlines())
    except Exception as e:
        utils.error("Error in extracting or reading {}: {}"
                    .format(taxdump_dir, e))
        sys.exit()

    utils.info('names.dmp and nodes.dmp successfully read\n')

    return (names_buf, nodes_buf)

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    do_extraction(config['extract'], verbose=config['extract']['verbose'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)
