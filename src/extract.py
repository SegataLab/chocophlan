# 20171006: This script takes as input the downloaded refseq/uniprot files
# and extracts information in such a manner that it can be conveniently
# accessed.

# Current plan: In order to extract sequences belonging to taxonomic units
# (i.e. species), we create a contig_ID to tax_ID mapping from the refseq
# catalogue file. We then read the bacterial genome files (data/refseq/genomes/*),
#  use the contig_ID to tax_ID mapping to write contig sequences to tax_ID files.


__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it)')
__version__ = '0.01'
__date__ = '22 Oct 2017'

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

    def __init__(self):
        pass

    def __init__(self, nodes_dmp, tax_ids_to_names=None):
        # Go through every line of Nodes file to construct tree. tmp_nodes will
        # be a dictionary pointing from the taxid to its clade
        tmp_nodes = {}
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

    # Set clade status values to "True" for sequence data and "final" or "draft" if it appears in accessions (taxid -> name, status, accessions)
            # if clade.tax_id in accessions:
            #    clade.sequence_data = True
            #    clade.status = clade.accession['status']

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

    # Recursively go through the tree, and at each node remove references to
    # child clades pertaining to plasmid DNA.
    def remove_plasmids(self, clade=None):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if 'plasmid' not in c.name]
        for c in clade.clades:
            self.remove_plasmids(c)

    def lookup_by_taxid(self):
        taxid={}
        for clade in self.tree.find_clades():
            if clade.tax_id:
                if clade.tax_id in taxid:
                    raise ValueError("Duplicate key: %s" % clade.tax_id)
                taxid[clade.tax_id] = clade
        return taxid

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
            if not clade.rank in reduced_tax_levels:
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
            cur_lev = reduced_tax_levels.index(clade.rank) if lev > 0 else 0

            jumps, to_add = [], ""
            for i, c in enumerate(clade.clades):
                if c.rank == 'taxon':
                    continue
                next_lev = reduced_tax_levels.index(c.rank)
                if next_lev == cur_lev + 1 or c.rank == 'superkingdom':
                    add_internal_missing_levels(c)
                    continue

                for i, l in enumerate(reduced_tax_levels[:-1]):
                    if clade.rank == l and c.rank != reduced_tax_levels[i + 1]:
                        jumps.append(c)
                        to_add = reduced_tax_levels[i + 1]
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
        remove_noranks(self.reduced_tree.root)
        #add_noranks( self.reduced_tree.root)
        utils.info('Adding taxon names to the taxonomyn\n')
        add_taxa(self.reduced_tree.root)
        utils.info('Adding nternal missing taxonomic levels\n')
        add_internal_missing_levels(self.reduced_tree.root, lev=-1)
        utils.info('Removing duplicated taxa\n')
        reduce_double_taxa(self.reduced_tree.root)

    # def save( self, out_file_name ):
    #    self.tree = self.tree.as_phyloxml()
    #    Phylo.write( self.tree, out_file_name, "phyloxml")


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
    catalogue_dir = config['download_base_dir'] + config['relpath_taxonomic_catalogue']
    genomes_dir = config['download_base_dir'] + config['relpath_bacterial_genomes']


    taxid_contigid_dict = map_taxid_to_contigid(catalogue_dir)
    pickle.dump(taxid_contigid_dict,
                open(config['download_base_dir'] + config['relpath_pickle_taxid_contigid'], "wb" ))

    contigid_to_filename = map_contigid_to_filename(genomes_dir)
    pickle.dump(contigid_to_filename,
                open(config['download_base_dir'] + config['relpath_pickle_contigid_filename'], "wb" ))

    testrun = extract_contigs_from_taxid(515619, taxid_contigid_dict, contigid_to_filename)

    names_buf, nodes_buf = read_taxdump(taxdump_dir)

    utils.info("Starting extraction of taxonomic names\n")
    names = Names(names_buf)
    utils.info("Finishing extraction of taxonomic names\n")

    utils.info("Creating taxonomy files from nodes.dmp\n")
    tax_tree = Nodes(nodes_buf, names.get_tax_ids_to_names())
    utils.info("Finished\n")

    utils.info("Refining tree\n")
    tax_tree.get_tree_with_reduced_taxonomy()
    utils.info('Finished postprocessing the taxonomy\n')

    utils.info("Pickling tree\n")
    pickle.dump(tax_tree, 
                open(config['download_base_dir'] + config['relpath_pickle_taxontree'], "wb" ))
    utils.info("Finished\n")

def map_taxid_to_contigid(catalogue_file_dir):
    utils.info("Reading and processing catalogue file\n")

    with gzip.open(catalogue_file_dir, 'r') as f:
        catalogue = f.readlines()

    # We need to use the decode method on the file object not only here, but also in the ported RepoPhlAn code
    # (took me 1 or 2 hours to figure this out). Otherwise, it will read in the string in byte format.
    # If you know a more elegant solution to this problem, please lemme know.
    catalogue = [x.decode('utf-8').strip().split("\t") for x in catalogue]
    utils.info("Finished reading catalogue file\n")
    taxids_contigids = [(x[0], x[2]) for x in catalogue]
    taxid_contigid_dict = defaultdict(list)
    utils.info("Creating taxid to contigid mapping\n")

    for key, value in [x for x in taxids_contigids]:
        taxid_contigid_dict[key].append(value)

    utils.info("Finished taxid to contigid mapping\n")

    return taxid_contigid_dict


def extract_contigs_from_taxid(taxid, taxid_contigid_mapping,
                               contigid_filename_mapping=None):
    # For now, this function simply searches ALL files. We can later refine this by providing the file names using the
    # files_to_search parameter. Once we created the contigid_to_filename mapping, we can easily query contigids for their
    # file membership and pass the results to this function.

    contigs_to_extract = taxid_contigid_mapping[taxid]

    if contigid_filename_mapping is None:
        # Search all files, as contigid_to_filename mapping hasn't been provided. Meant for sequenctial (instead of parallel)
        # computation.

        # get list of all filename in genome file directory
        # read them sequentially using BioPython
        # check whether fasta header contains contig id
        # if yes, extract fasta sequence

        # In the end, write the whole thing out.

        # For now, exit if contigid_filename_mapping is not provided
        utils.error("contig-ID to filename mapping not provided.")
        sys.exit()
    else:
        fasta_final = []
        for contig in contigs_to_extract:
            filename = contigid_filename_mapping[contig]

            with gzip.open(genomes_dir + "/" + filename, 'rt') as fasta_file:
                fasta_final.append([record for record in fasta_file if record.id == contig])


def map_contigid_to_filename(genomes_dir, verbose=False):
    # To Francesco:
    # I have not written this function yet.
    # I think we can later add this since the extract_contigs_from_taxid function has a parameter called "files_to_search".

    contigid_to_filename = defaultdict(list)
    utils.info("Building contig-ID to genome file name mapping\n")

    for genome_file in os.listdir(genomes_dir):
        with gzip.open(genomes_dir + "/" + genome_file, 'rt') as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                    contigid_to_filename[record.id] = genome_file

    utils.info("Finished building contig-ID to genome file name mapping\n")
    return contigid_to_filename


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
