# 20171006: This script takes as input the downloaded refseq/uniprot files 
# and extracts information in such a manner that it can be conveniently
# accessed. 

# Current plan: In order to extract sequences belonging to taxonomic units 
# (i.e. species), we create a contig_ID to tax_ID mapping from the refseq 
# catalogue file. We then read the bacterial genome files (data/refseq/genomes/*),
#  use the contig_ID to tax_ID mapping to write contig sequences to tax_ID files.



__author__ = ('Nicola Segata (nicola.segata@unitn.it), '
              'Francesco Beghini (@unitn.it)'
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

def extract_tax_and_contigs(config, verbose=False):

	utils.info('Reading the NCBI taxdump file from ' + config['refseq_taxdump'])
	try:
	    tarf = None
	    if os.path.isfile(config['refseq_taxdump']):
	        tarf = tarfile.open(config['refseq_taxdump'])
	    else:
	        loc = urllib2.urlopen(config['refseq_taxdump'])
	        compressedFile = StringIO.StringIO(loc.read())
	        tarf = tarfile.open(fileobj=compressedFile)
	    for m in tarf.getmembers():
	        if m.name == "names.dmp":
	            names_buf = (l.strip().split('\t')
	                         for l in tarf.extractfile(m))
	        if m.name == "nodes.dmp":
	            nodes_buf = (l.strip().split('\t')
	                         for l in tarf.extractfile(m))
	except Exception as e:
	    utils.error("Error in downloading, extracting, or reading " +
	                 config['refseq_taxdump'] + ": " + str(e))
	    sys.exit()

	utils.info(
	    'names.dmp and nodes.dmp successfully downloaded, extracted, and read')



if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)

    extract_tax_and_contigs(config['extract'], verbose=config['extract']['verbose'])

    t1 = time.time()

    utils.info('Total elapsed time {}s\n'.format(int(t1 - t0)), init_new_line=True)
    sys.exit(0)