# 20171006: This script takes as input the downloaded refseq/uniprot files 
# and extracts information in such a manner that it can be conveniently
# accessed. 

# Current plan: In order to extract sequences belonging to taxonomic units 
# (i.e. species), we create a contig_ID to tax_ID mapping from the refseq 
# catalogue file. We then read the bacterial genome files (data/refseq/genomes/*),
#  use the contig_ID to tax_ID mapping to write contig sequences to tax_ID files.
