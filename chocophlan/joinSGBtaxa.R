#!/shares/CIBIO-Storage/CM/mir/tools/R-3.5.1/bin/R Rscript
library(magrittr)
library(dplyr)
library(readr)
library(stringr)

sgbinfo <-
  read_csv2(
    '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/SupplementaryTable8-SGBsDescription_20190128.csv'
  )

gca2taxa <-
  read_tsv(
    '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/gca2taxa_07022019.txt',
    col_names = c('GCA', 'taxid', 'taxstr', 'todel') 
  ) %>% select(-todel)

gca2taxa2 <- read_tsv(
  '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/gca2taxon_filt.txt',
  col_names = c('GCA', 'taxid', 'taxstr')
)

gca2taxa %<>% bind_rows(gca2taxa2) %>% distinct()

gca2sgb <-
  read_tsv(
    '/home/francesco/Downloads/reconstructed_and_referencegenomes_ALL.txt',
    col_names = c('genome', 'sgb')
  ) %>% filter(str_detect(genome, 'GCA'))

ncbi_taxonomy <-
  read_tsv(
    '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/ncbi_taxonnomy',
    col_names = c('name', 'estimated_species_taxid', 'taxstr')
  ) %>% select(-taxstr)

res <-
  left_join(gca2sgb, gca2taxa, by = c('genome' = 'GCA')) %>% left_join(., sgbinfo, by = c('sgb' = 'SGB ID'))

missing <- read_delim(
  '/home/francesco/Downloads/assembly_result',
  col_names = c('genome', 'taxid'),
  delim = ' '
)

res[res$genome %in% missing$genome, 'taxid'] <- missing$taxid

res %<>%
  distinct(genome, .keep_all = TRUE)

res$estimated_species <-
  str_split(res$`Estimated taxonomy`,
            pattern = fixed("|"),
            simplify = TRUE)[, 7] %>% str_replace('s__', '')

res %<>% left_join(ncbi_taxonomy, by = c('estimated_species' = 'name'))

write_tsv(
  res,
  '/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/sgb_info.tsv'
)
