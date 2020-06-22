library(tidyverse)
library(taigr)

source(here::here('src', 'expt_meta_data.R'))

out_dir <- '~/CPDS/mix_seq_figshare'

all_CL_features <- read_rds(file.path(here::here('data', 'all_CL_features.rds')))
write_rds(all_CL_features, file.path(out_dir, 'all_CL_features.rds'))

for (cur_expt in sc_expts) {
  print(paste0('Processing ', cur_expt$expt_name))
  cur_dir <- file.path(out_dir, cur_expt$expt_name)
  cur_dir <- str_replace(cur_dir, 'expt10', 'expt4') #rename 10 to 4 for consistency with table S2
  if (!dir.exists(cur_dir)) {
    dir.create(cur_dir)
  }
  
  #load data files from taiga
  cur_expt_string <- cdsc:::get_permaname_from_expt_string(cur_expt$taiga_name)
  
  classifications <- load.from.taiga(data.name = cur_expt_string, data.version = cur_expt$taiga_version, data.file = 'classifications')
  genes <- load.from.taiga(data.name = cur_expt_string, data.version = cur_expt$taiga_version, data.file = 'genes')
  genes <- rbind(colnames(genes), genes) %>%
    magrittr::set_colnames(c('Ensembl_ID', 'Gene_Symbol'))
  barcodes <- load.from.taiga(data.name = cur_expt_string, data.version = cur_expt$taiga_version, data.file = 'barcodes')
  barcodes <- rbind(colnames(barcodes), barcodes) %>% magrittr::set_colnames('barcode')
  rMat <- Matrix::readMM(taigr::download.raw.from.taiga(data.name=cur_expt_string, data.file='matrix', data.version = cur_expt$taiga_version))
  
  #write files locally
  write_csv(classifications, file.path(cur_dir, 'classifications.csv'))
  write_tsv(genes, file.path(cur_dir, 'genes.tsv'), col_names = F)
  write_tsv(barcodes, file.path(cur_dir, 'barcodes.tsv'), col_names = F)
  Matrix::writeMM(rMat, file.path(cur_dir, 'matrix.mtx'))
}
