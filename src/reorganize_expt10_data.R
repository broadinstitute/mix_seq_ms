library(tidyverse)
library(taigr)
library(celllinemapr)
library(here)

source(here('src', 'MixSeq_helpers.R'))

data_dir <- '/Volumes/GoogleDrive/My Drive/Project Apollo/Apollo drugs /scRNAseq_010'

#get mapping between channels and treatments by cell line
pool_map <- read_csv(file.path(data_dir, 'Apollo_010_channel_treatment_metadata.csv'))
pool_map %<>% mutate(Treatment = plyr::revalue(Treatment, replace = c(AZD5991 = 'AZD5591'))) #fix typo in name of the MCL1 inh

#merge in list of cell lines in each pool
PRISM_meta <- read_csv('~/CPDS/data/Prism/CS1.4 OneshotPools_Metadata.csv')
pool_map %<>% left_join(PRISM_meta %>% 
                          mutate(PRISM_Pool_ID = str_match(Pool_ID, '^([:alnum:]+)\\.OS')[,2]) %>% 
                          dplyr::select(PRISM_Pool_ID, CCLE_Name))

#make sure each cell line has one entry for each drug
unique_cell_lines <- unique(pool_map$CCLE_Name)
length(unique_cell_lines)
for (cl in unique_cell_lines) {
  cur_df <- pool_map %>% filter(CCLE_Name == cl)
  stopifnot(nrow(cur_df) == 10)
  stopifnot(sum(duplicated(cur_df$Treatment)) == 0)
}

pool_map %<>% mutate(Treatment = str_replace_all(Treatment, '-', ''))
unique_drugs <- unique(pool_map$Treatment)
ddirs <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE) %>% 
  grep('^BP', ., value = T) #get list of data directories

#loop over drugs 
for (cur_drug in unique_drugs) {
  cur_dirs <- grep(cur_drug, ddirs, value = TRUE)
  stopifnot(length(cur_dirs) == 2)
  
  cur_pool_map <- pool_map %>% filter(Treatment == cur_drug)
  
  #from first channel
  cur_dir <- cur_dirs[[1]]
  cur_channel <- as.numeric(str_match(cur_dir, '^BP([0-9]+)_')[,2])
  cur_cell_lines <- pool_map %>% filter(Treatment == cur_drug, Channel == cur_channel) %>% .[['CCLE_Name']]
  
  rMat1 <- Matrix::readMM(file.path(data_dir, cur_dir, 'matrix.mtx'))
  genes1 <- read_tsv(file.path(data_dir, cur_dir, 'genes.tsv'), col_names = F) %>% 
    set_colnames(c('Ensembl_ID', 'Gene_Symbol'))
  barcodes1 <- read_tsv(file.path(data_dir, cur_dir, 'barcodes.tsv'), col_names = F) 
  classifications1 <- read_csv(file.path(data_dir, cur_dir, 'classifications.csv')) 
  stopifnot(nrow(classifications1) == nrow(barcodes1))
  # stopifnot(sum(!(classifications1$singlet_ID %in% unique_cell_lines)) == 0)
  classifications1 %<>% 
    filter(singlet_ID %in% cur_cell_lines)
  ucells <- match(classifications1$barcode, barcodes1$X1)
  rMat1 <- rMat1[, ucells]
  classifications1 %<>% 
    mutate(channel = cur_channel,
           barcode = paste0(barcode, '_ch', channel))
  barcodes1 <- classifications1$barcode
   
  
  #from second channel
  cur_dir <- cur_dirs[[2]]
  cur_channel <- as.numeric(str_match(cur_dir, '^BP([0-9]+)_')[,2])
  cur_cell_lines <- pool_map %>% filter(Treatment == cur_drug, Channel == cur_channel) %>% .[['CCLE_Name']]
  
  rMat2 <- Matrix::readMM(file.path(data_dir, cur_dir, 'matrix.mtx'))
  genes2 <- read_tsv(file.path(data_dir, cur_dir, 'genes.tsv'), col_names = F) %>% 
    set_colnames(c('Ensembl_ID', 'Gene_Symbol'))
  barcodes2 <- read_tsv(file.path(data_dir, cur_dir, 'barcodes.tsv'), col_names = F) 
  classifications2 <- read_csv(file.path(data_dir, cur_dir, 'classifications.csv')) 
  stopifnot(nrow(classifications2) == nrow(barcodes2))
  # stopifnot(sum(!(classifications2$singlet_ID %in% unique_cell_lines)) == 0)
  classifications2 %<>% 
    filter(singlet_ID %in% cur_cell_lines)
  ucells <- match(classifications2$barcode,barcodes2$X1)
  rMat2 <- rMat2[, ucells]
  classifications2 %<>% 
    mutate(channel = cur_channel,
           barcode = paste0(barcode, '_ch', channel))
  barcodes2 <- classifications2$barcode
  
  stopifnot(all.equal(genes1$Ensembl_ID, genes2$Ensembl_ID))
  
  #assemble new dataset
  barcodes <- c(barcodes1, barcodes2)
  classifications <- rbind(classifications1, classifications2)
  genes <- genes1
  rMat <- cbind(rMat1, rMat2)
  
  classifications %<>% mutate(singlet_ID = ccle.to.latest(singlet_ID)) #convert to updated cell line names (this fixes S117_SOFT_TISSUE, which has been renamed to S117_THYROID)
  
  #save to directory organized by treatment
  new_dir <- file.path(data_dir, cur_drug)
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }
  Matrix::writeMM(rMat, file.path(new_dir, 'matrix.mtx'))
  write_tsv(genes, file.path(new_dir, 'genes.tsv'), col_names = F) 
  write_lines(barcodes, file.path(new_dir, 'barcodes.tsv')) 
  write_csv(classifications, file.path(new_dir, 'classifications.csv')) 
}
