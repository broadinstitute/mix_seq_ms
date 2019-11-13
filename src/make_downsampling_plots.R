library(plyr)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)

mkae_downsampling_plots <- function(fig_dir) {

  # downsampling rates used in run_downsampling_analysis
  sampling_rates = c(0.9, 0.8, 0.7, 0.6,0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01)

  original_classifications <- readr::read_csv(here::here('data', 'classifications.csv'))
  out_of_pool_cls <- read_csv(here::here('data', 'bulk_reference_CLs.csv'))$CCLE_ID

  out_of_pool_percent <- numeric()
  mean_reads <- numeric()
  mean_SNPs <- numeric()
  error_min <- numeric()
  error_max <- numeric()
  number_of_SNPs <- numeric()
  mean_snps <- character()
  for(i in 1:length(sampling_rates)) {
    cur_result <- readr::read_csv(here::here('data', paste0('classifications_',sampling_rates[i],'.csv'))) %>% 
      as.data.frame()
    normal_cells <- which(original_classifications$cell_quality == 'normal')
    in_pool <- unique(original_classifications$singlet_ID)
    out_of_pool_percent <- c(out_of_pool_percent, length(which(!cur_result[normal_cells,'singlet_ID'] %in% in_pool))/length(normal_cells))
    
    mean_reads <- c(mean_reads, mean(cur_result$tot_reads[normal_cells], na.rm=T))
    mean_SNPs <- c(mean_SNPs, mean(cur_result$num_SNPs[normal_cells], na.rm=T))
    pt <- prop.test(length(which(!cur_result[normal_cells,'singlet_ID'] %in% in_pool)), length(normal_cells))
    
    error_min <- c(error_min, pt$conf.int[1])
    error_max <- c(error_max, pt$conf.int[2])
    
    mean_snps <- c(mean_snps, rep(round(mean(cur_result$num_SNPs[normal_cells], na.rm=T), 0), length(normal_cells)))
    number_of_SNPs <- c(number_of_SNPs, cur_result$num_SNPs[normal_cells])
    
  }
  
  ds_df <- cbind.data.frame(sampling_rates, out_of_pool_percent, mean_reads, mean_SNPs, error_min, error_max)
  in_pool_proportion <- length(unique(original_classifications$singlet_ID))/length(out_of_pool_cls)
  ds_df$error_rate  <- out_of_pool_percent/(1-in_pool_proportion)
  ds_df$error_rate_min  <- error_min/(1-in_pool_proportion)
  ds_df$error_rate_max  <- error_max/(1-in_pool_proportion)
  
  snp_dist <- cbind.data.frame(`number_of_snps`= dplyr::filter(original_classifications, cell_quality=='normal')$num_SNPs)
  
  log_scale <- c(0, round(emdbook::lseq(1, 1350, 14)))
  
  # plot distribution of number of SNP sites detected per cell
  snp_distribution <- ggplot2::ggplot(snp_dist, ggplot2::aes(number_of_snps)) +
    ggplot2::geom_histogram(fill="#666666", color='#333333') +
    ggplot2::scale_x_continuous(limits=c(8, 1350), trans='log2') +
    ggplot2::xlab("number of SNPs") +
    cdsr::theme_Publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=3), axis.title = ggplot2::element_text(size=4), 
                                      axis.text.y = ggplot2::element_text(size=3))
  
  # plot error rate as function of number of SNP sites detected
  downsampling_error_rate <- ggplot2::ggplot(ds_df, ggplot2::aes(mean_SNPs, error_rate)) +
    ggplot2::geom_point(size=0.2) + 
    ggplot2::scale_x_continuous(limits=c(8, 1350), trans='log2') +
    ggplot2::scale_y_continuous(breaks = round(seq(0.1, 1, by = 0.1),1)) +
    ggplot2::ylab("estimated error rate") + 
    ggplot2::xlab("mean number of SNPs") +
    #geom_errorbar(aes(ymin=error_rate_min, ymax=error_rate_max)) +
    ggplot2::geom_line(size=0.1) +
    cdsr::theme_Publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=3), axis.title = ggplot2::element_text(size=4), 
                                      axis.text.y = ggplot2::element_text(size=3))
  
  ggplot2::ggsave(filename="snp_distribution.png",  plot=snp_distribution, device='png', 
                  path=fig_dir, width=1.5, height=1, units='in')
  ggplot2::ggsave(filename="downsampling_error_rate.png",  plot=downsampling_error_rate, device='png', 
                  path=fig_dir, width=1.5, height=1, units='in')

}
