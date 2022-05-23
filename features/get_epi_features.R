library(RMySQL)
library(tidyverse)
library(hiAnnotator)
library(furrr)
library(gt23)
                    



if( !file.exists(file.path(.data_d,'intSites_full_ALL.rds'))) {
  dbConn  <- dbConnect(MySQL(), group='specimen_management')
  samples <- dbGetQuery(dbConn,'select * from gtsp where Trial="UPENN_CART19_ALL"')  
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 
                                    'specimen_management', 'intsites_miseq') %>%
    GenomicRanges::as.data.frame() %>% filter(refGenome == 'hg38') %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    stdIntSiteFragments(CPUs = .num_cores ) %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites(CPUs = .num_cores)
  saveRDS(intSites, file.path(.data_d,'intSites_full_ALL.rds'))
} else {
  intSites <- readRDS(file.path(.data_d,'intSites_full_ALL.rds'))
}
intSites <- intSites %>% as.data.frame() %>%
  mutate(GTPSposID=paste0(GTSP,posid  )) %>%
  group_by(GTSP) %>% mutate(nn=n()) %>%
  ungroup() %>% filter(nn>100) %>%
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

 test_sample <- intSites %>% as.data.frame() %>%
   filter(GTSP=='GTSP0567') %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)



epi_files <- list.files(.epigenetic_features_d)
names(epi_files) <- epi_files %>% str_remove(., ".rds")

## tested

all_names <- unique(intSites$GTSP)


all_epi <- lapply(epi_files,function(x){readRDS(file.path(.epigenetic_features_d, x))})
#kk_test <- getFeatureCounts(intSites, all_epi[[1]], 'test')

plan(sequential)
#plan(multisession, workers = .num_cores)
l_names <- length(all_names)
full_table2 <- lapply(all_names,function(c_gtsp){
  c_sample <- intSites %>% as.data.frame() %>%
    filter(GTSP==c_gtsp) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  c_num <- which(all_names==c_gtsp)
  print(paste0('starting work on ',c_gtsp,' ',c_num,'/',l_names,'--',format(Sys.time(), "%a %b %d %X %Y")))
  kk <- future_imap(all_epi, function(x,name){
    #print(paste0('- - - - starting with ',name))
    #epi_curr <- readRDS(file.path("epi_rds", x))
    p_kk <- getFeatureCounts(c_sample, x, name) %>%
      as.data.frame() %>%
      select(!colnames(as.data.frame(c_sample)))
    #print(paste0('--done with',name))
    return(p_kk)
  })
  #  plan(sequential)
  epi_field <- Reduce(cbind,kk)
  epi_field$GTPSposID <- c_sample$GTPSposID
  return(epi_field)
})

epi_final_data <- Reduce(rbind,full_table2)
final_final_data <- left_join(intSites %>% as.data.frame(),epi_final_data,by='GTPSposID')

saveRDS(epi_final_data,file=file.path(.features_d,'epi_data.rds'))
saveRDS(final_final_data,file=file.path(.features_d,'intSites_full_ALL_CLL_plus_epi.rds'))


