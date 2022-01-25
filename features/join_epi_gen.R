library(tidyverse)
library(RMySQL)
library(hiAnnotator)
#library(furrr)

gen_f <- readRDS(file.path(.features_d,'gen_data.rds'))
epi_f <- readRDS(file.path(.features_d,'epi_data.rds'))

### Fix an error on the aready computed epigenomic data
# epi_n <- epi_f %>% separate(GTPSposID,c('GTSP','posID'),'(?<=GTSP....)(?=chr.*$)',remove = FALSE) %>%
#   separate(posID,c('seqnames','start'),'[\\+\\-]') %>%
#   filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
#   select(-c(GTSP,seqnames,start))
# saveRDS(epi_n,file=file.path(.features_d,'epi_data.rds'))


if( !file.exists(file.path(.data_d,'intSites_full_ALL_CLL.rds'))) {
  dbConn  <- dbConnect(MySQL(), group='specimen_management')
  samples <- dbGetQuery(dbConn,'select * from gtsp where Trial="CART19_ALL" or Trial="CART19_CLL"')  
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 
                                    'specimen_management', 'intsites_miseq') %>%
    GenomicRanges::as.data.frame() %>% filter(refGenome == 'hg38') %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    stdIntSiteFragments(CPUs = numCores ) %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites(CPUs = numCores)
  saveRDS(intSites, file.path(.data_d,'intSites_full_ALL_CLL.rds'))
} else {
  intSites <- readRDS(file.path(.data_d,'intSites_full_ALL_CLL.rds'))
}
intSites_dd <- intSites %>% as.data.frame() %>%
  mutate(GTPSposID=paste0(GTSP,posid  )) %>%
  group_by(GTSP) %>% mutate(nn=n()) %>%
  ungroup() %>% filter(nn>100) %>%
  mutate(nn=NULL) %>%
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M")))

all_f <- left_join(gen_f,epi_f,by='GTPSposID')
intSite_features <- left_join(intSites_dd,all_f,by='GTPSposID')

saveRDS(intSite_features,file=file.path(.features_d,'intSite_feature.rds'))

gtsp_avg_epi_gen <- left_join(gen_f,epi_f,by='GTPSposID') %>%
  separate(GTPSposID,c('GTSP','posID'),'(?<=GTSP....)(?=chr.*$)') %>%
  select(-posID) %>%
  group_by(GTSP) %>%
  summarise_all(mean, na.rm = TRUE)
