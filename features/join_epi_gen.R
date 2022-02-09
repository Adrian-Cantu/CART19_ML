library(tidyverse)
library(RMySQL)
library(hiAnnotator)
library(vegan)
library(reldist)
#library(gintools)
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
  group_by(GTSP) %>%
  mutate(nn=n()) %>%
  ungroup() %>%
  filter(nn>100) %>%
  mutate(nn=NULL) %>%
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M")))

# intSites_dd %>%
#   filter(nn<=100) %>%
#   pull(var=nn) %>%
#   unique()


gtsp_avg_epi_gen <- left_join(gen_f,epi_f,by='GTPSposID') %>%
  separate(GTPSposID,c('GTSP','posID'),'(?<=GTSP....)(?=chr.*$)') %>%
  select(-posID) %>%
  group_by(GTSP) %>%
  summarise_all(mean, na.rm = TRUE)


#import response data

library(readxl)
CARTSite_Response <- read_excel("CARTSite.Response.BOR.01142022.xlsx", 
                     col_types = c("text", "numeric", "skip", 
                     "date", "skip", "skip", "skip", "date", 
                     "numeric", "text", "numeric", "text"))
# BOR= best overal response
colnames(CARTSite_Response) <- c('Trial','ID','inf_date','response_date','response_timepoint','BOR','BORc','note')
dbConn  <- dbConnect(MySQL(), group='specimen_management')

########
#samples <- dbGetQuery(dbConn,'select * from gtsp where Trial="CART19_ALL" or Trial="CART19_CLL"')
samples <- dbGetQuery(dbConn,'select * from gtsp')%>%
  filter(SpecimenAccNum %in% intSites$GTSP)


gtsp_to_pid <- samples %>% select(c("SpecimenAccNum","SamplePatientCode","Patient","Trial")) %>%
  mutate(SamplePatientCode = str_replace(SamplePatientCode,'^CHP','CHOP')) %>%
  separate(SamplePatientCode,c('Trial2',NA),'(?<=^(CHOP|UPCC))(?=\\d*(_|\\-))') %>%
  mutate(Trial2=ifelse(Trial2=='xxx','CHOP',Trial2)) %>%
  mutate(Patient=str_replace(Patient,'^\\D+','')) %>%
  mutate(Patient=str_replace(Patient,'\\-?[^0-9\\-].*','')) %>%
  mutate(PID=paste0(Trial2,Patient)) %>%
  mutate(GTSP=SpecimenAccNum) %>%
  select(c(GTSP,PID,Trial))

pid_to_bor <- CARTSite_Response %>%
  group_by(Trial,ID) %>%
  summarise(BORc=max(BORc), .groups = 'drop') %>%
  mutate(PID=paste0(Trial,'-',sprintf("%02d", ID))) %>%
  select(c(PID,BORc))
 
gtsp_to_bor <- left_join(gtsp_to_pid,pid_to_bor,by='PID') 

# merging features and responses class

features_responses <- left_join(gtsp_avg_epi_gen,gtsp_to_bor,by='GTSP')

# 
calculateUC50 <- function(abund){
  stopifnot(is.vector(abund) & is.numeric(abund))
  abund <- abund[order(abund)]
  accum <- sapply(1:length(abund), function(i){sum(abund[1:i])})
  length(accum[accum >= sum(abund)/2])
}


# Recreate columns used in published data.
pop_stats <- intSites_dd %>%
  mutate(timePoint=tolower(timePoint)) %>%
  group_by(GTSP,timePoint,cellType) %>%
#  mutate(within_gene=ifelse(nearest_geneDist == 0, TRUE, FALSE)) %>%
  summarise(
    "numUniqSites" = n(),
    "estAbund_avg" = mean(estAbund),
    "maxRelAbund" = max(relAbund),
    "numClones" = sum(estAbund), 
    "ShannonIndex" = diversity(estAbund),
    "GiniIndex" = gini(estAbund),
    "Chao1" = round(estimateR(estAbund, index='chao')[2], 0),
    "UC50" = calculateUC50(estAbund),
    "pctTxnUnit" = 100 * sum(inFeature, na.rm = TRUE)/n(),
    "pctSameOrt" = 100 * sum(inFeatureSameOrt, na.rm = TRUE) / 
      sum(inFeature, na.rm = TRUE),
    "pctNearTxnUn" = 100 * sum(
      inFeature == FALSE & abs(nearestFeatureDist) <= 5000, na.rm = TRUE) / 
      sum(inFeature == FALSE, na.rm = TRUE),
    "pctInOnco" = 100 * sum(abs(nearestOncoFeatureDist) <= 5000, na.rm = TRUE) / n(),
    "nearest_geneDist" = mean(nearestFeatureDist),
    "within_gene" = mean(ifelse(nearestFeatureDist == 0, 1, 0)),
    "same_ort" = mean(inFeatureSameOrt),
    .groups = 'drop'
  )



features_responses_pop <- left_join(features_responses,pop_stats,by='GTSP') %>%
  relocate(c(PID,Trial,BORc),.after=GTSP)

###
# add cluster

cluster_GTSP <- readRDS(file.path(.data_d,'cluster_GTSP.rds'))
old_f_names <- c('GTSP'          , 'timePoint', 'PID'    , 'cellType', 'estAbund_avg')
new_f_names <- c('SpecimenAccNum', 'Timepoint', 'Patient', 'CellType', 'estAbund')
features_responses_pop2 <- left_join(features_responses_pop,cluster_GTSP,by='GTSP') %>%
  rename_with(~ str_replace(.x,'Kb$','k')) %>%
  rename_with(~ new_f_names, all_of(old_f_names)) %>%
  mutate(Trial=ifelse(Trial=="Gill_CART19_18415","CART19_CLL",Trial)) %>%
  relocate(c(Patient,Trial,Timepoint,CellType,BORc),.after=SpecimenAccNum)

#save
saveRDS(features_responses_pop2,file=file.path(.features_d,'sample_features_20220209.rds'))

write.table(features_responses_pop2, file=file.path(.features_d,'sample_features_20220209.tsv'), quote=FALSE, sep='\t',row.names = FALSE)


wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "sample_features")
openxlsx::writeDataTable(wb,"sample_features",features_responses_pop2)
openxlsx::saveWorkbook(wb, file.path(.features_d,'sample_features_20220209.xlsx'), overwrite = TRUE)


#######################################
# 
# missGTSP <- features_responses_pop %>%
#   filter(is.na(BORc)) %>%
#   pull(var=GTSP)

# features_responses_pop2 %>%
#   dplyr::filter(if_any(.cols = everything(),is.na))

#missSample <- 
  
#  samples %>%
#  filter(SpecimenAccNum %in% missGTSP)

#samples_miss <- dbGetQuery(dbConn,'select * from gtsp where Trial="CART19_ALL" or Trial="CART19_CLL"')
#dbConn  <- dbConnect(MySQL(), group='specimen_management')
#arr <-paste0("SELECT * FROM gtsp WHERE SpecimenAccNum IN ('",paste(missGTSP,collapse = "','"),"')")
#samples_miss <- dbGetQuery(dbConn,arr)

  
#features_responses_pop$BORc
# library(readxl)
# John_CART <- read_excel("data/John_CART.xlsx")
# 
# topull <- 'H2BK5ac.10k'
# 
# for(topull in intersect(colnames(John_CART),colnames(features_responses_pop2))) {
# kk_j <- John_CART %>%
#   filter(!is.na(Chao1)) %>%
#   pull(var=all_of(topull),name=SpecimenAccNum)
# 
# kk_a <- features_responses_pop2 %>%
#   pull(var=all_of(topull),name=SpecimenAccNum)
# 
# blabla <- setdiff(names(kk_j),names(kk_a))
# int_blabla <- intersect(names(kk_j),names(kk_a))
# kk_aa <- kk_a[names(kk_a) %in% int_blabla]
# kk_jj <- kk_j[names(kk_j) %in% int_blabla]
# 
# print(topull)
# print(kk_aa[int_blabla][1:5])
# print(kk_jj[int_blabla][1:5])
# }
