library(RMySQL)
library(tidyverse)
library(hiAnnotator)
#library(furrr)

source('../.Rprofile')



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
intSites <- intSites %>% as.data.frame() %>%
  mutate(GTPSposID=paste0(GTSP,posid  )) %>%
  group_by(GTSP) %>% mutate(nn=n()) %>%
  ungroup() %>% filter(nn>100) %>%
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

test_sample <- intSites %>% as.data.frame() %>%
  filter(GTSP=='GTSP0567') %>% #GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M")))
## windows
window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_CpG_counts <- c("1k"=1e3, "10k"=1e4)
window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_GC <- c("100"=100, "1k"=1000, "10k"=1e4, "100k"=1e5, "1M"=1e6)
#window_size_GC <- c("100"=100, "1k"=1000)
window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_epi <- c("10k"=1e4)

#
refGenes <- readRDS(file.path(.genetic_features_d, "hg38.refSeq.rds"))
refGenes <- refGenes[seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))]
#
genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]
#
CpG_data <- cpg <- getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38") %>% 
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

#getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38")


CpG_islands <- GenomicRanges::GRanges(
  seqnames = CpG_data$chrom,
  ranges = IRanges::IRanges(
    start = CpG_data$chromStart, end = CpG_data$chromEnd
  ),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence)
)

mcols(CpG_islands) <- CpG_data
#
DNaseI_data <- getUCSCtable("wgEncodeRegDnaseClustered", "DNase Clusters", freeze = "hg38") %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

DNaseI <- GenomicRanges::GRanges(
  seqnames = DNaseI_data$chrom,
  ranges = IRanges::IRanges(start = DNaseI_data$chromStart, end = DNaseI_data$chromEnd),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence))

mcols(DNaseI) <- DNaseI_data
#
from_counts_to_density <- function(sites, column_prefix, window_size) {
  metadata <- mcols(sites)
  sapply(seq(window_size), function(i) {
    val <- window_size[i]
    name <- names(window_size)[i]
    column_name <- paste0(column_prefix, ".", name)
    metadata[[column_name]] <<- metadata[[column_name]]/val
  })
  mcols(sites) <- metadata
  sites
}
####

####
#all_names <- unique(intSites$GTSP)[100:104]
all_names <- unique(intSites$GTSP)
#plan(sequential)
#plan(multisession, workers = .num_cores)
l_names <- length(all_names)

hg38_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
hg38_seqlev <- seqlevels(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

full_table2 <- map(all_names, function(c_gtsp) {
    c_sample <- intSites %>% 
      as.data.frame() %>%
      filter(GTSP==c_gtsp) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    seqlevels(c_sample) <- hg38_seqlev
    seqinfo(c_sample) <- hg38_seqinfo
    
    c_num <- which(all_names==c_gtsp)
    print(paste0('starting work on ',c_gtsp,' ',c_num,'/',l_names,'--',format(Sys.time(), "%a %b %d %X %Y")))
  
    c_tab <- c_sample %>% #as.data.frame() %>%
      #GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo = hg38_seqinfo) %>% 
      hiAnnotator::getFeatureCounts(refGenes, "refSeq_counts", width = window_size_refSeq) %>%
      GCcontent::getGCpercentage("GC", window_size_GC, genome_sequence) %>%
      hiAnnotator::getFeatureCounts(CpG_islands, "CpG_counts", width = window_size_CpG_counts) %>%
      hiAnnotator::getFeatureCounts(CpG_islands, "CpG_density", width = window_size_CpG_density) %>%
      from_counts_to_density("CpG_density", window_size_CpG_density) %>%
      hiAnnotator::getFeatureCounts(DNaseI, "DNaseI_count", width = window_size_DNaseI) %>%
      as.data.frame()
    return(c_tab)
})

gen_final_data <- Reduce(rbind,full_table2) %>% select(-c(seqnames:nearestlymphomaFeatureStrand,nn))
final_final_data <- left_join(intSites %>% as.data.frame(),gen_final_data,by='GTPSposID')

saveRDS(gen_final_data,file=file.path(.features_d,'gen_data.rds'))
saveRDS(final_final_data,file=file.path(.features_d,'intSites_full_ALL_CLL_plus_gen.rds'))







