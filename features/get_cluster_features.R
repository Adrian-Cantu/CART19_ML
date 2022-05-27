#options(stringsAsFactors = FALSE)
#library(dplyr)
library(gt23)
#library(stringr)
#library(hiAnnotator)
#library(GCcontent)
library(BSgenome)
library(gintools)
#library(RMySQL)
library(geneRxCluster)
#library(ggplot2)
library(tidyverse)
numCores <- 2
#source('utils.R')
#source('supporting_functions.R')
#source('clusters.R')
#utilsDir <- 'utils'

get_top_sites <- function(sites, percent, rank_by){
  sites <- sites[order(mcols(sites)[,rank_by], decreasing = TRUE)]
  first_sites <- sites[!duplicated(sites$posid)]
  ranks <- rank(-(mcols(first_sites)[,rank_by]), ties.method = "min")
  cutoff <- round(length(first_sites)*percent/100)
  if(cutoff == 0){ cutoff <- 1 }
  first_sites[ranks <= cutoff]$posid
}

scan_sites_abundance <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    sum(sites[hits]$estAbund)
  })
}

scan_patients <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(sites[hits]$patient))
  })
}

scan_genes <- function(clus_gr, refGenes){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, refGenes))
    paste(unique(refGenes[hits]$name2), collapse = ", ")
  })
}

scan_orientation <- function(clus_gr, sites){

  sapply(1:length(clus_gr), function(i){
   # browser()
    clus <- clus_gr[i]
    site_hits <- findOverlaps(clus, sites)
    site_ort <- strand(sites[subjectHits(site_hits)])
    gene_ort <- sites[subjectHits(site_hits)]$in_geneOrt
    if(length(sites[subjectHits(site_hits)]) > 0){
      df <- data.frame(
        "posid" = generate_posid(sites[subjectHits(site_hits)]),
        "site_ort" = as.character(site_ort),
        "gene_ort" = as.character(gene_ort),
        "patient" = sites[subjectHits(site_hits)]$patient,
        stringsAsFactors = FALSE
      )
      df <- distinct(df)
      df$same_orientation <- df$site_ort == df$gene_ort
      score <- paste0(
        "T", length(grep("TRUE", df$same_orientation)), ":",
        "F", length(grep("FALSE", df$same_orientation)), ":",
        "N", length(grep("TRUE", is.na(df$same_orientation)))
      )
    }else{
      score <- "T0:F0:N0"
    }
    score
  })
}

scan_fisher_test_ort <- function(score1, score2){
  sapply(1:length(score1), function(i){
    grp1 <- unlist(strsplit(score1[i], ":"))
    grp2 <- unlist(strsplit(score2[i], ":"))
    x <- matrix(c(
      as.integer(substr(grp1[1], 2, 3)),
      as.integer(substr(grp1[2], 2, 3)),
      as.integer(substr(grp2[1], 2, 3)),
      as.integer(substr(grp2[2], 2, 3))),
      ncol = 2
    )
    fisher.test(x)$p.value
  })
}


scan_sites_count <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(hits))
  })
}

annotate_scan_clusters <- function(scanned_ranges, tdn_sites, 
                                   timepoint_sites, refGenes){
  
  #browser()
  scanned_ranges$n_sites_tdn <- scan_sites_count(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_sites_tp <- scan_sites_count(
    scanned_ranges, timepoint_sites)
  scanned_ranges$sum_abund_tdn <- scan_sites_abundance(
    scanned_ranges, tdn_sites)
  scanned_ranges$sum_abund_tp <- scan_sites_abundance(
    scanned_ranges, timepoint_sites)
  scanned_ranges$n_patients_tdn <- scan_patients(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_patients_tp <- scan_patients(
    scanned_ranges, timepoint_sites)
  scanned_ranges$genes_in_cluster <- scan_genes(
    scanned_ranges, refGenes)
  scanned_ranges$in_gene_ort_tdn <- scan_orientation(
    scanned_ranges, tdn_sites)
  scanned_ranges$in_gene_ort_tp <- scan_orientation(
    scanned_ranges, timepoint_sites)
  scanned_ranges$ort_fisher_test <- scan_fisher_test_ort(
    scanned_ranges$in_gene_ort_tdn, scanned_ranges$in_gene_ort_tp)
  scanned_ranges
}


createClusters <- function(tdn_sites,timePoint_sites,refGenes){
  # Transduction vs All timePoints -----------------------------------------------
  scan_sites <- gintools:::scan_format(tdn_sites, timePoint_sites, grouping = "patient")
  
  scanned_clus <- gRxCluster(
    object = scan_sites$chr,
    starts = scan_sites$pos,
    group = scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 7.5, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.35, na.rm = TRUE))
  
  scan_summary <- gRxSummary(scanned_clus)
  
  
  
  scanned_clus <- annotate_scan_clusters(scanned_clus, tdn_sites, timePoint_sites, refGenes)
  
  scanned_clus <- scanned_clus[scanned_clus$n_patients_tdn > 1 | scanned_clus$n_patients_tp > 1]
  
  enriched_tdn_clus <- scanned_clus[
    (scanned_clus$n_sites_tdn/length(tdn_sites)) >
      (scanned_clus$n_sites_tp/length(timePoint_sites))]
  
  enriched_timePoint_clus <- scanned_clus[
    (scanned_clus$n_sites_tdn/length(tdn_sites)) <
      (scanned_clus$n_sites_tp/length(timePoint_sites))]
  
  
  # Scans for clusters with orientation biases -----------------------------------
  pos_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "+"]
  pos_strand_timePoint_sites <- timePoint_sites[strand(timePoint_sites) == "+"]
  
  neg_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "-"]
  neg_strand_timePoint_sites <- timePoint_sites[strand(timePoint_sites) == "-"]
  
  pos_strand_scan_sites <- gintools:::scan_format(pos_strand_tdn_sites, pos_strand_timePoint_sites, grouping = "patient")
  neg_strand_scan_sites <- gintools:::scan_format(neg_strand_tdn_sites, neg_strand_timePoint_sites, grouping = "patient")
  
  pos_scanned_clus <- gRxCluster(
    object = pos_strand_scan_sites$chr,
    starts = pos_strand_scan_sites$pos,
    group = pos_strand_scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))
  
  pos_scan_summary <- gRxSummary(pos_scanned_clus)
  
  neg_scanned_clus <- gRxCluster(
    object = neg_strand_scan_sites$chr,
    starts = neg_strand_scan_sites$pos,
    group = neg_strand_scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))
  
  neg_scan_summary <- gRxSummary(neg_scanned_clus)
  
  ## Clusters are compared to all sites, not just the same strand
  pos_scanned_clus <- annotate_scan_clusters(
    pos_scanned_clus, tdn_sites[strand(tdn_sites) == "+"],
    timePoint_sites[strand(timePoint_sites) == "+"], refGenes)
  
  neg_scanned_clus <- annotate_scan_clusters(
    neg_scanned_clus, tdn_sites[strand(tdn_sites) == "-"],
    timePoint_sites[strand(timePoint_sites) == "-"], refGenes)
  
  pos_scanned_clus <- pos_scanned_clus[pos_scanned_clus$n_patients_tdn > 1 | pos_scanned_clus$n_patients_tp > 1]
  
  neg_scanned_clus <- neg_scanned_clus[neg_scanned_clus$n_patients_tdn > 1 | neg_scanned_clus$n_patients_tp > 1]
  
  pos_strand_timePoint_enriched_clus <- pos_scanned_clus[
    (pos_scanned_clus$n_sites_tp/length(pos_strand_timePoint_sites)) >
      (pos_scanned_clus$n_sites_tdn/length(pos_strand_tdn_sites))]
  
  neg_strand_timePoint_enriched_clus <- neg_scanned_clus[
    (neg_scanned_clus$n_sites_tp/length(neg_strand_timePoint_sites)) >
      (neg_scanned_clus$n_sites_tdn/length(neg_strand_tdn_sites))]
  
  # Scans for clusters in higher abundance ---------------------------------------
  # over lower abundance in other timePoints
  cutoff <- 2.5 # percent
  
  tp_high_abund_posids <- get_top_sites(timePoint_sites, cutoff, "estAbund")
  timePoint_sites$abund_status <- ifelse(
    timePoint_sites$posid %in% tp_high_abund_posids,
    "High Abundance", "Low Abundance")
  
  high_sites <- timePoint_sites[
    timePoint_sites$abund_status == "High Abundance"]
  low_sites <- timePoint_sites[
    timePoint_sites$abund_status == "Low Abundance"]
  
  high_low_sites <- gintools:::scan_format(
    low_sites, high_sites, grouping = "patient")
  
  high_low_clus <- gRxCluster(
    object = high_low_sites$chr,
    starts = high_low_sites$pos,
    group = high_low_sites$grp,
    kvals = c(10L:35L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 2, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.15, na.rm = TRUE))
  
  high_low_summary <- gRxSummary(high_low_clus)
  
  high_low_clus <- high_low_clus[width(high_low_clus) > 1]
  high_low_clus$n_sites_high <- scan_sites_count(high_low_clus, high_sites)
  high_low_clus$n_sites_low <- scan_sites_count(high_low_clus, low_sites)
  high_low_clus$n_patients_high <- scan_patients(high_low_clus, high_sites)
  high_low_clus$n_patients_low <- scan_patients(high_low_clus, low_sites)
  high_low_clus$genes_in_cluster <- scan_genes(high_low_clus, refGenes)
  high_low_clus$in_gene_ort_high <- scan_orientation(high_low_clus, high_sites)
  high_low_clus$in_gene_ort_low <- scan_orientation(high_low_clus, low_sites)
  high_low_clus$ort_fisher_test <- scan_fisher_test_ort(
    high_low_clus$in_gene_ort_high, high_low_clus$in_gene_ort_low)
  enriched_high_clus <- high_low_clus[
    (high_low_clus$n_sites_high/length(high_sites)) >
      (high_low_clus$n_sites_low/length(low_sites))]
  enriched_low_clus <- high_low_clus[
    (high_low_clus$n_sites_high/length(high_sites)) <
      (high_low_clus$n_sites_low/length(low_sites))]
  
  # Combine all clusters into the same frame
  # and reduce to find all unique clusters
  clusters <- list(
    "timePoint_clusters" = enriched_timePoint_clus,
    "positive_strand_clusters" = pos_strand_timePoint_enriched_clus,
    "negative_strand_clusters" = neg_strand_timePoint_enriched_clus,
    "high_abund_clusters" = enriched_high_clus)
  
  cart19_clusters <- unlist(GRangesList(lapply(
    1:4,
    function(i){
      clus_name <- c("timePoint", "Pos-Strand", "Neg-Strand", "Abundance")[i]
      cluster_group <- granges(clusters[[i]])
      cluster_group$clus.origin <- clus_name
      cluster_group$target.min <- clusters[[i]]$target.min
      cluster_group
    }
  )))
  
  red_clusters <- GenomicRanges::reduce(cart19_clusters, with.revmap = TRUE)
  
  red_clusters$cluster_origin <- sapply(red_clusters$revmap, function(x){
    paste(unique(cart19_clusters[x]$clus.origin), collapse = ", ")})
  
  red_clusters$revmap <- NULL
  
  red_clusters <- annotate_scan_clusters(
    red_clusters, tdn_sites, timePoint_sites, refGenes)
  
  red_df <- dplyr::select(as.data.frame(red_clusters), -strand)
  
  names(red_df) <- c(
    "Chr", "Start", "End", "Width", "Cluster Origin", "TDN: Num. Sites",
    "TP: Num. Sites", "TDN: Abundance Sum", "TP: Abundance Sum",
    "TDN: Num. Patients", "TP: Num. Patients", "Genes In Cluster",
    "TDN: Within Gene Ort.", "TP: Within Gene Ort.",
    "Orientation Fisher Test P-value")
  
  composite_clusters <-  c(clusters, list("all_clusters" = cart19_clusters, "red_clusters" = red_clusters))
  
  
  o <- list(enriched_tdn_clus, enriched_timePoint_clus, composite_clusters, clusters, cart19_clusters, red_clusters, red_df)
  
  return(o)
  
  
}


# real code ------------


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


# Recreate columns used in published data.
intSites$specimen <- intSites$GTSP
intSites$celltype <- intSites$cellType
intSites$timePoint <- tolower(intSites$timePoint) ## imp
intSites$timepoint <- intSites$timePoint
intSites$in_gene <- ifelse(intSites$nearestFeatureDist == 0, intSites$nearestFeature, FALSE)
intSites$in_geneOrt <- ifelse(intSites$nearestFeatureDist == 0, intSites$nearestFeatureStrand, NA) #imp
intSites$nearest_geneDist <- intSites$nearestFeatureDist
intSites$nearest_gene <- intSites$nearestFeature
intSites$nearest_geneOrt <- intSites$nearestFeatureStrand
intSites$gene_id_wo_annot <- intSites$nearestFeature

tdn_sites <- intSites[intSites$timePoint == "d0"]
### tdn_sites <- tdn_sites[tdn_sites$patient %in% std_clin_patients]
timePoint_sites <- intSites[intSites$timePoint != "d0"]
### timePoint_sites <- timePoint_sites[timePoint_sites$patient %in% std_clin_patients]
refGenes <- readRDS(file.path(.cluster_features_d,'hg38.refSeq.rds'))
refGenes <- refGenes[seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))]


# run clusters -----------
o <- createClusters(tdn_sites,timePoint_sites,refGenes)
enriched_tdn_clus <- o[[1]]
enriched_timePoint_clus <- o[[2]]
composite_clusters <- o[[3]]
clusters <- o[[4]]
cart19_clusters <- o[[5]]
red_clusters <- o[[6]]
red_df <- o[[7]]

# build object cluster -----------

uGTSP <- unique(intSites$GTSP)
# clustersRepresented
clustersRepresented_t <- sapply(uGTSP, function(x){
  sites <- intSites[intSites$specimen == x]
  hits <- findOverlaps(sites, red_clusters)
  length(unique(subjectHits(hits)))
})

numSitesInClusters_t <- sapply(uGTSP, function(x){
  sites <- intSites[intSites$specimen == x]
  hits <- findOverlaps(sites, red_clusters)
  length(unique(queryHits(hits)))
})

abundInClusters_t <- sapply(uGTSP, function(x){
  sites <- intSites[intSites$specimen == x]
  hits <- findOverlaps(sites, red_clusters)
  sum(sites[queryHits(hits)]$estAbund)
})

cluster_df <- data.frame(
  #GTSP=uGTSP,
  clustersRepresented=clustersRepresented_t,
  numSitesInClusters=numSitesInClusters_t,
  abundInClusters=abundInClusters_t
) %>% rownames_to_column(var='GTSP')

saveRDS(cluster_df, file = file.path(.features_d,'cluster_GTSP.rds'))
