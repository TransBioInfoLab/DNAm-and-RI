dir.base <- "."

dir.data.combp <- file.path(dir.base, "Mike")
dir_data <- file.path(dir.base, "data")
dir_results <- file.path(dir.base, "analysis_results")
dir_results_res <- file.path(dir_results, "resilience")
dir.data.aux <- file.path(dir.base,"../DATASETS/Aux_Sync/") 
dir.results.combp <- file.path(dir_results, "combp")

library(tidyverse)
library(rGREAT)
library(BiocParallel)

combp_results <- readxl::read_xlsx(
  file.path(dir.data.combp, "cnew.regions-p.bed.xlsx")
)
combp_results$`#chrom` <- ifelse(
  combp_results$`#chrom` == 0, "X", combp_results$`#chrom`
)
  

results <- read_csv(
  file.path(dir_results_res,
            "Resilience_rlm_model_results_bacon.csv")
)

# Auxillary function
## 1. Annotate results for EPICv2 only
anno_df2 <- readxl::read_xlsx(
  file.path(dir_data, "Supp/EPICv2_anno.xlsx"),
  skip = 1
)

AnnotateResults_EPICv2 <- function (lmmRes_df, nCores_int = 1L, UCSCinfo_df, ...) {
  stopifnot("data.frame" %in% class(lmmRes_df), all(c("chrom", 
                                                      "start", "end") %in% colnames(lmmRes_df)))
  lmmRes_df$start <- as.integer(lmmRes_df$start)
  lmmRes_df$end <- as.integer(lmmRes_df$end)

  locations_df <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations

  IslandsUCSCinfo_df <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Islands.UCSC
 
  locations_df <- as.data.frame(locations_df)
  locations_df$cpg <- gsub("_.*", "", row.names(locations_df))
  rownames(locations_df) <- NULL
  locations_df <- unique(locations_df)

  interestingColumns_char <-  c("geneNames", "RelationToGene", "distToTSS", "UCSC_RefGene_Accession")
  UCSCinfo_df <- UCSCinfo_df[, interestingColumns_char]

  IslandsUCSCinfo_df <- as.data.frame(IslandsUCSCinfo_df)
  IslandsUCSCinfo_df <- unique(IslandsUCSCinfo_df)
  cluster <- coMethDMR::CreateParallelWorkers(nCores_int, ...)
  resultsAnno_ls <- bplapply(seq_len(nrow(lmmRes_df)), function(row) {
    .AnnotateRow(row_df = lmmRes_df[row, ], loc_df = locations_df, 
                 info_df = UCSCinfo_df, island_df = IslandsUCSCinfo_df)
  }, BPPARAM = cluster)
  do.call(rbind, resultsAnno_ls)
}
.AnnotateRow <- function(row_df, loc_df, info_df, island_df){
  # browser()
  
  
  ###  Filter Data Frames  ###
  # Extract Row Region
  chr   <- row_df$chrom
  start <- row_df$start
  end   <- row_df$end
  
  # Find Probes in that Region
  chr_df  <- loc_df[loc_df$chr == chr, ]
  inRegion_lgl <- chr_df$pos >= start & chr_df$pos <= end
  out_df <- chr_df[inRegion_lgl, ]
  probes_char <- out_df$cpg
  
  # Find UCSC Annotation Information for those Probes
  infoOut_df <- info_df[probes_char, ]
  
  # Find UCSC Relation to Island Information for those Probes
  islandOut_df <- island_df[probes_char, ]
  
  
  ###  Wrangle UCSC Annotation  ###
  refGeneGroup_char <- .ExtractUCSCinfo(infoOut_df$RelationToGene)
  refGeneAcc_char   <- .ExtractUCSCinfo(infoOut_df$UCSC_RefGene_Accession)
  refGeneName_char  <- .ExtractUCSCinfo(infoOut_df$geneNames)
  refGenedist  <- .ExtractUCSCinfo(infoOut_df$distToTSS)
  
  refIslandRelation_char <- sort(unique(islandOut_df$Relation_to_Island))
  
  
  ###  Return Annotated 1-Row Data Frame  ###
  row_df$RelationToGene <-
    paste0(unique(refGeneGroup_char), collapse = ";")
  row_df$UCSC_RefGene_Accession <-
    paste0(unique(refGeneAcc_char), collapse = ";")
  row_df$UCSC_RefGene_Name <-
    paste0(unique(refGeneName_char), collapse = ";")
  row_df$Relation_to_Island <-
    paste0(unique(refIslandRelation_char), collapse = ";")
  row_df$distToTSS <-
    paste0(unique(refGenedist), collapse = ";")
  
  row_df
  
}
.ExtractUCSCinfo <- function(infoCol) {
  sort(
    unique(
      unlist(
        strsplit(infoCol, ";")
      )
    )
  )
}
## 2. Get CpGs in region for EPICv2 only
GetCpGsInRegion <- function(region, anno_df) {
  sp <- stringr::str_split(region,":|-")[[1]]
  region_chr <- sp[1]
  region_start <- as.numeric(sp[2])
  region_end <- as.numeric(sp[3])
  
  cpgs <- anno_df %>% 
    filter(seqnames %in% region_chr & start >= region_start & start <= region_end) %>%
    pull(Name)
  cpgs <- gsub("_.*", "",  cpgs) 
  cpgs <- unique(cpgs)
  
  cpgs
    
}
## 3. Add annotation
add_combp_annotation <- function(result, cpg, dir.data.aux, additional_info){
  
  # load auxilary data
  load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  
  data <- readr::read_tsv(
    file.path(dir.data.aux,"AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
  )
  CellType.selected <- readxl::read_xlsx(
    file.path(dir.data.aux,"Nassser study selected biosamples.xlsx"),col_names = FALSE
  ) %>% dplyr::pull(1)
  
  data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
    dplyr::filter(!isSelfPromoter)  %>% 
    dplyr::filter(class != "promoter")
  
  nasser.enhancer.gr <- data.filtered %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "chr",
    keep.extra.columns = TRUE
  )
  
  message("Annotating E073_15_coreMarks_segments")
  result$seqnames <- paste0("chr",result$`#chrom`) 
  result$region <- paste0(result$seqnames,":",result$start,"-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region,hits$region)]
  
  message("Annotating GREAT")
  
  job <- submitGreatJob(result.gr, species = "hg38")
  regionsToGenes_gr <- rGREAT::getRegionGeneAssociations(job)
  regionsToGenes <- as.data.frame(regionsToGenes_gr)
  
  GREAT_annotation <- lapply(
    seq_len(length(regionsToGenes$annotated_genes)),
    function(i) {
      g <- ifelse(regionsToGenes$dist_to_TSS[[i]] > 0,
                  paste0(regionsToGenes$annotated_genes[[i]],
                         " (+", regionsToGenes$dist_to_TSS[[i]], ")"),
                  paste0(regionsToGenes$annotated_genes[[i]],
                         " (", regionsToGenes$dist_to_TSS[[i]], ")"))
      paste0(g, collapse = ";")
    }
  )
  great <-  dplyr::select(
    regionsToGenes, "seqnames", "start", "end", "width")
  great <- data.frame(
    great, GREAT_annotation = unlist(GREAT_annotation))
  
  result <- dplyr::left_join(result, great,by = c("seqnames","start","end"))
  
  message("Annotating Island")
  
  result$chrom <- result$seqnames
  result <- AnnotateResults_EPICv2(result, nCores_int = 10, UCSCinfo_df = additional_info)
  result$chrom <- result$chr <- NULL
  
  message("Annotating enhancer")
  
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, nasser.enhancer.gr) %>% as.data.frame()
  result$nasser_is_enahncer <- FALSE
  result$nasser_is_enahncer[unique(hits$queryHits)] <- TRUE
  result$nasser_is_enahncer_cell_types <- NA
  result$nasser_is_enahncer_cell_types[unique(hits$queryHits)]  <- sapply(unique(hits$queryHits), function(x){nasser.enhancer.gr$CellType[hits$subjectHits[hits$queryHits %in% x]] %>% unique %>% paste(collapse = ",")})
  
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
    anno, start.field = "pos", end.field = "pos", keep.extra.columns = TRUE
  )
  anno_df <- as.data.frame(anno_gr)
  
  result$cpgs_in_region <- sapply(
    result$region, 
    function(x){
      GetCpGsInRegion(x, anno_df) %>% 
        intersect(cpg) %>% 
        paste(collapse = ",")
    })
  return(result)
  
}

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
  anno, start.field = "pos", end.field = "pos", keep.extra.columns = TRUE
)
anno_df <- as.data.frame(anno_gr)
anno_df <-  anno_df %>% 
  mutate(cpg = gsub("_.*", "",anno_df$Name)) %>%
  dplyr::select(cpg, seqnames, start, end, width, 
                UCSC_RefGene_Accession) %>% 
  unique()

anno_df2 <- anno_df2 %>% 
  dplyr::mutate(cpg = gsub("_.*", "", ProbeID)) %>%
  dplyr::select(cpg, geneNames, RelationToGene, distToTSS) %>% 
  unique() 
# Remove duplicates in cpgs without gene annotation
dup_cpgs <- anno_df2$cpg[duplicated(anno_df2$cpg)]
dup_cpgs_df <- anno_df2[anno_df2$cpg %in% dup_cpgs,] %>% 
  na.omit() %>%
  group_by(cpg) %>%
  dplyr::summarise(geneNames = paste(geneNames, collapse = ";"),
                   RelationToGene = paste(RelationToGene, collapse = ";"),
                   distToTSS = paste(distToTSS, collapse = ";"))
anno_df2 <- anno_df2 %>% filter(!cpg  %in% dup_cpgs)
anno_df2 <- rbind(anno_df2, dup_cpgs_df)

anno_df <- left_join(anno_df, 
                     anno_df2)
rownames(anno_df) <- anno_df$cpg

# Add annotation
results_anno <- add_combp_annotation(
  result = combp_results,
  cpg = results$cpg,
  dir.data.aux = dir.data.aux,
  anno_df
)

# Add direction
direction <- plyr::laply(
  str_split(results_anno$cpgs_in_region, ","),
  .fun = function(cpgs){
    est <- results %>% filter(cpg %in% cpgs) %>% pull(Estimate.bacon)
    paste0(ifelse(est < 0, "-", "+"), collapse = "")
  }
)
results_anno$direction <- direction
# Calculate percentage of direction
pct.direction <- plyr::laply(
  str_split(direction, ""),
  .fun = function(d){
    sum(d == "-")/length(d)
  }
)
results_anno$pct_direction <- pct.direction

write_csv(
  results_anno,
  file.path(dir.results.combp, "combp_resilience_rlm_results_annotated.csv")
)

results_filtered <- results_anno %>%
  filter(n_probes >= 3, z_sidak_p < 0.05, pct_direction %in% c(0,1), z_p < 1e-05)

write_csv(
  results_filtered,
  file.path(dir.results.combp, "combp_resilience_rlm_filtered_results_annotated.csv")
)
