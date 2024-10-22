# Manhattan plot
dir.base <- "."
dir.results <- file.path(dir.base, "analysis_results")
dir.results.res <- file.path(dir.results, "resilience")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync/") 
dir.results.combp <- file.path(dir.results, "combp")

source(file.path(dir.base, "code/utility/plot.R")) 

results <- read_csv(
  file.path(dir_results_res,
            "Resilience_rlm_model_results_bacon.csv"),
  show_col_types = F
)

cpgs.sig <- results %>% filter(pValue.bacon < 1e-05) %>%
  pull(cpg) 

sig.dmr.df <- read_csv(
  file.path(dir.results.combp, "combp_resilience_rlm_filtered_results_annotated.csv")
) 
sig.dmr <- str_split(sig.dmr.df$cpgs_in_region, ",") %>% unlist() ## 134
chr <- gsub("chr", "", results$seqnames)

cpgs_df <- data.frame(
  "CpG" = results$cpg,
  "pVal.final" = results$pValue.bacon,
  "fdr" = results$fdr.bacon,
  "pos" = results$start,
  "chr" = as.numeric(ifelse(chr == "X", 23, ifelse(chr == "Y", 24, chr))),
  "GREAT" = results$GREAT_annotation
) 
#cpgs_df$GREAT[cpgs_df$GREAT == "RBM14-RBM4"] <- "RBM14"

annotated_cpg <- results %>% 
  filter(cpg %in% cpgs.sig) %>%
  separate_longer_delim(GREAT_annotation,";") %>% 
  mutate(TSS = as.numeric(str_extract(GREAT_annotation, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2500 & abs(TSS) > 100) %>%
  pull(GREAT_annotation) %>%
  gsub("\\(.*| ", "",.)
#annotated_cpg[annotated_cpg == "RBM14-RBM4"]  <- "RBM14"
annotated_dmr <- results %>% 
  filter(cpg %in% sig.dmr & pValue.bacon < 1e-05) %>%
  separate_longer_delim(GREAT_annotation,";") %>% 
  mutate(TSS = as.numeric(str_extract(GREAT_annotation, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2500 & abs(TSS) > 100) %>%
  pull(GREAT_annotation) %>%
  gsub("\\(.*| ", "",.)

anno <- unique(c(annotated_cpg, annotated_dmr))
anno <- anno[!anno %in% "ENSG00000250644"]
anno <- c(anno, "CCDC32")
plot_manh(cpgs_df, annotated = anno, colored = cpgs.sig, thres = 1e-05)

ggplot2::ggsave(
  filename = file.path(dir.results.res,"manhattan_plot.pdf"),
  width = 10,
  height = 6
)
