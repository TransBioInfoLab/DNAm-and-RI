#######################################################################################################
# =================================================================================================== #
# Function for plot 
# =================================================================================================== #
#######################################################################################################
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
# =================================================================================================== 
# Manhattan plot
# =================================================================================================== 
plot_manh <- function(results, annotated, colored, thres = 0.05, top_n = 10) {
  
  don <- results %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(results, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each CpG
    arrange(chr, pos) %>%
    mutate( POScum=pos+tot) 

  annotated_df <- results %>% 
    separate_longer_delim(GREAT,";") %>% 
    mutate(GREAT = gsub("\\(.*| ", "", GREAT)) %>% 
    filter(GREAT %in% annotated) %>% 
    dplyr::select(pVal.final, GREAT, CpG) %>% 
    group_by(GREAT) %>% 
    slice_min(pVal.final, with_ties = F)
  
  don <- don %>%
    # Add annotation information
    mutate( is_annotate=ifelse(CpG %in% annotated_df$CpG & fdr < 0.05, "yes", "no"),
            is_red = ifelse(CpG %in% colored, "yes", "no")) 
  
  thres <- thres
    
  # Prepare X axis
  axisdf <- don %>% 
    group_by(chr) %>% 
    summarize(center=( max(POScum) + min(POScum) ) / 2 ) %>%
    mutate(chr = ifelse(chr %in% c(11,13,15,17,19,21), "", ifelse(
      chr == 23, "X", ifelse(chr == 24, "Y", chr)
    )))

  scaleFUN <- function(x) sprintf("%.0f", x)

  ggplot(don, aes(x=POScum, y=-log10(pVal.final))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), size=1) +
    scale_color_manual(values = rep(c("grey20", "blue2"), 24 )) +
    #geom_point(data = don %>% filter(is_red == "yes"), mapping = aes(x = POScum, y = -log10(pVal.final)), color = "red") +
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,-log10(1e-9)), labels=scaleFUN) +     # remove space between plot area and x axis
    
    # Add highlighted points
    geom_hline(aes(yintercept = -log10(thres)), color="red") +
    annotate("text", x = axisdf$center[2], y = -log10(1e-8), label = expression("pValue < "~"10"^-5)) +
    annotate("segment", x = 0.6 * axisdf$center[1], xend = 1.3*axisdf$center[1], y = -log10(1e-8), yend = -log10(1e-8),
             colour = "red") + 
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=subset(don %>% 
                                  separate_longer_delim(GREAT,";") %>% 
                                  mutate(GREAT = gsub("\\(.*| ", "", GREAT)), is_annotate=="yes"), 
                     aes(label=GREAT), 
                    max.overlaps = 20,
                     size=4) +
    xlab("Chromosome") +
    ylab(expression("-log"["10"]~"(p)")) + 
    # Custom the theme:
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
}
# =================================================================================================== 
# Histogram
# =================================================================================================== 
# Auxilary function
hist_chart <- function(data, Var1, Var2 = NULL, xmin = 65, xmax = 100, ymin = 110, ymax = 130,
                       save = T, dir.save = ".", prefix = "", width = 8, height = 6, ...){
  
  # Add summary table
  t0 <- data.frame("total",
                   data %>%
                     dplyr::summarise_at(Var1,  list(mean = ~round(mean(.x, na.rm = T), 2),
                                                     sd = ~round(sd(.x, na.rm = T), 2),
                                                     min = ~min(.x, na.rm = T),
                                                     max = ~max(.x, na.rm = T),
                                                     n = ~length(na.omit(.x)),
                                                     frac = ~paste0(round(length(na.omit(.x))/nrow(data),3) * 100, "%"))))
  
  
  # If stratified by Var2
  if(!is.null(Var2)){
    
    t <- data %>% 
      group_by(get(Var2)) %>%
      dplyr::summarise_at(Var1,  list(mean = ~round(mean(.x, na.rm = T), 2),
                                      sd = ~round(sd(.x, na.rm = T), 2),
                                      min = ~min(.x, na.rm = T),
                                      max = ~max(.x, na.rm = T),
                                      n = ~length(na.omit(.x)),
                                      frac = ~paste0(round(length(na.omit(.x))/nrow(data),3) * 100, "%")))
    
    colnames(t0)[1] <- Var2
    colnames(t)[1] <- Var2
    
    #t <- rbind(t, t0)
    
    data[[Var2]] <- as.factor(data[[Var2]])
    t <- na.omit(t)
    
    w <- rstatix::wilcox_test(data, as.formula(paste0(Var1,"~",Var2)))
    t$wilcox <- w$p
    
    
  } else {
    colnames(t0)[1] <- Var1
    t <- t0
  }
  
  
  # Plot histogram
  g1 <- ggpubr::gghistogram(
    data,
    x = Var1,
    fill = ifelse(is.null(Var2), "pink", Var2),
    palette = "jco",
    alpha = 0.5,
    bins = 30,
    position = position_dodge(),
    ...
  ) 
  t1 <- tableGrob(as.matrix(t), theme = ttheme_default(base_size = 6))
  grid.arrange(
    g1,t1,
    nrow = 2, heights = c(4.5, 1.5)
  ) %>% print()
  if(save){
    dev.new()
    pdf(
      file.path(dir.save, paste0(prefix, "_histogram.pdf")),
      width = width,
      height = height
    )
    grid.arrange(g1,t1,
                 nrow = 2, heights = c(4.5, 1.5)) 
    
    dev.off()
  }
}
