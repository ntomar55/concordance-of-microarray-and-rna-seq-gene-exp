#deseq analysis
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(gt)


#' Read in data as dataframe from input csv file.
#' Input csv file is a combined csv file consisting of all count data from 
#' all samples.
#' 
#' @param path
#' @return dataframe
#' path to data:/usr4/bf527/dlenci/Documents/project-3-saxophone-1/r_script_out/combined_counts.csv
read_data <- function(path){

  df <- read.csv(path, header = T)
  
  return(df)
  
}

#' Take all data files and generate a count matrix consisting of control
#' and non-control data. 
#' 
#' @param combined_counts: combined count matrix (excluding control data)
#' @param metadata: metadata df
#' @param control_df
#' 
#' @return count matrix
#' 
clean_data <- function(counts, meta, control){
  
  samples <- meta %>%
    dplyr::pull(Run)
  
  control <- control %>%
    dplyr::select_if(names(.) %in% samples)
  
  control_samples <- colnames(control)
  treated_samples <- colnames(counts)
  
  combined_counts <- bind_cols(counts, control)
  rownames(combined_counts) <- combined_counts[,1]
  combined_counts <- combined_counts[,-1]
  
  return(combined_counts)

}

#' Next we run the DE_seq analysis using the newly combined count data consisting
#' of both control and treated samples, and the corresponding meta data. Design 
#' for DESeq will be based off of the mode_of_action column in the meta data.
#' 
#' @param counts
#' @param meta_data
#' @param group: group comparing to controls (AhR, CAR/PXR, DNA_Damage)
#' 
#' @return DESeq results as dataframe
#' 

DE_analysis <- function(counts, meta, group){
  
  if(group=="AhR" | group=="CAR/PXR"){
    veh="CORN_OIL_100_%"
  } else{
    veh="SALINE_100_%"
  }
  
  #filter for genes with no expression
  #deseq recommends filtering out genes < 10 counts
  counts <- subset(counts, rowSums(counts) >= 10)
  #select only for the sample of interest.
  meta_filter <- meta %>%
    filter(mode_of_action=="Control" | mode_of_action==group) %>%
    filter(vehicle==veh)
  
  filtered_samples <- meta_filter %>%
    dplyr::pull(Run)
  
  subset_counts <- counts %>%
    select(filtered_samples)
  
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = meta_filter,
                                design = ~ mode_of_action)
  
  dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')
  dds <- DESeq(dds)
  res <- results(dds, contrast = c('mode_of_action', group, 'Control'))
  res <- lfcShrink(dds, coef=2)
  
  res <- as.data.frame(res) %>%
    arrange(pvalue) %>%
    filter(padj<0.05)
  
  return(res)
  
}

#' Generate a table to show the top 10 DE expressed genes from the 3 different
#' analyses. df_cleaner is helper function that takes a df an returns the top 10
#' DE genes by pvalue.
#' 
#' @param DESeq_results_df
#' @param DESeq_results_df
#' @param DESea_results_df
#' 
#' @return gt table
#' 
df_cleaner <- function(df){
  
  df <- rownames_to_column(df)
  df <- df %>%
    dplyr::rename(genes=rowname) %>%
    select(genes, pvalue, treatment) %>%
    dplyr::slice_head(n=10)
  df$pvalue <- formatC(df$pvalue, digits = 4)
  return(df)
}

top10_table <- function(AhR, CAR, DNA){
  
  AhR <- AhR %>%
    mutate(treatment="AhR")
  CAR <- CAR %>%
    mutate(treatment="CAR_PXR")
  DNA <- DNA %>%
    mutate(treatment="DNA_Damage")
  
  DESeq_res <- list(AhR, CAR, DNA)
  top_10 <- lapply(DESeq_res, df_cleaner)
  
  top_10 <- top_10 %>%
    bind_rows() %>%
    group_by(treatment) %>%
    mutate(row = row_number()) %>%
    mutate(pvalue=as.character(pvalue)) %>%
    pivot_wider(names_from = treatment,
                values_from = c(genes, pvalue)) %>%
    select(-row) %>%
    relocate(pvalue_AhR, .after = genes_AhR) %>%
    relocate(pvalue_CAR_PXR, .after = genes_CAR_PXR) 
  
  table <- top_10 %>%
    gt() %>%
    # group columns
    tab_spanner(label = "AhR vs Control",
                columns = c(genes_AhR, pvalue_AhR)) %>% 
    tab_spanner(label = "CAR/PXR vs Control",
                columns = c(genes_CAR_PXR, pvalue_CAR_PXR)) %>%
    tab_spanner(label = "DNA Damage vs Control",
                columns = c(genes_DNA_Damage, pvalue_DNA_Damage)) %>%
    # color to easily distinguish
    tab_style(style = list(cell_fill(color = "lightgrey")),
              locations = cells_body(columns = c(genes_CAR_PXR, pvalue_CAR_PXR))) %>%
    # add footnote
    tab_footnote(footnote = "AhR samples treated with leflunomide.",
                 locations = cells_column_spanners(spanners="AhR vs Control" )) %>%
    tab_footnote(footnote = "CAR/PXR samples treated with fluconazole.",
                 locations = cells_column_spanners(spanners="CAR/PXR vs Control")) %>%
    tab_footnote(footnote = "DNA Damaged samples treated with ifosfamide.",
                 locations = cells_column_spanners(spanners="DNA Damage vs Control")) %>%
    # change column labels
    cols_label(genes_AhR = "Genes",
               pvalue_AhR = "P-Value   ",
               genes_CAR_PXR = "Genes",
               pvalue_CAR_PXR = "P-Value   ",
               genes_DNA_Damage = "Genes",
               pvalue_DNA_Damage = "P-Value.  ") %>%
    # add title and subtitle
    tab_header(title = "Top 10 Differentially Expressed Genes",
               subtitle = "Generated By DESeq2")
  
  return(table)
  
}

#' Generate histograms showing the distribution of log2FC for significantly
#' DE genes for eahc DESeq Analysis.
#' 
#' @param DESeq_results_df
#' @param DESeq_results_df
#' @param DESea_results_df
#' 
#' @return plot
#' 
df_hist_cleaner <- function(df){
  df <- df %>%
    select(pvalue, log2FoldChange, treatment)
  return(df)
}


create_plots <- function(AhR,  CAR, DNA) {
  
  AhR <- AhR %>%
    mutate(treatment="AhR vs Control")
  CAR <- CAR %>%
    mutate(treatment="CAR/PXR vs Control")
  DNA <- DNA %>%
    mutate(treatment="DNA Damage vs Control")
  
  DESeq_res <- list(AhR, CAR, DNA)
  results <- lapply(DESeq_res, df_hist_cleaner)
  
  results <- results %>%
    bind_rows()
  
  labs <- c("AhR vs Control", "CAR/PXR vs Control", "DNA Damage vs Control")
  
  hist_plot <- results %>%
    ggplot(aes(x = log2FoldChange)) +
    geom_histogram(bins = 50, fill = "cyan", color = "black") +
    facet_wrap(~treatment, 
               nrow = 3,
               scales = "free_y") +
    theme_bw()
  
  scat_plot <- results %>%
    ggplot(aes(x = log2FoldChange, y = pvalue)) +
    geom_point() +
    facet_wrap(~treatment, 
               scales = "free_y") +
    theme_bw()
  
  plots <- list(scat_plot, hist_plot)
  
  return(plots)
}
