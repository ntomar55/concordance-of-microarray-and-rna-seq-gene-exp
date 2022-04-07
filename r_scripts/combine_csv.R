#'R-code for combining count data into single CSV file fro analysis.
#'Additionally, generates box plot of counts for each sample.
#'Outputs combined count data to csv and boxplot to png file.
#'libraries:
library(tidyverse)
library(ggplot2)

#'Helper function for read data function. Removes not needed data
#'and renames count column according to sample name.
#'
#'@param dataframe
#'@result a filtered dataframe consisting only of count data
#'

clean_df <- function(df){
  #get count column
  df <- df %>%
    dplyr::select(ncol(df))
  #rename column according to sample
  counts_name <- colnames(df)[[1]]
  dot <- "\\."
  counts_name <- str_split(counts_name, pattern = dot)[[1]][8]
  counts_name <- str_split(counts_name, pattern = '_')[[1]][1]
  colnames(df) <- counts_name
  
  return(df)
}

#'Reads in multiple count data files from a shared directory, 
#'and returns a single dataframe of the combined read counts
#'where rows are genes and cols are samples.
#'
#'@param directory
#'@return dataframe of combined data
#'
#'@example data <- read_data(path/to/conunt/data)
#'

read_data <- function(path){
  
  files <- list.files(path=path, pattern="*.counts") #get all file names in list
  files <- lapply(files, function(x){str_c(path, x, sep = "/")})
  data <- lapply(files, function(x){read.table(x, header=T)}) #read in each file
  
  genes <- as.data.frame(data[[1]])[,1]
  data <- lapply(data, function(x){clean_df(x)})
  
  combined_df <- data %>%
    bind_cols()
  
  rownames(combined_df) <- genes
  return(combined_df)
}


#'Generate box plot demonstrating the distribution of read counts
#'across all samples.
#'
#'@param dataframe
#'@return ggplot boxplot object
#'
#'@example plot <- box_plot(counts_df)
#'

box_plot <- function(count_df){
  
  long_count_df <- count_df %>%
    gather(key = sample, value = counts, factor_key = T)
  
  box_p <- long_count_df %>%
    ggplot(aes(x=sample,y=log(counts+1), fill = sample)) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, hjust=1)) 
  
  return(box_p)
}
