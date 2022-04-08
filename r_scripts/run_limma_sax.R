library(limma)

# Pre-requisite for probe conversions:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("rat2302.db")

library(annotate)
library(rat2302.db)

setwd('/projectnb2/bf528/users/saxophone/data_p3/limma_results') # navigate to data dir

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_3_mic_info.csv',
                    as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

# Compartamentalize analysis into function so it can be repeated for
# all three drugs. `drug` argument should be a string containing the name
# of the drug exactly as it appears in `samples$chemical`.
analyze_drug <- function(drug) {
  drug.idx <- samples$chemical == drug | samples$chemical == 'Control'
  
  # subset the full expression matrix to just those in this comparison
  rma.subset <- rma[paste0('X',samples$array_id[drug.idx])]
  
  # construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(
    ~factor(
      samples$chemical,
      levels=c('Control',drug)
    )
  )
  colnames(design) <- c('Intercept', drug)
  
  # run limma
  fit <- lmFit(rma.subset, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')
  symbols <- mapIds(rat2302.db, rownames(t), "SYMBOL", 'PROBEID',
                    multiVals = 'asNA') # Do not map ambiguous probes
  t$symbol <- symbols
  print(drug)
  print('Number of probes:')
  print(nrow(t))
  print('Number of significant probes:')
  print(sum(t$adj.P.Val < 0.05))
  print('Number of significant probes without unique mapping:')
  print(sum(is.na(t$symbol) & t$adj.P.Val < 0.05))
  
  t.no.na <- na.omit(t)
  t.no.duplicates <- t.no.na[!duplicated(t.no.na$symbol),]
  print('Number of significant genes:')
  print(sum(t.no.duplicates$adj.P.Val < 0.05))
  
  # write out the results to file
  write.csv(t, paste0(drug, '_limma_results.csv'))
  write.csv(t.no.duplicates, paste0(drug, '_limma_results_clean.csv'))
}

analyze_drug('LEFLUNOMIDE')
analyze_drug('FLUCONAZOLE')
analyze_drug('IFOSFAMIDE')
