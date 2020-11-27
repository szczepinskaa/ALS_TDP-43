library(limma)
library(edgeR)
library(dplyr)
library(R.utils)
library(EWCE)

mkdirs('Results')
mkdirs('Results/Tables')
mkdirs('Results/Figures')

data_5 <- read.delim('Data/GSE99353_5months_cortex_raw_counts.txt', sep='')
exp_matrix_5 <- data_5[,c(7,14:19, 26:33)]
data_20 <- read.delim('Data/GSE112575_20months_cortex_raw_counts.txt', sep='')
exp_matrix_20 <- data_20[,c(7,14:21, 32:41)]
exp_matrix <- merge(exp_matrix_5, exp_matrix_20, by='Feature')
exp_matrix <- as.matrix(exp_matrix)
rownames(exp_matrix) <- exp_matrix[,1]
exp_matrix <- exp_matrix[,-1]
class(exp_matrix) <- 'numeric'
colnames(exp_matrix)

condition <- c(rep('ctrl_5', 6), rep('Q331K_5', 8), rep('ctrl_20', 8), rep('Q331K_20', 10))
sample_annot <- data.frame(row.names = colnames(exp_matrix), condition)
design <- model.matrix(~0+sample_annot$condition)
colnames(design) <- gsub('sample_annot$condition','', colnames(design), fixed=TRUE)
contrasts <- makeContrasts(pre = Q331K_5 - ctrl_5,
                           post = Q331K_20 - ctrl_20,
                           levels=colnames(design))

counts <- DGEList(exp_matrix, samples=sample_annot)

keep <- filterByExpr(counts, design)
counts <- counts[keep,, keep.lib.sizes=FALSE]
counts <- calcNormFactors(counts)

voom <- voom(counts, design, plot=TRUE)
fit <- lmFit(voom, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)


contrasts_names <- colnames(contrasts)
for(i in contrasts_names){
  test <- topTable(fit, coef=i, adjust="BH", number=1000000, sort.by = 't')
  colnames(test)[1] <- 'MGI.symbol'
  write.csv(test, file=sprintf('Results/Tables/%s',  paste('DE', i, 'csv', sep='.')), row.names=FALSE)
}


# EWCE --------------------------------------------------------------------

load(file='ctd.rda')
t_EWCE <- list.files(path='Results/Tables/', pattern='DE.*.csv')

for (i in t_EWCE){
  t_data_EWCE <- read.csv(paste('Results/Tables/', i, sep=''))
  EWCE_results <- ewce_expression_data(sct_data=ctd,tt=t_data_EWCE, thresh=250,
                                       reps=1000000, annotLevel=1, ttSpecies="mouse",sctSpecies="mouse")
  EWCE_results$joint_results$p.adj = p.adjust(EWCE_results$joint_results$p)
  save(EWCE_results, file=sprintf('Results/Tables/%s', paste('EWCE_', gsub('.csv','', gsub('DE.','',i)), '.rda', sep='')))
  EWCE_results$joint_results[order(EWCE_results$joint_results$sd_from_mean),]
  write.csv(EWCE_results$joint_results, file=sprintf('Results/Tables/%s', paste(gsub('DE.','',i), sep='')))
  pdf(file=paste("Results/Figures/", 'EWCE_', gsub('.csv','',i), '.pdf', sep=''))
  print(ewce.plot(EWCE_results$joint_results))
  dev.off()
}

