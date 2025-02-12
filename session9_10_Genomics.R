if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

BiocManager::install("regioneR")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("karyoploteR")
BiocManager::install("biomaRt")

library("karyoploteR")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
library("regioneR")
library("biomaRt")

#Generate a random genotyping dataset----
set.seed(322)
createDataset <- function(num.snps = 20000, max.peaks = 5) {
 hg19.genome <- filterChromosomes(getGenome("hg19"))
 snps <- sort(createRandomRegions(nregions = num.snps, length.mean = 1, 
                                  length.sd = 0, genome=filterChromosomes(getGenome("hg19"))))
 names(snps) <- paste0("rs", seq_len(num.snps))
 snps$pval <- rnorm(n = num.snps, mean = 0.5, sd = 1)
 snps$pval[snps$pval < 0] <- -1*snps$pval[snps$pval < 0]
 #define the "significant peaks"
 peaks <- createRandomRegions(runif(1, 1, max.peaks), 8e6, 4e6)

 for(npeak in seq_along(peaks)) {
  snps.in.peak <- which(overlapsAny(snps, peaks[npeak]))
  snps$pval[snps.in.peak] <- runif(n = length(snps.in.peak), min = 0.1, max = runif(1, 6, 8))
 }
 snps$pval <- 10^(-1*snps$pval)
 return(list(peaks = peaks, snps = snps))
}


ds <- createDataset()
ds$snps

#Visualizing genotyping data----
#Create a "base" plot
?plotKaryotype
kp <- plotKaryotype(plot.type=4)
#Generate a Manhattan plot on top of the "base"
?kpPlotManhattan
kp <- kpPlotManhattan(kp, data = ds$snps)

#Re-create this plot...
kp <- plotKaryotype(plot.type=4)
#...with a manual region highlight
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = "chr3:1-30000000")

#Re-create this plot...
kp <- plotKaryotype(plot.type = 4)
#...with automatically highlighting significant variants
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, points.cex = 0.8)

#Different options for "base" visualization----
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set1", r0=autotrack(1,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "2blues", r0=autotrack(2,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "greengray", r0=autotrack(3,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "rainbow", r0=autotrack(4,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = c("orchid", "gold", "orange"), r0=autotrack(5,5))

#Creating a gradient for point visualization----
transf.pval <- -log10(ds$snps$pval)
points.col <- colByValue(transf.pval, colors=c("#BBBBBB", "orange"))
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = points.col)

#Modifying line and highlight colors----
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, 
                      highlight = ds$peaks, highlight.col = "orchid",
                      suggestive.col="orange", suggestive.lwd = 3,
                      genomewide.col = "red", genomewide.lwd = 6)

#Controlling axis look----
transf.pval <- -log10(ds$snps$pval)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, pval = transf.pval, logp = FALSE )
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
ymax <- kp$latest.plot$computed.values$ymax
ticks <- c(0, seq_len(floor(ymax)))
kpAxis(kp, ymin=0, ymax=ymax, tick.pos = ticks)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

#Labelling top significant SNP per chromosome----
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)


snps <- kp$latest.plot$computed.values$data
suggestive.thr <- kp$latest.plot$computed.values$suggestiveline
#Get the names of the top SNP per chr
top.snps <- tapply(seq_along(snps), seqnames(snps), function(x) {
 in.chr <- snps[x]
 top.snp <- in.chr[which.max(in.chr$y)]
 return(names(top.snp))
})
#Filter by suggestive line
top.snps <- top.snps[snps[top.snps]$y>suggestive.thr]
#And select all snp information based on the names
top.snps <- snps[top.snps]

top.snps


kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3)

#Fetch information about the significant variants----
snp.db <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start",                   
                      "chrom_strand","associated_gene",
                      "ensembl_gene_stable_id"),
                    values = names(top.snps),
                    mart = snp.db)

#Alternatively, downloading and parsing the database from dbSNP or its subset Clinvar is possible if biomaRt is slow:
#https://www.ncbi.nlm.nih.gov/snp/ --- https://www.ncbi.nlm.nih.gov/clinvar/

#Generate annotated SNP set----
#Genes taken from study at https://www.sciencedirect.com/science/article/pii/S1040842823001087
#Let's generate a processed large scale genetic screening dataset focusing on genes that are risk factors for 
#pancreatic cancer
genes_risk <- c("PDX1", "KLF5", "BCAR1", "HNF1B", "LINC00673", "GRP", "WNT2B", "NOC2L", "ZNRF3",
                "ETAA1", "TP63", "TERT", "TNS3", "SUGCT", "HNF4G", "PVT1")
genes_risk_rownames <- do.call(paste, c(expand.grid(genes_risk, c("low", "intermediate", "high")), sep = "_"))

prognosis_table <- matrix(nrow = length(genes_risk_rownames))
outcome_vector <- c()

set.seed(322)
for (patient in 1:1000) {
  outcome_seed <- sample(1:100, 1)
  if (outcome_seed <= 80) {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps <- sample(10:30, 1)
      snp_entry <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.6, 0.3, 0.1)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome <- c(gene_outcome, snp_entry)
    }
    outcome_vector <- c(outcome_vector, "low")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])

    
  } else if (outcome_seed >= 95) {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps <- sample(10:30, 1)
      snp_entry <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.2, 0.3, 0.5)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome <- c(gene_outcome, snp_entry)
    }
    outcome_vector <- c(outcome_vector, "high")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])
    
  } else {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps <- sample(10:30, 1)
      snp_entry <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.5, 0.3, 0.2)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome <- c(gene_outcome, snp_entry)
    }
    outcome_vector <- c(outcome_vector, "intermediate")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])
  }
}

#Clean up the results datasets
prognosis_table <- prognosis_table[, -1]
rownames(prognosis_table) <- genes_risk_rownames
colnames(prognosis_table) <- names(outcome_vector) <- paste0("patient", 1:1000)
prognosis_table[is.na(prognosis_table)] <- 0

prognosis_table[1:5, 1:5] # X (low high  intermediate )
outcome_vector[1:5] # patient (y-predicted )

#For added "realism" we can include missing data
missing_data <- sample(1:1000, 250)

for (patient in missing_data) {
  sim_missing <- prognosis_table[, patient]
  num_missing <- sample(genes_risk, sample(1:7, 1))
  sim_missing[grepl(paste(num_missing,collapse="|"), names(sim_missing))] <- NA
  prognosis_table[, patient] <- sim_missing
}

#If you decide to make a model in R, you could try employing "caret" package - https://topepo.github.io/caret/
#install.packages("caret")


library(caret)

# Split the dataset into training and testing sets
trainIndex <- createDataPartition(outcome_vector, p = 0.8, list = FALSE)

# Split the data and labels
train_data <- prognosis_table[, trainIndex]
test_data <- prognosis_table[, -trainIndex]
train_labels <- as.factor(outcome_vector[trainIndex])
test_labels <- as.factor(outcome_vector[-trainIndex])

# Impute missing data (replace NAs with column medians)
train_data[is.na(train_data)] <- apply(train_data, 1, function(x) median(x, na.rm = TRUE))
test_data[is.na(test_data)] <- apply(test_data, 1, function(x) median(x, na.rm = TRUE))

# Set up 5-fold cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Train the Random Forest model
set.seed(123)
model_rf <- train(x = t(train_data), y = train_labels, 
                  method = "rf",  # Random Forest
                  trControl = train_control)

# View the model details
print(model_rf)

# Make predictions on the test set
predictions <- predict(model_rf, t(test_data))

# View the first few predictions
head(predictions)

# Evaluate the model performance using a confusion matrix
confusionMatrix(predictions, test_labels)


