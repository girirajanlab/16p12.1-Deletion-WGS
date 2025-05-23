library(sva)
library(WGCNA)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(nlme)
library(org.Hs.eg.db)
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores())

# Input and output files
ABUNDANCE="/path/to/compiled/abundace/matrix.csv"
SAMPLE_INFO="/path/to/sample/information.csv"

OUTPUT_DIR="/path/to/output/directory/"

# Load data
count_matrix <- matrix(read.csv(ABUNDANCE))
samp_info <- read.csv(SAMPLE_INFO)

# Run batch correction with ComBat_seq
batch <- samp_info$Batch
tpm <- ComBat_seq(count_matrix, batch, group=NULL)

# Log transform TPM data
data <- log2(tpm+1)

# Keep the top 20,000 most variable genes
dat <- as.data.frame(t(data[order(apply(data,1,mad), decreasing = T)[1:20000],]))

# Unsigned hybrid analysis
# Define cut-off
R.sq_cutoff = 0.8 

# Identify optimal threshold
powers <- c(seq(1, 20, by = 1), seq(22, 30, by = 2))
sft <- pickSoftThreshold(dat,
                         networkType = "signed hybrid",
                         powerVector = powers,
                         RsquaredCut = R.sq_cutoff,
                         verbose = 5)
# Create diagnostic plot
pdf(paste(OUTPUT_DIR, 'WGCNA_soft_threshold_plot.pdf', sep=''))
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=R.sq_cutoff, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers,
     cex=cex1,
     col="red")
abline(h=100,col="red")
dev.off()

# Choose the optimal power
power <- sft$powerEstimate

# Build network
allowWGCNAThreads()
net <- blockwiseModules(dat,
                        power = power,
                        maxBlockSize = ncol(dat),
                        corType = "pearson", 
                        networkType = "signed hybrid",
                        TOMType = "unsigned",
                        minModuleSize = 100, 
                        mergeCutHeight = 0.25, 
                        numericLabels = TRUE,
                        saveTOMs = F,
                        verbose = 3
)

# Save network
saveRDS(net, paste(OUTPUT_DIR, 'WGCNA_network.rds', sep=''))

# Make a list of the genes in each module
moduleColors <- labels2colors(net$colors)
gene_module <- data.frame(Gene=colnames(dat),
                          Module=moduleColors)
write.csv(gene_module, file = paste(OUTPUT_DIR, "module_genes.csv", sep='') row.names = F)

# Calculate the Eigengenes for each module
datME <- moduleEigengenes(dat, moduleColors)$eigengene
datME <- orderMEs(datME)

# Save Eigengenes
write.csv(datME, file=paste(OUTPUT_DIR, "module_eigengenes.csv", sep=''))

# Relate Eigengene values to different sample traits
moddf <- cbind(samp_info, datME)
modules <- colnames(datME)
outdf <- data.frame()
for (m in modules) {
  mod <- lme(as.formula(paste(m, '~ Family + Carrier + Sex')), random= ~ 1 | Individual, data=moddf)
  
  # Parse results
  nfixed <- length(names(mod$coefficients$fixed))
  nrandom <- length(rownames(mod$coefficients$random$Individual))
  res <- data.frame(Module=rep(m, nfixed+nrandom),
                    Variable=c(names(mod$coefficients$fixed), rownames(mod$coefficients$random$Individual)),
                    Effect_Type=c(rep('Fixed', nfixed), rep('Random', nrandom)),
                    Effect=c(unlist(mod$coefficients$fixed), unlist(mod$coefficients$random$Individual)),
                    Std_Error=c(summary(mod)$tTable[, 'Std.Error'], rep(NA, nrandom)),
                    pvalue=c(summary(mod)$tTable[, 'p-value'], rep(NA, nrandom)))
  outdf <- rbind(outdf, res)
}

# Save
write.csv(outdf, paste(OUTPUT_DIR, 'WGCNA_LMM.csv', sep=''), row.names = F)
