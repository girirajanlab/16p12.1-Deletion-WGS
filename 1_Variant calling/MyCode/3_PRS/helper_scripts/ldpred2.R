
# Code adapted from: https://choishingwan.github.io/PRS-Tutorial/ldpred/
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(data.table)
library(magrittr)

args <- commandArgs(trailingOnly=T)

# Input and ouput files
INPUT_PATH <- args[1]
SUMMSTAT <- args[2]
PHENO <- args[3]
OUTPUT_PATH <- args[4]

HAPMAP3 <- "/path/to/HapMap3/variants/RData.rds"

print(args)

# Step 1: Load summary statistics file and HapMap dataset
info <- readRDS(HAPMAP3)
sumstats <- fread(SUMMSTAT, fill=T)

# Rename columns
dict <- c("CHR"="chr", "hg19chrc"="chr",
	"BP"="pos", "bp"="pos", "POS"="pos",
	"SNP"="rsid", "snpid"="rsid", "MarkerName"="rsid",
	"A1"="a1", "a1"="a1",
	"A2"="a0", "a2"="a0",
	"OR"="OR", "or"="OR",
	"Beta"="beta", "stdBeta"="beta",
	"SE"="SE", "se"="SE",
	"P"="p", "p"="p", "Pval"="p")

colnames(sumstats) <- dict[colnames(sumstats)]

# Ensure chromosomes are numeric (Skipping chrX SNPs anyway)
sumstats$chr <- gsub( "chr", "", as.character(sumstats$chr))
sumstats$chr <- as.numeric(sumstats$chr)

print("Converted chromosomes to numeric")

# Get betas for autism and schizophrenia (betas provided in intelligence and educational attainment summary stats)
if ("OR" %in% colnames(sumstats) & !("beta" %in% colnames(sumstats))) {
	sumstats$beta <- log(sumstats$OR) #Use log-transformed OR to generate betas. Betas are already provided in some files
}
sumstats$beta_se <- sumstats$SE #Use unmodified standard errors

# Ensure all OR are not NA
sumstats <- sumstats[!is.na(sumstats$beta), ]

print(sumstats)

print("Log-transformed OR")

sumstats <- data.frame(sumstats)
sumstats <- sumstats[sumstats$rsid %in% unlist(info["rsid"]),] #Filter out SNPs not in HapMap3 dataset

print(sumstats[1:10, 1:5])

print("Restricted to SNPs in both datasets")

# Calculate effective sample size from literature if not present in summary statistics. Formula is for binary traits: 4/(1/case + 1/control)
neff_dict <- c("autism"=(4 / (1/18381 + 1/27969)), "intelligence"=0, "educational_attainment"=766345, "schizophrenia"=(4 / (1/36989 + 1/113075)))
sumstats$n_eff <- neff_dict[pheno]
if (pheno=="intelligence") {
	sumstats$n_eff <- sumstats$N_analyzed
}

sumstats <- sumstats[, c("chr","rsid","pos","a0","a1","p","beta","beta_se","n_eff")]
print(sumstats[1:10, 1:5])

# Make sure alleles are upper case
sumstats$a0 <- unlist(lapply(sumstats$a0, toupper))
sumstats$a1 <- unlist(lapply(sumstats$a1, toupper))

# Step 2: Calculate LD matrix
# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL

# preprocess the bed file - this only need to be done once for each BED file, so check for processed RDS file
if (!file.exists(paste(INPUT_PATH, ".rds", sep=""))) {
	bed_loc=paste(INPUT_PATH, ".bed", sep="")
	snp_readBed(bed_loc)
}
# now attach the genotype object
rds_loc=paste(INPUT_PATH, ".rds", sep="")
obj.bigSNP <- snp_attach(rds_loc)
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]

names(map) <- c("chr", "snpid", "pos", "a1", "a0")
map <- subset(map,chr>0) #Remove SNPs w/o location information
map<-subset(map,chr<23) #Keep autosomal SNPs only (chr1-22)

print("Pre-SNP matching")

print(sumstats[1:10, 1:5])
print(map[1:10, ])

# perform SNP matching
info_snp <- snp_match(sumstats, map)

print("SNPs have been matched")

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes

# Rename the data structures
CHR <- map$chr
POS <- map$pos

# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".",ncores = NCORES)
# calculate LD
for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(genotype, ind.col = ind.chr2, ncores = NCORES, infos.pos = POS2[ind.chr2], size = 3 / 1000)
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,c("family.ID", "sample.ID"),c("FID", "IID"))

# Step 3: Perform LD score regression on betas
print('Peforming LD score regression')

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld,  length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]

# Step 4: Calculate PRS using auto model (formerly step 7)
print('Calculating PRS')

multi_auto <- snp_ldpred2_auto(corr,df_beta,h2_init = h2_est,vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),ncores = NCORES)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_auto <-big_prodMat(genotype,beta_auto,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)
# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-rowMeans(beta_auto[,abs(pred_scaled -median(pred_scaled)) < 3 * mad(pred_scaled)])
pred_auto <-big_prodVec(genotype,final_beta_auto,ind.row = ind.test,ind.col = info_snp$`_NUM_ID_`)

# Scale PRS values (mean=0, SD=1)
scores <- scale(pred_auto)

# Add sample IDs (samples go in the same order as in the FAM file)
fam_loc <- paste(INPUT_PATH, ".fam", sep="")
fam <- read.csv(fam_loc, sep=" ", header=F)
colnames(fam) <- c("FID", "IID", "Father", "Mother", "Sex", "Phenotype")

prs_col <- paste(PHENO, "PRS", sep="_")
fam[, prs_col] <- scores

# Save to file
outfile_name <- paste(OUTPUT_PATH, paste(PHENO, "PRS.csv", sep="_"), sep='/')
write.csv(fam[, c("FID", "IID", prs_col)], outfile_name, row.names=F, quote=F)

