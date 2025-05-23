#install.packages("data.table")
library("data.table")
#install.packages('curl')
library('curl')

print(paste("Started at:", Sys.time()))

# Get inputs
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
fam <- args[2]
output_dir <- args[3]
print(filename)

# Read in file
file <- fread(filename, sep = " ", header = FALSE, na.strings = "NA")
print(paste("Finished file upload at:", Sys.time())) 

# Read in file suffix
name_list <- as.character(unlist(strsplit(as.character(lapply(strsplit(filename, "/"), tail, 1)), "_")))
suffix <- name_list[3]

# Read in sample IDs
ukb_ids <- as.character(unlist(fread(fam, sep = " ", header = FALSE)[, 1]))

# Other file information
chrom <- name_list[2]
type <- name_list[1]

# Save each column as file into a separate folder for each person
for (i in 1:length(colnames(file))) {
	# Get column
	x <- file[ , ..i]
	# Output filename
	ukbid <- ukb_ids[i]
	outfile <- paste(output_dir, "/sample_split/", ukbid, "/", ukbid, "_", type, "_", chrom, "_", suffix, ".txt", sep = '')
	folder <- paste(output_dir, "/sample_split/", ukbid, "/", sep = '')
	print(outfile)

	# If folder does not already exist, make it
	if (!file.exists(folder)) {
		dir.create(folder)
	}
	fwrite(x, outfile, sep = " ", na = "NA", col.names = FALSE)
}
