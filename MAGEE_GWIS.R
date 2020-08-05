library(readr)
library(MAGEE)

Sys.setenv(MKL_NUM_THREADS=1)

# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure_names <- args[2]
gdsfile <- args[3]
groupfile <- args[4]
min_MAF <- as.numeric(args[5])
max_MAF <- as.numeric(args[6])
ncores <- as.integer(args[7])

# Read in null model object
null_model <- readRDS(null_modelfile)

# Parse exposures
exposures <- strsplit(exposure_names, split=" ")[[1]]

# Remove header from group file if necessary
if (grepl("group", readLines(groupfile, n=1))) {
  system(paste0("tail -n +2 ", groupfile, " > groupfile.tmp"))
  groupfile <- "groupfile.tmp"
}

# Run GWIS
if (groupfile == "") {
  glmm.gei(null_model, interaction=exposures, geno.file=gdsfile, ncores=ncores,
	   outfile="magee_res", MAF.range=c(min_MAF, max_MAF))
} else {
  res <- MAGEE(null_model, interaction=exposures, gdsfile, groupfile, 
  	       MAF.range=c(min_MAF, max_MAF), tests=c("JV", "JF", "JD"), 
	       ncores=ncores)
  write_delim(res, "magee_res", delim=" ")
}
