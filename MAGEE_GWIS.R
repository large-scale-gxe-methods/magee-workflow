library(readr)
library(MAGEE)

Sys.setenv(MKL_NUM_THREADS=1)

# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure_names <- args[2]
gdsfile <- args[3]
groupfile <- args[4]
ncores <- as.integer(args[5])

# Read in null model object
null_model <- readRDS(null_modelfile)

# Parse exposures
exposures <- strsplit(exposure_names, split=" ")[[1]]

# Run GWIS
res <- MAGEE(null_model, interaction=exposures, gdsfile, groupfile, 
	     tests=c("JV", "JF", "JD"), ncores=ncores)
write_delim(res, "magee_res", delim=" ")
