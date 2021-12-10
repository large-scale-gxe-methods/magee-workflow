library(readr)
library(MAGEE)

Sys.setenv(MKL_NUM_THREADS=1)

# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure_names <- args[2]
int_covar_names <- args[3]
gdsfile <- args[4]
min_MAF <- as.numeric(args[5])
max_MAF <- as.numeric(args[6])
ncores <- as.integer(args[7])
gds_filter <- args[8]

# Read in null model object
null_model <- readRDS(null_modelfile)

# Parse exposures
exposures <- strsplit(exposure_names, split=" ")[[1]]
int_covars <- if (int_covar_names == "") NULL else strsplit(int_covar_names, split=" ")[[1]]

# Apply gds filter if provided 
# (currently only string matches in "annotation/filter")
if (gds_filter == "") {  # No filter -> pass filename
  gds <- gdsfile
} else {  # Yes filter -> pass filtered GDS object
  keep_filter <- strsplit(gds_filter, split=" ")[[1]]
  gds <- SeqArray::seqOpen(gdsfile)
  keep_idx <- SeqArray::seqGetData(gds, "annotation/filter") %in% keep_filter
  SeqArray::seqSetFilter(gds, variant.sel=keep_idx)
}

# Run GWIS
glmm.gei(null_model, interaction=exposures, interaction.covariates=int_covars,
		 geno.file=gds, ncores=ncores,
	 outfile="magee_res", MAF.range=c(min_MAF, max_MAF), 
	 miss.cutoff=0.05, meta.output=TRUE)
