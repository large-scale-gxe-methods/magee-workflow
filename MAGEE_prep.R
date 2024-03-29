library(readr)
library(MAGEE)

Sys.setenv(MKL_NUM_THREADS=1)

# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure_names <- args[2]
int_covar_names <- args[3]
gdsfile <- args[4]
groupfile <- args[5]
gds_filter <- args[6]

# Read in null model object
null_model <- readRDS(null_modelfile)

# Parse exposures
exposures <- strsplit(exposure_names, split=" ")[[1]]
int_covars <- if (int_covar_names == "") NULL else strsplit(int_covar_names, split=" ")[[1]]

# Remove header from group file if necessary
if (groupfile != "") {
  if (grepl("group", readLines(groupfile, n=1))) {
    system(paste0("tail -n +2 ", groupfile, " > groupfile.tmp"))
    groupfile <- "groupfile.tmp"
  }
}

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

# Prep GWIS
prep <- MAGEE.prep(null_model, interaction=exposures, gds, groupfile,
				   interaction.covariates=int_covars)
saveRDS(prep, "magee_prep.rds")
