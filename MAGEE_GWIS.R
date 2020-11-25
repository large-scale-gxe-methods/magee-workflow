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
gds_filter <- args[8]

# Read in null model object
null_model <- readRDS(null_modelfile)

# Parse exposures
exposures <- strsplit(exposure_names, split=" ")[[1]]

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
  keep_idx <- SeqArray::seqGetData(gds, "annotation/info") %in% keep_filter
  SeqArray::seqSetFilter(gds, variant.sel=keep_idx)
}

# Run GWIS
if (groupfile == "") {
  glmm.gei(null_model, interaction=exposures, geno.file=gds, ncores=ncores,
	   outfile="magee_res", MAF.range=c(min_MAF, max_MAF), miss.cutoff=0.05)
} else {
  res <- MAGEE(null_model, interaction=exposures, gds, groupfile, 
  	       MAF.range=c(min_MAF, max_MAF), miss.cutoff=0.05, 
	       tests=c("JV", "JF", "JD"), ncores=ncores)
  write_delim(res, "magee_res", delim=" ")
}
