library(readr)
library(MAGEE)


# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure_names <- args[2]
int_covar_names <- args[3]
gdsfile <- args[4]
groupfile <- args[5]
gds_filter <- args[6]
meta_file_prefix <- args[7]
use.minor.allele<- args[8]
AF.strata.range<- args[9]
threads<- args[10]

Sys.setenv(MKL_NUM_THREADS=threads)

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

res=MAGEE(null.obj=null_model, interaction=exposures, geno.file=gdsfile, group.file=groupfile, group.file.sep = "\t",
      bgen.samplefile = NULL, interaction.covariates = int_covars, meta.file.prefix = meta_file_prefix,
      MAF.range = c(1e-7, 0.5), AF.strata.range = c(0, 1), MAF.weights.beta = c(1, 25),
      miss.cutoff = 1, missing.method = "impute2mean", method = "davies", tests = "JF",
      use.minor.allele = use.minor.allele, auto.flip = FALSE,
      Garbage.Collection = FALSE, is.dosage = FALSE, ncores = 1)
write_delim(res, "magee_res", delim=" ")
