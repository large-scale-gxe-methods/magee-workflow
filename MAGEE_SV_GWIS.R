library(readr)
library(MAGEE)

Sys.setenv(MKL_NUM_THREADS=1)

# Parse arguments
args <- commandArgs(T)

gdsfile <- args[1]
min_MAF <- as.numeric(args[2])
max_MAF <- as.numeric(args[3])
ncores <- as.integer(args[4])
gds_filter <- args[5]
meta_file_prefix <- args[6]
use.minor.allele<- args[7]
AF.strata.range<- args[8]

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

# Establish meta-analysis file prefix if there is one
mfp <- if (!(meta_file_prefix == "none")) meta_file_prefix else NULL

# Run GWIS
prep <- readRDS("magee_prep.rds")  # From MAGEE_prep.R script
prep$geno.file <- gds  # Hack for now: replace GDS object with subsetted version
res <- MAGEE.lowmem(prep, MAF.range=c(min_MAF, max_MAF), miss.cutoff=0.05,
       		    tests=c("JV", "JF", "JD"), ncores=ncores, meta.file.prefix=mfp,use.minor.allele=use.minor.allele,AF.strata.range=AF.strata.range)
write_delim(res, "magee_res", delim=" ")
