library(readr)
library(MAGEE)

# Parse arguments
args <- commandArgs(T)

null_modelfile <- args[1]
exposure <- args[2]
gdsfile <- args[3]
groupfile <- args[4]

# Read in null model object
null_model <- get(load(null_modelfile))

# Run GWIS
res <- MAGEE(null_model, interaction=exposure, gdsfile, groupfile, 
	     tests=c("JV", "JF", "JD"))
write_delim(res, "magee_res", delim=" ")
