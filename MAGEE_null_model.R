library(readr)
library(MAGEE)

# Parse arguments
args <- commandArgs(T)

phenofile <- args[1]
sample_id_header <- args[2]
outcome <- args[3]
binary_outcome <- as.logical(args[4])
exposure_names <- args[5]
covar_names <- args[6]
delimiter <- args[7]
missing <- args[8]
grmfile <- args[9]

# Read in phenotypes
phenos <- as.data.frame(read_delim(phenofile, delim=delimiter, na=missing))

# Null model formula
exposures <- strsplit(exposure_names, split=" ")[[1]]
covars <- strsplit(covar_names, split=" ")[[1]]
null_model_str <- paste0(outcome, "~", paste(c(exposures, covars), collapse="+"))

# Read in GRM
GRM <- readRDS(grmfile)

# Fit null model
family <- if (binary_outcome) binomial(link="logit") else "linear"
model0 <- GMMAT::glmmkin(as.formula(null_model_str), data=phenos, 
		  	 kins=GRM, id=sample_id_header, family=family)
saveRDS(model0, "null_model.rds")
