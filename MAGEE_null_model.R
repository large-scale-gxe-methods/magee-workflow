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
kinsfile <- args[9]

# Null model formula
exposures <- strsplit(exposure_names, split=" ")[[1]]
covars <- strsplit(covar_names, split=" ")[[1]]
null_model_str <- paste0(outcome, "~", paste(c(exposures, covars), collapse="+"))

# Read in phenotypes
phenos <- as.data.frame(read_delim(phenofile, delim=delimiter, na=missing))

# Read in kinship matrix
get_kinship_matrix <- function(kinsfile) {
  if (kinsfile == "") {
    NULL
  } else if (grepl('rds', kinsfile, ignore.case=T)) {
    readRDS(kinsfile)
  } else {
    as.matrix(read.csv(kinsfile, as.is=T, check.names=F, row.names=1))
  }
}
k_mat <- get_kinship_matrix(kinsfile)

# Choose regression family
family <- if (binary_outcome) binomial(link="logit") else gaussian(link="identity")

# Fit null model
model0 <- GMMAT::glmmkin(as.formula(null_model_str), data=phenos, 
		  	 kins=k_mat, id=sample_id_header, family=family)
saveRDS(model0, file="null_model.rds")
