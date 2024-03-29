library(readr)


# Parse arguments
args <- commandArgs(T)

phenofile <- args[1]
sample_id_header <- args[2]
outcome <- args[3]
binary_outcome <- as.logical(args[4])
exposure_names <- args[5]
int_covar_names <- args[6]
covar_names <- args[7]
delimiter <- args[8]
missing <- args[9]
kinsfiles <- args[10]
groups <- args[11]

# Null model formula
exposures <- strsplit(exposure_names, split=" ")[[1]]
int_covars <- strsplit(int_covar_names, split=" ")[[1]]
covars <- strsplit(covar_names, split=" ")[[1]]
null_model_str <- paste0(outcome, "~", paste(c(exposures, int_covars, covars), collapse="+"))

# Read in phenotypes
phenos <- as.data.frame(read_delim(phenofile, delim=delimiter, na=missing))

# Read in kinship matrix
get_kinship_matrix <- function(kinsfile) {
  if (kinsfile == "") {
    NULL
  } else if (grepl('rds', kinsfile, ignore.case=T)) {
    readRDS(kinsfile)
  } else if (grepl('RData', kinsfile, ignore.case=T)) {
    get(load(kinsfile))
  } else {
    as.matrix(read.csv(kinsfile, as.is=T, check.names=F, row.names=1))
  }
}
kinsfile_vec <- strsplit(kinsfiles, split=" ")[[1]]
k_mats <- if (length(kinsfile_vec) == 0) NULL else lapply(kinsfile_vec, get_kinship_matrix)

# Define grouping variable to allow heteroscedastic LMM
groups <- if (groups != "") groups else NULL

# Choose regression family
family <- if (binary_outcome) binomial(link="logit") else gaussian(link="identity")

# Fit null model
model0 <- GMMAT::glmmkin(as.formula(null_model_str), data=phenos, 
		  	 kins=k_mats, id=sample_id_header, family=family,
			 groups=groups)
saveRDS(model0, file="null_model.rds")
