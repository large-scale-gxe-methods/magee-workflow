FROM bioconductor/bioconductor_docker:devel
# Dependencies
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("SeqArray", "SeqVarTools", "foreach", "GMMAT", "CompQuadForm"))'
RUN Rscript -e 'install.packages(c("MAGEE", "doMC", "readr"), repos="http://cran.us.r-project.org")'


# Copy in R scripts
COPY MAGEE_null_model.R MAGEE_prep.R MAGEE_GWIS.R MAGEE_SV_GWIS.R MAGEE_META.R /

# Install tools for monitoring and resource tracking
RUN apt-get update && apt-get -y install dstat atop
