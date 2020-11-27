FROM quay.io/large-scale-gxe-methods/bioconductor_docker:RELEASE_3_11

# Dependencies
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("SeqArray", "SeqVarTools", "foreach", "GMMAT", "CompQuadForm"))'
RUN Rscript -e 'install.packages(c("readr", "doMC", "devtools"), repos="http://cran.us.r-project.org")'

# Install MAGEE R package
ADD https://api.github.com/repos/xwang21/MAGEE/git/refs/heads version.json
COPY MAGEE_dev.tar.gz /MAGEE_dev.tar.gz
RUN Rscript -e 'install.packages("MAGEE_dev.tar.gz", repos=NULL, type="source")'

# Copy in R scripts
COPY MAGEE_null_model.R MAGEE_GWIS.R /

# Install tools for monitoring and resource tracking
RUN apt-get update && apt-get -y install dstat atop
