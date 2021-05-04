FROM quay.io/large-scale-gxe-methods/bioconductor_docker:RELEASE_3_11

# Dependencies
RUN Rscript -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("SeqArray", "SeqVarTools", "foreach", "GMMAT", "CompQuadForm"))'
RUN Rscript -e 'install.packages(c("readr", "doMC", "devtools"), repos="http://cran.us.r-project.org")'

# Install MAGEE R package
ADD https://api.github.com/repos/xwang21/MAGEE/git/refs/heads version.json
RUN Rscript -e 'devtools::install_github("https://github.com/xwang21/MAGEE/tree/dev")'

# Copy in R scripts
COPY MAGEE_null_model.R MAGEE_prep.R MAGEE_GWIS.R MAGEE_SV_GWIS.R /

# Install tools for monitoring and resource tracking
RUN apt-get update && apt-get -y install dstat atop
