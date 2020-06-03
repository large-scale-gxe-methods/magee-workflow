FROM bioconductor/bioconductor_docker:RELEASE_3_11

# Dependencies
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c("SeqArray", "SeqVarTools", "foreach", "GMMAT", "CompQuadForm"))'
RUN Rscript -e 'install.packages("readr")'

# Install MAGEE R package
COPY MAGEE_0.1.1.tar.gz /MAGEE_0.1.1.tar.gz
RUN Rscript -e 'install.packages("MAGEE_0.1.1.tar.gz", repos=NULL, type="source")'

# Copy in R scripts
COPY MAGEE_null_model.R MAGEE_GWIS.R /
