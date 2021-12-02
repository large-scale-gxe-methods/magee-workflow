Implementation of MAGEE workflow for efficient, large-scale aggregate and single-variant tests of gene-environment interaction in related individuals.

This workflow implements the MAGEE (Mixed Model Association Test for GEne-Environment Interaction)  tool (https://github.com/large-scale-gxe-methods/MAGEE). MAGEE conducts genome-wide gene-environment interaction tests of single variants and aggregate tests of variant groups, while allowing for related individuals and inclusion of multiple exposures.

Author: Kenny Westerman (kewesterman@mgh.harvard.edu)

MAGEE tool information:

* Manuscript:  X Wang et al. Efficient gene–environment interaction tests for large biobank-scale sequencing studies. *Genetic Epidemiology*. 2020; 44(8): 908-923. https://doi.org/10.1002/gepi.22351.
* Source code: https://github.com/large-scale-gxe-methods/MAGEE

Workflow steps:

* Run null model
* Run MAGEE (single-variant or aggregate; scattered across an array of input files, usually chromosomes)
* Concatenate the outputs into a single summary statistics file

Inputs: 

See the "parameter_meta" section of the .wdl script.

Outputs:

* A summary statistics file containing estimates for genetic main effects and interaction effects as well as p-values for these along with a joint test of genetic main and interaction effects.
* A file containing additional score statistics and covariances necessary for meta-analysis.
