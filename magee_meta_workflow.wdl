workflow run_MAGEE_meta {
  
  Array[File] metafiles
  Boolean? use_minor_allele = false
  String? interaction
  String? out_prefix = "meta_output"
  File group_file
  String? group_file_sep = "NULL"
  String? E_match = "NULL"
  String? MAF_range = "1e-7 0.5"
  String? MAF_weights_beta = "1 25"
  Int? miss_cutoff = 1
  String? method = "davies"
  String? tests = "JF"
  String? cohort_group_idx = "NULL"
  Int? threads = 1
  Int? memory = 10
  Int? disk = 50
  Int? preemptible = 0
  Int? monitoring_freq = 1
  String prefix

  call run_meta {
    input:
      metafiles = metafiles,
      interaction = interaction,
      out_prefix= out_prefix,
      group_file = group_file,
      group_file_sep = group_file_sep,
      E_match = E_match,
      MAF_range = MAF_range,
      MAF_weights_beta = MAF_weights_beta,
      miss_cutoff = miss_cutoff,
      method = method,
      tests = tests,
      use_minor_allele = use_minor_allele,
      cohort_group_idx = cohort_group_idx,
      threads= threads,
      memory= memory,
      disk= disk,
      preemptible= preemptible,
      monitoring_freq= monitoring_freq,
      prefix=prefix
  }
  
  output {
    File? process_resource_usage = run_meta.process_resource_usage
  }
  
  parameter_meta {
    prefix: "String to use as the prefix for the proper summary statistics and covariance matrix for meta-analysis. If excluded, this file will not be created. Default is 'meta'"
    metafiles: "Prefixes to meta-analysis output files from the run_MAGEE workflow containing summary statistics and covariance matrix for meta-analysis. The names of the files are set by the 'meta_file_prefix' parameter."
    group_file: "For variant group test: A plain text file with 6 columns defining the test units for aggregate meta analysis. There should be no headers in the file, and the columns are group name, chromosome, position, reference allele, alternative allele and weight, respectively. This file should be matched to the genetic variants in both cohorts. If this field is not provided, single-variant meta-analysis will be performed."
    group_file_sep: "Delimiter of the group file."
    threads: "Optional number of compute cores to be requested and used for multi-threading during the genome-wide scan (default = 1)."
    cohort_group_idx: "For variant group test: String with the same length as metafiles, indicating which cohorts should be grouped together in the meta-analysis assuming homogeneous genetic effects. For example, 'a b a a b' means cohorts 1, 3, 4 are assumed to have homogeneous genetic effects, and cohorts 2, 5 are in another group with homogeneous genetic effects (but possibly heterogeneous with group 'a'). By default, all cohorts are in the same group."
    E_match: "For variant group test: Environmental factors that match the interactions (default = tab)"
    miss_cutoff: "Maximum missing rate allowed for a variant to be included (default = 1, including all variants). Filter applied to the combined samples."
    method: "For variant group test: Method to compute p-values for SKAT-type test statistics (default = 'davies'). 'davies' represents an exact method that computes a p-value by inverting the characteristic function of the mixture chisq distribution, with an accuracy of 1e-6. When 'davies' p-value is less than 1e-5, it defaults to method 'kuonen'. 'kuonen' represents a saddlepoint approximation method that computes the tail probabilities of the mixture chisq distribution. When 'kuonen' fails to compute a p-value, it defaults to method 'liu'. 'liu' is a moment-matching approximation method for the mixture chisq distribution."
    tests: "For variant group test: String vector indicating which MAGEE tests should be performed ('MV' for the main effect variance component test, MF for the main effect combined test of the burden and variance component tests using Fisher’s method, IV for the interaction variance component test, IF for the interaction combined test of the burden and variance component tests using Fisher’s method, JV for the joint variance component test for main effect and interaction, JF for the joint combined test of the burden and variance component tests for main effect and interaction using Fisher’s method, or JD for the joint combined test of the burden and variance component tests for main effect and interaction using double Fisher’s method.). The MV and IV test are automatically included when performing JV, and the MF and IF test are automatically included when performing JF or JD (default = JF)."
    use_minor_allele: "Boolean defining whether to use the alt allele as the coding allele (default = FALSE)."
    MAF_range: "Specifies the minimum and maximum MAFs for a variant to be included in the analysis. By default the minimum MAF is 1 × 10−7 and the maximum MAF is 0.5, meaning that only monomorphic markers in the sample will be excluded (if your sample size is no more than 5 million)."
    out_prefix: "Optional string to use as the file name prefix for writing the meta-analysis output file. Default is 'meta_output'"
    interaction: "For single-variant test: Numeric or a character vector indicating the environmental factors. If a numeric vector, it represents which indices in the order of covariates are the environmental factors; if a character vector, it represents the variable names of the environmental factors."
    memory: "Requested memory (in GB)."
    disk: "Requested disk space (in GB)."
    preemptible: "Optional number of attempts using a preemptible machine from Google Cloud prior to falling back to a standard machine (default = 0, i.e., don't use preemptible)."
    monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
  }
  
  meta {
    author: "Chris Bryan"
    email: "cjbryan@mgh.harvard.edu"
    description: "Run aggregate and single-variant meta-analysis using the MAGEE package and return a table of summary statistics for K-DF interaction and (K+1)-DF joint tests."
  }
}

task run_meta {
  Array[File]? metafiles
  String? out_prefix
  String? interaction
  File? group_file
  String? group_file_sep
  String? E_match
  String? MAF_range
  String? MAF_weights_beta
  Int? miss_cutoff
  String? method
  String? tests
  Boolean? use_minor_allele
  String? cohort_group_idx
  Int threads
  Int memory
  Int disk
  Int preemptible
  Int monitoring_freq
  String prefix
  
  command <<<
    dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
    atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &
    
    Rscript /MAGEE_META.R "${prefix}" "${interaction}" "${out_prefix}" "${threads}" "${group_file}" "${group_file_sep}" "${E_match}" "${MAF_range}" "${MAF_weights_beta}" "${miss_cutoff}" "${method}" "${tests}" "${use_minor_allele}" "${cohort_group_idx}" "${sep=' ' metafiles}"
  >>>
    
    runtime {
      docker: "getting-started"
      memory: "${memory} GB"
      cpu: "${threads}"
      disks: "local-disk ${disk} HDD"
      preemptible: "${preemptible}"
      maxRetries: 2
    }
  
  output {
    File system_resource_usage = "system_resource_usage.log"
    File process_resource_usage = "process_resource_usage.log"
  }
}
