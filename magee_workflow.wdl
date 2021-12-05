workflow run_MAGEE {

	File phenofile
	String? sample_id_header = "sampleID"
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String? delimiter = ","
	String? missing = "NA"
	Array[File]? kinsfiles
	String? var_group
	File? null_modelfile_input
	Array[File] gdsfiles
	String? gds_filter
	Array[File]? groupfiles
	Float? min_MAF = 0.0000001
	Float? max_MAF = 0.5
	String? meta_file_prefix = "none"
	Int? threads = 1
	Int? memory = 10
	Int? disk = 50
	Int? preemptible = 0
	Int? monitoring_freq = 1

	#Int null_memory = 2 * ceil(size(kinsfiles, "GB")) + 5

	if (!defined(null_modelfile_input)) {
		call run_null_model {
			input:
				phenofile = phenofile,
				sample_id_header = sample_id_header,
				outcome = outcome,
				binary_outcome = binary_outcome,
				exposure_names = exposure_names,
				covar_names = covar_names,
				delimiter = delimiter,
				missing = missing,
				kinsfiles = kinsfiles,
				var_group = var_group,
				memory = memory,
				disk = disk
		}
	}

	File? null_modelfile = if (defined(null_modelfile_input)) then null_modelfile_input else run_null_model.null_model
	String run_type = if defined(groupfiles) then "agg" else "sv"  # If groupfiles is defined, run aggregate tests, otherwise run single-variant tests

	if (run_type == "agg") {  
		scatter (i in range(length(gdsfiles))) {
			call run_gwis_agg {
				input:
					null_modelfile = null_modelfile,
					exposure_names = exposure_names,
					gdsfile = gdsfiles[i],
					gds_filter = gds_filter,
					groupfile = select_first([groupfiles])[i],  # select_first mechanism allows indexing of an optional array
					min_MAF = min_MAF,
					max_MAF = max_MAF,
					meta_file_prefix = meta_file_prefix,
					threads = threads,
					memory = memory,
					disk = disk,
					preemptible = preemptible,
					monitoring_freq = monitoring_freq
			}
		}
	}

	if (run_type == "sv") {  
		scatter (i in range(length(gdsfiles))) {
			call run_gwis_sv {
				input:
					null_modelfile = null_modelfile,
					exposure_names = exposure_names,
					gdsfile = gdsfiles[i],
					gds_filter = gds_filter,
					min_MAF = min_MAF,
					max_MAF = max_MAF,
					threads = threads,
					memory = memory,
					disk = disk,
					preemptible = preemptible,
					monitoring_freq = monitoring_freq
			}
		}
	}

	call cat_results {
		input:
			results_array = if run_type == "agg" then run_gwis_agg.res else run_gwis_sv.res
	}
	
	output {
		File? magee_null_model = null_modelfile
		File magee_results = cat_results.all_results
		Array[Array[File]]? meta_analysis_output = if (run_type == "agg") then run_gwis_agg.meta_file else [[]]
		Array[File]? system_resource_usage = if (run_type == "agg") then run_gwis_agg.system_resource_usage else run_gwis_sv.system_resource_usage
		Array[File]? process_resource_usage = if (run_type == "agg") then run_gwis_agg.process_resource_usage else run_gwis_sv.process_resource_usage
	}

	parameter_meta {
		phenofile: "Phenotype filepath."	
		sample_id_header: "Optional column header name of sample ID in phenotype file."
		outcome: "Column header name of phenotype data in phenotype file." 
		binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed."
		exposure_names: "Column header name(s) of the exposures for genotype interaction testing (space-delimited)."
		covar_names: "Column header name(s) of any covariates for which only main effects should be included selected covariates in the pheno data file (space-delimited). This set should not contain the exposure."
		delimiter: "Delimiter used in the phenotype file."
		missing: "Missing value key of phenotype file."
		kinsfiles: "Optional array of paths to files containing relationship matrices to be used as random effects in the null model. Sample IDs should be used as the row and column names. This input will usually be a single GRM/kinship matrix. Can be either a .rds file storing a matrix object or a .csv file. If excluded, a the null model will be fit as a GLM with no random effects."
		var_group: "Optional string allowing for a heteroscedastic null linear mixed model. If provided, the null model will be fit with differential residual variances for each value of this variable. Note: this is only a valid input for continuous outcomes."
		null_modelfile: "Optional path to file containing the pre-fitted null model in .rds format."
		gdsfiles: "Array of genotype filepaths in .gds format."
		gds_filter: "Optional space-delimited string of values defining variant filters to be retained for analysis. Currently, this is matched against only the annotation/filter field in the .gds file."
		groupfiles: "Optional array of variant group definition filepaths. Files should be tab-separated with the following fields: variant set, chromosome, position, reference allele, alternate allele, weight. If this field is not provided, single-variant tests will be performed."
		min_MAF: "Optional cutoff for the minimum minor allele frequency for variant inclusion (inclusive; default = 1e-7)."
		max_MAF: "Optional cutoff for the maximum minor allele frequency for variant inclusion (inclusive; default = 0.5)."
		meta_file_prefix: "Optional string to use as the prefix for writing the the proper summary statistics and covariance matrix for meta-analysis. If excluded, this file will not be created."
		threads: "Optional number of compute cores to be requested and used for multi-threading during the genome-wide scan (default = 1)."
		memory: "Requested memory (in GB)."
		disk: "Requested disk space (in GB)."
		preemptible: "Optional number of attempts using a preemptible machine from Google Cloud prior to falling back to a standard machine (default = 0, i.e., don't use preemptible)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

	meta {
		author: "Kenny Westerman"
		email: "kewesterman@mgh.harvard.edu"
		description: "Run aggregate and single-variant interaction tests using the MAGEE package and return a table of summary statistics for K-DF interaction and (K+1)-DF joint tests."
	}
}


task run_null_model {

	File phenofile
	String sample_id_header
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String delimiter
	String missing
	Array[File]? kinsfiles
	String? var_group
	Int memory
	Int disk

	command <<<
		Rscript /MAGEE_null_model.R ${phenofile} ${sample_id_header} ${outcome} ${binary_outcome} "${exposure_names}" "${covar_names}" ${delimiter} ${missing} "${sep=' ' kinsfiles}" "${var_group}"
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow@sha256:18d7be67fcc3df07f677ee29cfda9b9a9e45545c3b2598423e224e210a44459e"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		maxRetries: 2
	}

	output {
		File null_model = "null_model.rds"
	}
}


task run_gwis_agg {

	File null_modelfile
	String exposure_names
	File gdsfile
	String? gds_filter
	File groupfile
	Float min_MAF
	Float max_MAF
	String meta_file_prefix
	Int threads
	Int memory
	Int disk
	Int preemptible
	Int monitoring_freq
	
	command <<<
		dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
		atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &

		Rscript /MAGEE_prep.R ${null_modelfile} "${exposure_names}" ${gdsfile} ${groupfile} "${gds_filter}"
		Rscript /MAGEE_GWIS.R ${gdsfile} ${min_MAF} ${max_MAF} ${threads} "${gds_filter}" "${meta_file_prefix}"
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow@sha256:18d7be67fcc3df07f677ee29cfda9b9a9e45545c3b2598423e224e210a44459e"
		memory: "${memory} GB"
		cpu: "${threads}"
		disks: "local-disk ${disk} HDD"
		preemptible: "${preemptible}"
		maxRetries: 2
	}

	output {
		File res = "magee_res"
		Array[File] meta_file = glob("${meta_file_prefix}.*")
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
	}
}


task run_gwis_sv {

	File null_modelfile
	String exposure_names
	File gdsfile
	String? gds_filter
	Float min_MAF
	Float max_MAF
	Int threads
	Int memory
	Int disk
	Int preemptible
	Int monitoring_freq
	
	command <<<
		dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
		atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &

		Rscript /MAGEE_SV_GWIS.R ${null_modelfile} "${exposure_names}" ${gdsfile} ${min_MAF} ${max_MAF} ${threads} "${gds_filter}"
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow@sha256:18d7be67fcc3df07f677ee29cfda9b9a9e45545c3b2598423e224e210a44459e"
		memory: "${memory} GB"
		cpu: "${threads}"
		disks: "local-disk ${disk} HDD"
		preemptible: "${preemptible}"
		maxRetries: 2
	}

	output {
		File res = "magee_res"
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
	}
}


task cat_results {

	Array[File] results_array

	command <<<
		head -1 ${results_array[0]} > all_results.txt && \
						for res in ${sep=" " results_array}; do tail -n +2 $res >> all_results.txt; done
	>>>
	
	runtime {
		docker: "ubuntu:latest"
		disks: "local-disk 25 HDD"
	}

	output {
		File all_results = "all_results.txt"
	}
}
