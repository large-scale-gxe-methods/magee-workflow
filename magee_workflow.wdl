task run_null_model {

	File phenofile
	String sample_id_header
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String delimiter
	String missing
	File? kinsfile
	Int memory
	Int disk

	command <<<
		Rscript /MAGEE_null_model.R ${phenofile} ${sample_id_header} ${outcome} ${binary_outcome} "${exposure_names}" "${covar_names}" ${delimiter} ${missing} "${kinsfile}"
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
		File null_model = "null_model.rds"
	}
}

task run_gwis_agg {

	File null_modelfile
	String exposure_names
	File gdsfile
	File groupfile
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

		Rscript /MAGEE_GWIS.R ${null_modelfile} "${exposure_names}" ${gdsfile} ${groupfile} ${min_MAF} ${max_MAF} ${threads}
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow"
		memory: "${memory} GB"
		cpu: "${threads}"
		disks: "local-disk ${disk} HDD"
		preemptible: "${preemptible}"
	}

	output {
		File res = "magee_res"
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
	}
}


task run_gwis_sv {

	File null_modelfile
	String exposure_names
	File gdsfile
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

		Rscript /MAGEE_GWIS.R ${null_modelfile} "${exposure_names}" ${gdsfile} "" ${min_MAF} ${max_MAF} ${threads}
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow"
		memory: "${memory} GB"
		cpu: "${threads}"
		disks: "local-disk ${disk} HDD"
		preemptible: "${preemptible}"
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
		disks: "local-disk 10 HDD"
	}
	output {
		File all_results = "all_results.txt"
	}
}

workflow MAGEE {

	File phenofile
	String? sample_id_header = "sampleID"
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String? delimiter = ","
	String? missing = "NA"
	File? kinsfile
	File? null_modelfile_input
	Array[File] gdsfiles
	Array[File]? groupfiles
	Float? min_MAF = 0.0000007
	Float? max_MAF = 0.5
	Int? threads = 1
	Int? memory = 10
	Int? disk = 50
	Int? preemptible = 0
	Int? monitoring_freq = 1

	#Int null_memory = 2 * ceil(size(kinsfile, "GB")) + 5

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
				kinsfile = kinsfile,
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
					groupfile = select_first([groupfiles])[i],  # select_first mechanism allows indexing of an optional array
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

	if (run_type == "sv") {  
		scatter (i in range(length(gdsfiles))) {
			call run_gwis_sv {
				input:
					null_modelfile = null_modelfile,
					exposure_names = exposure_names,
					gdsfile = gdsfiles[i],
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
		kinsfile: "Optional path to file containing GRM/kinship matrix with sample IDs as the row and column names. Can be either a .rds file storing a matrix object or a .csv file. If excluded, a the null model will be fit as a GLM with no random effects."
		null_modelfile: "Optional path to file containing the pre-fitted null model in .rds format."
		gdsfiles: "Array of genotype filepaths in .gds format."
		groupfiles: "Optional array of variant group definition filepaths. Files should be tab-separated with the following fields: variant set, chromosome, position, reference allele, alternate allele, weight. If this field is not provided, single-variant tests will be performed."
		min_MAF: "Optional cutoff for the minimum minor allele frequency for variant inclusion (inclusive; default = 1e-7)."
		max_MAF: "Optional cutoff for the maximum minor allele frequency for variant inclusion (inclusive; default = 0.5)."
		threads: "Optional number of compute cores to be requested and used for multi-threading during the genome-wide scan (default = 1)."
		memory: "Requested memory (in GB)."
		disk: "Requested disk space (in GB)."
		preemptible: "Optional number of attempts using a preemptible machine from Google Cloud prior to falling back to a standard machine (default = 0, i.e., don't use preemptible)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

        meta {
                author: "Kenny Westerman"
                email: "kewesterman@mgh.harvard.edu"
                description: "Run aggregate and single-variant interaction tests using MAGEE and return a table of summary statistics for K-DF interaction and (K+1)-DF joint tests."
        }
}
