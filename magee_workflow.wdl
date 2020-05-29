task run_null_model {

	File phenofile
	String sample_id_header
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String delimiter
	String missing
	File grmfile
	Int memory
	Int disk

	command {
		Rscript /MAGEE_null_model.R ${phenofile} ${sample_id_header} ${outcome} ${binary_outcome} ${exposure_names} "${covar_names}" ${delimiter} ${missing} ${grmfile}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
		File null_model = "null_model.rds"
	}
}

task run_gwis {

	File null_modelfile
	String exposure
	File gdsfile
	File groupfile
	Int memory
	Int disk
	Int monitoring_freq
	
	command {
		dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
		atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &

		Rscript /MAGEE_GWIS.R ${null_modelfile} ${exposure} ${gdsfile} ${groupfile}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/magee-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
		File res = "magee_res"
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
	}
}

task cat_results {

	Array[File] results_array

	command {
		head -1 ${results_array[0]} > all_results.txt && \
			for res in ${sep=" " results_array}; do tail -n +2 $res >> all_results.txt; done
	}
	
	runtime {
		docker: "ubuntu:latest"
		disks: "local-disk 10 HDD"
	}
	output {
		File all_results = "all_results.txt"
	}
}

workflow MAGEE {

	Array[File] gdsfiles
	File phenofile
	String? sample_id_header = "sampleID"
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names
	String? delimiter = ","
	String? missing = "NA"
	File grmfile
	File groupfile
	Int? memory = 10
	Int? disk = 50
	Int? monitoring_freq = 1

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
			grmfile = grmfile,
			memory = memory,
			disk = disk
	}

	scatter (gdsfile in gdsfiles) {
		call run_gwis {
			input:
				null_modelfile = run_null_model.null_model,
				exposure = exposure_names,
				gdsfile = gdsfile,
				groupfile = groupfile,
				memory = memory,
				disk = disk,
				monitoring_freq = monitoring_freq
		}
	}

	call cat_results {
		input:
			results_array = run_gwis.res
	}

	output {
		File results = cat_results.all_results
		#Array[File] system_resource_usage = run_tests.system_resource_usage
		#Array[File] process_resource_usage = run_tests.process_resource_usage
	}

	parameter_meta {
		gdsfiles: "Array of genotype filepaths in .gds format."
		phenofile: "Phenotype filepath."	
		sample_id_header: "Optional column header name of sample ID in phenotype file."
		outcome: "Column header name of phenotype data in phenotype file."
                binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed."
		exposure_names: "Column header name(s) of the exposures for genotype interaction testing (space-delimited)."
		covar_names: "Column header name(s) of any covariates for which only main effects should be included selected covariates in the pheno data file (space-delimited). This set should not contain the exposure."
		delimiter: "Delimiter used in the phenotype file."
		missing: "Missing value key of phenotype file."
		grmfile: "Path to file containing kinship matrix (stored as binary .rds R object)."
		groupfile: "Path to variant group definition file. File should be tab-separated with the following fields: variant set, chromosome, position, reference allele, alternate allele, weight."
		memory: "Requested memory (in GB)."
		disk: "Requested disk space (in GB)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

        meta {
                author: "Kenny Westerman"
                email: "kewesterman@mgh.harvard.edu"
                description: "Run interaction tests using MAGEE and return a table of summary statistics for K-DF interaction and (K+1)-DF joint tests."
        }
}
