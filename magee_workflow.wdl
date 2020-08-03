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

	command {
		Rscript /MAGEE_null_model.R ${phenofile} ${sample_id_header} ${outcome} ${binary_outcome} "${exposure_names}" "${covar_names}" ${delimiter} ${missing} "${kinsfile}"
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
	String exposure_names
	File gdsfile
	File groupfile
	Int ncores
	Int memory
	Int disk
	Int monitoring_freq
	
	command {
		dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
		atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &

		Rscript /MAGEE_GWIS.R ${null_modelfile} "${exposure_names}" ${gdsfile} ${groupfile} ${ncores}
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
	File groupfile
	Int? ncores = 1
	Int? memory = 10
	Int? disk = 50
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

	scatter (gdsfile in gdsfiles) {
		call run_gwis {
			input:
				null_modelfile = null_modelfile,
				exposure_names = exposure_names,
				gdsfile = gdsfile,
				groupfile = groupfile,
				ncores = ncores,
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
		File? magee_null_model = null_modelfile
		File magee_results = cat_results.all_results
		Array[File] system_resource_usage = run_gwis.system_resource_usage
		Array[File] process_resource_usage = run_gwis.process_resource_usage
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
		groupfile: "Path to variant group definition file. File should be tab-separated with the following fields: variant set, chromosome, position, reference allele, alternate allele, weight."
		ncores: "Optional number of compute cores to be used for multi-threading during the genome-wide scan (default = 1)."
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
