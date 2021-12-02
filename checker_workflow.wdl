# version 1.0.1a

import "magee_workflow.wdl" as magee_wf

workflow checker {

	File phenofile
	String sample_id_header
	String outcome
	Boolean binary_outcome
	String exposure_names
	String covar_names
	Array[File] kinsfiles
	Array[File] gdsfiles
	Array[File] groupfiles
	Float max_MAF
	File expected_sumstats

  call magee_wf.run_MAGEE {
    input:
			phenofile = phenofile,
			sample_id_header = sample_id_header,
			outcome = outcome,
			binary_outcome = binary_outcome,
			exposure_names = exposure_names,
			covar_names = covar_names,
			kinsfiles = kinsfiles,
			gdsfiles = gdsfiles,
			groupfiles = groupfiles,
			max_MAF = max_MAF
  }

  call md5sum {
    input:
			sumstats = run_MAGEE.magee_results,
			expected_sumstats = expected_sumstats
  }

  meta {
    author: "Kenny Westerman"
    email: "kewesterman@mgh.harvard.edu"
  }

}


task md5sum {

	File sumstats
	File expected_sumstats

  command <<<

  md5sum ${sumstats} > sum.txt
  md5sum ${expected_sumstats} > expected_sum.txt

  # temporarily outputting to stderr for clarity's sake
  >&2 echo "Output checksum:"
  >&2 cat sum.txt
  >&2 echo "-=-=-=-=-=-=-=-=-=-"
  >&2 echo "Truth checksum:"
  >&2 cat expected_sum.txt
  >&2 echo "-=-=-=-=-=-=-=-=-=-"
  >&2 echo "Head of the output file:"
  >&2 head ${sumstats}
  >&2 echo "-=-=-=-=-=-=-=-=-=-"
  >&2 echo "Head of the truth file:"
  >&2 head ${expected_sumstats}

  echo "$(cut -f1 -d' ' expected_sum.txt) ${sumstats}" | md5sum --check

  >>>

  runtime {
    docker: "quay.io/large-scale-gxe-methods/ubuntu:focal-20210325"
  }

}


