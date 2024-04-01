#Load MAGEE CRAN library and set up the environment
library(MAGEE)

#Parse arguments
args <- commandArgs(T)

metafiles <- args[1]
interaction <- args[2]
out <- args[3]
threads <- args[4]
group_file <- args[5]
group_file_sep <- args[6]
E.match<- args[7]
MAF.range<- args[8]
MAF.weights.beta<- args[9]
miss.cutoff<- args[10]
method<- args[11]
tests<- args[12]
use.minor.allele<- args[13]
cohort.group.idx<- args[14]
metafile_prefix<- args[15]

Sys.setenv(MKL_NUM_THREADS = threads)

if(cohort.group.idx=="") {
  cohort.group.idx=NULL
} else {
  cohort.group.idx <- unlist(strsplit(cohort.group.idx, split=" "))
}

if(metafile_prefix=="") {
  metafile_prefix=NULL
} else {
  metafile_prefix <- unlist(strsplit(metafile_prefix, split=" "))
}

if(interaction=="") {
  interaction=NULL
}

if(E.match=="") {
  E.match=NULL
}

if(group_file=="") {
  group_file=NULL
}

if(group_file_sep=="") {
  group_file_sep="\t"
}

is_sv=is.null(group_file)

metafile_vec <- unlist(strsplit(metafiles, split=" "))
n.files=c()
if(!is_sv)
{
      for(i in 1:length(metafile_vec))
      {
        system(paste0('tar xvf ', metafile_vec[i]))
        n.files[i]=length(list.files(pattern=metafile_prefix[i]))/2
      }
  metafile_prefix <- unlist(strsplit(metafile_prefix, split=" "))
}

tests <- unlist(strsplit(tests, split=" "))
MAF.range <- as.numeric(unlist(strsplit(MAF.range, split=" ")))
MAF.weights.beta <- as.numeric(unlist(strsplit(MAF.weights.beta, split=" ")))
use.minor.allele=as.logical(use.minor.allele)

#Run meta analysis
if(is_sv)
{
for(file in metafile_vec)
{
  library(data.table)
  temp=fread(file)
  colnames(temp)=gsub('\\-','\\.',colnames(temp))
  fwrite(temp,file)
}
output=glmm.gei.meta(metafile_vec,
              interaction=interaction,
              outfile=out)
} else {
           output=MAGEE.meta(
           cohort.group.idx=cohort.group.idx,
           meta.files.prefix=metafile_prefix,
           group.file=group_file,
           MAF.range=MAF.range,
           MAF.weights.beta=MAF.weights.beta,
           miss.cutoff=miss.cutoff,
           method=method,
           tests=tests,
           use.minor.allele=use.minor.allele,
           group.file.sep=group_file_sep,
           E.match=E.match)
  write.table(output,out)
}
