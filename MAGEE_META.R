#Load MAGEE CRAN library and set up the environment
library(MAGEE)

#Parse arguments
args <- commandArgs(T)
print(args)
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
files<- args[15]
Sys.setenv(MKL_NUM_THREADS = threads)

if(cohort.group.idx=="NULL") {
  cohort.group.idx=NULL
} else {
  cohort.group.idx <- unlist(strsplit(cohort.group.idx, split=" "))
}

if(interaction=="") {
  interaction=NULL
}

if(E.match=="NULL") {
  E.match=NULL
} 
if(group_file_sep=="NULL") {
  group_file_sep="\t"
} 
metafile_vec <- unlist(strsplit(metafiles, split=" "))
files_vec <- unlist(strsplit(files, split=" "))
path=unlist(strsplit(files_vec[1],split='/'))
path=paste(path[-length(path)],collapse='/')
n.files=c()
for(i in 1:length(metafile_vec))
{
  n.files[i]=length(list.files(path=path,pattern=metafile_vec[i]))/2
}
metafile_vec=paste0(path,'/',metafile_vec)

  tests <- unlist(strsplit(tests, split=" "))
  MAF.range <- as.numeric(unlist(strsplit(MAF.range, split=" ")))
  MAF.weights.beta <- as.numeric(unlist(strsplit(MAF.weights.beta, split=" ")))
  use.minor.allele=as.logical(use.minor.allele)
  

#Run meta analysis
if(is.null(group_file))
{
output=glmm.gei.meta(metafile_vec,
              interaction=interaction,
              outfile=out)
write.table(output,out)
} else {
  

           output=MAGEE.meta(metafile_vec,
           cohort.group.idx=cohort.group.idx,
           group.file=group_file,
           MAF.range=MAF.range,
           MAF.weights.beta=MAF.weights.beta,
           miss.cutoff=miss.cutoff,
           method=method,
           tests=tests,
           use.minor.allele=use.minor.allele,
           group.file.sep=group_file_sep,
           E.match=E.match
           )
                   
} 
