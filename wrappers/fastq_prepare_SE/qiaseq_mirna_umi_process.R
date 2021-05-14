library(data.table)

run_all <- function(fastq_input,output_file){
  print("reading file")
  a <- Sys.time()
  fastq <- readLines(fastq_input,400000)
     
  print(Sys.time() - a)
  print("processing umi")
  split <- strsplit(fastq[seq(2,length(fastq),by = 4)],"AACTGTAGGCACCATCAAT")   
  split_lengths <- sapply(split,length)
  split <- unlist(split)
  seqs <- character(length(split_lengths))
  seqs[split_lengths == 1] <- split[cumsum(split_lengths)[split_lengths == 1]]
  seqs[split_lengths == 2] <- split[cumsum(split_lengths)[split_lengths == 2] - 1]
  umi <- character(length(split_lengths))
  umi[split_lengths == 2] <- stringi::stri_sub(split[cumsum(split_lengths)[split_lengths == 2]],1,12)
  cat(paste0(fastq[seq(1,length(fastq),by = 4)],"_",umi,"\n",seqs,"\n",
             fastq[seq(3,length(fastq),by = 4)],"\n",stringi::stri_sub(fastq[seq(4,length(fastq),by = 4)],1,nchar(seqs)),"\n"),file = output_file)
  print(Sys.time() - a)               
}

# develop and test
fastq_input <- "/mnt/ssd/ssd_1/snakemake/raw_fastq_temp/190213_2019-02-13-IVDplasmaMOII_AH3H3TBGX9/DV7511_S11_R1_001.fastq"
output_file <- "/mnt/ssd/ssd_1/snakemake/Ondrej_Slaby/sequencing_results/primary_data/190213_2019-02-13-IVDplasmaMOII/DV7511.fastq"
run_all(fastq_input,output_file)

# run as Rscript
# args <- commandArgs(trailingOnly = T)
# fastq_input <- args[1]
# output_file <- args[2]
# run_all(fastq_input,output_file)