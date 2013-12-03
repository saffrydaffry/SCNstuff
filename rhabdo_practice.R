#######
##  sequencing consensus sequence
#######
require(ggplot2)
require(stringr)
require(plyr)
require(seqinr)  ## contains codon table (DNA)
require(Biostrings)

##  Sequence Data
file.path <- "C:/Users/Mac/Dropbox/SCN-Thesis/Virus_Sequence/3_RhabroMuts.csv"
Rhabdo_dat <- read.csv(file.path)


##  Consensus Seq
genome.path <- "C:/Users/Mac/Dropbox/SCN-Thesis/Virus_Sequence/Rhabdo_seq.txt"
Consensus.seq <- read.table(genome.path,  colClasses = "character", stringsAsFactors = FALSE)[[1]][1]
Seq<- DNAString(Consensus.seq)

AA_string <-translate(Seq)

##  Identify ORFs
#### we know there are 5 in ScRV sequence
ORFs_num <- as.numeric(readline("Enter number of ORFs in the Sequence:"))
ORF_coord = matrix(rep(0), nrow = ORFs_num, ncol = 2)

##  Function to multiple inputs with readline
readlines <- function(...) {
  lapply(list(...), readline)
}

cat("One by one: Enter start and stop position for each ORF\n ")
for (i in 1: 5){
  ORF_coord[i,1] <- as.numeric(readline("Start Position: "))
  ORF_coord[i,2] <- as.numeric(readline("Stop Position: "))
}


ORF_start = c(115, 1915, 2886, 3669, 6067)
ORF_end = c(1848, 2799, 3497, 6032, 12651)
ORF_ranges <- IRanges(start = ORF_coord[,1], end = ORF_coord[,2])


####  Overlapping Bins

