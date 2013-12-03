## Reading in Excel file for Mutation Rate analysis in SCN Virus
## This is for the SCN Rhabdovirus (negative-sense ssRNA) which only has one ORF (one gene)

library(ggplot2)
library(stringr)
library(plyr)
library(seqinr)  ## contains codon table (DNA)
library(Biostrings)  ## sequence related algorithms
library(fastseg)     ## fast seqmentation
library(GenomicRanges)  ##GRanges class and method



Mut_rate_Rhabdo <- function(){
  
  require(ggplot2)
  require(stringr)
  require(plyr)
  require(seqinr)  ## contains codon table (DNA)
  require(Biostrings)
  
  ##  Include functions
  source("C:/Users/Mac/Dropbox/SCN-Thesis/Virus_Sequence/data_import.R")
  
  
  
#bin = 8  #divisible by 8
##  Grab data from two files and return a list of dataframe and character sequence

genome.dat <- data_import()
  
Rhabdo_dat = genome.dat$Mutation.data  ##
Seq = genome.dat$Consensus  ## Consensus sequence DNAString object
  genome.size <-length(Seq)
  
  
# not needed anymore: 
  #Consensus.vector <- laply(seq(1, nchar(Consensus.seq), 1), function(i) substr(Consensus.seq, i,i))  ## splits into individual bases: string array
  

## convert sequence to a Bioconductor Object
seq <- DNAString(Consensus.seq)
## Amino Acid Sequence
Consensus.codon <- laply(seq(1, nchar(Consensus.seq), 3), function(i) substr(Consensus.seq, i,i+2))  ## splits into 3 base segments 
AAstring_all <- translate(seq)

##---Reverse funciton----##
#strReverse <- function(x){
#  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")}
##-----------------------##

  
  
##  Bins
##  Set bin steps in bin.width
bin <- as.numeric(readline(" Enter number of bins per ORF or enter '0' for analysis at every codon: "))
if(bin!=0){
  if(genome.size %% as.integer(bin) ==0){
    binwidth = genome.size/bin
  } else{ 
    binwidth = floor(genome.size/bin) ## will potentially cut off the very last base
  }
} else if(bin == 0){   ## bins are coerced to integers to compare, need to be factor of genome.size otherwise need to round 
  ##  create a function to separate by codon and generate peptide sequence
  bin = genome.size/3
  break
} 

## Base Changes
#Base_from = Rhabdo_dat$Reference
#Base_to = Rhabdo_dat$Allele


## Amino Acid Changes
n = length(Rhabdo_dat$Amino.acid.change)
Rhabdo_dat$Amino.acid.change = as.character(Rhabdo_dat$Amino.acid.change)

#  Set memory for vectors
Amino.acid.change = rep(0,n)
AA_from = rep(0,n)
AA_to =rep(0,n)

##  Split the AA's, mutate from and mutate to:
for(i in 1:n){
    Amino.acid.change[i] <- strsplit(Rhabdo_dat$Amino.acid.change[i], "[.]")[[1]][2]
  
    if (is.na(Amino.acid.change[i])){
      AA_from[i] <-NA
      AA_to[i] <-NA
    }
    else{
      AA_from[i] <-strsplit(as.character(Amino.acid.change[i]), "[[:digit:]]{1,}")[[1]][1]
      AA_to[i] <-strsplit(as.character(Amino.acid.change[i]), "[[:digit:]]{1,}")[[1]][2]
    }
}

Rhabdo_sub <- subset(Rhabdo_dat, select = c(Reference.Position, Type,Reference, Allele, Frequency, Coding.region.change))

## Add the amino acid changes
Rhabdo_sub$AA_from = AA_from
Rhabdo_sub$AA_to = AA_to


## Calculate Ka/Ks per bin
#  first locate nucleotide positions with substitution mutations
Rhabdo_SNV <- subset(Rhabdo_sub, Rhabdo_sub$Type == "SNV")
Rhabdo_SNV$Type <- NULL

## Ka = nonsynonomous substitutions, Ks = synonomous
Ka_dat <- subset( Rhabdo_SNV, select= c(Reference.Position,Frequency, AA_to), is.na(Rhabdo_SNV$AA_to) == FALSE) # Note: this includes "*" which i'm not sure the meaning of
Ks_dat <- subset( Rhabdo_SNV, select= c(Reference.Position, Frequency, AA_to), is.na(Rhabdo_SNV$AA_to) == TRUE)

##  Take the cumulative frequencies at each bin
bin.break = rep(1,bin+1)
  
## Piecing together reference positions and bin.breaks for Ka/Ks calculation
while(bin.break[bin+1]==1){
#for(i in 1:bin){
  # First bin
 # for(k in 1:nrow(Ka_dat)){
  #  if(Ka_dat$Reference.Position[k] <= (genome.size/bin)){  
   #    bin.break[i+1] =which.min(Ka_dat$Reference.Position <= (genome.size/bin))  #not sure why which.min works and not which.max
       ##  which.mins cut off at one step higher than the end of the which vector
       ##  Thus, the index we get will be the starting point of the next bin as opposed to the ending point of the current bin.
    #}
  #}
#}
  ##  remaining bins
 for(i in 1:bin){
   for(k in bin.break[i+1]:nrow(Ka_dat)){
    if(Ka_dat$Reference.Position[k] > (genome.size/bin*i)){                                               
       #bin.break[i+1] = which.min(xor(Ka_dat$Reference.Position >= (genome.size/bin*i), Ka_dat$Reference.Position < genome.size)
       bin.break[i+1] = which.max(Ka_dat$Reference.Position >= (genome.size/bin*i))
     }
    } 
    ##got. here bin.break is up to 1 1, 8 and 9
   }
 
   ## got here

## Check for empty bins--will show up as 1's toward the end
 for(i in 2:(bin+1)){  ## Last bins if "empty"
     ##  got here
   if(bin.break[i]==1){
     split = (nrow(Ka_dat)-bin.break[i-1])/(bin+2-i)
     got.here = 1  ## troubleshoot
     end = i
     for(n in end:bin){
       bin.break[n] = bin.break[n-1]+split
       bin.break[bin+1] = nrow(Ka_dat)
       got.here = 1 # troubleshoot
     }
     break
   }
  }     
}

## Ka/Ks per bin
## Set up empty vectors
Ka.Ks = rep(0,bin)
Ka_bin = rep(0,bin)
Ks_bin = rep(0, bin)
##  Calculation Ka/Ks per bin
  for(i in 1:bin){
    freqs_a = subset(Ka_dat, Ka_dat$Reference.Position >= bin.break[i] & Ka_dat$Reference.Position < bin.break[i+1])
    freqs_s = subset(Ks_dat, Ks_dat$Reference.Position >= bin.break[i] & Ks_dat$Reference.Position < bin.break[i+1])
    Ka_bin[i]= sum(freqs_a$Frequency)
    Ks_bin[i]= sum(freqs_s$Frequency)
    Ka.Ks[i] = sum(freqs_a$Frequency)/sum(freqs_s$Frequency)

  }  ## 0/0 become NaNs

##  Rewire some of the uncalculable values
for(i in 1:bin){ 
  if(Ka.Ks[i]== "NaN"){## NaNs will be turned to 0
    Ka.Ks[i] = 0
  }else if (Ka.Ks[i] ==Inf){##  Infs will turn to NaN
    Ka.Ks[i]=NaN
  }
}
## Coordinate x positions
genome.position= seq(genome.size/bin,genome.size, by=genome.size/bin )
## Plot bar graph of Ka/Ks
ggplot()+geom_bar(aes(x= genome.position, y = Ka.Ks), stat = "identity") + xlab("Nucleotide Position")+ylab("Ka/Ks")

  
## Connection to RNAfold program (only on Unix machines)
#system("./RNAfold -p -d2 -n -noLP<sequence.fna")
## ouptput file as dot.ps  
}  ## end of function