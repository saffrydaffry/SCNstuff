####  FUnction for Illumina csv data import 
#Rhabdo_dat <- read.csv("C:/Users/Mac/Dropbox/SCN-Thesis/Virus_Sequence/3_RhabroMuts.csv")
#genome.path <- "C:/Users/Mac/Dropbox/SCN-Thesis/Virus_Sequence/Rhabdo_seq.txt"

data_import <- function(){

  require(Biostrings)
  
  file.path <- readline("\n Enter the file path of the spreadsheet containing mutation data,
                     \n e.g. 'C:/Users/username/folder/folder.ext' 
                      \n -->")
  
  file.ext <- strsplit(file.path, "[.]")[[1]][2]
  
  
  if(file.ext =="csv"){
    #
    Rhabdo_dat <- read.csv(file.path)
  } else if(file.ext=="xlsx"|file.ext=="xls"){  ## Ideally is a csv file, can just read
    file.path <- strsplit(file.path,"[.]")[[1]][1]
    file.path <- paste(file.path, "", sep = ".csv")
    Rhabdo_dat <- read.csv(file.path)
  }  else if(file.ext=="txt"){  ##  if excel file, convert to csv and then read
    Rhabdo_dat <- read.table(file.path)
  } else{
    file.path <-readline("\n Incompatible file format, try again
                       \n e.g. 'C:/Users/username/folder/folder.ext' \n --> ")
  }   ## if none of the above, re-eneter the file name
  
  ##  consensus sequence
  genome.path <- readline("Enter the text file path of the document containing the consensus sequence  
                         \n -->")
  Consensus.seq <- DNAString(read.table(genome.path,  colClasses = "character", stringsAsFactors = FALSE)[[1]][1])
  
  ##Return list with consensus.seq and rhabdo data
  
  result <- list(Mutation.data = Rhabdo_dat, Consensus = Consensus.seq)
  return(result)
  
}