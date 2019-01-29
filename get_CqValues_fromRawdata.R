#!/usr/bin/env Rscript
#
# Author: Miiko
# This script processes raw qPCR data and obtains "cycle threschold" or Ct/Cq values
# for each sample, using qpcR package that uses second derivative cpD2 to estimate cycle threshold
# This method is ought to be more robust in estimating the Ct/Cq value.
#
# Input: Comma separated Rn-file from qpcr program (contains raw fluorescence values for each cycle).
# The file has following columns:
# Well, Detector, Reporter Dye and Cycle 1:Cycle 40.
# 
# This script transforms the Rn.csv data and results in following columns:
# cycleNr and well numbers from 1 to whatever has been in use. Well numbers
# correspond to rows in 96-well plate in such a way that row A is 1-12,
# row B is 13-24, C is 25-36, D is 37-48, E is 49-60, F is 61-72, G is 73-84
# and H is 85- 96.
#
# Pictures of Fluorescence values vs. cycle number are saved as tiff files for
# each of the detector group.
#
# The raw data is then processed with qpcR package to yield cpD2 values for each
# detector. If flag duplicate samples is on (default on), duplicate Cq values are
# presented in another table on their own columns, with additional columns for
# mean and range.
#
# Outputs (in two directories, output_Cq and output_flGraphs):
# Cq.csv:
#   Contains Cq cycle threshold values obtained using qpcR
#   package by Andrej-Nikolai Spiess. Cq values are defined using LM-Algorithm.
# table_Cq.csv:
#   If replicate flag is set, the replicate Cq values are stored in their respective
#   columns (rep1, rep2), with their means and ranges.
# <detector>.tif:
#   If "skip plots" flag is not set, tiff picture for each detector is generated.
#   Contains fluorescence vs. cycle Nr data for visual inspection of how qPCR progressed.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(qpcR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))

# FUNCTIONS

# Helper function to load data in R and transform it into the correct form
transformData <- function(file.data){
  df.data <- read.csv(file.data, header = TRUE, sep = ",")
  wellID <- df.data$Well # extract well id
  wellID.t <- t(wellID) # transform well id and store it in a variable
  df.data <- df.data[c(4:43)] # select only cycle data columns
  df.data <- t(df.data) # transform cycle data
  colnames(df.data) <- wellID.t # add wellID information as column names
  cycleNr <- seq.int(nrow(df.data))
  row.names(df.data) <- cycleNr # change row names
  df.data <- cbind(cycleNr,df.data)
  df.data <- data.frame(df.data, stringsAsFactors = TRUE) #this was needed to make pcrfit functional
  return(df.data)
}

# Helper function to subset fluorescence data for plotting.
subsetFluoData <- function( file.data, df.names ){
  # Set vectors based on well ID, first row 1:12 etc. until 96.
  # Each row must contain data from one detector.
  # Get the data and subset into dataframes based on rows
  df.fldata <- read.csv(file.data, header = TRUE, sep = ",")
  df.fldata <- df.fldata[,c(1, 4:43)]
  df.fldata <- data.frame(df.names, df.fldata)
  df.fldata <- df.fldata[,c(4,2,3,5:44)]
  colnames(df.fldata) <- c("Well", "detector", "sample", 1:40) # change the column names to numbers
  
  df.fldata.A <- subset(df.fldata, df.fldata$Well %in% 1:12)
  df.fldata.B <- subset(df.fldata, df.fldata$Well %in% 13:24)
  df.fldata.C <- subset(df.fldata, df.fldata$Well %in% 25:36)
  df.fldata.D <- subset(df.fldata, df.fldata$Well %in% 37:48)
  df.fldata.E <- subset(df.fldata, df.fldata$Well %in% 49:60)
  df.fldata.F <- subset(df.fldata, df.fldata$Well %in% 61:72)
  df.fldata.G <- subset(df.fldata, df.fldata$Well %in% 73:84)
  df.fldata.H <- subset(df.fldata, df.fldata$Well %in% 85:96)
  
  list.flSubsets <- list(df.fldata.A, df.fldata.B,
                         df.fldata.C, df.fldata.D,
                         df.fldata.E, df.fldata.F,
                         df.fldata.G, df.fldata.H)
  return(list.flSubsets)
}

# Main function to get the Cq values.
# Depends on the function transformData
getCqValues <- function(file.data){
  df.data <- transformData(file.data) #get the data in through another function
  df.cpD2 <- data.frame(sample = character(), cpD2 = numeric()) #set the empty table
  vectorSample <- c()
  vectorDP2 <- c()
  ncol.df <- NCOL(df.data)
  while (ncol.df > 1){
    df.sample <- df.data[,c(1:2)] # get data from one sample
    m1 <- pcrfit(df.sample, cyc = 1, fluo = 2, model = l4)   #build model
    m2 <- efficiency(m1, plot = FALSE) #pass "m1" into efficiency and bind to "m2"
    cpD2 <- m2$cpD2 #extract m2$cpD2 value and store with the rows into a new table
    col.name <- colnames(df.sample)
    col.name <- col.name[2]
    vectorSample <- append(vectorSample, col.name)
    vectorDP2 <- append(vectorDP2, cpD2)
    df.data <- df.data[-c(2)]
    ncol.df <- NCOL(df.data)
  }
  # Bind the results of the while loop into df.cpD2
  df.cpD2 <- data.frame(vectorSample, vectorDP2)
  colnames(df.cpD2) <- c("well", "cpD2")
  df.cpD2 <- cbind(df.cpD2, df.names)
  
  # The following lines are required to remove an extra 'X' in the well column.
  # The 'X's are added by the last line in transformData, I believe. This is the
  # result of checking names to be syntactically correct (see script above).
  # Check also https://www.rforexcelusers.com/remove-currency-dollar-sign-r/
  well <- as.numeric(gsub("X", "", df.cpD2$well))
  df.cpD2 <- df.cpD2[,c("wellNo", "detector", "sample", "cpD2")]
  df.cpD2 <- cbind(df.cpD2, well)
  
  if ( df.cpD2$well == df.cpD2$wellNo ){ # check that the well number information matches
    df.cpD2 <- df.cpD2[,c("well", "detector", "sample", "cpD2")]
    write.table(df.cpD2, file.path("output_Cq", "Cq.csv"),
                sep =',', row.names = FALSE) # write results in a comma separated file
  } else {
    stop("Well IDs do not match between sample file and raw fluorescence data file.")
  }
  
  # If replicate flag is on, add an if loop to transform replicated data
  # to a new file, where the replicate samples are in their own columns,
  # and add new columns for mean and range.
  if ( opt$replicate == TRUE){
    # Define empty vectors for loop
    vector.rep1 <- as.numeric(c())
    vector.rep2 <- as.numeric(c())
    vector.detector <- as.character(c())
    vector.sample <- as.character(c())
    for ( row in 1:(nrow(df.cpD2)-1) ){
      nextrow <- row + 1
      if ( df.cpD2$detector[row] == df.cpD2$detector[nextrow] &
           df.cpD2$sample[row] == df.cpD2$sample[nextrow] ){
        vector.rep1 <- append(vector.rep1, df.cpD2$cpD2[row])
        vector.rep2 <- append(vector.rep2, df.cpD2$cpD2[nextrow])
        vector.detector <- append(vector.detector, as.character(df.cpD2$detector[row]))
        vector.sample <- append(vector.sample, as.character(df.cpD2$sample[row]))
      }
    }
    df.Cq <- data.frame("detector" = vector.detector, "sample" = vector.sample,
                        "rep1" = round(vector.rep1, 2),
                        "rep2" = round(vector.rep2, 2))
    # Write the results in table
    df.Cq <- mutate(df.Cq, "mean" = round((rep1+rep2)/2, 2),
                    "range" = round(abs(rep1-rep2), 2))
    write.table(df.Cq, file.path("output_Cq", "table_Cq.csv"),
                sep =',', row.names = FALSE)
  }
}

# Main function to get automatically/manually generated Ct values from a file
# obtained from qPCR machine. Called if the Ct file
# is provided as an argument.
getCtValues <- function(file.manualCt){
  df.Ctdata <- read.csv(file.manualCt, header = TRUE, sep = ",",
                        colClasses =c("integer", "NULL", "NULL",
                                      "NULL", "numeric"), skip = 2)
  colnames(df.Ctdata) <- c("well", "Ct")
  df.Ctdata <- cbind(df.Ctdata, df.names) # Add names from supplied file
  
  # Check that the well number information matches and write results
  # in a comma separated file.
  if ( df.Ctdata$well == df.Ctdata$wellNo ){
    df.Ctdata <- df.Ctdata[,c("well", "detector", "sample", "Ct")]
    write.table(df.Ctdata, file.path("output_Cq", "Ct.csv"),
                sep =',', row.names = FALSE)
  } else {
    stop("Well IDs do not match between sample file and raw fluorescence data file.")
  }
  
  # If replicate flag is on, add an if loop to transform replicated data
  # to a new file, where the replicate samples are in their own columns,
  # and add new columns for mean and range
  if ( opt$replicate == TRUE ){
    # Define empty vectors for loop
    vector.rep1 <- as.numeric(c())
    vector.rep2 <- as.numeric(c())
    vector.detector <- as.character(c())
    vector.sample <- as.character(c())
    for ( row in 1:(nrow(df.Ctdata)-1) ){
      nextrow <- row + 1
      if ( df.Ctdata$detector[row] == df.Ctdata$detector[nextrow] &
           df.Ctdata$sample[row] == df.Ctdata$sample[nextrow] ){
        vector.rep1 <- append(vector.rep1, df.Ctdata$Ct[row])
        vector.rep2 <- append(vector.rep2, df.Ctdata$Ct[nextrow])
        vector.detector <- append(vector.detector, as.character(df.Ctdata$detector[row]))
        vector.sample <- append(vector.sample, as.character(df.Ctdata$sample[row]))
      }
    }
    df.Ct <- data.frame("detector" = vector.detector, "sample" = vector.sample,
                        "rep1" = round(vector.rep1, 2),
                        "rep2" = round(vector.rep2, 2))
    # Write the results in table
    df.Ct <- mutate(df.Ct, "mean" = round((rep1+rep2)/2, 2),
                    "range" = round(abs(rep1-rep2), 2))
    write.table(df.Ct, file.path("output_Cq", "table_Ct.csv"),
                sep =',', row.names = FALSE)
  }
}

# Main function to get fluorescence plots.
# SubsetFluoData needs to be executed first.
getPlots <- function(flSubset){
  detector <- flSubset[1,"detector"]
  df.flmelted <- melt(flSubset, id.vars = 1:3)
  colnames(df.flmelted) <- c("well", "detector", "sample", "cycleNr", "flValue")
  
  if ( nrow(flSubset) != 0 ){
    plot <- ggplot(df.flmelted, aes(cycleNr, flValue, col = sample, group = well)) +
      geom_line() +
      theme_classic(14)
    ggsave(filename = paste(detector, ".tif", sep=""),
           device = "tiff", path = "output_flGraphs",
           dpi = 100, width = 12, height = 10, units = "in")
  }
}


# COMMAND LINE ARGUMENTS
list.options <- list( 
  make_option(c("-i", "--input"), type="character", 
              help="Input file name. Directly obtained raw fluorescence values from qPCR machine as a csv file.",
              metavar="filename.csv"),
  make_option(c("-s", "--sampleNames"), type="character",
              help="File name that contains sample list in column format.
              Must contain following columns: wellNo, detector, sample.
              Must be comma separated file.",
              metavar="samplefilename.csv"),
  make_option(c("-r", "--replicate"), action = "store_true", default = FALSE,
              help="Prepares another table based on assumption that the adjacent
              two rows contain replicate samples. The table contains mean and
              range columns for the replicate samples",
              metavar="default False"),
  make_option(c("-m", "--manualCt"), type="character", default = c(),
              help="A csv file obtained from qPCR machine with manually
              obtained Ct values. Must be comma separated file.",
              metavar="ct-file.csv"),
  make_option(c("-p", "--skipPlots"), action = "store_true", default = FALSE,
              help="If set, does not create plots from the raw fluorescence data",
              metavar="default False")
)

opt <- parse_args(OptionParser(option_list=list.options))
file.data <- opt$input
file.samplelist <- opt$sampleNames
file.manualCt <- opt$manualCt

# EXECUTE
dir.working <- getwd()
setwd(dir.working)
dir.create("output_Cq", showWarnings = FALSE) # output directory for Cq/Ct tables
df.names <- read.csv(file.samplelist, header = TRUE)
if ( is.null(file.manualCt) == FALSE ){ getCtValues(file.manualCt) }
if ( opt$skipPlots != TRUE ) {
  dir.create("output_flGraphs") # output directory for fluorescence graphs
  flSubset <- subsetFluoData(file.data, df.names)
  mapply(getPlots, flSubset) # get fluorescence graphs
}
getCqValues(file.data) #get the Cq value table

# FOR TROUBLESHOOTING
#setwd("/Users/Miiko/Desktop/kumar-ct")
#file.data <- "StagedLarve-EcR-31July2018_rn.csv"
#file.samplelist <- "samplelist.csv"
# file.manualCt <- "manualfile"
  
  
  
  