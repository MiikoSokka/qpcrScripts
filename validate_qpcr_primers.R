#!/usr/bin/env Rscript

# qPCR primer validation by Miiko (some lines of code by John Urban)
# 
# These lines of code will output a table that shows the slope,
# reaction efficiency and R-square for each primer pair.
# 
# Input: clean-data.tsv file that has the following columns:
# Detector, Sample, Ct.
#
# Output: a primerValidationResults.csv table with reaction efficiency,
# slope and R-square values for each tested primer.
#
# Acceptable range = 0.9 - 1.1 for reaction efficiency
#
# that the Cts reliably match with sample identity.
#
# TODO: Detector and ct info may be enough, but find out how to make sure
# TODO: How this code deals with NA samples?

library(tidyverse)
library(dplyr)



setwd("/Users/Miiko/Documents/projects/ns-seq_method-dev/experiments_yeast-ns-seq/spikeIn/validation_qpcr/2019-05-28_qpcrValidation/2019-05-28_qpcrValidation_0528_qpcr")
input.file <- "clean-data.tsv"
data <- read.table(input.file, header=T, colClasses=c("factor", "factor", "numeric"))
data <- filter(data, Sample != "dil7")
data <- filter(data, Sample != "dil8")
data
#data <- data %>% filter(grepl("ssDNA", Sample)) #filter only subset of data, either ssDNA series or dsDNA (change as appropriate)

primerPairs <- as.factor(levels(data$Detector))
primerPairs


# create a dilSer, which is a vector of template concentrations in ng
# create a log of dilSer
startingAmount <- 1
dilutionFactor <- 10
numTechRep <- 1
numPoints <- length(unique(data$Sample))
dilSer <- c()
for(i in 1:numPoints){
  nextDil <- startingAmount/(dilutionFactor^(i-1))
  dilSer <- c(dilSer, nextDil)
  #    for (j in 1:numTechRep){
  #      dilSer <- c(dilSer, nextDil)
  #    }
}
logNg <- log10(dilSer)
logNg
dilSer

# create a function to calculate reaction efficiency
reactionEfficiency <- function(slope){
  ## reaction efficiency is 10^(-1/slope) -1
  eff <- (10^(-1/slope) - 1)
  return(as.numeric(eff))
}

# iterate over all primer pairs in the dataset to produce a table of reaction efficiency, slope, and R-square.
# produces also Ct ~ dilSer plots with a linear model.
primerStats <- data.frame()
for (primer in primerPairs){
  primerSubset <- subset(data, Detector %in% primer, select = c(Detector, Sample, Ct))
  print(primerSubset)
  dataSubset <- data.frame(primerSubset, dilSer)
  print(dataSubset)
  # Create vector of concentrations for standard curve 
  Cts <- dataSubset$Ct
  dilSerLog <- logNg[!is.na(logNg)]#previously was: logNg[!isNA], but this was wrong
  linMod <- lm(Cts ~ dilSerLog)
  slope <- as.numeric(linMod$coeff[2])
  Rsqr <- (cor(dilSerLog, Cts))^2
  efficiency <- reactionEfficiency(slope)
  
  # bind the results to primerStats
  tempDF <- data.frame(primerPair = primer, slope = sprintf("%.3f", round(slope, 2)), efficiency = sprintf("%.3f", round(efficiency, 3)), Rsqr = sprintf("%.3f", round(Rsqr, 3)))
  primerStats <- rbind(primerStats, tempDF)
  # write results to a file
  dir.create("output_validation", showWarnings = FALSE)
  write.table(primerStats, "output_validation/primerValidationResults.csv", sep =',', row.names = FALSE, quote = FALSE)

  # plot dilution series and linear model
  dilSer
  plot <-
    ggplot(dataSubset, mapping = aes(dilSer, Ct)) +
    geom_point() +
    scale_x_log10() +
    theme_classic(14) +
    stat_smooth(method = "lm", col = "red", size = 0.3)
  ggsave(filename = paste(primer,".tif", sep = ""),
         path = "output_validation",
         device = "tiff",
         dpi = 300)
  }
View(primerStats)

