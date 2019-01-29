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
# that the Cts reliably match with sample identity.
#
# TODO: Detector and ct info may be enough, but find out how to make sure
# TODO: How this code deals with NA samples?

library(tidyverse)



setwd("/Users/Miiko/Documents/laboratory/qPCR_yeast-ORIs/2nd_primerSet/validation/2018-06-15/analysis")
input.file <- "clean-data.tsv"
data <- read.table(input.file, header=T, colClasses=c("factor", "factor", "factor", "numeric"))
# View(data)
primerPairs <- as.factor(levels(data$Detector))
# View(primerPairs)

# create a dilSer, which is a vector of template concentrations in ng
# create a log of dilSer
startingAmount <- 0.5
dilutionFactor <- 10
numPoints <- 5
dilSer <- c()
for(i in 1:numPoints){
    nextDil <- startingAmount/(dilutionFactor^(i-1))
    dilSer <- c(dilSer, nextDil)
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
# produces also Ct ~ dilSer plots with the linear model.
primerStats <- data.frame()
for (primer in primerPairs){
  primerSubset <- subset(data, Detector %in% primer, select = c(Detector, Sample, Ct))
  #print(primerSubset)
  dataSubset <- data.frame(primerSubset, dilSer)
  #print(dataSubset)
  # Create vector of concentrations for standard curve 
  Cts <- dataSubset$Ct
  dilSerLog <- logNg[!isNA]
  linMod <- lm(Cts ~ dilSerLog)
  slope <- as.numeric(linMod$coeff[2])
  Rsqr <- (cor(dilSerLog, Cts))^2
  efficiency <- reactionEfficiency(slope)
  
  # bind the results to primerStats
  tempDF <- data.frame(primerPair = primer, slope = sprintf("%.3f", round(slope, 2)), efficiency = sprintf("%.3f", round(efficiency, 3)), Rsqr = sprintf("%.3f", round(Rsqr, 3)))
  primerStats <- rbind(primerStats, tempDF)
  # write results to a file
  write.table(primerStats, "primerValidationResults.csv", sep =',', row.names = FALSE, quote = FALSE)

  # plot dilution series and linear model
  dilSer
  plot <-
    ggplot(dataSubset, mapping = aes(dilSer, Ct)) +
    geom_point() +
    scale_x_log10() +
    theme_classic(14) +
    stat_smooth(method = "lm", col = "red", size = 0.3)
  ggsave(input.file = paste(primer,".tif", sep = ""),
         device = "tiff",
         dpi = 300)
  }
View(primerStats)

