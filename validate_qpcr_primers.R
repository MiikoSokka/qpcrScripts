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


setwd("/Users/Miiko/Documents/projects/ns-seq_method-dev/experiments_yeast-ns-seq/spikeIn/validation_qpcr/2019-05-30_0523-0528_combined")
input.file <- "clean-data.tsv"
data <- read.table(input.file, header=T, colClasses=c("factor", "numeric", "numeric"))
#data <- filter(data, Sample != "dil7")
#data <- filter(data, Sample != "dil8")
data
#data <- data %>% filter(grepl("ssDNA", Sample)) #filter only subset of data, either ssDNA series or dsDNA (change as appropriate)

primerPairs <- as.factor(levels(data$Detector))
primerPairs


# create a dilSer, which is a vector of template concentrations in ng
# create a log of dilSer
startingAmount <- max(data$Sample)
startingAmount
dilutionFactor <- 10
numTechRep <- 2 # allows two values; either 1 or 2
numPoints <- length(unique(data$Sample))
dilSer <- c()
for(i in 1:numPoints){
  nextDil <- startingAmount/(dilutionFactor^(i-1))
  if (numTechRep == 1) {
  dilSer <- c(dilSer, nextDil)
  } else {
      for (j in 1:numTechRep){
        dilSer <- c(dilSer, nextDil)
      }
    }
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
  tick.labels <- c(startingAmount/1e7, startingAmount/1e6, startingAmount/1e5, startingAmount/1e4, startingAmount/1e3, startingAmount/100, startingAmount/10, startingAmount)
  plot <-
    ggplot(dataSubset, mapping = aes(dilSer, Ct)) +
    labs(title = primer, x = "Copy number", y = "Cycle threshold, Cq") +
    geom_point() +
    annotate("text", x = startingAmount/10, y = 30, label = paste("slope = ",round(slope,digits = 2))) +
    annotate("text", x = startingAmount/10, y = 28, label = paste("efficiency = ",round(efficiency,digits = 2))) +
    annotate("text", x = startingAmount/10, y = 26, label = paste("Rsqr = ",round(Rsqr,digits = 3))) +
    ylim(0, 40) +
    scale_x_log10(labels = tick.labels, breaks = tick.labels) +
    annotation_logticks(sides = "b") +
    theme_classic(14) +
    stat_smooth(method = "lm", col = "red", size = 0.3)
  ggsave(filename = paste(primer,".tif", sep = ""),
         path = "output_validation",
         device = "tiff",
         dpi = 300)
  }
View(primerStats)

