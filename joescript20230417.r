# install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("flowViz")
BiocManager::install("flowVS")
BiocManager::install("flowAI")
BiocManager::install("flowAI")
BiocManager::install("PeacoQC")
BiocManager::install("CATALYST")
BiocManager::install("SingleCellExperiment")

if(!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")

devtools::install_github('saeyslab/CytoNorm')
install.packages("uwot")
install.packages("knitr")
install.packages("xlsx")

# Load necessary libraries
set.seed(123)
library(flowCore)
library(flowViz)
library(flowVS)
library(flowAI)
library(PeacoQC)
library(CATALYST)
library(CytoNorm)
library(SingleCellExperiment)
library(uwot)
library(knitr)
library(xlsx)

# Load data from *
