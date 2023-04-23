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

#
# for now, working on this with just concatenated live small CD3+ T cells
#

# Load data from *FCS files
pathname <- "/Volumes/jwreda/013 flow files/LN Tfh/concat_1_small T cells - 013 LN Tfh - 20230420 - unmixed1.fcs"

# Read in the data using flowCore
fcs_data <- read.FCS(pathname,
                    transformation = FALSE,
                    truncate_max_range = FALSE)

# rewrite the above but as reading a flowSet
fcs_data <- read.flowSet(pathname,
                         transformation = FALSE,
                         truncate_max_range = FALSE)

fcs_colname <- colnames(fcs_data)

# fcs_fluors <- c("SampleID",
#               "FSC-A","FSC-H",
#               "SSC-A", "SSC-B-A",
#               "SSC-B-H", "SSC-H", "CD3",
#               "PD-1", "CD62L", "FoxP3",
#               "CD4", "CD8",
#               "CD25", "SiglecF", "CD19", "LD",
#               "CD40L", "Bcl6", "Gata3",
#               "CD44", "Time")

antigen <- pData(parameters(fcs_data[[1]]))$desc

# head(fcs_colname) #nolint

marker_class <- rep("none",
                    ncol(fcs_data[[1]])
                    )



# head(marker_class) #nolint

# identify index of markers for typing
marker_state <- c(9, 18)
marker_type <- c(8, 10:17, 19:21)

# assign marker class
marker_class[marker_state] <- "state"
marker_class[marker_type] <- "type"

# define panel table
panel <- data.frame(fcs_colname,
                    antigen,
                    marker_class,
                    row.names = NULL
                    )

# save panel as excel file
write.xlsx(panel,
          file = "panel_A.xlsx",
          sheetName = "Panel_A"
          )

# use kable to view panel in console
kable(panel)

### Transforming spectral flow cytometry data
markerstotransform <- panel$fcs_colname[c(8:21)]


Downsampling_FlowSet <- function(x, samplesize, replace = TRUE, prob = NULL) {
  if(missing(samplesize))
    samplesize <- min(flowCore::fsApply(x,nrow))
  flowCore::fsApply(x, function(ff){
    i <- fcs_data(nrow(ff), size = samplesize, replace = replace, prob)
    ff[i,]
  })
}

fcs_data_small <- Downsampling_FlowSet(x = fcs_data, samplesize = 2000) 
            # samplesize is the number of cells included, you can include more cells.

cofactors <- estParamFlowVS(fcs_data_small, channels = markerstotransform)
cofactordata <- data.frame(markerstotransform, cofactors)