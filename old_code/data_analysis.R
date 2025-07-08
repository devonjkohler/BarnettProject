
library(tidyverse)
library(MSstats)
library(data.table)

source("data_processing.R")

# Load and prep data -----------------------------------------------------------
diann_file = fread("Data/all_search.tsv")
annotation = fread("Data/meta.csv")

ub_ligase_list = fread("Data/Ubiquitin_ligases.csv")
ub_ligase_list$gene_match = paste0(ub_ligase_list$`Gene Symbol`, "_MOUSE")

# Prep for MSstats converter
diann_file$File.Name = diann_file$Run
annotation$BioReplicate = substring(str_split_i(annotation$Run, "_", 4), 2, 3)
# Include all covariates for differential analysis
annotation$Condition = paste(annotation$Tissue, 
                             annotation$Age, 
                             annotation$Gender, 
                             sep="_")
annotation = annotation[, c("Run", "Time", "BioReplicate", "Condition")]

# Half life data ---------------------------------------------------------------
half_life_data = convert_LH_data(diann_file, annotation)
# TODO: why weird inf values
half_life_data = half_life_data %>% filter(is.finite(Intensity))

# TODO: dataProcess might not (probably) wont work well here
summarized_hl_data = dataProcess(half_life_data, 
                                 normalization=FALSE,
                                 MBimpute = FALSE,
                                 remove50missing=FALSE,
                                 maxQuantileforCensored=1)

summarized_hl_data = half_life_data %>% 
  group_by(ProteinName, Run) %>% 
  summarize(LogIntensities = median(Intensity, na.rm=TRUE))
  
summarized_hl_data %>% 
  write.csv("Data/half_life_data.csv", row.names = FALSE)


# Protein intensity data -------------------------------------------------------
# Converter
msstats_input = MSstatsConvert::DIANNtoMSstatsFormat(diann_file, annotation)
save(msstats_input, file="Data/msstats_input.rds")

# Data processing
summarized_data = dataProcess(msstats_input, 
                              normalization=FALSE,
                              featureSubset="topN", 
                              n_top_feature=20,
                              numberOfCores=8)
save(summarized_data, file="Data/summarized_data.rds")
load(file="Data/summarized_data.rds")
summarized_data$ProteinLevelData %>% 
  write.csv("Data/quant_data.csv", row.names = FALSE)

dataProcessPlots(summarized_data, type="ProfilePlot",
                 featureName = "NA",
                 which.Protein = "AFF4_MOUSE", address=FALSE)

summarized_data$ProteinLevelData$Tissue = str_split_i(
  summarized_data$ProteinLevelData$GROUP, "_", 1)

# Ub Ligase per tissue
summarized_data$ProteinLevelData %>% 
  filter(Protein %in% ub_ligase_list$gene_match) %>% 
  group_by(Tissue) %>% summarize(n_distinct(Protein))

# All measured Ub Ligase
summarized_data$ProteinLevelData %>% 
  filter(Protein %in% ub_ligase_list$gene_match) %>% distinct(Protein) %>% 
  write_csv("Data/measured_ub_ligase.csv")



# Differential analysis (ignore) -----------------------------------------------
head(summarized_data$FeatureLevelData)
levels(summarized_data$ProteinLevelData$GROUP)
comparison <- matrix(c(-1/4,-1/4,-1/4,-1/4,
                       1/4,1/4,1/4,1/4,
                       0,0,0,0,
                       0,0,0,0),nrow=1)
row.names(comparison) <- "BM-Brain"
groups = levels(summarized_data$ProteinLevelData$GROUP)
colnames(comparison) <- groups[order(as.numeric(groups))]

model = groupComparison(comparison, summarized_data)
save(model, file="Data/model.rds")

