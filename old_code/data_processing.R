
library(data.table)
library(tidyverse)

#' Function to convert diann heavy/light data into MSstats format
convert_LH_data = function(diann_report, annotation){
  
  data = as.data.table(diann_report)
  
  data$PSM = paste(data$Stripped.Sequence,
                   data$Precursor.Charge, sep="_")
  
  # Only keep peptides cleaved at K aa
  data$cleavage = str_sub(data$Stripped.Sequence, -1)
  data = data %>% filter(cleavage == 'K')
  
  # Q val filtering
  data = data %>% filter(Translated.Q.Value < .01)
  
  # Seperate into heavy and light buckets
  data$label = str_sub(data$Precursor.Id ,-3, -3)
  light_data = data %>% filter(label == 'L')
  heavy_data = data %>% filter(label == 'H')
  
  # Keep only PSMs with both labels
  matching_psms = intersect(light_data$PSM, heavy_data$PSM)
  
  light_data = light_data %>% filter(PSM %in% matching_psms)
  heavy_data = heavy_data %>% filter(PSM %in% matching_psms)
  
  # Join and calc half life
  data = merge(light_data %>% select(Protein.Names, Run, PSM, Ms1.Area),
               heavy_data %>% select(Protein.Names, Run, PSM, Ms1.Area),
               by=c("Protein.Names", "Run", "PSM"),
               suffixes = c(".light", ".heavy"))

  data$half_life = log(data$Ms1.Area.heavy/data$Ms1.Area.light + 1)
  
  # Adjust by time
  data = merge(data, 
               annotation %>% select(Run, Time, BioReplicate, Condition),
               by = c("Run"), all.x=TRUE, all.y=FALSE)
  data$half_life = data$half_life / data$Time
  
  # Random stuff MSstats needs
  data$PeptideSequence = str_split_i(data$PSM, "_", 1)
  data$PrecursorCharge = str_split_i(data$PSM, "_", 2)
  data$FragmentIon = NA
  data$ProductCharge = NA
  data$IsotopeLabelType = "L"
  data$Fraction = 1
  half_life_data$half_life = 2**half_life_data$half_life
  
  data = data %>% select(Protein.Names, PeptideSequence, PrecursorCharge, 
                         FragmentIon, ProductCharge, IsotopeLabelType, 
                         Condition, BioReplicate, Run, Fraction, half_life)
  
  setnames(data, c("Protein.Names", "half_life"), c("ProteinName","Intensity"))
  
  return(data)
  
}