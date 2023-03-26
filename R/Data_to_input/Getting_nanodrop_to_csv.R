
# install libraries
library(dplyr)
library(XML)
library(methods)
library(plyr)
library(readr)

#read tsv table

results <- read_tsv("Rodents_26032023.tsv")

# remove ffirst column
results <-  results[,-c(1, 3)]

#change the column names
results <- results %>%
    dplyr::rename("Sample" = "V2", "Date" = "V4", "Time" = "V5", "Nucleic_Acid" = "V6", 
                  "Unit" = "V7", "A260" = "V8", "A280" = "V9", "260/280" = "V10", 
                  "260/230" = "V11", "Sample_type" =  "V12", "Factor" = "V13" )
write.csv(results, 
          "~/GitHub/Namibia_project/Data/Nanodrop_measurements/CSV/Rodents_26032023.csv",
          row.names = FALSE)

