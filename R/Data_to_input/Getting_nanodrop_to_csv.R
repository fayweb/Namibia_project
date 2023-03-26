
# install libraries
library(dplyr)
library(XML)
library(methods)
library(plyr)
library(readr)

#read tsv table
results <- read_tsv("CSV/Düppel_26022023.csv")

# remove ffirst column
results <-  results[,-1]

#change the column names

write.csv(results, 
          "~/GitHub/Namibia_project/Data/Nanodrop_measurements/CSV/Düppel_26022023.csv",
          row.names = FALSE)


