####Sorting and Cleaning Data from master excel sheet 
#### Started by Jay Phadke 9-06-2025

###Useful resources: 
## Useful AI : https://www.phind.com/
## Use cytoscape post this

## Before beginning, save the required sheet in the drive folder as Torpor_genes.csv if absent/not updated

##Install all packages before reading in the libraries
##Reading in the Libraries####
library(dplyr)
library(stringr)


#setting the working directory
here()

##Reading in the Data####
tgdat <- read.csv("Torpor_genes.csv")
##Making extar columns####
##make a combined column for genus and species
tgdat <- tgdat|> mutate(Genus_Species=paste(Genus,Species,sep="_"))
##convert columns with discrete options into factors
tgdat <- tgdat|> mutate(across(c("Genus_Species","Dormancy_type","DataType","Datacode"),as.factor))

##Filtering out columns####
##making a list of all useful columns
sortcol <- c("GeneID","GENEID_CAPS","Genus_Species","Tissue","Fold_change","Log2_Fold_Change","P_value","Dormancy_type","DataType","Datacode","Significant")
##Keep only required columns
tgdat_colsort<- tgdat |> select(sortcol)

##Filtering out Data####
##filtering all results with significant p-value (Significant is 1 if p<0.05 and 0 otherwise)
tgdat_valsort <- tgdat_colsort |> filter (Significant==1) |>
##filtering based on fold change 
  filter(abs(Fold_change)>2) |>
##filtering only torpor vs normothermy. Refer to database key for more info  
  filter(str_detect(Datacode,"A")) |>
##filtering only transcriptomic data
  filter(DataType=="T") |>
##filtering only hibernation data  
  filter(Dormancy_type=="H")
##Removing more than one gene per Tissue (can add Genus_species instead)
tgdat_fil <- tgdat_valsort |> group_by(Tissue,GENEID_CAPS) |>
  filter(n()==1) |>
  ungroup() 

##Outputting clean and filtered Data####
##outputting as a csv file
write.csv(tgdat_fil, "filtered_genelist.csv", row.names = FALSE)

