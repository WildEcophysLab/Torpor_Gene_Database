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
library(here)


#setting the working directory
setwd(here())

##Reading in the Data####
tgdat <- read.csv("Torpor_genes.csv")
##Making extra columns####
##make a combined column for genus and species
tgdat <- tgdat|> mutate(Genus_Species=paste(Genus,Species,sep="_"))
##convert columns with discrete options into factors
tgdat <- tgdat|> mutate(across(c("Genus_Species","Dormancy_type","DataType","Datacode"),as.factor))

##Filtering out columns####
##making a list of all useful columns
sortcol <- c("GeneID","GENEID_CAPS","Genus_Species","Tissue","Fold_change","Log2_Fold_Change","P_value","Dormancy_type","DataType","Datacode","Significant")
##Keep only required columns
tgdat_colsort<- tgdat |> select(all_of(sortcol))

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

##Filtering out species that have more than 2000 entries
tgdat_fil <- tgdat_fil %>%
  group_by(Genus_Species)

##lower threshold for classyfying large networks (for running commandspost vs commandsrun)
lim <- 1500

# Extract large groups (≥ 2000 rows)
tgdat_fil_large <- tgdat_fil %>%
  filter(n() >= lim) %>%
  ungroup()

# Extract small groups (< 2000 rows)
tgdat_fil <- tgdat_fil %>%
  filter(n() < lim) %>%
  ungroup()

##Outputting clean and filtered Data####
##outputting as a csv file
write.csv(tgdat_fil, "filtered_genelist.csv", row.names = FALSE)
write.csv(tgdat_fil_large, "filtered_genelist_large.csv", row.names = FALSE)

##Creating a list of all the species for cytoscape analysis ####
species_list <- tgdat_fil %>%
  mutate(species = as.character(Genus_Species)) %>%
  distinct(species, .keep_all = FALSE) %>%
  rename(Species_list = species) %>%
  mutate(ref_species_id = NA)

species_list_large <- tgdat_fil_large %>%
  mutate(species = as.character(Genus_Species)) %>%
  distinct(species, .keep_all = FALSE) %>%
  rename(Species_list = species) %>%
  mutate(ref_species_id = NA)

##Writing out species list as csv####
write.csv(species_list, "species_list.csv", row.names = FALSE)
write.csv(species_list_large, "species_list_large.csv", row.names = FALSE)

print("END")



