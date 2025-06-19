####Automating Cytoscape Workflow
#### Started by Jay Phadke 9-06-2025

###Useful resources: 
##http://127.0.0.1:1234/v1/swaggerUI/swagger-ui/index.html?url=http://127.0.0.1:1234/v1/commands/swagger.json#/
## Useful AI : https://www.phind.com/
## STRINGapp API: https://github.com/RBVI/stringApp/blob/master/doc_automation.md#retrieve-functional-enrichment
## Code for RCy3
## Maybe figure out the GUI if it exists for cytoscape

##Install all packages before reading in the libraries
##Open cytoscape in background
##make sure STRINGapp is installed on cytoscape

##Reading in the library####
library(RCy3)
library(dplyr)
library(stringr)
library(here)
cytoscapePing()

##Input####
#setting the working directory
wd <- here()
setwd(wd)

##Before you begin, make sure to have the following:
# 1. A file in the drive folder saved as Torpor_genes.csv (containing the sheet downladed from google sheets)
# 2. A list of all the species you are working with
# 3. A list of the species ID for the reference species

#Readng in the data (After running Sorting_Gene_List.R)
tgdat_filtered <- read.csv("filtered_genelist_large.csv")
tgdat_filtered <- tgdat_filtered|>group_by(Genus_Species)
species_list <- read.csv("species_list_large.csv")

##Loop for different species####
## Makes sure the file species_list contains the reference species IDs added manually
for (i in 1:nrow(species_list)) {
  
  ##Initialisations####
  species <- species_list$Species_list[i]
  output_path <- here("cytoscape_results",species)
  species_id <- species_list$ref_species_id[i]
  cluster_network_name <- paste0(species,"--clustered")
  ##go into the outputs folder
  setwd(here::here("./cytoscape_results/"))
  
  ##create new folder for every species, if it does not already
  if (!dir.exists(species)) {
    dir.create(species)
  }
  
  ##Extracting species specific data####
  gene_list <- tgdat_filtered |> filter(Genus_Species==species)
  
  ##creating a list of all genes within the species
  string_query <- paste(gene_list$GeneID, collapse = ",")
  
  ##Creating STRING network####
  commandsPOST(paste0(
    'string protein query query="', string_query,
    '" TaxonID=',species_id,' cutoff=0.5 limit=0 newNetName=',species))
  renameNetwork(species)
  
  #Pause system in case processing has not happened
  Sys.sleep(27)
  
  ##Run MCL clustering####
  commandsRun('cluster mcl inflation_parameter=2')
  #visualise the cluster seperately
  commandsRun('clusterviz clusterview attribute="cluster"')
  #Get data for all nodes(genes), their SUID(specific ID) and cluster number
  node_table <- getTableColumns("node", c("name", "__mclCluster","SUID"))
  
  #Pause system in case processing has not happened
  Sys.sleep(27)
  
  # Calculate cluster sizes
  cluster_sizes <- node_table %>%
    group_by(`__mclCluster`) %>%
    summarise(size = n())
  
  # Set lower threshold for cluster sizes and collect node IDs for all nodes that do not belong to these clusters
  small_clusters <- cluster_sizes %>%
    filter(size < 4) %>%
    pull(`__mclCluster`)
  
  #create a list of all nodes to be deleted
  nodes_to_remove <- node_table %>%
    filter(is.na(`__mclCluster`) | `__mclCluster` %in% small_clusters)
  #select all nodes to be deleted
  selectNodes(nodes_to_remove$SUID)
  #delete all selected nodes
  deleteSelectedNodes()
  
  #select the clustered network
  setCurrentNetwork(cluster_network_name)
  
  #Get node data for all nodes post thresholding clusters
  cluster_node_table <- getTableColumns("node", c("name","__mclCluster","SUID"))
  #create a list of all unique clusters
  clusters <- unique(cluster_node_table$`__mclCluster`)
  
  #Pause system in case processing has not happened
  Sys.sleep(27)
  
  
  #Loop for processing each cluster individually
  for (cluster_id in clusters) {
    # Get genes in this cluster
    setCurrentNetwork(cluster_network_name)
    genes <- cluster_node_table %>% 
      filter(`__mclCluster` == cluster_id)
    
    # Select these nodes in Cytoscape
    selectNodes(genes$SUID,preserve.current.selection = FALSE,network=cluster_node_table)
    
    #naming the cluster based on its mcl number
    cluster_net_name <- paste0("cluster_", cluster_id)
    
    ##Create a Cluster specific subnetwork####
    createSubnetwork(nodes=genes$SUID, subnetwork.name = cluster_net_name)
    
    #Pause system in case processing has not happened
    Sys.sleep(27)
    
    # Set it as the current network
    setCurrentNetwork(cluster_net_name)
    #Reassign table_suid to null
    
    
    ## Run and export STRING enrichment ####
    ##capture output to get SUID of enrichment network
    output_str <- capture.output({commandsRun("string retrieve enrichment selectedNodesOnly=FALSE")})
    ##making the enrichment visible on cytoscape
    commandsRun('string show enrichment')
    commandsRun('string show enrichment')
    
    ##Outputting the file in csv format####
    #remove this step to prevent overwriting of files
    # First, check if file exists and remove it if necessary
    setwd(here::here("./cytoscape_results/"))
    #Set output directory
    output_file_path <- normalizePath(here("cytoscape_results", species, paste0(cluster_net_name, ".csv")), winslash = "/", mustWork = FALSE)
    
    if (file.exists(output_file_path)) {
      file.remove(output_file_path)
    }
    
    ##clean the output to obtain exact number
    table_suid <- str_extract(output_str, "(?<=: )\\d+")
    
    # If no valid SUID was found, skip export
    if (is.na(table_suid)) {
      warning(paste("No valid SUID found for cluster", cluster_net_name, "in species", species, "- skipping export."))
      next
    }
    ## create directory for output
    dir.create(dirname(output_file_path), recursive = TRUE, showWarnings = FALSE)
    ## command for exporting the enrichment table
    expt_cmd <- paste0('table export outputFile="',output_file_path,'" table=',"SUID:",table_suid)
    commandsRun(expt_cmd) 
    
    setCurrentNetwork(species)      
  }
  setwd(here())
}
print("END")


