####Automating Cytoscape Workflow
#### Started by Jay Phadke 9-06-2025

###Useful resources: 
##http://127.0.0.1:1234/v1/swaggerUI/swagger-ui/index.html?url=http://127.0.0.1:1234/v1/commands/swagger.json#/
## Useful AI : https://www.phind.com/
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
here()

#Readng in the data (After running Sorting_Gene_List.R)
tgdat_filtered <- read.csv("filtered_genelist.csv")
tgdat_filtered <- tgdat_filtered|>group_by(Genus_Species)

##Initialisations####
species <- "Urocitellus_parryii"
output_path <- here("cytoscape_results",Species)
species_id <- "10090"
cluster_network_name <- paste0(species,"--clustered")

##Extracting species specific data####
trial <- tgdat_filtered |> filter(Genus_Species==species)

#creating a list of all genes within the species
string_query <- paste(trial$GeneID, collapse = ",")

##Creating STRING network####
commandsRun(paste0(
  'string protein query query="', string_query,
  '" TaxonID=10090 cutoff=0.5 limit=0 newNetName=',species))
renameNetwork(species)


##Run MCL clustering####
commandsRun('cluster mcl inflation_parameter=3')
#visualise the cluster seperately
commandsRun('clusterviz clusterview attribute="cluster"')
#Get data for all nodes(genes), their SUID(specific ID) and cluster number
node_table <- getTableColumns("node", c("name", "__mclCluster","SUID"))

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
cluster_node_table <- getTableColumns("node", c("name", "__mclCluster","SUID"))
#create a list of all unique clusters
clusters <- unique(cluster_node_table$`__mclCluster`)

#Loop for processing each cluster individually
for (cluster_id in clusters) {
  
  #initialisations (cluster specific)
  output_file <- paste0(output_path,species,"/cluster",cluster_id,".csv")
  
  # Get genes in this cluster
  setCurrentNetwork(cluster_network_name)
  genes <- cluster_node_table %>% 
    filter(`__mclCluster` == cluster_id)
  
  # Select these nodes in Cytoscape
  selectNodes(genes$SUID,preserve.current.selection = FALSE,network=cluster_node_table)
  
  #naming the cluster basedon its mcl number
  cluster_net_name <- paste0("cluster_", cluster_id)
  
  ##Create a Cluster specific subnetwork####
  createSubnetwork(nodes=genes$SUID, subnetwork.name = cluster_net_name)
  
  # Set it as the current network
  setCurrentNetwork(cluster_net_name)
  
  ## Run STRING enrichment ####
  #commandsRun('string retrieve enrichment')
  #commandsRun('string show enrichment')
  
  setCurrentNetwork(species)
}
