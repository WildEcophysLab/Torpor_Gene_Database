cluster_id <- 1

# Get genes in this cluster
setCurrentNetwork(cluster_network_name)
genes <- cluster_node_table %>% 
  filter(`__mclCluster` == cluster_id)

# Select these nodes in Cytoscape
selectNodes(genes$SUID,preserve.current.selection = FALSE,network=cluster_node_table)

output_file <- paste0(output_path,species,"/cluster",cluster_id,".csv")
# Create a subnetwork for this cluster

cluster_net_name <- paste0("cluster_", cluster_id)
createSubnetwork(nodes=genes$SUID, subnetwork.name = cluster_net_name)

# Set it as the current network
setCurrentNetwork(cluster_net_name)

# Run STRING enrichment (must have stringApp installed)
commandsRun('string show enrichment')

setCurrentNetwork(network_name)


temp_sink <- NULL
cmd <- paste0('table export table=',final_suid_no,'outputFile="', output_file, '"')
print(cmd)
commandsRun(cmd)


table_suid_no <- str_extract(temp_sink, "\\d+")
final_suid_no <- paste0('"SUID:',table_suid_no,'"')
print(final_suid_no)


export table STRING Enrichment: All to C:\Users\Jay Phadke\OneDrive\Desktop\WEL\u_p_c1.csv







# Create a single folder
dir.create("my_new_folder")

# Create multiple folders at once
dir.create(c("folder1", "folder2", "folder3"))

# Create nested folders (parent/child structure)
dir.create("parent/child", recursive = TRUE)

