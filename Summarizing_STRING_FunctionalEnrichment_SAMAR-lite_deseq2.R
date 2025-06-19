#####################SAMAR_lite################
############ DESEQ2 #################
############# Biological Processes (BP) ###################
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ComplexUpset)
library(patchwork)
library(ggVennDiagram)
library(ggplot2)
library(purrr)
library(here)

# Set base directory with all species subdirectories (e.g., Gut1, Gut2, etc.)
base_dir <- here("./cytoscape_results")

# Get all immediate species subdirectories
species_dirs <- list.dirs(base_dir, recursive = FALSE)

# Initialize list to store GO BP descriptions by species
go_descriptions_by_species <- list()

# Loop through species directories
for (species_crnt_dir in species_dirs) {
  # Get species name from the folder path
  species_name <- basename(species_crnt_dir)
  
  # List all CSV files in the species directory
  file_list <- list.files(path = species_crnt_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a vector to store all GO BP terms for this species
  species_go_terms <- c()
  
  for (file_path in file_list) {
    df <- read_delim(file_path, col_types = cols(.default = "c"))
    
    if ("category" %in% names(df) && "description" %in% names(df)) {
      terms <- df %>%
        filter(category == "GO Biological Process") %>%
        pull(description)
      
      species_go_terms <- c(species_go_terms, terms)
    } else {
      warning(paste("Skipping file due to missing columns:", file_path))
    }
  }
  
  # Store unique sorted terms for this species
  go_descriptions_by_species[[species_name]] <- unique(species_go_terms)
}

go_long <- enframe(go_descriptions_by_species, name = "species", value = "GO_Term") %>%
  unnest(GO_Term)

# Create wide binary matrix: 1 if term appears in species, else 0
go_binary_matrix <- go_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = species, values_from = value, values_fill = 0)

# Preview the binary matrix
print(head(go_binary_matrix))

go_binary_filtered <- go_binary_matrix %>%
  rowwise() %>%
  mutate(species_count = sum(c_across(-GO_Term))) %>%
  ungroup() %>%
  filter(species_count >= 5)



SAMAR_lite_DESEQ2_StringEnrichment<-go_binary_filtered
write.csv(SAMAR_lite_DESEQ2_StringEnrichment,here("./Network_Analysis/String_FunctionalEnrichment/SAMAR_lite_DESEQ2_StringEnrichment_BP.csv"))


# Upset plot
a<-ComplexUpset::upset(
  go_binary_filtered %>% select(-species_count),
  intersect = names(go_descriptions_by_species),
  name = "GO BP Terms_speciess_DESEq2 (≥5 speciess)"
)

# Venn Diagram
go_binary_filtered_1<-go_binary_filtered[,-7]

go_list <- lapply(names(go_binary_filtered_1)[-1], function(species) {
  go_binary_filtered_1$GO_Term[go_binary_filtered_1[[species]] == 1]
})
names(go_list) <- names(go_binary_filtered_1)[-1]


b<-ggVennDiagram(go_list, 
                 label_alpha = 0, 
                 label = "count", 
                 edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "firebrick2") +
  theme(legend.position = "none")

#### Find out what GO terms these are? 
term_counts <- table(unlist(go_descriptions_by_species))

# Filter for terms that appear in at least 4 speciess
common_terms_3_or_more <- names(term_counts[term_counts >= 3])
common_terms_4_or_more <- names(term_counts[term_counts >= 4])
common_terms_5_or_more <- names(term_counts[term_counts >= 5])

common_terms_3_or_more
common_terms_4_or_more
common_terms_5_or_more

# Optionally reorder GO terms by species_count
c<- go_binary_filtered %>%
  ggplot(aes(x = reorder(GO_Term, species_count), y = species_count)) +
  geom_col(fill = "#5ab4ac") +
  coord_flip() +  # Flip for readability
  labs(
    title = "GO BP Terms Shared Across speciess",
    x = "",
    y = "Number of speciess"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

library(patchwork)
layout <- patchwork::area(
  t = c(1, 2, 1),  # top row for a and c, 2nd row for b
  l = c(1, 1, 2),  # left col for a and b, right col for c
  b = c(1, 2, 2),
  r = c(1, 1, 2)
)

# Assign plots to layout
layout_design <- patchwork::wrap_plots(
  a, b, c,
  design = layout
)

layout_design
ggsave(here("./Network_Analysis/SAMAR_lite_DESEQ2_BP_FunctionalEnrichment_AcrossspeciesTypes.jpg"), plot = last_plot(), device = NULL, path = NULL, width = 21, height = 15, dpi = 300)

#####################SAMAR_lite################
############ DESEQ2 #################
################ Cellular Component (CC) ###################
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ComplexUpset)
library(patchwork)
library(ggVennDiagram)
library(ggplot2)
library(purrr)
library(here)



# Set base directory with all species subdirectories (e.g., Gut1, Gut2, etc.)
base_dir <- here("./Network_Analysis/speciess_DESEq2")

# Get all immediate species subdirectories
species_dirs <- list.dirs(base_dir, recursive = FALSE)

# Initialize list to store GO BP descriptions by species
go_descriptions_by_species <- list()

# Loop through species directories
for (species_crnt_dir in species_dirs) {
  # Get species name from the folder path
  species_name <- basename(species_crnt_dir)
  
  # List all CSV files in the species directory
  file_list <- list.files(path = species_crnt_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a vector to store all GO BP terms for this species
  species_go_terms <- c()
  
  for (file_path in file_list) {
    df <- read_delim(file_path, col_types = cols(.default = "c"))
    
    if ("category" %in% names(df) && "description" %in% names(df)) {
      terms <- df %>%
        filter(category == "GO Cellular Component") %>%
        pull(description)
      
      species_go_terms <- c(species_go_terms, terms)
    } else {
      warning(paste("Skipping file due to missing columns:", file_path))
    }
  }
  
  # Store unique sorted terms for this species
  go_descriptions_by_species[[species_name]] <- unique(species_go_terms)
}

go_long <- enframe(go_descriptions_by_species, name = "species", value = "GO_Term") %>%
  unnest(GO_Term)

# Create wide binary matrix: 1 if term appears in species, else 0
go_binary_matrix <- go_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = species, values_from = value, values_fill = 0)

# Preview the binary matrix
print(head(go_binary_matrix))

go_binary_filtered <- go_binary_matrix %>%
  rowwise() %>%
  mutate(species_count = sum(c_across(-GO_Term))) %>%
  ungroup() %>%
  filter(species_count >= 2)

SAMAR_lite_DESEQ2_StringEnrichment<-go_binary_filtered
write.csv(SAMAR_lite_DESEQ2_StringEnrichment,here("./Network_Analysis/String_FunctionalEnrichment/SAMAR_lite_DESEQ2_StringEnrichment_CC.csv"))

# Upset plot
a<-ComplexUpset::upset(
  go_binary_filtered %>% select(-species_count),
  intersect = names(go_descriptions_by_species),
  name = "GO CC Terms_speciess_DESEq2 (≥2 speciess)"
)

# Venn Diagram
go_binary_filtered_1<-go_binary_filtered[,-9]

go_list <- lapply(names(go_binary_filtered_1)[-1], function(species) {
  go_binary_filtered_1$GO_Term[go_binary_filtered_1[[species]] == 1]
})
names(go_list) <- names(go_binary_filtered_1)[-1]


b<-ggVennDiagram(go_list, 
                 label_alpha = 0, 
                 label = "count", 
                 edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "firebrick2") +
  theme(legend.position = "none")

#### Find out what GO terms these are? 
term_counts <- table(unlist(go_descriptions_by_species))

# Filter for terms that appear in at least 4 speciess
common_terms_3_or_more <- names(term_counts[term_counts >= 3])
common_terms_4_or_more <- names(term_counts[term_counts >= 4])
common_terms_5_or_more <- names(term_counts[term_counts >= 5])

common_terms_3_or_more
common_terms_4_or_more
common_terms_5_or_more

# Optionally reorder GO terms by species_count
c<- go_binary_filtered %>%
  ggplot(aes(x = reorder(GO_Term, species_count), y = species_count)) +
  geom_col(fill = "#5ab4ac") +
  coord_flip() +  # Flip for readability
  labs(
    title = "GO CC Terms Shared Across speciess",
    x = "",
    y = "Number of speciess"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

library(patchwork)
layout <- patchwork::area(
  t = c(1, 2, 1),  # top row for a and c, 2nd row for b
  l = c(1, 1, 2),  # left col for a and b, right col for c
  b = c(1, 2, 2),
  r = c(1, 1, 2)
)

# Assign plots to layout
layout_design <- patchwork::wrap_plots(
  a, b, c,
  design = layout
)

layout_design
ggsave(here("./Network_Analysis/SAMAR_lite_DESEQ2_CC_FunctionalEnrichment_AcrossspeciesTypes.jpg"), plot = last_plot(), device = NULL, path = NULL, width = 21, height = 15, dpi = 300)

#####################SAMAR_lite################
############ DESEQ2 #################
################ Molecular Function (MF) ###################
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ComplexUpset)
library(patchwork)
library(ggVennDiagram)
library(ggplot2)
library(purrr)
library(here)

# Set base directory with all species subdirectories (e.g., Gut1, Gut2, etc.)
base_dir <- here("./Network_Analysis/speciess_DESEq2")

# Get all immediate species subdirectories
species_dirs <- list.dirs(base_dir, recursive = FALSE)

# Initialize list to store GO BP descriptions by species
go_descriptions_by_species <- list()

# Loop through species directories
for (species_crnt_dir in species_dirs) {
  # Get species name from the folder path
  species_name <- basename(species_crnt_dir)
  
  # List all CSV files in the species directory
  file_list <- list.files(path = species_crnt_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a vector to store all GO BP terms for this species
  species_go_terms <- c()
  
  for (file_path in file_list) {
    df <- read_delim(file_path, col_types = cols(.default = "c"))
    
    if ("category" %in% names(df) && "description" %in% names(df)) {
      terms <- df %>%
        filter(category == "GO Molecular Function") %>%
        pull(description)
      
      species_go_terms <- c(species_go_terms, terms)
    } else {
      warning(paste("Skipping file due to missing columns:", file_path))
    }
  }
  
  # Store unique sorted terms for this species
  go_descriptions_by_species[[species_name]] <- unique(species_go_terms)
}

go_long <- enframe(go_descriptions_by_species, name = "species", value = "GO_Term") %>%
  unnest(GO_Term)

# Create wide binary matrix: 1 if term appears in species, else 0
go_binary_matrix <- go_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = species, values_from = value, values_fill = 0)

# Preview the binary matrix
print(head(go_binary_matrix))

go_binary_filtered <- go_binary_matrix %>%
  rowwise() %>%
  mutate(species_count = sum(c_across(-GO_Term))) %>%
  ungroup() %>%
  filter(species_count >= 2)

SAMAR_lite_DESEQ2_StringEnrichment<-go_binary_filtered
write.csv(SAMAR_lite_DESEQ2_StringEnrichment,here("./Network_Analysis/String_FunctionalEnrichment/SAMAR_lite_DESEQ2_StringEnrichment_MF.csv"))

# Upset plot
a<-ComplexUpset::upset(
  go_binary_filtered %>% select(-species_count),
  intersect = names(go_descriptions_by_species),
  name = "GO MF Terms_speciess_DESEq2 (≥2 speciess)"
)

# Venn Diagram
go_binary_filtered_1<-go_binary_filtered[,-9]

go_list <- lapply(names(go_binary_filtered_1)[-1], function(species) {
  go_binary_filtered_1$GO_Term[go_binary_filtered_1[[species]] == 1]
})
names(go_list) <- names(go_binary_filtered_1)[-1]


b<-ggVennDiagram(go_list, 
                 label_alpha = 0, 
                 label = "count", 
                 edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "firebrick2") +
  theme(legend.position = "none")

#### Find out what GO terms these are? 
term_counts <- table(unlist(go_descriptions_by_species))

# Filter for terms that appear in at least 4 speciess
common_terms_3_or_more <- names(term_counts[term_counts >= 3])
common_terms_4_or_more <- names(term_counts[term_counts >= 4])
common_terms_5_or_more <- names(term_counts[term_counts >= 5])

common_terms_3_or_more
common_terms_4_or_more
common_terms_5_or_more

# Optionally reorder GO terms by species_count
c<- go_binary_filtered %>%
  ggplot(aes(x = reorder(GO_Term, species_count), y = species_count)) +
  geom_col(fill = "#5ab4ac") +
  coord_flip() +  # Flip for readability
  labs(
    title = "GO MF Terms Shared Across speciess",
    x = "",
    y = "Number of speciess"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

library(patchwork)
layout <- patchwork::area(
  t = c(1, 2, 1),  # top row for a and c, 2nd row for b
  l = c(1, 1, 2),  # left col for a and b, right col for c
  b = c(1, 2, 2),
  r = c(1, 1, 2)
)

# Assign plots to layout
layout_design <- patchwork::wrap_plots(
  a, b, c,
  design = layout
)

layout_design
ggsave(here("./rsults"), plot = last_plot(), device = NULL, path = NULL, width = 21, height = 15, dpi = 300)


#####################SAMAR_lite################
############ DESEQ2 #################
################ KEGG Pathways ###################
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ComplexUpset)
library(patchwork)
library(ggVennDiagram)
library(ggplot2)
library(purrr)
library(here)

# Set base directory with all species subdirectories (e.g., Gut1, Gut2, etc.)
base_dir <- here("./Network_Analysis/speciess_DESEq2")

# Get all immediate species subdirectories
species_dirs <- list.dirs(base_dir, recursive = FALSE)

# Initialize list to store KEGG descriptions by species
go_descriptions_by_species <- list()

# Loop through species directories
for (species_crnt_dir in species_dirs) {
  # Get species name from the folder path
  species_name <- basename(species_crnt_dir)
  
  # List all CSV files in the species directory
  file_list <- list.files(path = species_crnt_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize a vector to store all GO BP terms for this species
  species_go_terms <- c()
  
  for (file_path in file_list) {
    df <- read_delim(file_path, col_types = cols(.default = "c"))
    
    if ("category" %in% names(df) && "description" %in% names(df)) {
      terms <- df %>%
        filter(category == "KEGG Pathways") %>%
        pull(description)
      
      species_go_terms <- c(species_go_terms, terms)
    } else {
      warning(paste("Skipping file due to missing columns:", file_path))
    }
  }
  
  # Store unique sorted terms for this species
  go_descriptions_by_species[[species_name]] <- unique(species_go_terms)
}

go_long <- enframe(go_descriptions_by_species, name = "species", value = "GO_Term") %>%
  unnest(GO_Term)

# Create wide binary matrix: 1 if term appears in species, else 0
go_binary_matrix <- go_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = species, values_from = value, values_fill = 0)

# Preview the binary matrix
print(head(go_binary_matrix))

go_binary_filtered <- go_binary_matrix %>%
  rowwise() %>%
  mutate(species_count = sum(c_across(-GO_Term))) %>%
  ungroup() %>%
  filter(species_count >= 2)

SAMAR_lite_DESEQ2_StringEnrichment<-go_binary_filtered
write.csv(SAMAR_lite_DESEQ2_StringEnrichment,here("./Network_Analysis/String_FunctionalEnrichment/SAMAR_lite_DESEQ2_StringEnrichment_KEGG.csv"))

# Upset plot
a<-ComplexUpset::upset(
  go_binary_filtered %>% select(-species_count),
  intersect = names(go_descriptions_by_species),
  name = "KEGG Terms_speciess_DESEq2 (≥2 speciess)"
)

# Venn Diagram
go_binary_filtered_1<-go_binary_filtered[,-9]

go_list <- lapply(names(go_binary_filtered_1)[-1], function(species) {
  go_binary_filtered_1$GO_Term[go_binary_filtered_1[[species]] == 1]
})
names(go_list) <- names(go_binary_filtered_1)[-1]


b<-ggVennDiagram(go_list, 
                 label_alpha = 0, 
                 label = "count", 
                 edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "firebrick2") +
  theme(legend.position = "none")

#### Find out what GO terms these are? 
term_counts <- table(unlist(go_descriptions_by_species))

# Filter for terms that appear in at least 4 speciess
common_terms_3_or_more <- names(term_counts[term_counts >= 3])
common_terms_4_or_more <- names(term_counts[term_counts >= 4])
common_terms_5_or_more <- names(term_counts[term_counts >= 5])

common_terms_3_or_more
common_terms_4_or_more
common_terms_5_or_more

# Optionally reorder GO terms by species_count
c<- go_binary_filtered %>%
  ggplot(aes(x = reorder(GO_Term, species_count), y = species_count)) +
  geom_col(fill = "#5ab4ac") +
  coord_flip() +  # Flip for readability
  labs(
    title = "GO KEGG Terms Shared Across speciess",
    x = "",
    y = "Number of speciess"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

library(patchwork)
layout <- patchwork::area(
  t = c(1, 2, 1),  # top row for a and c, 2nd row for b
  l = c(1, 1, 2),  # left col for a and b, right col for c
  b = c(1, 2, 2),
  r = c(1, 1, 2)
)

# Assign plots to layout
layout_design <- patchwork::wrap_plots(
  a, b, c,
  design = layout
)

layout_design
ggsave(here("./Network_Analysis/SAMAR_lite_DESEQ2_KEGG_FunctionalEnrichment_AcrossspeciesTypes.jpg"), plot = last_plot(), device = NULL, path = NULL, width = 21, height = 15, dpi = 300)


