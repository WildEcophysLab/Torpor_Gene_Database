####Querying on GEO database for downloading metadata
#### Started by Jay Phadke 29-06-2025

###mainly used Chatgpt for code. look up rentrez API for further details
##Install all packages before reading in the libraries

## Initializing Libraries ####
library(rentrez)
library(dplyr)
library(here)

## Setting the Search query ####
query <- '((torpor[MeSH Terms] OR torpor[All Fields]) OR (hibernation[MeSH Terms] OR hibernation[All Fields])) AND "mammals"[porgn]'
## Searching GEO database ####
search_result <- entrez_search(db = "gds", term = query, retmax = 2000)
## Fetching list of all search result GEO IDs
geo_ids <- search_result$ids
cat("Found", length(geo_ids), "GDS entries.\n")

## Split IDs into batches to avoid URI Too Long errors
batch_ids <- function(id_vector, batch_size = 200) {
  split(id_vector, ceiling(seq_along(id_vector) / batch_size))
}
geo_id_batches <- batch_ids(geo_ids, batch_size = 200)

## Fetch summaries ####
summaries <- unlist(
  lapply(geo_id_batches, function(batch) {
    Sys.sleep(5)  # polite pause
    tryCatch(
      entrez_summary(db = "gds", id = batch),
      error = function(e) {
        warning("Batch failed: ", conditionMessage(e))
        return(NULL)
      }
    )
  }),
  recursive = FALSE
)

## Convert summary into a dataframe, withoot mismatching rows and columns (flattening) ####
flatten_summary <- function(s) {
  flat <- lapply(s, function(x) {
    if (is.null(x)) {
      return(NA_character_)
    } else if (is.atomic(x) && length(x) == 1) {
      return(as.character(x))
    } else if (is.atomic(x)) {
      return(paste(as.character(x), collapse = "; "))
    } else if (is.data.frame(x)) {
      # Flatten a data frame to "row1 | row2" format
      df_flat <- apply(x, 1, function(row) paste(row, collapse = ":"))
      return(paste(df_flat, collapse = " | "))
    } else if (is.list(x)) {
      return(paste(sapply(x, as.character), collapse = "; "))
    } else {
      return(as.character(x))
    }
  })
  
  # Convert the flattened list into a 1-row data frame
  as.data.frame(flat, stringsAsFactors = FALSE)
}

## Flatten all summaries into a data frame ####
summary_df_list <- lapply(summaries, flatten_summary)
summary_df_clean <- Filter(Negate(is.null), summary_df_list)

## Combine all rows into a full metadata frame
metadata_df <- bind_rows(summary_df_clean)

## Set Working directory to where you want to save the file
setwd(here())
## Write to CSV ####
write.csv(metadata_df, "metadata_output.csv", row.names = FALSE)
cat("Saved", nrow(metadata_df), "summaries to CSV ✅\n")

print("END")