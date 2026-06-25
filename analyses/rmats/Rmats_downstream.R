# Get Command Line Inputs
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("You need to provide at least three arguments: input_directory, output_directory and Comparison")
}

input_dir <- args[1]
output_dir <- args[2]
comparison <- args[3]

cat("Input Directory:", input_dir, "\n")
cat("Output Directory:", output_dir, "\n")
cat("Analysis Being run for comparison:",comparison, "\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# ==============================================================================
# Function Definitions
# ==============================================================================

# Function to read and annotate JCEC files
annotate_file <- function(event_type, analysis_dir){
  master_table <- fread(file.path(analysis_dir,paste0(event_type,".MATS.JCEC.txt")), header = T) %>% 
    select(!all_of(which(duplicated(names(.)))))
  
  novel_junction_ids <- fread(file.path(analysis_dir,paste0("fromGTF.novelJunction.",event_type,".txt")), header = T) %>%
    select(ID) %>% 
    unique() %>% 
    unlist()
  
  novel_splicesite_ids <- fread(file.path(analysis_dir,paste0("fromGTF.novelSpliceSite.",event_type,".txt")), header = T) %>% 
    select(ID) %>% 
    unique() %>% 
    unlist()
  
  junction_counts_ids <- fread(file.path(analysis_dir,paste0(event_type,".MATS.JC.txt")), header = T,) %>% 
    select(!all_of(which(duplicated(names(.))))) %>% 
    select(ID) %>% 
    unique() %>% 
    unlist(use.names = F) 
  # The above code is required to make the dataframe into a list and then remove the name of the vector.
  
  master_table_annotated <- master_table %>% 
    mutate(annotation = case_when(ID %in% novel_junction_ids ~ "Novel Junction",
                                  ID %in% novel_splicesite_ids ~ "Novel Splice Site",
                                  TRUE ~ "GTF Annotated" # All other cases will be labeled as this
    )) %>% 
    mutate(readType = ifelse(ID %in% junction_counts_ids, "JC", "JCEC")) %>% 
    mutate(eventType = event_type) %>% 
    rowwise() %>%   # Apply functions row-wise
    mutate(
      SAMPLE_1_COUNTS = numeric_from_comma_string(IJC_SAMPLE_1) + numeric_from_comma_string(SJC_SAMPLE_1),
      SAMPLE_2_COUNTS = numeric_from_comma_string(IJC_SAMPLE_2) + numeric_from_comma_string(SJC_SAMPLE_2),
      IJC_DIFF = numeric_from_comma_string(IJC_SAMPLE_1) - numeric_from_comma_string(IJC_SAMPLE_2),
      SJC_DIFF = numeric_from_comma_string(SJC_SAMPLE_1) - numeric_from_comma_string(SJC_SAMPLE_2)
    ) %>% 
    ungroup()  # Ungroup to return to a standard data frame
  
  return(master_table_annotated)
}

# Function to seperate Comma Seperated Rows in the Table
numeric_from_comma_string <- function(comma_string, return_mean = TRUE) {
  numeric_values <- as.numeric(unlist(strsplit(comma_string, split = ',')))
  
  if (return_mean) {
    return(mean(numeric_values))
  } else {
    return(numeric_values)
  }
}

# Function to find min PSI for an event
min_psi_from_row <- function(row) {
  sample_1_psi_values <- numeric_from_comma_string(row["IncLevel1"], return_mean = FALSE)
  sample_2_psi_values <- numeric_from_comma_string(row["IncLevel2"], return_mean = FALSE)
  min_sample_1 <- min(sample_1_psi_values)
  min_sample_2 <- min(sample_2_psi_values)
  
  return(min(min_sample_1, min_sample_2))
}

# Function to find Max PSI for an event
max_psi_from_row <- function(row) {
  sample_1_psi_values <- numeric_from_comma_string(row["IncLevel1"], return_mean = FALSE)
  sample_2_psi_values <- numeric_from_comma_string(row["IncLevel2"], return_mean = FALSE)
  max_sample_1 <- max(sample_1_psi_values)
  max_sample_2 <- max(sample_2_psi_values)
  
  return(max(max_sample_1, max_sample_2))
}

filtering <- function(table){
  table_filtered <- table %>% 
    filter(FDR <= 0.05) %>% 
    filter(IncLevelDifference <= 0) %>% 
    filter(SAMPLE_1_COUNTS > 50 & SAMPLE_2_COUNTS > 50) %>% 
    arrange(desc(IncLevelDifference))
  
  return(table_filtered)
}

# ==============================================================================
# Variable Definitions
# ==============================================================================
splicing_events <- c("SE","RI","MXE","A3SS","A5SS")

output_unfiltered <- list()
output_filtered <- list()

# ==============================================================================
# Script
# ==============================================================================
for (event in splicing_events) {
  table <- annotate_file(event_type = event, analysis_dir = input_dir) 
  output_unfiltered[[event]] <- table
  output_filtered[[event]] <- table %>% filtering()
  
}

output_unfiltered_result <- do.call(rbind, output_unfiltered)
output_filtered_result <- do.call(rbind, output_filtered)

output_unfiltered_path <- file.path(output_dir,paste0(comparison,"_table_allevents.csv"))
output_filtered_path <- file.path(output_dir,paste0(comparison,"_filtered_table_allevents.csv"))

write.csv(output_unfiltered_result, output_unfiltered_path, row.names = FALSE)
write.csv(output_filtered_result, output_filtered_path)

