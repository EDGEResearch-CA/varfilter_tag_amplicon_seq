# Script Summary ----
## Function call_tag_amp_mutations to apply a custom mutation confidence scoring system. 
## Application of thresholds to determine the confidence of each variant being a real mutation.
## IT considers formalin-induced artifacts, which can falsely look like mutations.


## Summary of rules ----
## Minimum of 2 supporting amplicons
## VAF (Variant Allele Frequency) >= 2%
## For 2 supporting amplicons: Read support for amplicon >= 20  
### except when the alteration is C>T/G>A then it must be >= 50
## For 3 supporting amplicons: Read support for amplicon >= 3  
### except when the alteration is C>T/G>A then it must be >= 20

# Load required packages ----
library(tidyverse)     # includes dplyr, ggplot2, etc. for data manipulation and visualization
library(here)          # helps create platform-independent file paths
library(readxl)        # used to read Excel files (.xlsx)
library(writexl)       # used to write Excel files (.xlsx)

## Define a function to assess whether a mutation is real ----
## Includes a custom system of rules
assess_mutation <- function(single_row, vaf_cutoff = 2) { #single_row represents a single variant for one single case (1 row)
  
  # Extract values from each row of the data frame - correspondent to one single variant for one case
  vaf <- as.numeric(single_row[["VAF"]]) #as.numeric to ensure value is treated as a number 
  # assumes that VAF is already converted to % (if not, add * 100)
  
  n_support_amplicons <- as.numeric(single_row[["SupportingAmplicons"]])  
  # how many amplicons support the mutation (3 maximum possible by current strategy)
  
  # Extract read counts from the 3 amplicon columns 
  amplicon_reads_vector <- c(                       # creates a vector containg the read os support from the 3 amplicons
    as.numeric(single_row[["AmpliconProbeHits_1"]]),
    as.numeric(single_row[["AmpliconProbeHits_2"]]),
    as.numeric(single_row[["AmpliconProbeHits_3"]]) 
  )
  
  ### Deamination alteration that is a signature of formalin-induced artifacts ----
  # Check if this mutation is a C>T or G>A substitution
  # Formalin can chemically changes in DNA in FFPE samples
  # Cytosine deamination can be caused by formalin leading to false positives in mutation calling
  is_ffpe_deamination <- str_detect(single_row[["VariantDescription"]], "C>T|G>A")
  
  # Remove NAs and zeroes to reflect only real amplicon support.
  nonzero_reads <- amplicon_reads_vector[amplicon_reads_vector > 0]
  
  
  # SAFETY CHECK: Number of valid (non-zero) reads must match number of supporting amplicons
  if (length(nonzero_reads) !=  n_support_amplicons) {
    stop(paste0(
      "❌ Fatal Error: mismatch in amplicon support!\n",
      "  Sample: ", single_row[["Sample"]], "\n",
      "  Variant: ", single_row[["VariantDescription"]], "\n",
      "  Number of supporting amplicons: ", n_support_amplicons, "\n",
      "  Number of Non-zero reads: ", length(nonzero_reads), "\n"
    ))
    
  }
  
  # Early exit: If not enough amplicon support or VAF is too low, return FALSE (variant is not a real mutation)
  if (n_support_amplicons < 2 || vaf < vaf_cutoff) {
    return(FALSE)
  }
  
  # Use all non-zero read values as reads_to_check
  reads_to_check <- nonzero_reads
  
  ### Apply confidence thresholds based on rules ----
  if (n_support_amplicons == 2) {
    if (is_ffpe_deamination) {
      return(all(reads_to_check >= 50)) # stricter threshold for C>T or G>A
    } else {
      return(all(reads_to_check >= 20))
    }
  } else if (n_support_amplicons >= 3) {
    if (is_ffpe_deamination) {
      return(all(reads_to_check >= 20)) # stricter threshold for C>T or G>A
    } else {
      return(all(reads_to_check >= 3))
    }
  }
  
  # Fallback (should not occur under normal circumstances)
  warning(
    paste0(
      "Undefined mutation status given current system of rules!\n",
      "  Sample: ", single_row[["Sample"]], "\n",
      "  Variant: ", single_row[["Variant"]], "\n",
      "  VariantDescription: ", single_row[["VariantDescription"]], "\n"
    )
  )
  return(NA)  # returns FALSE if none of the above conditions are met
}

# Define a reusable function that takes a data set as input (pre-processed data from variant calling) ----
# This function uses the previous 'assess_mutation' function defined above
call_tag_amp_mutations <- function(
    input_filename, #Change to input_filename = relative path of the input file you want to process
    output_prefix, #edit as needed to personalize output file names
    output_folder = "output",
    vaf_cutoff = 2
) {

  ## Start files and folders ----
  ## Create input and output folders if they don't exist ----
  dir.create(here::here(output_folder), showWarnings = FALSE) #suppresses “folder already exists” messages
  
  ## Read the Excel input file containing pre-processed variant data ----
  variant_support_data <- read_excel( # it reads and store as an R data frame
    input_filename) #input_filename defined in the main function
  
  
  ## Loop through each row and assess the mutation status ----
  mutation_calls <- logical(nrow(variant_support_data))  # Initialize an empty logical vector to store results
  
  for (i in 1:nrow(variant_support_data)) {
    single_row <- variant_support_data[i, ]  # Get one row (one variant) from the data frame
    mutation_calls[i] <- assess_mutation(single_row)  # Call the assess_mutation function and store the result
  }
  
  # Add the results as a new column to your original data ----
  variant_support_data$IsRealMutation <- mutation_calls  # Add a new column to the data frame with mutation status
  
  
  # Exporting final results ----
  ## Write the full processed data to a new Excel file 
  write_xlsx(
    variant_support_data,  
    here::here(output_folder, paste0("mutation_calls_tag_amp_seq_", output_prefix, ".xlsx"))
  )
  
  ## Save a spreadsheet with true mutations only 
  variant_support_data_true_mutations <- variant_support_data %>% 
    dplyr::filter(IsRealMutation)
  
  write_xlsx(
    variant_support_data_true_mutations,
    here::here(output_folder, paste0("true_mutation_calls_tag_amp_seq_", output_prefix, ".xlsx"))
  )
}

# Calling function ----
## Input bianca_data_2025
call_tag_amp_mutations(
    input_filename = here::here("input", "DA-1987_variant_support_final.xlsx"),
    output_prefix = "bianca",
    output_folder = here::here("output"),
    vaf_cutoff = 2
)
