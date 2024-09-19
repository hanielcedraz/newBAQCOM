#' @export
#' @name loadSamplesFile
#' @title Loading sample file
#' @author Haniel Cedraz
#' @details December 2021
#' @usage
#' loadSamplesFile(file, reads_folder, column = "SAMPLE_ID", libraryType = "pairEnd")
#' @description
#' Function to load the sample file
#'
#' @param file
#' \code{Character.} The filename of the sample file. Default samples.txt.
#' @param reads_folder
#' \code{Character.} Directory where the raw sequence data is stored. Default 00-Fastq.
#' @param column
#' \code{Character.} Column name from the sample sheet to use as read folder names. Default SAMPLE_ID
#' @param libraryType
#' \code{Character.} The library type to use. Available: 'pairEnd' or 'singleEnd'. Default pairEnd
#' @importFrom tools file_ext
#' @importFrom dplyr mutate %>%
#' @importFrom glue glue
#' @importFrom readr read_table
#' @importFrom purrr modify_if
#' @export



## loadSampleFile
loadSamplesFile <- function(file = "samples.txt", column = "SAMPLE_ID", libraryType = "pairEnd"){
  ## debug
  #file = opt$samplesFile; reads_folder = opt$Raw_Folder; column = opt$samplesColumn
  ##

  aceptedLibraryType <- c("pairEnd", "singleEnd")
  if (!libraryType %in% aceptedLibraryType) {
    cli_abort(c("x" = "Library type ({libraryType}) not found, please provide one of 'pairEnd' or 'singleEnd'"))
    
  }

  if (!file.exists(file)) {
    cli_abort(c("x" = "Sample file {file} does not exist"))

  }
  
  
  ### column SAMPLE_ID should be the sample name
  ### rows can be commented out with #
  if (libraryType == "singleEnd") {
    targets <- read_table(file, col_types = list("c", "c"), comment = "#") %>%
      modify_if(~is.double(.), ~as.character(.)) %>%
      as.data.frame()

  } else if (libraryType == "pairEnd") {
    targets <- read_table(file, col_types = list("c", "c", "c"), comment = "#") %>%
      modify_if(~is.double(.), ~as.character(.)) %>%
      as.data.frame()
  }


  if (libraryType == "pairEnd") {
    if (!all(c(column, "Read_1", "Read_2") %in% colnames(targets))) {
      
      cli_abort(c("x" = "Expecting three columns SAMPLE_ID, Read_1 and Read_2 in samples file (tab-delimited)"))

    }
  }
  
  cat("\n\n");cli_alert_success("{file} contains {nrow(targets)} samples to process");cat("\n\n")
  
  return(targets)
}

