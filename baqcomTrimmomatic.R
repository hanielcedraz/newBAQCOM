#!/usr/bin/env Rscript

#trimmomatic_dir <- paste('XXX/Trimmomatic-0.39/')
trimmomatic_dir <- paste('../BAQCOM/Trimmomatic-0.39/')
trimmomatic <- paste(trimmomatic_dir, 'trimmomatic-0.39.jar', sep = "")

########################################
### LOADING PACKAGES
########################################
if (suppressPackageStartupMessages(!require(pacman))) {
    suppressPackageStartupMessages(install.packages("pacman", repos = "http://cran.us.r-project.org", dependencies = TRUE, lib = "~/R/lib"))
}
suppressPackageStartupMessages(
    p_load(
        tools, 
        parallel, 
        optparse, 
        stringr,
        #baqcomPackage, 
        dplyr, 
        data.table,
        cli,
        glue,
        purrr,
        readr,
        rlang
    )
)




########################################
### SETING PARAMETERS
########################################
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
    make_option(
        opt_str = c("-f", "--file"), 
        type = "character",
        default = "samples.txt",
        help = "The filename of the sample file [default %default]",
        dest = "samplesFile"
    ),
    make_option(
        opt_str = c("-d", "--directory"), 
        type = "character",
        default = "00-Fastq",
        help = "Directory where the raw sequence data is stored [default %default]",
        dest = "Raw_Folder"
    ),
    make_option(
        opt_str = c("-c", "--column"), 
        type = "character", 
        default = "SAMPLE_ID",
        help = "Column name from the sample sheet to use as read folder names [default %default]",
        dest = "samplesColumn"
    ),
    make_option(
        opt_str = c("-l", "--fastqc"), 
        action = 'store_true',
        type = "logical",
        default = FALSE,
        help = "Use this option if you want to run FastQC software  [default %default]",
        dest = "fastqc"
    ),
    make_option(
        opt_str = c("-r", "--multiqc"),
        action = 'store_true',
        type = "logical",
        default = FALSE,
        help = "Use this option if you want to run multiqc software  [default %default]",
        dest = "multiqc"
    ),
    make_option(
        opt_str = c("-o", "--output"), 
        type = "character", 
        default = "01-CleanedReads",
        help = "Output folder [default %default]",
        dest = "output"
    ),
    make_option(
        opt_str = c("-p", "--processors"), 
        type = "integer", 
        default = 8,
        help = "Number of processors to use [default %default]",
        dest = "procs"
    ),
    make_option(
        opt_str = c("-q", "--sampleprocs"),
        type = "integer",
        default = 2,
        help = "Number of samples to process at each time [default %default]",
        dest = "sampleToprocs"
    ),
    make_option(
        opt_str = c("-a", "--adapters"),
        type  = 'character',
        default = 'Data/TruSeq3-PE-2.fa',
        help = "fastaWithAdaptersEtc: specifies the whole path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used. Options: TruSeq3-PE-2.fa, TruSeq3-SE.fa, NexteraPE-PE.fa, TruSeq3-PE.fa, TruSeq2-PE.fa, TruSeq2-SE.fa [default %default]",
        dest = "adapters"
    ),
    make_option(
        opt_str = c("-M", "--seedMismatches"),
        type  = 'integer',
        default = 2,
        help = "SeedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed. [default %default]",
        dest = "seedMis"
    ),
    make_option(
        opt_str = c("-P", "--palindromeClipThreshold"),
        type  = 'integer',
        default = 30,
        help = "PalindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment. [default %default]",
        dest = "palinClipThresh"
    ),
    make_option(
        opt_str = c("-n", "--simpleClipThreshold"),
        type  = 'character',
        default = 10,
        help = "simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read. [default %default]",
        dest = "simClipThresh"
    ),
    make_option(
        opt_str = c('-s', '--sliding'),
        type = 'integer',
        default = 20,
        help = 'Quality sliding to use during trimming. Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. [default %default]',
        dest = 'qual'
    ),
    make_option(
        opt_str = c('-w', '--window'), 
        type = 'integer', 
        default = 10,
        help = 'Quality window to use during trimming. Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. [default %default]',
        dest = 'window'
    ),
    make_option(
        opt_str = c('-L', '--leading'), 
        type = 'integer', 
        default = 3,
        help = 'Remove leading low quality or N bases. Cut bases off the start of a read, if below a threshold quality [default %default]',
        dest = 'leading'
    ),
    make_option(
        opt_str = c('-t', '--trailing'), 
        type = 'integer', 
        default = 3,
        help = 'Remove trailing low quality or N bases. Cut bases off the end of a read, if below a threshold quality [default %default]',
        dest = 'trailing'
    ),
    make_option(
        opt_str = c("-z", "--libraryType"),
        type  = 'character', 
        default = "pairEnd",
        help = "The library type to use. Available: 'pairEnd' or 'singleEnd'. [ default %default]",
        dest = "libraryType"
    ),
    make_option(
        opt_str = c("-m", "--miniumumLength"),
        type = "integer",
        default = 50,
        help = "Discard reads less then minimum length [default %default]",
        dest = "minL"
    ),
    make_option(
        opt_str = c("-C", "--showCommand"),
        action = 'store_true',
        type = "logical",
        default = FALSE,
        help = "Show command line - For debug only [default %default]",
        dest = "showCommand"
    )
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, description =  paste('Authors: OLIVEIRA, H.C. & CANTAO, M.E.', 'Version: 1.0', 'E-mail: hanielcedraz@gmail.com', sep = "\n", collapse = '\n')))



cat("\n\n");cli_rule(center = col_blue("Starting {opt$libraryType} Quality Control"));cat("\n\n")


functionsToSource <- list.files("src", full.names = TRUE)

cat("\n");cli_h1("Loading functions")
for (i in functionsToSource) {
    
    source(i)
    cli_alert_success("Loading script: {basename(i)}")
}



########################################
### PERFORMING QC ANALYSIS
########################################


multiqc <- system('which multiqc 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)
if (opt$multiqc) {
    if (multiqc != 0) {
        cat("\n\n");cli_abort(c("x" = "Multiqc is not installed. Remove -r parameter or install it and run again"))
        
    }
}


# verify if sample_file exist
if (!file.exists(opt$samplesFile)) {
    cat("\n\n");cli_abort(c("x" = "Sample file {opt$samplesFile} does not exist"))
    

}



# samples <- baqcomPackage::loadSamplesFile(file = "samplesPair.txt", reads_folder = "00-FastqPair//", column = opt$samplesColumn, libraryType = "pairEnd")
# 
# baqcomPackage::createSampleList(samples = samples, reads_folder = "00-FastqPair/", column = opt$samplesColumn, fileType = "fastq.gz", libraryType = "pairEnd", step = "QualityControl")


samples <- loadSamplesFile(file = opt$samplesFile, column = opt$samplesColumn, libraryType = opt$libraryType)
cat("\n\n");cli_h1("Samples")
print(as_tibble(samples))

cat("\n\n");cli_h1("Number of processors to use")
procs <- prepareCore(nThreads = opt$procs)


cat("\n\n");cli_h1("Preparing jobs")
qcquery <- createQueryFastqc(samples = samples, readsFolder = opt$Raw_Folder, column = opt$samplesColumn, libraryType = opt$libraryType)

#print(qcquery)


# create report folders
reports <- '02-Reports'
if (!file.exists(file.path(reports))) dir.create(file.path(reports), recursive = TRUE, showWarnings = FALSE)

report_folder <- 'reportBaqcomQC'
if (!file.exists(file.path(paste0(reports,'/',report_folder)))) dir.create(file.path(paste0(reports,'/',report_folder)), recursive = TRUE, showWarnings = FALSE)



#Creating Fastqc plots before Quality Control
# if (opt$fastqc) {
#     beforeQC <- 'FastQCBefore'
#     cat("\n\n");cli_rule(center = col_blue("Start fastqc - Raw fastq"));cat("\n\n")
#     fastq.defore <- runFastqc(
#         inputFolder = opt$Raw_Folder, 
#         reportsPath = glue("{reports}/{beforeQC}"), 
#         procs = opt$procs, 
#         nSamples = opt$sampleToprocs
#     )
# }






## create output folder
output_Folder <- opt$output
if (!file.exists(file.path(output_Folder))) dir.create(file.path(output_Folder), recursive = TRUE, showWarnings = FALSE)



# Trimmomatic analysis function
# pigz <- system('which pigz 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE)
## Running trimmomatic. The command line is wraped in the funcion runTrimmomatic saved in the src folder.
if (opt$libraryType == "pairEnd") {
    cat("\n\n");cli_rule(center = col_blue("START Trimmomatic PE"));cat("\n\n")
    reporstsFolder <- glue("{reports}/{report_folder}/")
    trimmomatic <- "../BAQCOM/Trimmomatic-0.39/trimmomatic-0.39.jar"
    
    
    trimmomatic.pair <- runTrimmomatic(
        query = qcquery, 
        trimmomatic = trimmomatic,
        libType = "PE", 
        cores = opt$procs, 
        nSamples = opt$sampleToprocs, 
        inputFolder = opt$Raw_Folder, 
        outputFolder = opt$output, 
        trimAdapter = opt$adapters, 
        leading = opt$leading,
        trailing = opt$trailing, 
        window = opt$window, 
        qual = opt$qual, 
        minLength = opt$minL, 
        reportsPath = reporstsFolder,
        showCommand = opt$showCommand
    )
    
    if (!all(sapply(trimmomatic.pair , "==", 0L))) {
        
        
        out <- list.files(reporstsFolder, pattern = ".log", full.names = TRUE)
        
        errorLog <- lapply(out, function(x) {
            errorLog <- grep("^Exception in thread|^sh:", readLines(x), value = TRUE)
            
            
            if (str_detect(errorLog, "Exception in thread")) {
                errorLog <- str_remove(errorLog, "Exception in thread \"main\" java.io.FileNotFoundException: ")
            } 
            
            return(errorLog)
            
        }) %>% unique()
        
        
        
        quit_cli_abort(errorLog)
        
    } else {
        
        out <- list.files(reporstsFolder, pattern = ".log", full.names = TRUE)
        
        cat("\n\n");cli_alert_success("TrimmomaticPE Completed successfully");cat("\n")
        
        successLog <- lapply(out, function(x) {
            successLog <- grep("^Input Read Pairs:", readLines(x), value = TRUE)
        
            cli_alert_success("{successLog}")   
            
        }) %>% unique()
        
        
        
    }
    
    
    
} else if (opt$libraryType == "singleEnd") {
    cat("\n\n");cli_rule(center = col_blue("START Trimmomatic SE"));cat("\n\n")
    reporstsFolder <- glue("{reports}/{report_folder}/")
    trimmomatic <- "../BAQCOM/Trimmomatic-0.39/trimmomatic-0.39.jar"
    
    trimmomatic.single <- runTrimmomatic(
        query = qcquery, 
        trimmomatic = trimmomatic,
        libType = "SE", 
        cores = opt$procs, 
        nSamples = opt$sampleToprocs, 
        inputFolder = opt$Raw_Folder, 
        outputFolder = opt$output, 
        trimAdapter = opt$adapters, 
        leading = opt$leading,
        trailing = opt$trailing, 
        window = opt$window, 
        qual = opt$qual, 
        minLength = opt$minL, 
        reportsPath = reporstsFolder
    )
    
    
    if (!all(sapply(trimmomatic.pair , "==", 0L))) {
        
        
        out <- list.files(reporstsFolder, pattern = ".log", full.names = TRUE)
        
        
        
        errorLog <- lapply(out, function(x) {
            errorLog <- grep("^Exception in thread", readLines(x), value = TRUE)
            
            str_remove(errorLog, "Exception in thread \"main\" java.io.FileNotFoundException: ")
            
        }) %>% unique()
        
        
        cli_abort(c("x" = str_remove(errorLog, "Exception in thread \"main\" java.io.FileNotFoundException: ")))
        
    } else {
        cat("\n\n");cli_alert_success("Trimmomatic SE Completed successfully");cat("\n\n")
    }
    
}





# if (opt$fastqc) {
#     
#     afterQC <- 'FastQCAfter'
#     cat("\n\n");cli_rule(center = col_blue("Start fastqc - Cleaned fastq"));cat("\n\n")
#     fastq.after <- runFastqc(
#         inputFolder = opt$output, 
#         reportsPath = glue("{reports}/{afterQC}"), 
#         procs = opt$procs, 
#         nSamples = opt$sampleToprocs
#     )
#     
# }


# if (opt$multiqc) {
#     cat("\n\n");cli_rule(center = col_blue("Running Multiqc"));cat("\n\n")
#     system(glue("multiqc {reports} -o {reports}"))
# }


# Creating samples report
cat("\n\n");cli_rule(center = col_blue("Creating samples report"));cat("\n\n")
getSummaryTable <- function(x) {
    if (opt$libraryType == "pairEnd") {
        final <- data.frame(
            "SampleName" = x[1,3],
            'InputReadPairs' = x[1,2],
            'SurvivedReadsPairs' = x[2,2],
            'SurvivingReadsPairPercent' = x[3,2],
            'ForwardOnlySurvivingReads' = x[4,2],
            'ForwardOnlySurvivingReadPercent' = x[5,2],
            'ReverseOnlySurvivingReads' = x[6,2],
            'ReverseOnlySurvivingReadPercent' = x[7,2],
            'DroppedReads' = x[8,2],
            'DroppedReadPercent' = x[9,2])
    } else if (opt$libraryType == "singleEnd") {
        final <- data.frame(
            'InputReadPairs' = x[1,2],
            'SurvivedReadsPairs' = x[2,2],
            'SurvivingReadPercent' = x[3,2],
            'DroppedReads' = x[4,2],
            'DroppedReadPercent' = x[5,2]
        )
    }
    return(final)
}




report_sample <- list()
trimReportFolder <- glue("{reports}/{report_folder}")
trimStats <- list.files(trimReportFolder, pattern = "statsSummaryFile", full.names = TRUE)
x <- "02-Reports/reportBaqcomQC/HE21-100K_statsSummaryFile.txt"

finalSumDF <- lapply(trimStats, function(x) {
    read.table(x, header = F, as.is = T, fill = TRUE, sep = ':', text = TRUE) %>% 
        mutate(SampleName = str_remove_all(basename(x), "_statsSummaryFile.txt")) %>% 
        getSummaryTable()   
}) %>% 
    bind_rows()


finalSumDF %>% 
    select(SampleName, 
           InputReadPairs, 
           SurvivedReadsPairs,
           DroppedReads,
           DroppedReadPercent) %>% 
    as_tibble()

cat("\n");cli_alert_info("See the file {reports}/QualityControlReportSummary.txt for the whole summary")


write.table(finalSumDF, file = paste0(reports, '/', 'QualityControlReportSummary.txt'), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
# 
# system2('cat', paste0(reports,'/','QualityControlReportSummary.txt'))
# readr::read_table(paste0(reports,'/','QualityControlReportSummary.txt'))
cat("\n\n")

write(glue("How to cite: Please, visit https://github.com/hanielcedraz/BAQCOM/blob/master/how_to_cite.txt or see the file how_to_cite.txt"), stderr())
