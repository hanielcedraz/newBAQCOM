



createQueryFastqc <- function(libraryType = "pairEnd", samples = "samples.txt", readsFolder = "00-Fastq", column = "SAMPLE_ID") {
    
    fastqList <- list()
    reads <- dir(path = file.path(readsFolder), pattern = "fastq.gz$", full.names = TRUE)
    
    if (libraryType == "pairEnd") {
        
        if (class(samples) == "character") {
            samples <- read.table(samples, header = TRUE)
        }
        
        
        for (i in 1:nrow(samples)) {
            
            if (all(str_detect(samples$Read_1, "PE1"))) {
                map <- lapply(c("_PE1", "_PE2"), grep, x = reads, value = TRUE)
                names(map) <- c("PE1", "PE2")
                map$sampleName <-  samples[i,column]
                map$PE1 <- map$PE1[i]
                map$PE2 <- map$PE2[i]
                
            } else if (all(str_detect(samples$Read_1, "R1"))) {
                map <- lapply(c("_R1", "_R2"), grep, x = reads, value = TRUE)
                names(map) <- c("R1", "R2")
                map$sampleName <-  samples[i,column]
                map$R1 <- map$R1[i]
                map$R2 <- map$R2[i]
                
            } else {
                #cli_alert_danger("Filenames out of patern: {samples$Read_1}")

                cli_abort("Filenames out of patern: {samples$Read_1}")
            }
            
            fastqList[[paste(map$sampleName)]] <- map
            fastqList[[paste(map$sampleName, sep = "_")]]
        }

        
        
    } else if (libraryType == "singleEnd") {
        
        for (i in 1:nrow(samples)) {
            map <- lapply(c("_SE"), grep, x = reads, value = TRUE)
            names(map) <- c("SE")
            map$sampleName <-  samples[i,column]
            map$SE <- map$SE[i]
            fastqList[[paste(map$sampleName)]] <- map
        }
        
        
    }
    
    cat("\n\n");cli_alert_success("Setting up {length(fastqList)} jobs");cat("\n\n")
    return(fastqList)
}




runFastqc <- function(inputFolder, reportsPath, procs, nSamples) {
    
    if (!file.exists(reportsPath)) {
        dir.create(reportsPath, recursive = TRUE, showWarnings = FALSE)
    }
    query <- list.files(inputFolder)
    
    fastq.defore <- mclapply(query, function(index){
        cli_process_start("Running sample {index}")
        fastqcRun <- try({
            system(
                glue("fastqc",
                      "{inputFolder}/{index}",
                      "-o {reportsPath}",
                      "-t {procs}",
                     "-q",
                     .sep = " "
                )
            )
        })
        cli_process_done()
        return(fastqcRun)
    }, mc.cores = nSamples)
    
    return(fastq.defore)
}
