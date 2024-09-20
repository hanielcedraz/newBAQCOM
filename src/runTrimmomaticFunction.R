

runTrimmomatic <- function(query, libType = "PE", cores, nSamples, outputFolder, trimAdapter, leading, trailing, window, qual, minLength, reportsPath, showCommand = FALSE) {

    
    trins <- mclapply(query, function(index){
        if (libType == "PE") {

            if (showCommand) {
                
                write(glue("trimmomatic {libType} -threads {cores} {index$R1} {index$R2} {outputFolder}/{index$sampleName}_trimmed_R1.fastq.gz {outputFolder}/{index$sampleName}_trimmed_SE1.fastq.gz {outputFolder}/{index$sampleName}_trimmed_R2.fastq.gz {outputFolder}/{index$sampleName}_trimmed_SE2.fastq.gz -summary {reportsPath}/{index$sampleName}_statsSummaryFile.txt ILLUMINACLIP:{trimAdapter}:2:30:10 LEADING:{leading} TRAILING:{trailing} SLIDINGWINDOW:{window}:{qual} MINLEN:{minLength} 2> {reportsPath}/{index$sampleName}_trimmomatic_out.log \n\n"), stdout())
            }
            
            
            cli_process_start("Running sample {index$sampleName}")
            tryingTrim <- try({
                system(
                    glue("trimmomatic {libType} -threads {cores}",
                         "{index$R1}",
                         "{index$R2}",
                         "{outputFolder}/{index$sampleName}_trimmed_R1.fastq.gz",
                         "{outputFolder}/{index$sampleName}_trimmed_SE1.fastq.gz",
                         "{outputFolder}/{index$sampleName}_trimmed_R2.fastq.gz",
                         "{outputFolder}/{index$sampleName}_trimmed_SE2.fastq.gz",
                         "-summary {reportsPath}/{index$sampleName}_statsSummaryFile.txt",
                         "ILLUMINACLIP:{trimAdapter}:2:30:10",
                         "LEADING:{leading}",
                         "TRAILING:{trailing}",
                         "SLIDINGWINDOW:{window}:{qual}",
                         "MINLEN:{minLength}",
                         "2> {reportsPath}/{index$sampleName}_trimmomatic_out.log",
                         .sep = " "
                    )
                )
                
            })
            
            
            if (tryingTrim == 0) {
                cli_process_done()
            } else {
                cli_process_failed()
            }
            
            
        } else {
            cli_process_start("Running sample {index$sampleName}")
            
            tryingTrim <- try({
                system(
                    glue("trimmomatic {libType} -threads {cores}",
                         "{index$SE}",
                         "{outputFolder}/{index$sampleName}_trimmed_SE.fastq.gz",
                         "-summary {reportsPath}/{index$sampleName}_statsSummaryFile.txt",
                         "ILLUMINACLIP:{trimAdapter}:2:30:10",
                         "LEADING:{leading}",
                         "TRAILING:{trailing}",
                         "SLIDINGWINDOW:{window}:{qual}",
                         "MINLEN:{minLength}",
                         "2> {reportsPath}/{index$sampleName}_trimmomatic_out.log",
                         .sep = " "
                    )
                )
                
            })
            
            if (class(tryingTrim) == "try-error") {
                cli_process_done()
            } else {
                cli_process_failed()
            }
            
        }
        
        
        
        
        
        
        return(tryingTrim)
    }, mc.cores = nSamples)
    
    
    
    return(trins)
}

