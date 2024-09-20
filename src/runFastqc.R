




runFastqc <- function(query, cores, nSamples, outputFolder, libType) {
    
    
    if (!dir.exists(outputFolder)) {
        dir.create(outputFolder, recursive = TRUE)
    }
    
    if (!dir.exists("Logs")) {
        dir.create("Logs", recursive = TRUE)
    }
    
    
    
    trins <- mclapply(query, function(index){
        if (libType == "PE") {
            
            
            cli_process_start("Running sample {index$sampleName}")
            tryingTrim <- try({
                system(
                    glue("fastqc -t {cores} ",
                         "{index$R1}",
                         "{index$R2}",
                         "-o {outputFolder}",
                         "1> Logs/{index$sampleName}.log",
                         "2> Logs/{index$sampleName}.err",
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
                    glue("fastqc -t {cores} ",
                         "{index$SE}",
                         "-o {outputFolder}",
                         "1> Logs/{index$sampleName}.log",
                         "2> Logs/{index$sampleName}.err",
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

