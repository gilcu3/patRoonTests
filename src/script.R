script <- "process.R"
source(script)

exps <- c(
		"small"
		)
exptypes <- c(
		"pos"
		#"neg"
		)
		
dbids <- c(
		"PubChem_EXPOSOMICS"
		)

path <- "../"
setwd(path)

for(exp in exps){
	for(exptype in exptypes){
		for(dbid in dbids){
			
			
			
				
			timestamp = as.integer(Sys.time())
			output <- paste("results/report", timestamp, sep="_")
			dir.create(output)
			zstdout <- file(paste(output,"stdout.log", sep = "/"), open="wt")
			zstderr <- file(paste(output,"stderr.log", sep = "/"), open="wt")
			sink(zstderr, type="message")
			sink(zstdout, type="output")
			
			try(process_f(exp, exptype, dbid, output))
			
			sink(type="message")
			sink(type="output")
			
			
		}
	}
}







