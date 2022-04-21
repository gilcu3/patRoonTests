process_f <- function(data, datatype, dbid, output){
	
	
	
	myLocalDatabasePath <- "exp/PubChemLite_31Oct2020_exposomics.csv"
	
	analysis_input <- "exp/analysis.csv"
	
	cache_path = paste(output, "cache.sqlite", sep="/")

	#Set all options
	options(patRoon.path.MetFragCL = "opt/metfrag/MetFragCommandLine-2.4.7.jar")
	options(patRoon.MP.logPath = paste(output, "log", sep = "/"))
	options(patRoon.cache.fileName = cache_path)
	# options(patRoon.MP.maxProcs = 1)
	
	library(patRoon, warn.conflicts=FALSE)
	
	print(paste("Executing experiment: ", data, datatype, myLocalDatabasePath, output, "Number of cores:", getOption("patRoon.MP.maxProcs")))
	
	
	
	
	
	anaInfo <- read.csv(analysis_input)
	
	
	# This line was executed on Windows to obtain the mzml files
	# convertMSFiles("raw", "mzml", dirs = TRUE, from = "thermo", to = "mzML", algorithm="pwiz", centroid = "vendor", filters = "msLevel 1")
	
	

	library(xcms, warn.conflicts=FALSE)

	fList <- findFeatures(anaInfo, "xcms3", param = xcms::CentWaveParam(ppm = 10, 
						  peakwidth = c(10, 30), mzdiff = -0.001, mzCenterFun =
						  "wMean", integrate = 1L, snthresh = 10, prefilter = 
						  c(3, 100), noise = 500, fitgauss = FALSE, firstBaselineCheck = TRUE), 
						  verbose = TRUE)
	
	print(paste("fList", length(fList)))
	
	# Selects centerSample as the first mzml that is a sample
	for(i in 1:nrow(anaInfo["analysis"])){
		if(startsWith(anaInfo[i, "group"], "sample")){
			centerSample <- i
			break
		}
	}
	print(paste("Using centerSample =", centerSample))
	
	# BiocParallel::register(BiocParallel::SerialParam(), default = TRUE) #IMPORTANT TO AVOID BIOCPARALLEL ERRORS
	fGroups <- groupFeatures(fList, "xcms3", rtalign = TRUE,
							 groupParam = xcms::PeakDensityParam(sampleGroups = 
							 analysisInfo(fList)$group, bw = 10, minFraction = 0.5, 
							 minSamples = 1, binSize = 0.01, maxFeatures = 50), 
							 retAlignParam = xcms::ObiwarpParam(binSize = 1,
							 centerSample = centerSample, response = 100, distFun = "cor_opt",
							 gapInit = 0.3, gapExtend = 2.4, factorDiag = 2, 
							 factorGap = 1, localAlignment = FALSE, 
							 initPenalty = 0))
	
	print(paste("fGroups", length(fGroups)))
	
	fGroups <- filter(fGroups, absMinIntensity = 100,
					  relMinReplicateAbundance = NULL, maxReplicateIntRSD = NULL,
					  blankThreshold = 3, removeBlanks = TRUE,
					  retentionRange = NULL, mzRange = NULL)
	
	print(paste("fGroups", length(fGroups)))
	
	if(datatype == "pos"){
		adduct <- "[M+H]+"
		elements <- "CHNOClBrSP"
	}
	else{
		adduct <- "[M-H]-"
		elements <- "CHNOCl"
	}
	
	
	# -------------------------
	# annotation
	# -------------------------


	# Retrieve MS peak lists
	avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
	mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 15, 
								   precursorMzWindow = 0.5, avgFeatParams = 
								   avgPListParams, avgFGroupParams =
								   avgPListParams)
	print(paste("mslists", length(mslists)))

	mslists <- filter(mslists, withMSMS = TRUE, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
						relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
						deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)
	print(paste("mslists", length(mslists)))
	
	

	# Calculate formulae candidates
	formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 5,
								 adduct = adduct, elements = elements,
								 calculateFeatures = FALSE, featThreshold = 0.75,
								 # timeout = 120
								 )
	print(paste("formulas", length(formulas)))
	# Find compound structure candidates
	compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
								   dbRelMzDev = 5, 
								   fragRelMzDev = 5, 
								   fragAbsMzDev = 0.002, 
								   adduct = adduct,
								   database = "csv" , extraOpts = 
									list(LocalDatabasePath = myLocalDatabasePath), scoreTypes = c("fragScore","metFusionScore","score", "individualMoNAScore"),
								   # maxCandidatesToStop = 500, errorRetries = 200, timeoutRetries = 20
								   )
	
	print(paste("compounds", length(compounds)))
	
	#formulascoring                              
	
	compounds <- addFormulaScoring(compounds, formulas, TRUE) 
	
	# -------------------------
	# reporting
	# -------------------------
	
	
	# Fixes bug "NULL" -> NA in the compounds table, to avoid error later in reportPDF
#	cTable <- compoundTable(compounds)
#	for(i in 1:length(cTable)){
#		if("XlogP" %in% colnames(cTable[[i]])){
#			# when this code is executed some warnings will be shown
#			# like "In process_f(exp, exptype, dbid) : NAs introduced by coercion"
#			cTable[[i]]$XlogP <- as.double(cTable[[i]]$XlogP)
#		}
#	}
#	compounds@compounds <- cTable
	
	
	
	# Summary of MetFrag Results in a a Single Table 
	MFsummary <- as.data.table(compounds)
	outputSummary <- paste(output, "MFsummary.csv", sep = "/")
	
	write.csv(MFsummary, outputSummary)
	

	reportCSV(fGroups, path = output, reportFeatures = FALSE, formulas = formulas,
			  compounds = compounds, compoundsNormalizeScores = "max")
			  

	reportPDF(fGroups, path = output, reportFGroups = TRUE, formulas = formulas,
			  reportFormulaSpectra = TRUE, compounds = compounds, 
			  compoundsNormalizeScores = "max",
			  MSPeakLists = mslists)
	
	unlink(cache_path)
}


