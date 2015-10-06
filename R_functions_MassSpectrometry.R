###INSTALL THE REQUIRED PACKAGES
requiredPackages <- c ("parallel", "MALDIquant", "MALDIquantForeign", "caret", "tcltk", "beepr")
installedPackages <- installed.packages () [,1]
missingPackages <- character()
for (p in 1:length(requiredPackages)) {
	if ((requiredPackages[p] %in% installedPackages) == FALSE) {
		missingPackages <- append (missingPackages, requiredPackages[p])
	}
}
if (length(missingPackages) > 0) {
	install.packages (missingPackages)
}

### LOAD THE REQUIRED PACKAGES
library ("parallel")
library ("MALDIquant")
library ("MALDIquantForeign")
library ("tcltk")
library ("caret")
#library ("e1071")
library ("beepr")
#library ("gWidgets")
#library ("gWidgetstcltk")


##### MEMORY LIMIT
# Set the memory (RAM) limit to 100GB (100000MB)
memory.limit (100000)



######### PARALLELISATION
# Detect the number of CPU threads
#CPUthreadNumber <- detectCores (logical=TRUE)
# Turn it into the number of actual cores
#CPUcoreNumber <- CPUthreadNumber/2

CPUcoreNumber <- 2


################################################ INSTALL REQUIRED PACKAGES
installAndLoadRequiredPackages <- function (requiredPackages) {
installedPackages <- installed.packages () [,1]
missingPackages <- character()
for (p in 1:length(requiredPackages)) {
	if ((requiredPackages[p] %in% installedPackages) == FALSE) {
		missingPackages <- append (missingPackages, requiredPackages[p])
	}
}
if (length(missingPackages) > 0) {
	install.packages (missingPackages)
}
for (p in requiredPackages) {
	library (p)
}
}





################################################ SPECTRA QUALITY CONTROL
spectraQualityControl <- function (spectra) {
## Empty spectra
print (paste ("Are there any empty spectra?", (any (sapply (spectra, isEmpty)))))
## Same number of data points
dataPointsTable <- table (sapply (spectra, length))
print ("Same number of data points?")
print (dataPointsTable)
## Same distance between data points?
print (paste("Same distance between data points?", (all (sapply (spectra, isRegular)))))
## Flat spectra
flatSpectra <- findEmptyMassObjects (spectra)
print (paste ("Flat spectra in the dataset:", length (flatSpectra)))
}





######################################### REPLACE THE metaData$file WITH THE SAMPLE NAME (SPECTRA or PEAKS)
replaceSampleName <- function (spectra, folder, fileFormat="imzml") {
setwd (folder)
# Read the files in the folder
folderFiles <- readSpectraFiles (folder, fileFormat=fileFormat, fullPath=FALSE)
if (fileFormat == "imzml" | fileFormat == "imzML") {
	for (s in 1:length(folderFiles)) {
		for (i in 1:length(spectra)) {
			if (length(grep(folderFiles[s],spectra[[i]]@metaData$file, fixed=TRUE)) != 0) {
				spectra[[i]]@metaData$file <- folderFiles[s]
				spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(".imzML")))[1]
			}
		}
	}
}
#if (fileFormat == "imzml" | fileFormat == "imzML") {
#	for (i in 1:length(spectra)) {
#		spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(folder)))[2]
#		spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(".imzML")))[1]
 #   }
#}
return (spectra)
}





#########################################################################





######################################### REPLACE THE metaData$file WITH THE CLASS NAME (SPECTRA or PEAKS)
replaceClassName <- function (spectra, classList, fileFormat="imzml") {
classList <- sort (classList)
if (fileFormat == "imzml" | fileFormat == "imzML") {
	for (w in 1:length(classList)) {
		for (i in 1:length(spectra)) {
			if (length(grep(classList[w],spectra[[i]]@metaData$file)) == 1) {
				spectra[[i]]@metaData$file <- classList [w]
			}
		}
	}
}
if (fileFormat == "brukerflex") {
	for (w in 1:length(classList)) {
		for (i in 1:length(spectra)) {
			if (length(grep(classList[w],spectra[[i]]@metaData$sampleName)) == 1) {
				spectra[[i]]@metaData$sampleName <- classList [w]
			}
		}
	}
}
return (spectra)
}





########################################################################





######################################### SPECTRA PRE-PROCESSING
# Preprocess spectra
preprocessSpectra <- function (spectra, tofMode="linear", smoothingStrength="medium", processInPackagesOf=length(spectra)) {
if (processInPackagesOf <= 0 || preprocessInPackagesOf > length(spectra)) {
	processInPackagesOf <- length(spectra)
}
# Create the list containing the processed spectra
spectraProcessed <- list()
# Calculate the kFold
kFold <- floor (length(spectra) / processInPackagesOf)
# Determine the spectra to be randomly allocated into processing folds
index <- createFolds (y=seq(1:length(spectra)), k=kFold)
for (i in 1:length(index)) {
	# Make the cluster (one for each core/thread)
	cl <- makeCluster (CPUcoreNumber)
	clusterEvalQ (cl, {library(MALDIquant)})
	# Allocate the spectra (in the package) to be processed temporarily
	spectraTemp <- spectra [index[[i]]]
	if (tofMode == "linear") {
		## Remove flat spectra
		#spectra <- removeEmptyMassObjects (spectra)
		if (smoothingStrength == "small") {
			## Smoothing (Savitzky-Golay filter, with window size 5, 11 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=5))
		}
		if (smoothingStrength == "medium") {
			## Smoothing (Savitzky-Golay filter, with window size 10, 21 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=10))
		}
		if (smoothingStrength == "strong") {
			## Smoothing (Savitzky-Golay filter, with window size 20, 41 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=20))
		}
		if (smoothingStrength == "veryStrong") {
			## Smoothing (Savitzky-Golay filter, with window size 30, 61 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=30))
		}
		## Baseline correction
		spectraTemp <- parLapply (cl, spectraTemp, fun= function (spectra) removeBaseline (spectra, method="TopHat"))
		## Normalisation (TIC)
		spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) calibrateIntensity (spectra, method="TIC"))
	}
	if (tofMode == "reflectron" || tofMode == "reflector") {
		## Remove flat spectra
		#spectra <- removeEmptyMassObjects (spectra)
		if (smoothingStrength == "small") {
			## Smoothing (Savitzky-Golay filter, with window size 2, 5 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=1))
		}
		if (smoothingStrength == "medium") {
			## Smoothing (Savitzky-Golay filter, with window size 6, 13 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=3))
		}
		if (smoothingStrength == "strong") {
			## Smoothing (Savitzky-Golay filter, with window size 12, 25 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=6))
		}
		if (smoothingStrength == "veryStrong") {
			## Smoothing (Savitzky-Golay filter, with window size 18, 37 points)
			spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=9))
		}
		## Baseline correction
		spectraTemp <- parLapply (cl, spectraTemp, fun= function (spectra) removeBaseline (spectra, method="TopHat"))
		## Normalisation (TIC)
		spectraTemp <- parLapply (cl, spectraTemp, fun=function (spectra) calibrateIntensity (spectra, method="TIC"))
	}
	# Close the processes
	stopCluster(cl)
	# Append the processed spectra package to the global list
	spectraProcessed <- append (spectraProcessed, spectraTemp)
}
return (spectraProcessed)
}





########################################################################





######################################### AVERAGE SIGNAL NUMBER
# Average number of signals with that S/N in the peaklist set
averageSignalNumber <- function (spectra, SNR=5) {
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
avgSignalNumber <- vector (length=0)
signalNumberVector <- vector (length=0)
# For each peaklist
for (p in 1:length(peaks)) {
	# Determine the number of peaks
	signalNumber <- length(peaks[[p]]@mass)
	signalNumberVector <- append (signalNumberVector, signalNumber)
}
avgSignalNumber <- sum (signalNumberVector) / length (signalNumberVector)
return (avgSignalNumber)
}





#########################################################################





######################################### REPRESENTATIVE SPECTRA (one imzML per patient)
# Number of spectra with more than a defined number of signals with that SNR
representativeSpectraNumber <- function (spectra, SNR=15, signalThreshold=20) {
counter <- 0
## Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
# For each peaklist
for (p in 1:length(peaks)) {
	# Increase the counter if the number of peaks is more than the defined threshold
	if (length(peaks[[p]]@mass) >= signalThreshold) {
		counter <- counter + 1
	}
}
return (counter)
}





########################################################################





######################################### REPRESENTATIVE SPECTRA DATASET (one imzML per patient)
# Average number of spectra with more than a defined number of signals with that SNR
averageRepresentativeSpectraNumber <- function (filepath, SNR=15, signalThreshold=20, tofMode="linear",
	fileFormat="imzml") {
counterVector <- vector (length=0)
for (f in 1:length(filepath)) {
	# Import Spectra
	if (fileFormat == "imzml" | fileFormat == "imzML") {
		spectra <- importImzMl (filepath[f])
	}
	if (fileFormat == "brukerflex") {
		spectra <- importBrukerFlex (filepath[f])
	}
	# Pre-processing
	spectra <- preprocessSpectra (spectra, tofMode=tofMode)
	# Counter
	counter <- 0
	## Peak picking
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
	# For each peaklist
	for (p in 1:length(peaks)) {
		# Increase the counter if the number of peaks is more than the defined threshold
		if (length(peaks[[p]]@mass) >= signalThreshold) {
			counter <- counter + 1
		}
	}
	counterVector <- append (counterVector, counter)
}
averageCounter <- mean (counterVector)
return (counterVector)
}





########################################################################





########################################## SPECTRA FILTERING
# Remove all the bad spectra, keep only the best ones (Number of signals with a certain SNR)
filterSpectraSignalNumber <- function (spectra, SNRfilter=15, signalThreshold=15) {
	# Make the cluster (one for each core/thread)
	cl <- makeCluster (CPUcoreNumber)
	clusterEvalQ (cl, {library(MALDIquant)})
##################################################### FILTERING FUNCTION
spectraFilteringFunction <- function (spectra) {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNRfilter)
	####### Select only the desired peaklists, based upon the threshold characteristics
		if (length(peaks@mass) < signalThreshold) {
			# Remove the bad spectrum from the list
			spectra <- NULL
		}
	return (spectra)
}
########################################################################
	spectraPurified <- parLapply (cl, spectra, fun=function(spectra) spectraFilteringFunction(spectra))
	# Close the processes
	stopCluster(cl)
	### Keep only the elements that are different from NULL
	spectraPurified <- spectraPurified [!sapply(spectraPurified, is.null)]
	return (spectraPurified)
}






#########################################################################





########################################## SPECTRA GROUPING (PATIENTS) (one imzML file per patient)
# Obtain a certain number of spectra per patient -> Representative spectra
groupSpectra <- function (spectra, spectraPerPatient=1, fileFormat="imzml", tofMode="linear",
	seed=0, algorithm="random", discardedNodes=1, balanced=TRUE, method="mean") {
### Create the file Vector
if (fileFormat == "imzml" | fileFormat == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		fileVector <- append (fileVector, spectra[[i]]@metaData$file)
	}
	patientVector <- unique (fileVector)
}
if (fileFormat == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		fileVector <- append (fileVector, spectra[[i]]@metaData$sampleName)
	}
	patientVector <- unique (fileVector)
}
######################################## Balancing
############ The spectraPerPatient get adjusted based on the shortest spectra list
if (balanced == TRUE) {
	########## All the patients should have the same number of spectra, always
	### Determine the shortest patient's spectra number
	#### Set the lowest as the row number of the first patient
	# First patient
	spectraFirstPatient <- list()
	for (s in 1:length(spectra)) {
		if (spectra[[s]]@metaData$file == patientVector[[1]]) {
			spectraFirstPatient <- append (spectraFirstPatient, spectra[[s]])
		}
	}
	lowestObservationNumber <- length (spectraFirstPatient)
	# Scroll the other patients (if there are any)
	if (length(patientVector) > 1) {
		for (p in 2:length(patientVector)) {
			# Isolate the spectra of that patient
			spectraPatient <- list()
			for (s in 1:length(spectra)) {
				if (spectra[[s]]@metaData$file == patientVector[[p]]) {
					spectraPatient <- append (spectraPatient, spectra[[s]])
				}
			}
			# If its lower than the reference
			if (length(spectraPatient) < lowestObservationNumber) {
				lowestObservationNumber <- length(spectraPatient)
			}
		}
	}
	#### Rows per patient should not be more than the minimum number of rows
	if (spectraPerPatient >= lowestObservationNumber) {
		spectraPerPatient <- lowestObservationNumber
	}
}
# If the spectra per patient is equal to zero or one, do the normal averaging
if (spectraPerPatient == 0 || spectraPerPatient == 1) {
	if (method == "mean") {
		patientSpectraFinal <- averageMassSpectra (spectra, labels=fileVector, method="mean")
	}
	if (method == "skyline") {
		patientSpectraFinal <- groupSpectraSkyline (spectra, fileFormat=fileFormat)
	}
}
# Run this script if there has to be two or more representative spectra per patient
if (spectraPerPatient > 1) {
	# If there is only one spectrum, use it
	if (length(spectra) == 1) {
		patientSpectraFinal <- spectra
	}
	# If there are more spectra...
	if (length(spectra) > 1) {
		# Make the randomness reproducible
		if (seed != 0) {
			set.seed (seed)
		}
		############################# OUTPUTS
		# Generate the final list of spectra
		patientSpectraFinal <- list()
		discardedSpectra <- list()
		discardedSpectraAVG <- list()
		# Create a new list of spectra for plotting purposes (the intensities will be replaced)
		spectraForPlotting <- list()
		# List of additional plots
		plots <- list()
		# Generate the final list of MSI images
		msiPlots <- list()
		##############################
		# For each patient (p)
		for (p in 1:length(patientVector)) {
			# Add the single patient spectra to a list
			patientSpectra <- list()
			for (s in 1:length(spectra)) {
				if (spectra[[s]]@metaData$file == patientVector[p]) {
					patientSpectra <- append (patientSpectra, spectra[[s]])
				}
			}
			############################################## RANDOMNESS
			if (algorithm == "random") {
				# Do this if the spectra length is more than two, otherwise it does not make any sense
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
				index <- createFolds (fileVector[fileVector==patientVector[p]], k = spectraPerPatient)
				# For each fold
				for (k in 1:length(index)) {
					# Generate a temporary list, where the patient spectra of one fold will be stored (they will be averaged)
					temporaryPatientSpectra <- list()
					# Create a new list of spectra for plotting purposes (the intensities will be replacd)
					temporaryPatientSpectraForPlotting <- list()
					# Take the corresponding spectra subset (of the selected fold)
					for (i in 1:length(index[[k]])) {
						foldIndex <- index[[k]] [i]
						temporaryPatientSpectra <- append (temporaryPatientSpectra, patientSpectra [[foldIndex]])
						temporaryPatientSpectraForPlotting <- append (temporaryPatientSpectraForPlotting, patientSpectra [[foldIndex]])
					}
					# Replace the intensities with the K number for plotting purposes
					for (n in 1:length(temporaryPatientSpectraForPlotting)) {
						temporaryPatientSpectraForPlotting[[n]]@intensity <- rep (k, length(temporaryPatientSpectraForPlotting [[n]]@intensity))
					}
					# Add these modified spectra to the final list of spectra for plotting purposes
					spectraForPlotting <- append (spectraForPlotting, temporaryPatientSpectraForPlotting)
					# Average them
					if (method == "mean") {
						temporaryPatientAvg <- averageMassSpectra (temporaryPatientSpectra, method="mean")
					}
					if (method == "skyline") {
						temporaryPatientAvg <- generateSkylineSpectrum (temporaryPatientSpectra)
					}
					# Store them in the final spectra list
					patientSpectraFinal <- append (patientSpectraFinal, temporaryPatientAvg)
				}
				# Store the plot into the list
				plotMsiSlice (spectraForPlotting, center=spectraForPlotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msiPlots [[p]] <- recordPlot ()
			}
			######################################### SIMILARITY
			############################ HCA
			if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectraPerPatient <- spectraPerPatient + discardedNodes
				# Detect and align peaks
				if (tofMode=="linear") {
					peaks <- detectPeaks (patientSpectra, method="MAD", SNR=3, halfWindowSize=20)
					peaks <- alignAndFilterPeaks (peaks, tolerancePPM=2000, peaksFiltering=FALSE, frequencyThreshold=0.25)
				}
				if (tofMode=="reflectron" | tofMode=="reflector") {
					peaks <- detectPeaks (patientSpectra, method="MAD", SNR=3, halfWindowSize=5)
					peaks <- alignAndFilterPeaks (peaks, tolerancePPM=200, peaksFiltering=FALSE, frequencyThreshold=0.25)
				}
				# Generate the peaklist matrix
				peaklist <- intensityMatrix (peaks, patientSpectra)
				# Compute the distance matrix
				distanceMatrix <- dist (peaklist, method="euclidean")
				hca <- hclust (distanceMatrix)
				# Store the plot
				plot (hca, xlab="MSI elements", main="Hierarchical clustering analysis", sub="")
				legendText <- paste ("Cluster method:", hca$method, "\nDistance method: ", hca$dist.method, "\n")
				legend ("topright", title="Hierarchical clustering", legend=legendText, xjust=0.5, border="black")
				plots[[p]] <- recordPlot()
				# Associate to each row/spectrum the subgroup of the tree to which it belongs
				if (spectraPerPatient > 56) {
					spectraPerPatient <- 56
				}
				hcaGroups <- cutree (hca, k=spectraPerPatient)
				if (discardedNodes != 0) {
				# For each subgroup to be isolated...
					for (d in 1:discardedNodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which (hcaGroups == d)
						spectraHCA <- patientSpectra[index]
						spectraHCAForPlotting <- patientSpectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectraHCAForPlotting)) {
							spectraHCAForPlotting[[n]]@intensity <- rep (0, length(spectraHCAForPlotting[[n]]@intensity))
						}
						spectraForPlotting <- append (spectraForPlotting, spectraHCAForPlotting)
						# They do not get added to the final list of patient spectra
						# Add them to the discarded spectra
						discardedSpectra <- append (discardedSpectra, spectraHCA)
						# Average the spectra
						if (method == "mean") {
							spectraAVG <- averageMassSpectra (spectraHCA, method="mean")
						}
						if (method == "skyline") {
							spectraAVG <- generateSkylineSpectrum (spectraHCA)
						}
						# Add the average to the final list of discarded spectra AVG
						discardedSpectraAVG <- append (discardedSpectraAVG, spectraAVG)
					}
					# For each subgroup to be isolated...
					for (d in (discardedNodes+1):spectraPerPatient) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which (hcaGroups == d)
						spectraHCA <- patientSpectra[index]
						spectraHCAForPlotting <- patientSpectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectraHCAForPlotting)) {
							spectraHCAForPlotting[[n]]@intensity <- rep (d, length(spectraHCAForPlotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectraForPlotting <- append (spectraForPlotting, spectraHCAForPlotting)
						# Average the spectra in this HCA subgroup
						if (method == "mean") {
							spectraAVG <- averageMassSpectra (spectraHCA, method="mean")
						}
						if (method == "skyline") {
							spectraAVG <- generateSkylineSpectrum (spectraHCA)
						}
						# Add the average to the final list
						patientSpectraFinal <- append (patientSpectraFinal, spectraAVG)
					}
				} else {
				# For each subgroup to be isolated...
				for (s in 1:spectraPerPatient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which (hcaGroups == s)
					spectraHCA <- patientSpectra[index]
					spectraHCAForPlotting <- patientSpectra[index]
					# Replace the intensities with the S number for plotting purposes
					for (n in 1:length(spectraHCAForPlotting)) {
						spectraHCAForPlotting[[n]]@intensity <- rep (s, length(spectraHCAForPlotting[[n]]@intensity))
					}
					# Add these modified spectra to the final list of spectra for plotting purposes
					spectraForPlotting <- append (spectraForPlotting, spectraHCAForPlotting)
					# Average the spectra in this HCA subgroup
					if (method == "mean") {
						spectraAVG <- averageMassSpectra (spectraHCA, method="mean")
					}
					if (method == "skyline") {
						spectraAVG <- generateSkylineSpectrum (spectraHCA)
					}
					# Add the average to the final list
					patientSpectraFinal <- append (patientSpectraFinal, spectraAVG)
				}
				}
				# Store the plot into the list
				plotMsiSlice (spectraForPlotting, center=spectraForPlotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msiPlots [[p]] <- recordPlot ()
			}
			######################### K-MEANS
			if (algorithm == "kMeans" | algorithm == "kmeans" | algorithm == "k-Means") {
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectraPerPatient <- spectraPerPatient + discardedNodes
				# Detect and align peaks
				if (tofMode=="linear") {
					peaks <- detectPeaks (patientSpectra, method="MAD", SNR=3, halfWindowSize=20)
					peaks <- alignAndFilterPeaks (peaks, tolerancePPM=2000, peaksFiltering=FALSE, frequencyThreshold=0.25)
				}
				if (tofMode=="reflectron" | tofMode=="reflector") {
					peaks <- detectPeaks (patientSpectra, method="MAD", SNR=3, halfWindowSize=5)
					peaks <- alignAndFilterPeaks (peaks, tolerancePPM=200, peaksFiltering=FALSE, frequencyThreshold=0.25)
				}
				# Generate the peaklist matrix
				peaklist <- intensityMatrix (peaks, patientSpectra)
				# Compute the k-Means clustering
				kMeans <- kmeans (peaklist, centers=spectraPerPatient)
				# Associate to each row/spectrum the subgroup/cluster to which it belongs
				kMeansGroups <- kMeans$cluster
				if (discardedNodes != 0) {
				# For each subgroup to be isolated (to discard)...
					for (d in 1:discardedNodes) {
						# Index the spectra under in the selected subgroup
						index <- which (kMeansGroups == d)
						spectraKMeans <- patientSpectra[index]
						spectraKMeansForPlotting <- patientSpectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectraKMeansForPlotting)) {
							spectraKMeansForPlotting[[n]]@intensity <- rep (0, length(spectraKMeansForPlotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectraForPlotting <- append (spectraForPlotting, spectraKMeansForPlotting)
						# They do not get averaged and added to the final list of patient spectra
						# Add them to the discarded spectra
						discardedSpectra <- append (discardedSpectra, spectraKMeans)
						# Average the spectra
						if (method == "mean") {
							spectraAVG <- averageMassSpectra (spectraKMeans, method="mean")
						}
						if (method == "skyline") {
							spectraAVG <- generateSkylineSpectrum (spectraKMeans)
						}
						# Add the average to the final list of discarded spectra AVG
						discardedSpectraAVG <- append (discardedSpectraAVG, spectraAVG)
					}
					# For each subgroup to be isolated (to keep)...
					for (d in (discardedNodes+1):spectraPerPatient) {
						# Index the spectra under in the selected subgroup
						index <- which (kMeansGroups == d)
						spectraKMeans <- patientSpectra[index]
						spectraKMeansForPlotting <- patientSpectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectraKMeansForPlotting)) {
							spectraKMeansForPlotting[[n]]@intensity <- rep (d, length(spectraKMeansForPlotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectraForPlotting <- append (spectraForPlotting, spectraKMeansForPlotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectraAVG <- averageMassSpectra (spectraKMeans, method="mean")
						}
						if (method == "skyline") {
							spectraAVG <- generateSkylineSpectrum (spectraKMeans, method="mean")
						}
						# Add the average to the final list
						patientSpectraFinal <- append (patientSpectraFinal, spectraAVG)
					}
				} else {
					# For each subgroup to be isolated...
					for (s in 1:spectraPerPatient) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which (kMeansGroups == s)
						spectraKMeans <- patientSpectra[index]
						spectraKMeansForPlotting <- patientSpectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectraKMeansForPlotting)) {
							spectraKMeansForPlotting[[n]]@intensity <- rep (s, length(spectraKMeansForPlotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectraForPlotting <- append (spectraForPlotting, spectraKMeansForPlotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectraAVG <- averageMassSpectra (spectraKMeans, method="mean")
						}
						if (method == "skyline") {
							spectraAVG <- generateSkylineSpectrum (spectraKMeans, method="mean")
						}
						# Add the average to the final list
						patientSpectraFinal <- append (patientSpectraFinal, spectraAVG)
					}
					}
				# Store the plot into the list
				plotMsiSlice (spectraForPlotting, center=spectraForPlotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msiPlots [[p]] <- recordPlot ()
			}
			#####################
		}
	}
}
# Processing before returning
patientSpectraFinal <- removeBaseline (patientSpectraFinal, method="TopHat")
patientSpectraFinal <- calibrateIntensity (patientSpectraFinal, method="TIC")
#
return (list(spectra=patientSpectraFinal, msImages=msiPlots, discarded=discardedSpectra, discardedAvg=discardedSpectraAVG, plots=plots))
}
################################################################################





########################################## SPECTRA GROUPING (PATIENTS) - SKYLINE
# Obtain one spectrum per patient (average)
groupSpectraSkyline <- function (spectra, fileFormat="imzml") {
if (fileFormat == "imzml" | fileFormat == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		fileVector <- append (fileVector, spectra[[i]]@metaData$file)
	}
	# Put the class names into a vector
	pathSpectra <- unique (fileVector)
	# Create the empty avgList
	spectraGrouped <- list ()
	# For each filepath...
	for (p in 1:length (pathSpectra)) {
		# Create the empty list in which spectra will be allocated
		spectraList <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length (spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$file == pathSpectra [p]) {
				spectraList <- append (spectraList, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same path)
		skylineSpectrum <- generateSkylineSpectrum (spectraList)
		# Add the average spectrum to another final list of average spectra
		spectraGrouped <- append (spectraGrouped, skylinepectrum)
	}
}
if (fileFormat == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		fileVector <- append (fileVector, spectra[[i]]@metaData$sampleName)
	}
	# Put the file names into a vector
	pathSpectra <- unique (fileVector)
	# Create the empty avgList
	spectraGrouped <- list ()
	# For each filepath...
	for (p in 1:length (pathSpectra)) {
		# Create the empty list in which spectra will be allocated
		spectraList <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length (spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$sampleName == pathSpectra [p]) {
				spectraList <- append (spectraList, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same path)
		skylineSpectrum <- generateSkylineSpectrum (spectraList)
		# Add the average spectrum to another final list of average spectra
		spectraGrouped <- append (spectraGrouped, skylineSpectrum)
	}
}
return (spectraGrouped)
}





########################################################################





########################################## SPECTRA GROUPING (CLASSES)
# Obtain one spectrum per class (average)
groupSpectraClass <- function (spectra, classList, fileFormat="imzml") {
classList <- sort (classList)
if (fileFormat == "imzml" | fileFormat == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replaceClassName (spectra, classList=classList, fileFormat="imzml")
	# Put the filenames/classes in a vector
	# Create the empty vector
	classVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		classVector <- append (classVector, spectra[[i]]@metaData$file)
	}
	classSpectraGrouped <- averageMassSpectra (spectra, labels=classVector, method="mean")
}
if (fileFormat == "brukerflex") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replaceClassName (spectra, classList=classList, fileFormat="brukerflex")
	# Put the filenames/classes in a vector
	# Create the empty vector
	classVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		classVector <- append (classVector, spectra[[i]]@metaData$sampleName)
	}
	classSpectraGrouped <- averageMassSpectra (spectra, labels=classVector, method="mean")
}
return (classSpectraGrouped)
}





########################################################################





########################################## SPECTRA GROUPING (CLASSES) - SKYLINE
# Obtain one spectrum per class (average)
groupSpectraClassSkyline <- function (spectra, classList, fileFormat="imzml") {
classList <- sort (classList)
if (fileFormat == "imzml" | fileFormat == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replaceClassName (spectra, classList=classList, fileFormat="imzml")
	# Put the filenames/classes in a vector
	# Create the empty vector
	classVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		classVector <- append (classVector, spectra[[i]]@metaData$file)
	}
	# Put the class names into a vector
	classSpectra <- unique (classVector)
	# Create the empty avgList
	classSpectraGrouped <- list ()
	# For each class...
	for (p in 1:length (classSpectra)) {
		# Create the empty list in which spectra will be allocated
		spectraList <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length (spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$file == classSpectra [p]) {
				spectraList <- append (spectraList, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same class)
		skylineSpectrum <- generateSkylineSpectrum (spectraList)
		# Add the average spectrum to another final list of average spectra
		classSpectraGrouped <- append (classSpectraGrouped, skylineSpectrum)
	}
}
if (fileFormat == "brukerflex") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replaceClassName (spectra, classList=classList, fileFormat="brukerflex")
	# Put the filenames/classes in a vector
	# Create the empty vector
	classVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		classVector <- append (classVector, spectra[[i]]@metaData$sampleName)
	}
	# Put the class names into a vector
	classSpectra <- unique (classVector)
	# Create the empty avgList
	classSpectraGrouped <- list ()
	# For each class...
	for (p in 1:length (classSpectra)) {
		# Create the empty list in which spectra will be allocated
		spectraList <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length (spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$sampleName == classSpectra [p]) {
				spectraList <- append (spectraList, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same class)
		skylineSpectrum <- generateSkylineSpectrum (spectraList)
		# Add the average spectrum to another final list of average spectra
		classSpectraGrouped <- append (classSpectraGrouped, skylineSpectrum)
	}
}
return (classSpectraGrouped)
}





################################################################################





######################################### ADD THE CLASS AND THE SAMPLE NAME TO THE MATRIX
matrixAddClassAndSample <- function (signalMatrix, peaks=list(), classList=list(), fileFormat="imzml", sampleOutput=TRUE, classOutput=TRUE) {
signalMatrix <- as.matrix (signalMatrix)
spectraNumber <- length (peaks)
### The name of the rows will be either the sample name or the class name (depending on the function parameter)
# If the rows are named according to the sample name, an additional column for the class is added
# Otherwise, use the class names as the row names
if ((classOutput==FALSE && sampleOutput==TRUE) || (classOutput==TRUE && length(classList)==0 && sampleOutput==TRUE)) {
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length (if averaged, it contains the name of all the single spectra)
		if (fileFormat == "imzml" | fileFormat == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file[1])
			}
		}
		if (fileFormat == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	# Create the sample matrix column and append it to the global matrix
	sampleColumn <- matrix (0, ncol=1, nrow=spectraNumber)
	colnames (sampleColumn) <- "Sample"
	sampleColumn [,1] <- cbind (fileVector)
	signalMatrix <- cbind (signalMatrix, sampleColumn)
}
if (classOutput==TRUE && length(classList)>=1 && sampleOutput==TRUE) {
	### Put the sample names in the first column of the matrix
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length (if averaged, it contains the name of all the single spectra)
		if (fileFormat == "imzml" | fileFormat == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file[1])
			}
		}
		if (fileFormat == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	# Create the sample matrix column and append it to the global matrix
	sampleColumn <- matrix (0, ncol=1, nrow=spectraNumber)
	colnames (sampleColumn) <- "Sample"
	sampleColumn [,1] <- cbind (fileVector)
	signalMatrix <- cbind (signalMatrix, sampleColumn)
	### Add the class column
	classList <- sort (classList)
	classColumn <- matrix (0, ncol=1, nrow=spectraNumber)
	colnames (classColumn) <- "Class"
	# Rename the classes according to the classList vector
	classVector <- fileVector
	for (p in 1:length(classVector)) {
		for (w in 1:length(classList)) {
			if (length(grep(classList[w],classVector[p], ignore.case=TRUE)) == 1) {
				classVector[p] <- classList [w]
			}
		}
	}
	# Fill in the matrix column with the fileVector classes and samples
	classColumn [,1] <- cbind (classVector)
	signalMatrix <- cbind (signalMatrix, classColumn)
}
if (classOutput==TRUE && length(classList)>=1 && sampleOutput==FALSE) {
	classList <- sort (classList)
	### Put the sample names in the first column of the matrix
	# Create the empty vector
	fileVector <- vector (length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length (if averaged, it contains the name of all the single spectra)
		if (fileFormat == "imzml" | fileFormat == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$file[1])
			}
		}
		if (fileFormat == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	### Add the class column
	classColumn <- matrix (0, ncol=1, nrow=spectraNumber)
	colnames (classColumn) <- "Class"
	# Rename the classes according to the classList vector
	classVector <- vector (length=0)
	for (w in 1:length(classList)) {
		for (p in 1:length(fileVector)) {
			if (length(grep(classList[w],fileVector[p], ignore.case=TRUE)) !=0) {
				classVector[p] <- classList [w]
			}
		}
	}
	# Fill in the matrix column with the fileVector classes and samples
	classColumn [,1] <- cbind (classVector)
	signalMatrix <- cbind (signalMatrix, classColumn)
}
### Add these matrix columns to the peaklist matrix
return (signalMatrix)
}





########################################################################





###################################### NO HEMOGLOBIN
noHemoglobin <- function (spectra, hemoglobinMass=15400, tolerancePPM=2000, averageSpectrumEvaluation=FALSE) {
# Make the cluster (one for each core/thread)
cl <- makeCluster (CPUcoreNumber)
clusterEvalQ (cl, {library(MALDIquant)})
#################################################### FILTERING FUNCTION
hemoglobinRemovalFunction <- function (spectra, hemoglobinMass, tolerancePPM) {
## Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=3)
hemoglobinBin <- numeric()
# Do this if there are peaks
if (length(peaks@mass) > 0) {
for (s in 1:length(peaks@mass)) {
	# Select the Hemoglobin peaks (based on the difference from the exact mass)
	if ((abs(peaks@mass[s] - hemoglobinMass)*10^6/hemoglobinMass) <= tolerancePPM) {
		# That is the peak of interest, add its intensity to a bin
		hemoglobinBin <- append (hemoglobinBin, peaks@intensity[s])
	}
}
# Evaluate the peaks in the hemoglobin Bin (if there are any)...
if (length(hemoglobinBin) > 0) {
	# Start assuming that the spectrum is good
	isSpectrumToBeRemoved <- FALSE
	# If one of this is the base peak...
	for (h in 1:length(hemoglobinBin)) {
		if (hemoglobinBin[h] == max(peaks@intensity)) {
			# The spectrum is bad and to be removed
			isSpectrumToBeRemoved <- TRUE
		}
	}
	# Remove the spectrum if it is to be removed
	if (isSpectrumToBeRemoved == TRUE) {
		spectra <- NULL
	}
}
}
return (spectra)
}
########################################################################
spectraNoHemoglobin <- parLapply (cl, spectra, fun=function (spectra) hemoglobinRemovalFunction (spectra, hemoglobinMass, tolerancePPM))
# Close the processes
stopCluster(cl)
### Keep only the elements that are different from NULL
spectraNoHemoglobin <- spectraNoHemoglobin [!sapply(spectraNoHemoglobin, is.null)]
##### Now we have discarded all the spectra with Hemoglobin, but the average can still be bad
# Evaluate the average spectrum
if (averageSpectrumEvaluation == TRUE) {
	# Compute the average spectrum of this filtered dataset
	averageSpectrum <- averageMassSpectra (spectraNoHemoglobin, method="mean")
	## Peak picking on the AVG
	peaksAvg <- detectPeaks (averageSpectrum, method="MAD", SNR=3)
	hemoglobinBin <- numeric()
	# Do this if there are peaks
	if (length(peaksAvg@mass) > 0) {
	for (s in 1:length(peaksAvg@mass)) {
		# Select the Hemoglobin peaks (based on the difference from the exact mass)
		if ((abs(peaksAvg@mass[s] - hemoglobinMass)*10^6/hemoglobinMass) <= tolerancePPM) {
			# That is the peak of interest, add its intensity to a bin
			hemoglobinBin <- append (hemoglobinBin, peaksAvg@intensity[s])
		}
	}
	# Evaluate the peaks in the hemoglobin Bin (if there are any)...
	if (length(hemoglobinBin) > 0) {
	# Start assuming that the spectrum is good
	isSpectrumToBeRemoved <- FALSE
	# If one of this is the base peak...
	for (h in 1:length(hemoglobinBin)) {
		if (hemoglobinBin[h] == max(peaksAvg@intensity)) {
			# The spectrum is bad and to be removed
			isSpectrumToBeRemoved <- TRUE
		}
	}
	# Remove the spectra if they are to be removed
	if (isSpectrumToBeRemoved == TRUE) {
		spectraNoHemoglobin <- list()
	}
	}
	}
}
return (spectraNoHemoglobin)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE
biotyperPerformance <- function (sampleVector, classList, diseaseName="Diseased", healthyName="Controls") {
classList <- sort (classList)
# Declare the parameters
truePositive <- 0
trueNegative <- 0
falsePositive <- 0
falseNegative <- 0
# For each tested sample
for (s in 1:sampleNumber) {
	# For each class taken into account
	for (w in 1:length(classList)) {
		# If the sample and the class match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) & (length(grep(classList[w],diseaseName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true positive number increases
			truePositive <- truePositive + 1
		}
		# If the sample and the class match, the class is HP and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) & (length(grep(classList[w],healthyName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true negative number increases
			trueNegative <- trueNegative + 1
		}
		# If the sample and the class do not match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) & (length(grep(classList[w],diseaseName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false positive number increases
			falsePositive <- falsePositive + 1
		}
		# If the sample and the class do not match, the class is HP and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) & (length(grep(classList[w],healthyName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false negative number increases
			falseNegative <- falseNegative + 1
		}
    }
}
# Calculate the performance values
sensitivity <- (truePositive / (truePositive + falseNegative))*100
specificity <- (trueNegative / (falsePositive + trueNegative))*100
ppv <- truePositive / (truePositive + falsePositive)
npv <- trueNegative / (trueNegative + falseNegative)
fdr <- 1 - ppv
forValue <- 1 - npv
# Create the output matrix
testPerformance <- matrix (0, nrow=6, ncol=1)
testPerformanceRows <- c("Sensitivity", "Specificity", "PPV", "NPV", "FDR", "FOR")
rownames (testPerformance) <- testPerformanceRows
colnames (testPerformance) <- "Biotyper"
testPerformance [1,1] <- paste (sensitivity,"%")
testPerformance [2,1] <- paste (specificity, "%")
testPerformance [3,1] <- ppv
testPerformance [4,1] <- npv
testPerformance [5,1] <- fdr
testPerformance [6,1] <- forValue
return (testPerformance)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE 2
biotyperPerformance2 <- function (sampleVector, classList, diseaseName="Diseased", healthyName="Controls") {
classList <- sort (classList)
# Declare the parameters
truePositive <- 0
trueNegative <- 0
falsePositive <- 0
falseNegative <- 0
# For each tested sample
for (s in 1:sampleNumber) {
	# For each class taken into account
	for (w in 1:length(classList)) {
		# If the sample and the class match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) & (length(grep(classList[w],diseaseName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true positive number increases
			truePositive <- truePositive + 1
		}
		# If the sample and the class match, the class is HP and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) & (length(grep(classList[w],healthyName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true negative number increases
			trueNegative <- trueNegative + 1
		}
		# If the sample and the class do not match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) & (length(grep(classList[w],diseaseName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false positive number increases
			falsePositive <- falsePositive + 1
		}
		# If the sample and the class do not match, the class is HP and it is a YES (score >=2)
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) & (length(grep(classList[w],healthyName, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false negative number increases
			falseNegative <- falseNegative + 1
		}
	}
}
# Calculate the performance values
sensitivity <- (truePositive / (truePositive + falseNegative))*100
specificity <- (trueNegative / (falsePositive + trueNegative))*100
ppv <- truePositive / (truePositive + falsePositive)
npv <- trueNegative / (trueNegative + falseNegative)
fdr <- 1 - ppv
forValue <- 1 - npv
# Create the output matrix
testPerformance <- matrix (0, nrow=3, ncol=3)
rownames (testPerformance) <- c("Test outcome positive", "Test outcome negative", "")
colnames (testPerformance) <- c("Condition positive", "Condition negative", "")
testPerformance [1,1] <- paste ("True positive =", truePositive)
testPerformance [1,2] <- paste ("False positive =", falsePositive)
testPerformance [2,1] <- paste ("False negative =", falseNegative)
testPerformance [2,2] <- paste ("True negative =", trueNegative)
testPerformance [3,1] <- paste (sensitivity,"%")
testPerformance [3,2] <- paste (specificity, "%")
testPerformance [1,3] <- paste ("Positive predictive value =", ppv)
testPerformance [2,3] <- paste ("Negative predictive value =", npv)

return (testPerformance)
}





########################################################################





################# PLOT THE MEAN SPECTRUM WITH THE SD BARS ON THE AVERAGE
avgSpectrumBars <- function (spectra, SNR=5, tolerancePPM=2000, massRangePlot=c(4000,15000),
	graphTitle="Spectrum", avgSpectrumColour="black", peakPoints="yes", pointsColour="red",
	barWidth=40, barColour="blue") {
# Generate the average spectrum
averageSpectrum <- averageMassSpectra (spectra, method="mean")
averageSpectrum <- removeBaseline (averageSpectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaksAvg <- detectPeaks (averageSpectrum, method="MAD", SNR=SNR)
# Peak picking on the dataset
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
peaks <- removeEmptyMassObjects (peaks)
# Alignment: merge the two lists and align them all
peaksGlobal <- list()
peaksGlobal <- peaks
peaksGlobal <- append (peaksGlobal, peaksAvg)
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
# Empty the lists
peaksAvg <- list()
peaks <- list()
# Re-fill the lists with the aligned peaklists
for (i in 1:(length(peaksGlobal) - 1)) {
	peaks <- append (peaks, peaksGlobal[[i]])
}
peaksAvg <- peaksGlobal[[length(peaksGlobal)]]
######## Generate a vector with the standard deviations of the peaks in the average peaklist
stDevIntensity <- vector (length=0)
# For each peak in the average peaklist
for (a in 1:length(peaksAvg@mass)) {
	intensityVector <- vector (length=0)
	# Search for it in the dataset peaklist
	for (p in 1:length(peaks)) {
		for (i in 1:length(peaks[[p]]@mass)) {
			if (abs(peaks[[p]]@mass[i] == peaksAvg@mass[a])) {
				intensityVector <- append (intensityVector, peaks[[p]]@intensity[i])
			}
		}
	}
	stDevIntensity <- append (stDevIntensity, sd(intensityVector))
}
####### Plot
# Spectrum
plot (averageSpectrum, main=graphTitle, col.main=avgSpectrumColour, xlab="m/z", ylab="Intensity (a.i.)", xlim=massRangePlot, col=avgSpectrumColour)
# Peaks
if (peakPoints == "yes") {
	points (peaksAvg, pch=4, col=pointsColour)
}
# Bars
# epsilon: length of the horizontal segment
epsilon = barWidth
for (i in 1:length(peaksAvg@mass)) {
	# Define the upper and the lower limit of the vertical segment
	up <- peaksAvg@intensity[i] + stDevIntensity [i]
	low <- peaksAvg@intensity[i] - stDevIntensity [i]
	# Vertical bar (x,y x,y)
	segments (peaksAvg@mass[i], low, peaksAvg@mass[i], up, col=barColour)
	# Horizontal segments (x,y , x,y)
	segments (peaksAvg@mass[i]-epsilon, low, peaksAvg@mass[i]+epsilon, low, col=barColour)
	segments (peaksAvg@mass[i]-epsilon, up, peaksAvg@mass[i]+epsilon, up, col=barColour)
}
}





########################################################################





################# PLOT THE SIGNALS OF INTEREST WITH THE SD BARS ON THE AVERAGE
avgSpectrumBarsSignalsOfInterest <- function (spectra, SNR=5, signalsOfInterest=peaksAvg@mass,
	tolerancePPM=2000, halfWindowPlot=1000, graphTitle="Spectrum", avgSpectrumColour="black",
	peakPoints=TRUE, pointsColour="red", barWidth=40, barColour="blue") {
# Generate the average spectrum
averageSpectrum <- averageMassSpectra (spectra, method="mean")
averageSpectrum <- removeBaseline (averageSpectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaksAvg <- detectPeaks (averageSpectrum, method="MAD", SNR=SNR)
# Peak picking on the dataset
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
peaks <- removeEmptyMassObjects (peaks)
######## Generate a vector with the standard deviations of the peaks of interest
# For each peak in the average peaklist
for (a in 1:length(signalsOfInterest)) {
	# Search for the custom mass value in the average spectrum, in order to assign the real value
	for (j in 1:length(peaksAvg@mass)) {
		if ((abs(signalsOfInterest[a] - peaksAvg@mass[j])*10^6/peaksAvg@mass[j]) <= tolerancePPM) {
			signalsOfInterest[a] <- peaksAvg@mass[j]
			# Position of the peak of interest in the avg peak list
			peakIndex <- j
		}
	}
	# Generate the intensity vector
	intensityVector <- vector (length=0)
	# Search for it in the dataset peaklist
	for (p in 1:length(peaks)) {
		for (i in 1:length(peaks[[p]]@mass)) {
			if ((abs((signalsOfInterest[a] - peaks[[p]]@mass[i])/peaks[[p]]@mass[i])*10^6) <= tolerancePPM) {
				intensityVector <- append (intensityVector, peaks[[p]]@intensity[i])
			}
		}
	}
	# Calculate the standard deviation
	stDevIntensity <- sd(intensityVector)
	####### Plot
	# Define the limits on the x-axis
	lowMassPlot <- signalsOfInterest[a] - halfWindowPlot
	upMassPlot <- signalsOfInterest[a] + halfWindowPlot
	# Spectrum
	X11() #Opens a new graphic instance
	plot (averageSpectrum, main=graphTitle, col.main=avgSpectrumColour, xlab="m/z", ylab="Intensity (a.i.)", xlim=c(lowMassPlot,upMassPlot), col=avgSpectrumColour)
	# Peaks
	if (peakPoints == TRUE) {
		points (peaksAvg[peakIndex], pch=4, col=pointsColour)
	}
	# Bars
	# epsilon: length of the horizontal segment
	epsilon = barWidth
	# Define the upper and the lower limit of the vertical segment
	up <- peaksAvg@intensity[peakIndex] + stDevIntensity
	low <- peaksAvg@intensity[peakIndex] - stDevIntensity
	# Vertical bar (x,y x,y)
	segments (signalsOfInterest[a], low, signalsOfInterest[a], up, col=barColour)
	# Horizontal segments (x,y , x,y)
	segments (signalsOfInterest[a]-epsilon, low, signalsOfInterest[a]+epsilon, low, col=barColour)
	segments (signalsOfInterest[a]-epsilon, up, signalsOfInterest[a]+epsilon, up, col=barColour)
}
}





########################################################################





################################### MEMORY EFFICIENT IMPORTING
################## FOR EACH PATIENT: IMPORT, PREPROCESS, FILTER, AVERAGE, REPRESENTATIVE SPECTRA
memoryEfficientImport <- function (folder, tofMode="linear", ticPurification=FALSE, absoluteTICthreshold=0,
	smoothingStrength="medium", massRange=c(0,0), preProcessSpectra=TRUE, processInPackagesOf=length(spectra),
	generateRepresentativeSpectra=FALSE, spectraPerPatient=5, representativeAlgorithm="hca", discardedNodes=1,
	averagePatient=FALSE, skyline=FALSE, alignSpectra=FALSE, alignmentTolerancePPM=2000, fileFormat="imzml") {
setwd (folder)
folderFiles <- readSpectraFiles (folder, fileFormat=fileFormat, fullPath=FALSE)
spectraDataset <- list ()
# For each imzML file (patient)...
for (i in 1:length(folderFiles)) {
	### Load the spectra
	if (fileFormat == "imzml" | fileFormat == "imzML") {
		spectra <- importImzMl (folderFiles[i])
	}
	if (fileFormat == "brukerflex") {
		spectra <- importBrukerFlex (folderFiles[i])
	}
	if (ticPurification == TRUE) {
		spectra <- spectraTICpurification (spectra, tofMode=tofMode, absoluteTICthreshold=absoluteTICthreshold)
	}
	if (averagePatient == FALSE) {
	### Preprocessing
		if (preProcessSpectra == TRUE) {
			spectra <- preprocessSpectra (spectra, tofMode=tofMode, smoothingStrength=smoothingStrength, processInPackagesOf=processInPackagesOf)
		}
	}
	### Average the spectra
	if (averagePatient == TRUE) {
		# At least two spectra for the averaging
		if (length(spectra) > 1) {
			if (tofMode == "linear") {
				spectra <- averageMassSpectra (spectra, method="mean")
				spectra <- smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=10)
				spectra <- removeBaseline (spectra, method="TopHat")
				spectra <- calibrateIntensity (spectra, method="TIC")
			}
			if (tofMode == "reflector" || tofMode == "reflectron") {
				spectra <- averageMassSpectra (spectra, method="mean")
				spectra <- smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=2.5)
				spectra <- removeBaseline (spectra, method="TopHat")
				spectra <- calibrateIntensity (spectra, method="TIC")
			}
		}
		# If there is only one spectrum left, use it
		if (length(spectra) == 1) {
			spectra <- spectra
		}
	}
	if (skyline == TRUE) {
		# At least two spectra for the skyline
		if (length(spectra) > 1) {
			spectra <- generateSkylineSpectrum (spectra)
			spectra <- removeBaseline (spectra, method="TopHat")
			spectra <- calibrateIntensity (spectra, method="TIC")
		}
		# If there is only one spectrum left, use it
		if (length(spectra) == 1) {
			spectra <- spectra
		}
	}
	### Add this purified spectra list to a global list (only if there is something to actually add)
	if (length(spectra) >= 1) {
		spectraDataset <- append (spectraDataset, spectra)
	}
	### Free the memory
	rm (spectra)
	gc()
}
### Trimming the entire dataset
	if (massRange[1] == 0 && massRange[2] == 0) {
		spectraDataset <- trim (spectraDataset)
	}
	if ((massRange[1] != 0 || massRange[2] == 0) && massRange[2] != 0) {
		spectraDataset <- trim (spectraDataset, range=massRange)
	}
	if (massRange[1] != 0 && massRange[2] == 0) {
		upperCroppingBoundary <- Inf
		spectraDataset <- trim (spectraDataset, range=massRange)
	}
	spectraDataset <- calibrateIntensity (spectraDataset, method="TIC")
### Spectral alignment
if (alignSpectra==TRUE) {
	# Determine the warping function
	peaks <- detectPeaks (spectraDataset, SNR=5)
	referenceForAlignment <- referencePeaks (peaks, method="relaxed", tolerance=(alignmentTolerancePPM/10^6), minFrequency=0.9)
	if (length(referenceForAlignment@mass) > 0) {
		if (tofMode == "linear") {
			spectraDataset <- alignSpectra (spectraDataset, halfWindowSize=20, noiseMethod="MAD", SNR=3,
					reference=referenceForAlignment, tolerance=(alignmentTolerancePPM/10^6), warpingMethod="lowess")
		}
		if (tofMode == "reflector" | tofMode == "reflectron") {
			spectraDataset <- alignSpectra (spectraDataset, halfWindowSize=5, noiseMethod="MAD", SNR=3,
					reference=referenceForAlignment, tolerance=(alignmentTolerancePPM/10^6), warpingMethod="lowess")
		}
	}
}
### Generate representative spectra
if (generateRepresentativeSpectra == TRUE) {
	if (length(spectra) > 1) {
		spectra <- groupSpectra (spectra, spectraPerPatient=spectraPerPatient, fileFormat=fileFormat, tofMode=tofMode, seed=0, algorithm=representativeAlgorithm, discardedNodes=discardedNodes)
	}
	# If there is only one spectrum left, use it
	if (length(spectra) == 1) {
		spectra <- spectra
	}
}
return (spectraDataset)
}





########################################################################





#################### SIGNAL COUNTER IN THE AVERAGE SPECTRUM
##### Count signals in a mass range
countSignalsRangeAverageSpectrum <- function (spectra, SNR=5, lower=4000, intervalWidth=1000, max=20000) {
avg <- averageMassSpectra (spectra, method="mean")
avg <- removeBaseline (avg, method="TopHat")
# Peak picking on the avg
peaksAvg <- detectPeaks (avg, method="MAD", SNR=SNR)
# Vector with signals
peakVector <- peaksAvg@mass
#Create the resultMatrix
resultMatrix <- matrix (ncol=1, nrow=1)
resultMatrix [1,1] <- "Signal counter"
colnames (resultMatrix) <- "Signals"
# Count the signals in the specified mass range, with a counter
counter <- 0
upper <- lower + intervalWidth
for (n in 1:((max-lower)/intervalWidth)) {
	for (i in 1:length(peakVector)) {
		if (peakVector [i] >= lower & peakVector [i] <= upper & upper <= max) {
			counter <- counter + 1
		}
	}
	# Put the counter into a matrix column
	signalCounter <- matrix (ncol=1, nrow=1)
	colnames (signalCounter) <- paste (lower, "-", upper)
	signalCounter [1,1] <- counter
	# Add the signalCounter column to the resultMatrix, each time
	resultMatrix <- cbind (resultMatrix, signalCounter)
	# Move to the next interval
	lower <- lower + intervalWidth
	upper <- upper + intervalWidth
	counter <- 0
}
return (resultMatrix)
}





########################################################################





#################### SIGNAL COUNTER IN THE AVERAGE SPECTRUM
##### Count signals in a mass range
countSignalsRangeSpectra <- function (spectra, SNR=5, lower=4000, intervalWidth=1000, max=20000) {
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
#Create the resultMatrix
resultMatrix <- matrix ("Number of signals", ncol=(((max-lower)/intervalWidth)+1), nrow=1)
# Count the signals in the specified mass range, with a counter
counter <- 0
# The first upper value would be
upper <- lower + intervalWidth
# Then it will increase during the loop
# Scroll the "peaks" list and create a vector each time containing the peaklist of that spectrum
for (s in 1:length(peaks)) {
	peakVector <- peaks[[s]]@mass
	#Create the resultMatrixRow (this would be referred to each spectrum)
	resultMatrixRow <- matrix (ncol=1, nrow=1)
	# Reiterate for every interval (from 1 to the number of intervals)
	for (n in 1:((max-lower)/intervalWidth)) {
		# Check every value of the peakVector and adjust the counter based on the condition
		for (i in 1:length(peakVector)) {
			if (peakVector [i] >= lower & peakVector [i] <= upper & upper <= max) {
				counter <- counter + 1
			}
		}
		# Put the counter into a matrix column
		signalCounter <- matrix (ncol=1, nrow=1)
		colnames (signalCounter) <- paste (lower, "-", upper)
		signalCounter [1,1] <- counter
		# Add the signalCounter column to the resultMatrix, each time
		resultMatrixRow [1,1] <- s
		resultMatrixRow <- cbind (resultMatrixRow, signalCounter)
		# Move to the next interval
		lower <- lower + intervalWidth
		upper <- upper + intervalWidth
		counter <- 0
	}
	# Add this row to the resultMatrix
	resultMatrix <- rbind (resultMatrix, resultMatrixRow)
	# Restore the initial values, in order to reiterate for the next spectrum
	lower <- 4000
	counter <- 0
	upper <- lower + intervalWidth
}
resultMatrix [1,1] <- "Spectrum number"
return (resultMatrix)
}





########################################################################





################################## GENERATE A SKYLINE SPECTRUM
generateSkylineSpectrum <- function (spectra) {
# Generate the matrix containing all the intensities of the data points for each spectrum
spectraMatrix <- matrix (0, nrow=(length (spectra)+1), ncol=length(spectra[[1]]@mass))
colnames (spectraMatrix) <- spectra[[1]]@mass
# Fill in the matrix rows with the intensity values of the data points for each spectrum
for (s in 1:length(spectra)) {
	spectraMatrix [s,] <- spectra[[s]]@intensity
}
# In the last row, calculate the maximum intensity value of all the data points
for (i in 1:length(spectraMatrix [(length(spectra)+1),])) {
	spectraMatrix [(length(spectra)+1),i] <- max (spectraMatrix[,i])
}
# Generate the mean spectrum with the R function (it will be replaced by the skyline)
skylineSpectrum <- averageMassSpectra (spectra, method="mean")
# Replace the mean intensities with the maximum intensities
skylineSpectrum@intensity <- spectraMatrix [(length(spectra)+1),]
# Return the function
return (skylineSpectrum)
}





########################################################################





####################### LOWEST NUMBER OF DATAPOINTS
estimateLowestDataPoints <- function (spectra) {
lowestDataPoints <- length (spectra[[1]])
# Compare it with all the others
for (s in 2:length(spectra)) {
	if (length(spectra[[s]]@mass) < lowestDataPoints) {
		lowestDataPoints <- length(spectra[[s]]@mass)
	}
}
return (lowestDataPoints)
}





########################################################################





###################### SPECTRA BINNING
resampleSpectra <- function (spectra, finalDataPoints=lowestDataPoints, method="sum") {
####################################################### BINNING FUNCTION
binningFunction <- function (spectra, finalDataPoints, method) {
	# Create the new spectraBinned list
	spectraBinned <- spectra
	# Empty the mass and intensity values
	for (s in 1:length(spectraBinned)) {
		spectraBinned@mass <- numeric()
		spectraBinned@intensity <- numeric()
	}
	# Calculate the number of datapoints per bin
	dataPointsPerBin <- length (spectra@mass) / finalDataPoints
	dataPointsPerBin <- floor (dataPointsPerBin)
	# Define the indexes
	index1 <- 1
	index2 <- dataPointsPerBin
	# For each bin...
	for (d in 1:finalDataPoints) {
		# Create the (temporary) bin vectors (mass and intensity), where the data points will be stored for the binning
		binMass <- numeric()
		binIntensity <- numeric()
		# Scroll the data points, grouping them by bins
		for (i in index1:index2) {
			binMass <- append (binMass, spectra@mass[i])
			binIntensity <- append (binIntensity, spectra@intensity[i])
		}
		# Calculate the value of the new data point
		dataPointMass <- mean (binMass)
		if (method == "sum") {
			dataPointIntensity <- sum (binIntensity)
		}
		if (method == "mean") {
			dataPointIntensity <- mean (binIntensity)
		}
		if (method == "median") {
			dataPointIntensity <- median (binIntensity)
		}
		if (method == "RMS" | method == "rms") {
			dataPointIntensity <- sqrt(sum(binIntensity^2))
		}
		# Store it in the new spectraBinned list
		spectraBinned@mass <- append (spectraBinned@mass, dataPointMass)
		spectraBinned@intensity <- append (spectraBinned@intensity, dataPointIntensity)
		# Increase the indexes
		index1 <- index1 + dataPointsPerBin
		index2 <- index2 + dataPointsPerBin
	}
	return (spectraBinned)
}
#######################################################################
# Calculate the lowest amount of data points, that corresponds to the maximum
# number of data points that can be used for the binning
# Use this value as default if finalDataPoints is not specified (Function)
lowestDataPoints <- estimateLowestDataPoints (spectra)
# Check the qualities of the spectral dataset
dataPointsTable <- table (sapply (spectra, length))
datapointsDataset <- as.numeric (names (dataPointsTable))
equalityDataPoints <- length (dataPointsTable)
######## Do not bin if all the spectra are of the same length and have the same number of datapoints as defined
cl <- makeCluster (CPUcoreNumber)
if ((equalityDataPoints == 1) && (datapointsDataset == finalDataPoints)) {
	spectraBinned <- spectra
}
######## Perform the binning if the spectra are not of the same length or
# they are of the same length but with a different number of datapoints than
# the desired one
if (!(equalityDataPoints == 1) || ((equalityDataPoints == 1) && !(datapointsDataset == finalDataPoints))) {
	# Do the binning only if the number of final data points is lower than the
	# lowest number of original data points
	if (finalDataPoints > lowestDataPoints) {
		finalDatapoints <- lowestDataPoints
		print ("Binning at this sample rate is not possible, the highest number of data points possible will be used")
		spectraBinned <- parLapply (cl, spectra, fun=function (spectra) binningFunction (spectra, finalDataPoints, method))
	}
	if (finalDataPoints <= lowestDataPoints) {
		spectraBinned <- parLapply (cl, spectra, fun=function (spectra) binningFunction (spectra, finalDataPoints, method))
	}
}
# Close the processes
stopCluster(cl)
return (spectraBinned)
print (table (sapply (spectraBinned, length)))
print (paste ("Equal distance between datapoints", (all (sapply (spectraBinned, isRegular)))))
}





########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N MOST INTENSE
mostIntenseSignals <- function (spectra, signalsToTake=20) {
####################################################### PICKING FUNCTION
pickingFunction <- function (peaks, signalsToTake) {
	# Create a dataframe with mass and intensity
	peaksDf <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Sort the dataframe according to the SNR
	peaksDf <- peaksDf [order(-peaksDf$snr),]
	# Select only the first most intense signals
	selectedSignals <- peaksDf [1:signalsToTake,]
	# Sort the dataframe back according to mass
	selectedSignals <- selectedSignals [order(selectedSignals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selectedSignals$mass
	peaks@intensity <- selectedSignals$intensity
	peaks@snr <- selectedSignals$snr
	return (peaks)
}
########################################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=3)
peaksMostIntense <- mclapply (peaks, FUN= function (peaks) pickingFunction (peaks, signalsToTake=signalsToTake))
return (peaksMostIntense)
}






########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N LESS INTENSE
lessIntenseSignals <- function (spectra, baseSNR, signalsToTake) {
############################################################
pickingFunction <- function (peaks, signalsToTake) {
	# Create a dataframe with mass and intensity
	peaksDf <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Sort the dataframe according to the SNR
	peaksDf <- peaksDf [order(peaksDf$snr),]
	# Select only the first most intense signals
	selectedSignals <- peaksDf [1:signalsToTake,]
	# Sort the dataframe back according to mass
	selectedSignals <- selectedSignals [order(selectedSignals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selectedSignals$mass
	peaks@intensity <- selectedSignals$intensity
	peaks@snr <- selectedSignals$snr
	return (peaks)
}
###########################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=baseSNR)
peaksLessIntense <- mclapply (peaks, FUN= function (peaks) pickingFunction (peaks, signalsToTake=signalsToTake))
return (peaksLessIntense)
}





########################################################################
####### LESS AND MOST INTENSE SIGNALS
highestAndLowestIntensitySignals <- function (spectra, baseSNR=3, signalsToTake) {
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=baseSNR)
# For each peaklist...
for (p in 1:length(peaks)) {
	# Create a dataframe with mass and intensity
	peaksDf <- data.frame (mass = peaks[[p]]@mass, intensity = peaks[[p]]@intensity, snr = peaks[[p]]@snr)
	# Sort the dataframe according to the SNR
	peaksDfMostIntense <- peaksDf [order(-peaksDf$snr),]
	peaksDfLessIntense <- peaksDf [order(peaksDf$snr),]
	# Select only the first most intense signals
	selectedSignalsMostIntense <- peaksDfMostIntense [1:(signalsToTake/2),]
	selectedSignalsLessIntense <- peaksDfLessIntense [1:(signalsToTake/2),]
	# Sort the dataframe back according to mass
	selectedSignalsMostIntense <- selectedSignalsMostIntense [order(selectedSignalsMostIntense$mass),]
	selectedSignalsLessIntense <- selectedSignalsLessIntense [order(selectedSignalsLessIntense$mass),]
	# Put these signals back into the peaklist
	peaks[[p]]@mass <- selectedSignalsLessIntense$mass
	peaks[[p]]@intensity <- selectedSignalsLessIntense$intensity
	peaks[[p]]@snr <- selectedSignalsLessIntense$snr
	peaks[[p]]@mass <- append (peaks[[p]]@mass, selectedSignalsMostIntense$mass)
	peaks[[p]]@intensity <- append (peaks[[p]]@intensity, selectedSignalsMostIntense$intensity)
	peaks[[p]]@snr <- append (peaks[[p]]@snr, selectedSignalsMostIntense$snr)
}
return (peaks)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE
biotyperPerformanceWide <- function (biotyperScore) {
classList <- colnames (biotyperScore)
sampleVector <- rownames (biotyperScore)
sampleNumber <- length (sampleVector)
# Declare the parameters
correctlyIdentified <- 0
falselyIdentified <- 0
notIdentified <- 0
# For each tested sample
for (s in 1:sampleNumber) {
	# For each class taken into account
	for (w in 1:length(classList)) {
		# If the sample and the class match and it is a YES
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) && (length(grep("YES", score[s,w], ignore.case=TRUE)) == 1)) {
			# The correctly identified number increases
			correctlyIdentified <- correctlyIdentified + 1
		}
		# If the sample and the class match and it is a NI
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) && (length(grep("NI", 	score[s,w], ignore.case=TRUE)) == 1)) {
			# The correctly identified number increases
			correctlyIdentified <- correctlyIdentified + 0.5
		}
		# If the sample and the class do not match and it is a YES
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) && (length(grep("YES", score[s,w], ignore.case=TRUE)) == 1)) {
			# The falsely identified number increases
			falselyIdentified <- falselyIdentified + 1
		}
		# If the sample and the class do not match and it is a NI
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 0) && (length(grep("NI", score[s,w], ignore.case=TRUE)) == 1)) {
			# The falsely identified number increases
			falselyIdentified <- falselyIdentified + 0.5
		}
		# If the sample and the class match and it is a NO
		if ((length(grep(classList[w],sampleVector[s], ignore.case=TRUE)) == 1) & (length(grep("NO", score[s,w], ignore.case=TRUE)) == 1)) {
			# The not identified number increases
			notIdentified <- notIdentified + 1
		}
	}
}
# Calculate the performance values
sensitivity <- (correctlyIdentified / sampleNumber)*100
falsePositiveRate <- (falselyIdentified / sampleNumber)*100
missedIdentifications <- (notIdentified / sampleNumber)*100
# Create the output matrix
testPerformance <- matrix (0, nrow=3, ncol=1)
testPerformanceRows <- c("Sensitivity", "False positive rate", "Missed Identifications")
rownames (testPerformance) <- testPerformanceRows
colnames (testPerformance) <- "Biotyper"
testPerformance [1,1] <- paste (sensitivity,"%")
testPerformance [2,1] <- paste (falsePositiveRate, "%")
testPerformance [3,1] <- paste (missedIdentifications, "%")
return (testPerformance)
}





########################################################################





####################### SPECTRA LIST CLASS
spectraListClass <- function (spectra, className, fileFormat="imzml") {
spectraList <- list()
if (fileFormat == "imzml" | fileFormat == "imzML") {
	for (s in 1:length(spectra)) {
		if (length(grep(className,spectra[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
			spectraList <- append (spectraList, spectra[[s]])
		}
	}
}
if (fileFormat == "brukerflex") {
	for (s in 1:length(spectra)) {
		if (length(grep(className,spectra[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
			spectraList <- append (spectraList, spectra[[s]])
		}
	}
}
return (spectraList)
}




########################################################################





########################### ALIGN THE DATAPOINTS AFTER TRUNCATION OR BINNING
alignDataPointsAfterTruncation <- function (spectra, tofMode="linear") {
# Make the cluster (one for each core/thread)
cl <- makeCluster (CPUcoreNumber)
########################################################################
dataPointsAlignmentFunction <- function (spectra, referenceSpectrum, tofMode) {
dataPointsDifference <- spectra@mass[1] - referenceSpectrum@mass[1]
if (tofMode == "linear") {
	if (dataPointsDifference <= 4) {
		for (m in 1:length(spectra@mass)) {
			spectra@mass[m] <- spectra@mass[m] - dataPointsDifference
		}
	}
}
if (tofMode == "reflector") {
	if (dataPointsDifference <= 1) {
		for (m in 1:length(spectra@mass)) {
			spectra@mass[m] <- spectra@mass[m] - dataPointsDifference
		}
	}
}
return (spectra)
}
########################################################################
#### Align the data points after the truncation
referenceSpectrum <- spectra[[1]]
spectra <- parLapply (cl, spectra, fun=function(spectra) dataPointsAlignmentFunction (spectra, referenceSpectrum, tofMode=tofMode))
# Close the processes
stopCluster(cl)
return (spectra)
}





##################################### BIOTYPER SCORE
biotyperScore <- function (filepathSamples, peaksTest, peaksLibrary, classList, tolerancePPM=2000, intensityTolerancePercent=50, fileFormat="brukerflex", spectraPathOutput=TRUE, replicatesAverage=FALSE, scoreOnly=TRUE, standardDeviationEvaluation=FALSE, nStDev=1) {
sampleNumber <- length (peaksTest)
librarySize <- length (peaksLibrary)
# Create the sample vector
sampleVector <- vector(length=sampleNumber)
if (fileFormat == "brukerflex") {
	if (replicatesAverage == FALSE) {
		for (s in 1:sampleNumber) {
			sampleVector[s] <- peaksTest[[s]]@metaData$sampleName
		}
	}
	# If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	if (replicatesAverage == TRUE) {
		for (s in 1:sampleNumber) {
			sampleVector[s] <- peaksTest[[s]]@metaData$sampleName[1]
		}
	}
}
if (fileFormat == "imzml" | fileFormat == "imzML") {
	for (s in 1:sampleNumber) {
		sampleVector[s] <- peaksTest[[s]]@metaData$file
	}
}
# Read the list of spectra folders in the sample mother folder
spectraPathVector <- readSpectraFiles (filepathSamples, fileFormat="brukerflex")
################# PEAK ALIGNMENT
# Create a global list for the alignment (with library and tests)
peaksGlobal <- peaksLibrary
peaksGlobal <- append (peaksGlobal, peaksTest)
# Alignment (R function)
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksGlobal <- binPeaks (peaksGlobal, method="relaxed", tolerance=(tolerancePPM/10^6))
# Empty the two original lists
peaksLibrary <- list()
peaksTest <- list()
# Re-fill the peaklists with the aligned peaks
for (i in 1:librarySize) {
	peaksLibrary <- append (peaksLibrary, peaksGlobal[[i]])
}
for (i in (librarySize+1):length(peaksGlobal)) {
	peaksTest <- append (peaksTest, peaksGlobal[[i]])
}
############################################################ SCORE (FRI)
if (standardDeviationEvaluation == FALSE) {
	# Compare the peaks in the single sample peaklists with the peaks in each peaklist in the library
	# COUNTER 1 - FIT
	# Create a counter, symmetrical to the database Peaklist
	counter1 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each sample
	for (s in 1:sampleNumber) {
		# For each peaklist in the Library
		for (l in 1:librarySize) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaksTest[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaksLibrary[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
					counter1 [s,l] <- counter1 [s,l] + 1
					}
				}
			}
		}
	}
	for (s in 1:sampleNumber) {
		for (l in 1:librarySize) {
			counter1[s,l] <- counter1[s,l] / length(peaksTest[[s]]@mass)
		}
	}
	# Compare the peaks in the library peaklists with the peaks in each sample
	# COUNTER 2 - RETRO FIT
	# Create a counter, symmetrical to the database Peaklist
	counter2 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each peaklist in the library
	for (l in 1:librarySize) {
		# For each sample
		for (s in 1:sampleNumber) {
			# Scroll the peaks in the library
			for (j in 1:length(peaksLibrary[[l]]@mass)) {
				# Scroll the peaklist in the sample
				for (i in 1:length(peaksTest[[s]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
						counter2 [s,l] <- counter2 [s,l] + 1
					}
				}
			}
		}
	}
	for (l in 1:librarySize) {
		for (s in 1:sampleNumber) {
			counter2[s,l] <- counter2[s,l] / length(peaksLibrary[[l]]@mass)
		}
	}
	# COUNTER 3
	# Symmetry -> comparison between intensities
	# Create a counter, symmetrical to the database Peaklist
	counter3 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each sample
	for (s in 1:sampleNumber) {
		# For each peaklist in the Library
		for (l in 1:librarySize) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaksTest[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaksLibrary[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
						# Evaluate the difference in intensity
						if ((abs(peaksTest[[s]]@intensity[i] - peaksLibrary[[l]]@intensity[j])*100/peaksLibrary[[l]]@intensity[j]) < intensityTolerancePercent) {
							counter3[s,l] <- counter3[s,l] + 1
						}
					}
				}
			}
		}
	}
	for (s in 1:sampleNumber) {
		for (l in 1:librarySize) {
			counter3[s,l] <- counter3[s,l] / length(peaksTest[[s]]@mass)
		}
	}
	### Score calculation
	score <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	for (i in 1:sampleNumber) {
		for (j in 1:length(peaksLibrary)) {
			score[i,j] <- log10(counter1[i,j]*counter2[i,j]*counter3[i,j]*1000)
		}
	}
	#### Output the classification
	output <- matrix ("NO", nrow=sampleNumber, ncol=librarySize)
	colnames (output) <- sort(classList)
	if ((spectraPathOutput == FALSE) || (replicatesAverage == TRUE)) {
		rownames (output) <- sampleVector
	}
	if ((spectraPathOutput == TRUE) && (replicatesAverage == FALSE)) {
		rownames (output) <- spectraPathVector
	}
	if (scoreOnly == TRUE) {
		for (r in 1:sampleNumber) {
			for (w in 1:length(peaksLibrary)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste ("YES","(", round (score[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste ("NO", "(", round (score[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste ("NI","(", round (score[r,w], digits=3), ")")
				}
			}
		}
	}
	if (scoreOnly == FALSE) {
		for (r in 1:sampleNumber) {
			for (w in 1:length(peaksLibrary)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste ("YES","(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste ("NO", "(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste ("NI","(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
			}
		}
	}
}
####################################################### SCORE (FR-STDEV)
############################ peaksLibrary SHULD HAVE SNR REPLACED BY SD!
if (standardDeviationEvaluation == TRUE) {
	# Compare the peaks in the single sample peaklists with the peaks in each peaklist in the library
	# COUNTER 1 - FIT
	# Create a counter, symmetrical to the database Peaklist
	counter1 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each sample
	for (s in 1:sampleNumber) {
		# For each peaklist in the Library
		for (l in 1:librarySize) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaksTest[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaksLibrary[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
					counter1 [s,l] <- counter1 [s,l] + 1
					}
				}
			}
		}
	}
	for (s in 1:sampleNumber) {
		for (l in 1:librarySize) {
			counter1[s,l] <- counter1[s,l] / length(peaksTest[[s]]@mass)
		}
	}
	# Compare the peaks in the library peaklists with the peaks in each sample
	# COUNTER 2 - RETRO FIT
	# Create a counter, symmetrical to the database Peaklist
	counter2 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each peaklist in the library
	for (l in 1:librarySize) {
		# For each sample
		for (s in 1:sampleNumber) {
			# Scroll the peaks in the library
			for (j in 1:length(peaksLibrary[[l]]@mass)) {
				# Scroll the peaklist in the sample
				for (i in 1:length(peaksTest[[s]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
						counter2 [s,l] <- counter2 [s,l] + 1
					}
				}
			}
		}
	}
	for (l in 1:librarySize) {
		for (s in 1:sampleNumber) {
			counter2[s,l] <- counter2[s,l] / length(peaksLibrary[[l]]@mass)
		}
	}
	# COUNTER 3
	# Symmetry -> comparison between intensities
	# Create a counter, symmetrical to the database Peaklist
	counter3 <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	# For each sample
	for (s in 1:sampleNumber) {
		# For each peaklist in the Library
		for (l in 1:librarySize) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaksTest[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaksLibrary[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaksTest[[s]]@mass[i] == peaksLibrary[[l]]@mass[j]) {
						# Evaluate the difference in intensity
						if ((peaksTest[[s]]@intensity[i] <= (peaksLibrary[[l]]@intensity[j] + (nStDev * peaksLibrary[[l]]@snr[j]))) && (peaksTest[[s]]@intensity[i] >= (peaksLibrary[[l]]@intensity[j] - (nStDev * peaksLibrary[[l]]@snr[j])))) {
							counter3[s,l] <- counter3[s,l] + 1
						}
					}
				}
			}
		}
	}
	for (s in 1:sampleNumber) {
		for (l in 1:librarySize) {
			counter3[s,l] <- counter3[s,l] / length(peaksTest[[s]]@mass)
		}
	}
	### Score calculation
	score <- matrix (0, nrow=sampleNumber, ncol=librarySize)
	for (i in 1:sampleNumber) {
		for (j in 1:length(peaksLibrary)) {
			score[i,j] <- log10(counter1[i,j]*counter2[i,j]*counter3[i,j]*1000)
		}
	}
	#### Output the classification
	output <- matrix ("NO", nrow=sampleNumber, ncol=librarySize)
	colnames (output) <- sort(classList)
	if ((spectraPathOutput == FALSE) || (replicatesAverage == TRUE)) {
		rownames (output) <- sampleVector
	}
	if ((spectraPathOutput == TRUE) && (replicatesAverage == FALSE)) {
		rownames (output) <- spectraPathVector
	}
	if (scoreOnly == TRUE) {
		for (r in 1:sampleNumber) {
			for (w in 1:length(peaksLibrary)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste ("YES","(", round (score[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste ("NO", "(", round (score[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste ("NI","(", round (score[r,w], digits=3), ")")
				}
			}
		}
	}
	if (scoreOnly == FALSE) {
		for (r in 1:sampleNumber) {
			for (w in 1:length(peaksLibrary)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste ("YES","(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste ("NO", "(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste ("NI","(Score:", round (score[r,w], digits=3), "), ","(Fit:", round (counter1[r,w], digits=3), "), ","(Retrofit:", round (counter2[r,w], digits=3), "), ","(Intensity:", round (counter3[r,w], digits=3), ")")
				}
			}
		}
	}
}
########################################################################
return (output)
}





########################################################################





###################### SIGNAL FOLLOWER
signalFollowerStatistics <- function (filepath, signalList, massLabels=list(), SNR=5, fileFormat="imzml") {
# Define the column and row headers (if mass labels are provided or not)
if (length(massLabels) != 0) {
	columnNamesStDev <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesStDev[i] <- paste (signalList[i], "-", massLabels[i], "- StDev")
	}
	columnNamesCoeffVar <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesCoeffVar[i] <- paste (signalList[i], "-", massLabels[i], "- CV")
	}
	columnNamesMeanAbsDev <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesMeanAbsDev[i] <- paste (signalList[i], "-", massLabels[i], "- MeanAbsDev")
	}
}
if (length(massLabels) == 0) {
	columnNamesStDev <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesStDev[i] <- paste (signalList[i], "- StDev")
	}
	columnNamesCoeffVar <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesCoeffVar[i] <- paste (signalList[i], "- CV")
	}
	columnNamesMeanAbsDev <- vector (length=length(signalList))
	for (i in 1:length(signalList)) {
		columnNamesMeanAbsDev[i] <- paste (signalList[i], "- MeanAbsDev")
	}
}
sampleNames <- filepath
### Create the result matrices
stDevResultMatrix <- matrix (ncol=length(signalList), nrow=length(filepath))
colnames (stDevResultMatrix) <- columnNamesStDev
rownames (stDevResultMatrix) <- sampleNames
coeffVarResultMatrix <- matrix (ncol=length(signalList), nrow=length(filepath))
colnames (coeffVarResultMatrix) <- columnNamesCoeffVar
rownames (coeffVarResultMatrix) <- sampleNames
meanAbsDevResultMatrix <- matrix (ncol=length(signalList), nrow=length(filepath))
colnames (meanAbsDevResultMatrix) <- columnNamesMeanAbsDev
rownames (meanAbsDevResultMatrix) <- sampleNames
# The script is run for every imzML file
for (lib in 1:length(filepath)) {
	# Spectra import
	spectra <- importImzMl (filepath[lib])
	### Smoothing (Savitzky-Golay filter, with window size 10, 21 points)
	spectra <- smoothIntensity (spectra, method="SavitzkyGolay")
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
	### Normalisation (TIC)
	spectra <- calibrateIntensity (spectra, method="TIC")
	### Peak Picking and alignment
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
	#### Create the partial result matrices
	stDevMatrix <- matrix (ncol=length(signalList), nrow=1)
	colnames (stDevMatrix) <- columnNamesStDev
	rownames (stDevMatrix) <- spectra[[1]]@metaData$file
	coeffVarMatrix <- matrix (0, ncol=length(signalList), nrow=1)
	colnames (coeffVarMatrix) <- columnNamesCoeffVar
	rownames (coeffVarMatrix) <- spectra[[1]]@metaData$file
	meanAbsDevMatrix <- matrix (0, ncol=length(signalList), nrow=1)
	colnames (meanAbsDevMatrix) <- columnNamesMeanAbsDev
	rownames (meanAbsDevMatrix) <- spectra[[1]]@metaData$file
	# For each signal in the signalList
	for (s in 1:length(signalList)) {
		intensityVector <- vector(length=0)
		# For each peaklist
		for (p in 1:length(peaks)) {
			# Scroll the peaklist
			for (m in 1:length(peaks[[p]]@mass)) {
				# If there is a match
				if (abs(signalList[s] - peaks[[p]]@mass[m]) <= 4) {
					# Put the intensity of that peak in the peaklist in a vector
					intensityVector <- append (intensityVector, peaks[[p]]@intensity[m])
				}
			}
		}
		# Calculate the standard deviation and the coefficient of variation
		stDevInt <- sd (intensityVector, na.rm=TRUE)
		meanAbsDevInt <- mad (intensityVector, na.rm=TRUE)
		meanInt <- mean (intensityVector, na.rm=TRUE)
		medianInt <- median (intensityVector, na.rm=FALSE)
		coeffVarInt <- (stDevInt / meanInt)*100
		# Fill in the partial result matrices with the values
		stDevMatrix [1,s] <- paste (meanInt, "+/-", stDevInt)
		coeffVarMatrix [1,s] <- paste (coeffVarInt, "%", "( mean =", meanInt, ")")
		meanAbsDevMatrix [1,s] <- paste (medianInt, "+/-", meanAbsDevInt)
	}
	# Put the partial result matrices together in one matrix
	stDevResultMatrix [lib,] <- stDevMatrix
	coeffVarResultMatrix [lib,] <- coeffVarMatrix
	meanAbsDevResultMatrix [lib,] <- meanAbsDevMatrix
}
resultMatrix <- stDevResultMatrix
resultMatrix <- cbind (resultMatrix, coeffVarResultMatrix)
resultMatrix <- cbind (resultMatrix, meanAbsDevResultMatrix)
return (resultMatrix)
}





########################################################################





################# READ THE LIST OF imzML FILES FROM A FOLDER
readSpectraFiles <- function (folder, fileFormat="imzml", fullPath=FALSE) {
if (fileFormat=="imzml") {
	# Read all the files
	folderAllFiles <- list.files (folder, full.names=fullPath)
	# Create the empty vector in which only the imzML files will be listed
	spectraFiles <- character()
	# Put the imzML files in the new vector discarding the ibd files
	for (l in 1:length (folderAllFiles)) {
		if (length(grep(".imzML", folderAllFiles [l], fixed=TRUE)) == 1) {
			spectraFiles <- append (spectraFiles, folderAllFiles [l])
		}
	}
}
if (fileFormat=="brukerflex") {
	# Read all the folders
	sampleDirectoryListAllFolders <- dir (folder, ignore.case=TRUE, full.names=fullPath, recursive=TRUE, include.dirs=TRUE)
	# Take only the fid paths (spectra)
	spectraFiles <- character()
	for (i in 1:length (sampleDirectoryListAllFolders)) {
		if ((length(grep("fid", sampleDirectoryListAllFolders[i], ignore.case=TRUE)) == 1)) {
			spectraFiles <- append (spectraFiles, sampleDirectoryListAllFolders[i])
		}
	}
}
return (spectraFiles)
}





########################################################################





################### AVERAGE THE REPLICATES (BY NAME)
averageReplicatesByName <- function (spectra, patternForReplicates="nameOfTheFile", numberOfCharactersForSimilarity=10, fileFormat="brukerflex") {
## Put the filenames in a vector
# Create the empty vector
sampleVector <- vector (length=0)
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(spectra)) {
	sampleVector <- append (sampleVector, spectra[[i]]@metaData$sampleName)
}
if (patternForReplicates == "nameOfTheFile") {
	for (s in 1: length(sampleVector)) {
		refSample <- sampleVector [s]
		for (j in s: length (sampleVector)) {
			# Determine the matching with each other sample
			matching <- pmatch (unlist(strsplit(refSample,"")), unlist(strsplit(sampleVector[j],"")))
			# If the FIRST characters match but the LAST character does not match, it is a replicate
			trueMatching <- seq (from=1, to=numberOfCharactersForSimilarity)
			if ((length(intersect(trueMatching, matching)) == length(trueMatching)) && is.na (matching[length(matching)])) {
				sampleVector [j] <- refSample
			}
		}
	}
}
if ((patternForReplicates != "nameOfTheFile") || (length(patternForReplicates) != 1)) {
	for (p in 1:length(patternForReplicates)) {
		for (s in 1: length(sampleVector)) {
			refSample <- sampleVector [s]
			for (j in s: length (sampleVector)) {
				# Determine the matching with each other sample
				matching <- pmatch (unlist(strsplit(refSample,"")), unlist(strsplit(sampleVector[j],"")))
				# If the FIRST characters match but the LAST character does not match, it is a replicate
				trueMatching <- seq (from=1, to=numberOfCharactersForSimilarity)
				if ((length(intersect(trueMatching, matching)) == length(trueMatching)) && is.na (matching[length(matching)]) && length(grep(patternForReplicates[p], refSample, ignore.case=TRUE)) == 1 && length(grep(patternForReplicates[p], sampleVector [j], ignore.case=TRUE))) {
					sampleVector [j] <- refSample
				}
			}
		}
	}
}
spectraReplicatesAveraged <- averageMassSpectra (spectra, labels=sampleVector)
return (spectraReplicatesAveraged)
}





###########################################################################





########################################## PEAKLIST FILTERING
# Remove all the bad peaklists, keep only the best ones
filterPeaklist <- function (peaks, SNRfilter=15, signalThreshold=15) {
if (SNRfilter != 0 && signalThreshold != 0) {
	# Create the clusters
	cl <- makeCluster (CPUcoreNumber)
	############################################# FILTERING FUNCTION
	peaksFilteringFunction <- function (peaks, SNRfilter, signalThreshold) {
		### Create a dataframe with the snr
		peaksDf <- data.frame (snr = peaks@snr)
		### elect only the snr more than the threshold
		peaksDf <- peaksDf [peaksDf$snr >= SNRfilter,]
		### Select only the desired peaklists, based upon the threshold characteristics
			if (length(peaksDf) < signalThreshold) {
				# Remove the bad spectrum from the list
				peaks <- NULL
			}
		return (peaks)
	}
	################################################################
	peaksPurified <- parLapply (cl, peaks, fun=function(peaks) peaksFilteringFunction(peaks, SNRfilter=SNRfilter, signalThreshold=signalThreshold))
	### Keep only the elements that are different from NULL
	peaksPurified <- peaksPurified [!sapply(peaksPurified, is.null)]
	# Stop the clusters
	stopCluster (cl)
	return (peaksPurified)
}
}





########################################################################





#################################### SPECTRA GROUPING (PEAKLIST MERGING)
mergePeaklist <- function (peaks, fileFormat="imzml", tofMode="linear") {
### Alignment
if (tofMode == "linear") {
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
}
if (tofMode == "reflectron" || tofMode == "reflector") {
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
}
if (fileFormat == "imzml" | fileFormat == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		fileVector <- append (fileVector, peaks[[i]]@metaData$file)
	}
	peaksMerged <- mergeMassPeaks (peaks, labels=fileVector, method="mean", ignore.na=FALSE)
}
if (fileFormat == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
	}
	peaksMerged <- mergeMassPeaks (peaks, labels=fileVector, method="mean", ignore.na=FALSE)
}
return (peaksMerged)
}





########################################################################





############## SPECTRA GROUPING CLASS (PEAKLIST MERGING)
mergePeaklistClass <- function (peaks, classList, fileFormat="imzml", tofMode="linear") {
### Alignment
if (tofMode == "linear") {
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.002)
}
if (tofMode == "reflectron" || tofMode == "reflector") {
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks (peaks, method="relaxed", tolerance=0.0002)
}
if (fileFormat == "imzml" | fileFormat == "imzML") {
	# Replace the name with the class name
	peaks <- replaceClassName (peaks, classList=classList, fileFormat="imzml")
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		fileVector <- append (fileVector, peaks[[i]]@metaData$file)
	}
	peaksMerged <- mergeMassPeaks (peaks, labels=fileVector, method="mean", ignore.na=FALSE)
}
if (fileFormat == "brukerflex") {
	# Replace the name with the class name
	peaks <- replaceClassName (peaks, classList=classList, fileFormat="brukerflex")
	# Put the filenames in a vector
	# Create the empty vector
	fileVector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
	}
	peaksMerged <- mergeMassPeaks (peaks, labels=fileVector, method="mean", ignore.na=FALSE)
}
return (peaksMerged)
}






########################################################################





######################### PEAK INTENSITY PLOT
peakIntensityPlot <- function (peaks, peakOfInterest, classList=list(), tolerancePPM=2000) {
if (length(classList) > 2) {
	print ("This algorithm can be computed only with no classes or one class")
}
# If there are one or no classes, compute the plot on the entire dataset
if (length(classList) == 0 || length(classList) == 1) {
	intensityVector <- numeric()
	# Scroll the peaklists and add the corresponding intensities (of the peak of interese) to a vector
	for (s in 1:length(peaks)) {
		for (p in 1:length(peaks[[s]]@mass)) {
			if (abs(peaks[[s]]@mass[p] - peakOfInterest)*10^6/peakOfInterest <= tolerancePPM) {
				intensityVector <- append (intensityVector, peaks[[s]]@intensity[p])
			}
		}
	}
	# Plot the results
	peakIntensityPlot <- plot (intensityVector)
}
# If there are one or no classes, compute the plot on the entire dataset
if (length(classList) == 2) {
	classNames <- classList
	intensityVector <- list(class1=numeric(), class2=numeric())
	for (w in 1:length(classList)) {
		# Scroll the peaklists and add the corresponding intensities (of the peak of interese) to a vector
		for (s in 1:length(peaks)) {
			# Check if the peaklist belongs to the class
			if (length(grep(classList[w], peaks[[s]]@metaData$file)) == 1) {
				for (p in 1:length(peaks[[s]]@mass)) {
					if (abs(peaks[[s]]@mass[p] - peakOfInterest)*10^6/peakOfInterest <= tolerancePPM) {
						intensityVector[[w]] <- append (intensityVector[[w]], peaks[[s]]@intensity[p])
					}
				}
			}
		}
	}
	# Plot the results
	peakIntensityPlot <- plot (x=intensityVector[[1]], xlab=classList[1], y=intensityVector[[2]], ylab=classList[2])
}
return (peakIntensityPlot)
}





########################################################################





##################### AVERAGE THE REPLICATES (BY FOLDER)

averageReplicatesByFolder <- function (spectra, folder, fileFormat="brukerflex") {
# List the spectra files
folderFiles <- readSpectraFiles (folder, fileFormat=fileFormat, fullPath=TRUE)
# Split the path into individual folders
folderFilesSplitted <- list()
# Split the paths into folders
for (f in 1:length(folderFiles)) {
	folderFilesSplitted[f] <- strsplit (folderFiles[f], "/")
}
# Check for replicates in the same folder
for (f in 1:(length(folderFilesSplitted)-1)) {
	# If the spectrum folder name is different between two samples but the folder above is the same, they are replicates (based on the brukerflex folder structure)
	if ((folderFilesSplitted[[f]] [length(folderFilesSplitted[[f]]) - 4] != folderFilesSplitted[[f+1]] [length(folderFilesSplitted[[f+1]]) - 4]) && (folderFilesSplitted[[f]] [length(folderFilesSplitted[[f]]) - 5] == folderFilesSplitted[[f+1]] [length(folderFilesSplitted[[f+1]]) - 5])) {
		folderFilesSplitted[[f+1]] [length(folderFilesSplitted[[f+1]]) - 4] <- folderFilesSplitted[[f]] [length(folderFilesSplitted[[f]]) - 4]
	}
}
# Recreate the sample Vector
sampleVector <- character(length=length(folderFilesSplitted))
for (f in 1:length(folderFilesSplitted)) {
	sampleVector[f] <- folderFilesSplitted[[f]] [length(folderFilesSplitted[[f]]) -4]
}
spectraReplicatesAveraged <- averageMassSpectra (spectra, labels=sampleVector)
return (spectraReplicatesAveraged)
}





########################################################################





######################################## LIBRARY CREATION (for Biotyper)
libraryCreation <- function (filepathLibrary, classList=list(), classGrouping=TRUE, massRange=c(3000,15000), replicatesAverage=FALSE, patientsAverage=FALSE, SNR=5, mostIntensePeaks=TRUE, signalsToTake=20, standardDeviation=FALSE, coeffOfVar=FALSE, lowIntensityPeakRemoval=FALSE, intensityThresholdPercent=0.1, tofMode="linear", fileFormat="brukerflex") {
if (tofMode == "linear") {
	tolerancePPM <- 2000
}
if (tofMode == "reflectron" || tofMode == "reflector") {
	tolerancePPM <- 200
}
if (fileFormat == "brukerflex") {
	### Load the spectra
	spectra <- importBrukerFlex (filepathLibrary)
}
if (fileFormat == "imzml" | fileFormat == "imzML") {
	### Load the spectra
	spectra <- importImzMl (filepathLibrary)
}
### Truncation
spectra <- trim (spectra, range = massRange)
### Preprocessing
spectra <- preprocessSpectra (spectra, tofMode=tofMode)
### Average the replicates
if (replicatesAverage == TRUE) {
	spectra <- averageReplicatesByFolder (spectra, filepathLibrary, fileFormat=fileFormat)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
if (patientsAverage == TRUE) {
	spectra <- groupSpectra (spectra, spectraPerPatient=1, fileFormat=fileFormat, tofMode=tofMode)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
### Peak picking on the individual spectra
if (classGrouping == TRUE) {
	### Replace the sample name with the class name
	spectra <- replaceClassName (spectra, classList=classList, fileFormat=fileFormat)
	### Spectra grouping (class)
	spectra <- groupSpectraClass (spectra, classList, fileFormat=fileFormat)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
if (classGrouping == FALSE) {
	spectra <- spectra
}
##########################
### Peak picking
if (mostIntensePeaks == TRUE) {
	peaksLibrary <- mostIntenseSignals (spectra, signalsToTake=signalsToTake)
}
if (mostIntensePeaks == FALSE) {
	peaksLibrary <- detectPeaks (spectra, method="MAD", SNR=SNR)
}
###########################
if (standardDeviation == TRUE && coeffOfVar == FALSE) {
	# For each peaklist in the library, compute the SD
	for (j in 1:length(peaksLibrary)) {
		peaksLibrary[[j]] <- replaceSNRwithSDinPeaklist (peaksLibrary[[j]], peaks, tolerancePPM=tolerancePPM, fileFormat=fileFormat)
	}
}
if (standardDeviation == FALSE && coeffOfVar == TRUE) {
	# For each peaklist in the library, compute the CV
	for (j in 1:length(peaksLibrary)) {
		peaksLibrary[[j]] <- replaceSNRwithCVinPeaklist (peaksLibrary[[j]], peaks, tolerancePPM=tolerancePPM, fileFormat=fileFormat)
	}
}
if (standardDeviation == TRUE && coeffOfVar == TRUE) {
	# For each peaklist in the library, compute the SD
	for (j in 1:length(peaksLibrary)) {
		peaksLibrary[[j]] <- replaceSNRwithSDinPeaklist (peaksLibrary[[j]], peaks, tolerancePPM=tolerancePPM, fileFormat=fileFormat)
	}
}
############
if (lowIntensityPeakRemoval == TRUE) {
	peaksLibrary <- removeLowIntensityPeaks (peaksLibrary, intensityThresholdPercent=intensityThresholdPercent)
}
####
libraryList <- list (spectra = spectra, peaks = peaksLibrary)
return (libraryList)
}





########################################################################





######################## PEAK ALIGNMENT
alignAndFilterPeaks <- function (peaks, tolerancePPM=2000, peaksFiltering=TRUE, frequencyThreshold=0.25, lowIntensityPeaksRemoval=FALSE, intensityThresholdPercent=0.1) {
peaksAligned <- binPeaks (peaks, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksAligned <- binPeaks (peaks, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksAligned <- binPeaks (peaks, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksAligned <- binPeaks (peaks, method="relaxed", tolerance=(tolerancePPM/10^6))
peaksAligned <- binPeaks (peaks, method="relaxed", tolerance=(tolerancePPM/10^6))
if (peaksFiltering == TRUE) {
	peaksAligned <- filterPeaks (peaksAligned, minFrequency=frequencyThreshold)
}
if (lowIntensityPeaksRemoval == TRUE) {
	peaksAligned <- removeLowIntensityPeaks (peaksAligned, intensityThresholdPercent=intensityThresholdPercent)
}
return (peaksAligned)
}





########################################################################





############################ PEAK FILTERING
discardPeaks <- function (peaks, listOfPeaks=list(), massType="monoisotopic", clusterSize=4, tofMode="reflectron") {
if (massType == "monoisotopic") {
	# For each peak in the peaklist of interest
	for (p in 1:length(listOfPeaks)) {
		# Add the peaks of the cluster
		for (i in 1:(clusterSize-1)) {
			listOfPeaks <- append (listOfPeaks, listOfPeaks[p]+i)
		}
	}
}
if (massType == "average") {
	# For each peak in the peaklist of interest
	for (p in 1:length(listOfPeaks)) {
		# Add the peaks of the cluster
		for (i in 1:(clusterSize/2)) {
			listOfPeaks <- append (listOfPeaks, listOfPeaks[p]-i)
			listOfPeaks <- append (listOfPeaks, listOfPeaks[p]+i)
		}
	}
}
if (tofMode == "linear") {
	tolerancePPM <- 2000
}
if (tofMode == "reflectron" || tofMode == "reflector") {
	tolerancePPM <- 1000
}
###################################################### ERASING FUNCTION
peaksErasingFunction <- function (peaks, listOfPeaks) {
	# Create a dataframe with mass and intensity
	peaksDf <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Select the peaks of interest and remove them from the peaklist
	for (p in 1: length (listOfPeaks)) {
		for (s in 1:length(peaksDf$mass)) {
			if (abs(peaksDf$mass[s] - listOfPeaks[p])*10^6/listOfPeaks[p] <= tolerancePPM) {
				peaksDf$mass[s] <- 0
			}
		}
	}
	# Keep only the signals that are not zero
	peaksDf <- peaksDf [peaksDf$mass != 0, ]
	# Put these signals back into the peaklist
	peaks@mass <- peaksDf$mass
	peaks@intensity <- peaksDf$intensity
	peaks@snr <- peaksDf$snr
	return (peaks)
}
#######################################################################
if (length(listOfPeaks) == 0) {
	peaks <- peaks
}
if (length(listOfPeaks) > 0) {
	cl <- makeCluster (CPUcoreNumber)
	peaks <- parLapply (cl, peaks, fun=function(peaks) peaksErasingFunction (peaks, listOfPeaks))
	stopCluster (cl)
}
return (peaks)
}





########################################################################





########################## FILTER SPECTRA
filterSpectraNoise <- function (spectra, SNRfilter=20, percentageAboveThreshold=10) {
cl <- makeCluster (CPUcoreNumber)
clusterEvalQ (cl, {library(MALDIquant)})
################################################## NOISY SPECTRA REMOVAL
noisySpectraRemovalFunction <- function (spectra, SNRfilter, percentageAboveThreshold) {
	# Peak picking
	peaks <- detectPeaks (spectra, method="MAD", SNR=3)
	# Keep only the spectra with a certain percentage of signals above the threshold
	peaksNew <- detectPeaks (spectra, method="MAD", SNR=SNRfilter)
	if ((length(peaksNew@mass) / length(peaks@mass))*100 < percentageAboveThreshold) {
		spectra <- NULL
	}
	return (spectra)
}
########################################################################
spectra <- parLapply (cl, spectra, fun=function (spectra) noisySpectraRemovalFunction (spectra, SNRfilter, percentageAboveThreshold))
stopCluster (cl)
### Keep only the elements that are different from NULL
spectra <- spectra [!sapply(spectra, is.null)]
return (spectra)
}





########################################################################





######################### REPLACE THE SNR WITH THE STDEV IN THE PEAKS
### Compute the standard deviation of each peak of a average spectrum peaklist,
# by replacing the existing SNR slot with the SD (each peak of the average peaklist
# is searched across the dataset
replaceSNRwithSDinPeaklist <- function (peaksAvg, peaks, tolerancePPM=2000, fileFormat="imzml") {
	# Scroll the peaks of the average peaklist
	for (p in 1:length(peaksAvg@mass)) {
		# Create an empty vector where to allocate the intensity of this peak into the dataset
		intensityVector <- numeric()
		# Pass every peaklist of the dataset (one for each spectrum)
		for (s in 1:length(peaks)) {
			# Do the operation only for the peaklists that matches the AVG
			if (fileFormat == "brukerflex") {
				if (length(grep(peaksAvg@metaData$sampleName, peaks[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaksAvg@mass[p])*10^6/peaksAvg@mass[p]) <= tolerancePPM) {
							# Record its intensity into a vector
							intensityVector <- append (intensityVector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
			if (fileFormat == "imzml" | fileFormat == "imzML") {
				if (length(grep(peaksAvg@metaData$file, peaks[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaksAvg@mass[p])*10^6/peaksAvg@mass[p]) <= tolerancePPM) {
							# Record its intensity into a vector
							intensityVector <- append (intensityVector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
		}
		# Calculate the standard deviation of the peak (whose intensities are in the vector)
		stDev <- sd (intensityVector)
		# Replace the snr slot with the stDev (for the peak)
		peaksAvg@snr[p] <- stDev
	}
return (peaksAvg)
}




################################################################################





######################### REPLACE THE SNR WITH THE CV IN THE PEAKS
### Compute the coefficient of variation of each peak of a average spectrum peaklist,
# by replacing the existing SNR slot with the CV (each peak of the average peaklist
# is searched across the dataset
replaceSNRwithCVinPeaklist <- function (peaksAvg, peaks, tolerancePPM=2000, fileFormat="imzml") {
	# Scroll the peaks of the average peaklist
	for (p in 1:length(peaksAvg@mass)) {
		# Create an empty vector where to allocate the intensity of this peak into the dataset
		intensityVector <- numeric()
		# Pass every peaklist of the dataset (one for each spectrum)
		for (s in 1:length(peaks)) {
			# Do the operation only for the peaklists that matches the AVG
			if (fileFormat == "brukerflex") {
				if (length(grep(peaksAvg@metaData$sampleName, peaks[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaksAvg@mass[p])*10^6/peaksAvg@mass[p]) <= tolerancePPM) {
							# Record its intensity into a vector
							intensityVector <- append (intensityVector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
			if (fileFormat == "imzml" | fileFormat == "imzML") {
				if (length(grep(peaksAvg@metaData$file, peaks[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaksAvg@mass[p])*10^6/peaksAvg@mass[p]) <= tolerancePPM) {
							# Record its intensity into a vector
							intensityVector <- append (intensityVector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
		}
		# Calculate the standard deviation and the mean of the peak (whose intensities are in the vector)
		stDev <- sd (intensityVector)
		meanIntensity <- mean (intensityVector)
		# Calculate the coefficient of variation
		coeffVarInt <- (stDev / meanIntensity) * 100
		# Replace the snr slot with the CV (for the peak)
		peaksAvg@snr[p] <- coeffVarInt
	}
return (peaksAvg)
}





################################################################################





######## RICHEST PEAKLIST AS REFERENCE FOR THE SPECTRA ALIGNMENT/WARPING
richestPeaklistForAlignment <- function (spectra, SNR=5, tolerancePPM=2000) {
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
# Define the indexes
index <- 0
maxLength <- 0
# Scroll the peaklist list
for (i in 1:length(peaks)) {
# Determine the longest peaklist and record its position
	if (length(peaks[[i]]@mass) >= maxLength) {
	maxLength <- length(peaks[[i]]@mass)
	index <- i
	}
}
# Create the vector containing the richest peaklist
referencePeaklist <- peaks[[index]]
# SPECTRA WARPING / ALIGNMENT
spectra <- alignSpectra (spectra, reference=referencePeaklist, tolerance=(tolerancePPM/10^6), noiseMethod="MAD", SNR=SNR, warpingMethod="quadratic")
return (spectra)
}





################################################################################





############################# MATRIX ZERO FILLING
# Replace the NA values with a zero
matrixZeroFilling <- function (inputMatrix) {
cl <- makeCluster (CPUcoreNumber)
######################################################## NA REPLACEMENT FUNCTION
naReplacementFunction <- function (x) {
	if (is.na(x) == TRUE) {
		x <- 0
	}
return (x)
}
################################################################################
inputMatrix <- parApply (cl, inputMatrix, c(1,2), FUN=function(inputMatrix) naReplacementFunction (inputMatrix))
stopCluster(cl)
return (inputMatrix)
}





###############################################################################





######################################### ADD THE THY CLASS TO THE MATRIX (THYROID)
matrixAddThy <- function (signalMatrix, peaks, fileFormat="imzml") {
spectraNumber <- length (peaks)
# Create the empty vector
fileVector <- character()
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(peaks)) {
	if (fileFormat == "imzml" | fileFormat == "imzML") {
		fileVector <- append (fileVector, peaks[[i]]@metaData$file)
	}
	if (fileFormat == "brukerflex") {
		fileVector <- append (fileVector, peaks[[i]]@metaData$sampleName)
	}
}
# Create the thy vector, symmetrical to the file vector
thyVector <- fileVector
# Find the "THY" in the file name and store its value
for (t in 1:length(thyVector)) {
	############## THY
	# THY1
	if (length(grep("THY1", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 1
	}
	# THY2
	if (length(grep("THY2", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 2
	}
	# THY3
	if (length(grep("THY3", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 3
	}
	# THY4
	if (length(grep("THY4", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 4
	}
	# THY5
	if (length(grep("THY5", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 5
	}
	############## TIR
	# TIR1
	if (length(grep("TIR1", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 1
	}
	# TIR2
	if (length(grep("TIR2", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 2
	}
	# TIR3
	if (length(grep("TIR3", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 3
	}
	# TIR4
	if (length(grep("TIR4", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 4
	}
	# TIR5
	if (length(grep("TIR5", thyVector[t], ignore.case=TRUE)) == 1) {
		thyVector[t] <- 5
	}
}
# Create the sample matrix column and append it to the global matrix
thyColumn <- matrix (0, ncol=1, nrow=spectraNumber)
colnames (thyColumn) <- "THY"
# Fill in the matrix thy column with the thyVector and attach it to the matrix
thyColumn [,1] <- cbind (thyVector)
signalMatrix <- cbind (signalMatrix, thyColumn)
#
return (signalMatrix)
}





################################################################################





############################## PEAK STATISTICS (on processed Spectra)
peakStatistics <- function (spectra, SNR=3, classList=list(),
	tofMode="linear", fileFormat="imzml", tolerancePPM=2000, peaksFiltering=TRUE, frequencyThreshold=0.25, removeOutliers=TRUE) {
# Determine the number of classes
if (length(classList) == 0 | length(classList) == 1) {
numberOfClasses <- 1
}
if (length(classList) > 1) {
	numberOfClasses <- length(classList)
}
# Detect and Align Peaks
if (tofMode == "linear") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR, halfWindowSize=20)
}
if (tofMode == "reflector" | tofMode == "reflectron") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR, halfWindowSize=20)
}
peaks <- alignAndFilterPeaks (peaks, tolerancePPM=tolerancePPM, peaksFiltering=peaksFiltering, frequencyThreshold=frequencyThreshold)
# Generate the matrix (and convert it into a data frame)
signalMatrix <- intensityMatrix (peaks, spectra)
# Peak vector
#peakVector <- as.numeric (names(signalMatrix))
############################################################## ONE CLASS
if (numberOfClasses == 1) {
	# Fix the signalMatrix (Add the sample column)
	signalMatrix <- matrixAddClassAndSample (signalMatrix, peaks=peaks,
		fileFormat=fileFormat, sampleOutput=TRUE, classOutput=FALSE)
	# Output matrix
	peakStatMatrix <- matrix (0, nrow=(ncol(signalMatrix)-1), ncol=9)
	rownames (peakStatMatrix) <- as.numeric (colnames(signalMatrix)[1:(ncol(signalMatrix)-1)])
	colnames (peakStatMatrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter", "Sample")
	# For each peak
	for (p in 1:(ncol(signalMatrix)-1)) {
		intensityVector <- as.numeric (signalMatrix[,p])
		if (removeOutliers == TRUE) {
			intensityVector <- outliersRemoval (intensityVector)
			intensityVector <- intensityVector$vector
		}
		# Calculate the statistical parameters on the intensity values in the vector
		# Normality
		if (length(intensityVector) >= 3 & length(intensityVector) <= 5000) {
		shapiroTest <- shapiro.test (intensityVector)
			if (shapiroTest$p.value < 0.05) {
			distributionType <- "Non-normal"
			}
			if (shapiroTest$p.value >= 0.05) {
			distributionType <- "Normal"
			}
		}
		if (length(intensityVector) < 3) {
		distributionType <- "Not determinable, number of samples too low"
		}
		if (length(intensityVector) > 5000) {
		distributionType <- "Number of samples too high, assume it is normal"
		}
		# Other parameters
		spectraNames <- levels (as.factor(signalMatrix[,ncol(signalMatrix)]))
		stDevIntensity <- sd (intensityVector)
		summaryIntensityVector <- summary (intensityVector)
		meanIntensity <- summaryIntensityVector [4]
		coeffVariation <- (stDevIntensity / meanIntensity) *100
		medianIntensity <- summaryIntensityVector [3]
		firstQuartile <- summaryIntensityVector [2]
		thirdQuartile <- summaryIntensityVector [5]
		interQuartileRange <- thirdQuartile - firstQuartile
		spectraCounter <- length (intensityVector)
		# Fill the matrix with the values
		peakStatMatrix [p,1] <- distributionType
		peakStatMatrix [p,2] <- meanIntensity
		peakStatMatrix [p,3] <- stDevIntensity
		peakStatMatrix [p,4] <- coeffVariation
		peakStatMatrix [p,5] <- medianIntensity
		peakStatMatrix [p,6] <- interQuartileRange
		peakStatMatrix [p,7] <- paste ("1st quartile", firstQuartile, "; 3rd quartile", thirdQuartile)
		peakStatMatrix [p,8] <- spectraCounter
		peakStatMatrix [p,9] <- spectraNames
	}
}
############################################################ TWO OR MORE CLASSES
# Every variable now is a list, each element of which corresponds to a certain value from a class
# So every variable is a list with the same length of the class list (each element of the list
# is referred to a class
if (numberOfClasses > 1) {
	# Fix the signalMatrix (Add the sample column)
	signalMatrix <- matrixAddClassAndSample (signalMatrix, peaks=peaks, classList=classList,
		fileFormat=fileFormat, sampleOutput=TRUE, classOutput=TRUE)
	# Output matrix
	peakStatMatrix <- matrix (0, nrow=(ncol(signalMatrix)-2), ncol=14)
	rownames (peakStatMatrix) <- as.numeric (names(signalMatrix)[1:(ncol(signalMatrix)-2)])
	colnames (peakStatMatrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %",
		"Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)",
		"t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
	# For each peak
	for (p in 1:(ncol(signalMatrix)-2)) {
		# Put the intensity of that peak into one vector per class (in a global list)
		intensityVector <- list()
		# Scroll the peaklists and Add the peak intensity to a vector (one for each class)
		for (l in 1:length(classList)) {
			# Allocate in the intensity vector the rows for that peak belonging to the certain class
			intensityVector [[l]] <- as.numeric (signalMatrix [signalMatrix[,ncol(signalMatrix)] == classList[l],p])
		}
		if (removeOutliers == TRUE) {
			for (i in 1:length(intensityVector)) {
				intensityVector[[l]] <- outliersRemoval (intensityVector[[l]])
				intensityVector[[l]] <- intensityVector[[l]]$vector
			}
		}
		######################## STATISTICAL PARAMETERS
		############################################### Normality for each class
		shapiroTest <- list()
		distributionType <- list()
		for (l in 1:length(classList)) {
			if (length(intensityVector[[l]]) >= 3 && length(intensityVector[[l]]) <=5000) {
				shapiroTest[[l]] <- shapiro.test (intensityVector[[l]])
				if (shapiroTest[[l]]$p.value < 0.05) {
				distributionType[[l]] <- "Non-normal"
				}
				if (shapiroTest[[l]]$p.value >= 0.05) {
				distributionType[[l]] <- "Normal"
				}
			}
			if (length(intensityVector[[l]]) < 3) {
			distributionType[[l]] <- "Not determinable, number of samples too low"
			}
			if (length(intensityVector) > 5000) {
			distributionType[[l]] <- "Number of samples too high, assume it is normal"
			}
		}
		##################################################### Homoscedasticity
		if (length(classList) == 2) {
			varianceTestPar <- var.test (intensityVector[[1]], intensityVector[[2]])
		}
		if (length(classList) >= 2) {
			varianceTestNonPar <- bartlett.test (as.numeric(signalMatrix[,p]), g=as.factor(signalMatrix[,ncol(signalMatrix)]))
		}
		########################################### Other parameters (per class)
		stDevIntensity <- list()
		summaryIntensityVector <- list()
		meanIntensity <- list()
		coeffVariation <- list()
		medianIntensity <- list()
		firstQuartile <- list()
		thirdQuartile <- list()
		interQuartileRange <- list()
		spectraCounter <- list()
		variance <- list()
		for (l in 1:length(classList)) {
			stDevIntensity[[l]] <- sd (intensityVector[[l]])
			summaryIntensityVector [[l]] <- summary (intensityVector[[l]])
			meanIntensity[[l]] <- summaryIntensityVector[[l]] [4]
			coeffVariation[[l]] <- (stDevIntensity[[l]] / meanIntensity[[l]]) *100
			medianIntensity[[l]] <- summaryIntensityVector[[l]] [3]
			firstQuartile[[l]] <- summaryIntensityVector[[l]] [2]
			thirdQuartile[[l]] <- summaryIntensityVector[[l]] [5]
			interQuartileRange[[l]] <- thirdQuartile[[l]] - firstQuartile[[l]]
			spectraCounter[[l]] <- length (intensityVector[[l]])
			variance[[l]] <- var (intensityVector[[l]])
		}
		############################################# Parameters between classes
		# T-test
		if (length(classList) == 2) {
			tTest <- t.test (intensityVector[[1]], intensityVector[[2]])
		}
		# ANOVA TEST
		if (length(classList) >= 2) {
		anovaTest <- aov (signalMatrix[,p] ~ signalMatrix[,ncol(signalMatrix)])
		}
		# WILCOXON - MANN-WHITNEY TEST
		if (length(classList) == 2) {
			WilcoxonTest <- wilcox.test (intensityVector[[1]], intensityVector[[2]])
		}
		# KRUSKAL-WALLIS TEST
		if (length(classList) >= 2) {
			KruskalWallisTest <- kruskal.test (signalMatrix[,p], g=as.factor(signalMatrix[,ncol(signalMatrix)]))
		}
		######################################## Fill the matrix with the values
		# Distribution Type
		distributionTypeName <- character()
		for (l in length(classList)) {
			distributionTypeName <- paste (distributionTypeName, " ", distributionType[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,1] <- paste (distributionTypeName)
		# Mean
		meanIntensityName <- character()
		for (l in length(classList)) {
			meanIntensityName <- paste (meanIntensityName, " ", meanIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,2] <- meanIntensityName
		# Standard Deviation
		StDevName <- character()
		for (l in length(classList)) {
			StDevName <- paste (StDevName, " ", stDevIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,3] <- StDevName
		# Coefficient of Variation
		coeffVariationName <- character()
		for (l in length(classList)) {
			coeffVariationName <- paste (coeffVariationName, " ", coeffVariation[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,4] <- coeffVariationName
		# Median
		medianIntensityName <- character()
		for (l in length(classList)) {
			medianIntensityName <- paste (medianIntensityName, " ", medianIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,5] <- medianIntensityName
		# Interquartile Range (IQR)
		interQuartileRangeName <- character()
		for (l in length(classList)) {
			interQuartileRangeName <- paste (interQuartileRangeName, " ", interQuartileRange[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,6] <- interQuartileRangeName
		# Spectra counter
		spectraCounterName <- character()
		for (l in length(classList)) {
			spectraCounterName <- paste (interQuartileRangeName, " ", spectraCounter[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,7] <- spectraCounterName
		# Class
		className <- character()
		for (l in length(classList)) {
			className <- paste (className, " ", classList[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,8] <- className
		# Homoscedasticity (Parametric)
		if (varianceTestPar$p.value < 0.05) {
		homoscedasticityP <- paste ("Non homoscedastic data", "(p-value:", varianceTestPar$p.value, ")")
		}
		if (varianceTestPar$p.value >= 0.05) {
		homoscedasticityP <- paste ("Homoscedastic data", "(p-value:", varianceTestPar$p.value, ")")
		}
		if (varianceTestNonPar$p.value < 0.05) {
		homoscedasticityNP <- paste ("Non homoscedastic data", "(p-value:", varianceTestNonPar$p.value, ")")
		}
		if (varianceTestNonPar$p.value >= 0.05) {
		homoscedasticityNP <- paste ("Homoscedastic data", "(p-value:", varianceTestNonPar$p.value, ")")
		}
		peakStatMatrix [p,9] <- homoscedasticityP
		peakStatMatrix [p,10] <- homoscedasticityNP
		# t-Test
		peakStatMatrix [p,11] <- tTest$p.value
		# ANOVA
		peakStatMatrix [p,12] <- summary(anovaTest)[[1]]$"Pr(>F)"[1]
		# Wilcoxon / Mann-Whitney test
		peakStatMatrix [p,13] <- WilcoxonTest$p.value
		# Kruskal-Wallis test
		peakStatMatrix [p,14] <- KruskalWallisTest$p.value
	}
}
return (peakStatMatrix)
}



################################################################################





########################## REMOVE LOW INTENSITY PEAKS
removeLowIntensityPeaks <- function (peaks, intensityThresholdPercent=0.1) {
cl <- makeCluster (CPUcoreNumber)
########################################### INTENSITY FILTERING FUNCTION
intensityFilteringFunction <- function (peaks, intensityThresholdPercent) {
	# Filter out the peaks whose intensity is below a certain threshold
	for (p in 1:length(peaks@mass)) {
		if ((abs(peaks@intensity[p] - max(peaks@intensity,na.rm=TRUE))*100/max(peaks@intensity,na.rm=TRUE)) < intensityThresholdPercent) {
			peaks@intensity[p] <- NA
			peaks@mass[p] <- NA
		}
	}
	# Discard the NA values
	peaks@mass <- peaks@mass [!is.na (peaks@mass)]
	peaks@intensity <- peaks@intensity [!is.na (peaks@intensity)]
	#
	return (peaks)
}
########################################################################
peaksFiltered <- parLapply (cl, peaks, fun = function (peaks) intensityFilteringFunction (peaks, intensityThresholdPercent=intensityThresholdPercent))
return (peaksFiltered)
stopCluster (cl)
}





################################################################################





###################### SPECTRA PURIFICATION BASED ON THE TIC (Before processing)
spectraTICpurification <- function (spectra, tofMode="reflectron", absoluteTICthreshold=0) {
# If not specified, set it based on the TOF mode
if (absoluteTICthreshold == 0) {
	if (tofMode == "linear") {
		absoluteTICthreshold <- 2000
	}
	if (tofMode == "reflectron" | tofMode == "reflector") {
		absoluteTICthreshold <- 200
	}
}	else {absoluteTICthreshold <- absoluteTICthreshold}
# Before preprocessing (and thus normalisation), evaluate the TIC
ticVector <- numeric()
for (s in 1:length(spectra)) {
	# Calculate the spectrum TIC
	ticSpectrum <- sum (spectra[[s]]@intensity)
	# Add this value to the global TIC vector
	ticVector <- append (ticVector, ticSpectrum)
}
# Remove the spectra with an absolute TIC value of less than ...
lowTICposition <- which (ticVector <= absoluteTICthreshold)
if (length(lowTICposition) > 0) {
	# Remove the bad spectra both from the spectra list and from the ticVector
	spectra <- spectra [-lowTICposition]
	ticVector <- ticVector [-lowTICposition]
}	else {spectra <- spectra}
###### Outliers detection
# Output the summary of the vector (quartiles, mean, median, max and min)
summaryTicVector <- summary (ticVector)
# Calculate the interquartile range
interQuartileRange <- summaryTicVector[5] - summaryTicVector[2]
# Calculate the fences, beyond which the spectrum is an outlier
IQRfences <- c((summaryTicVector[2] - 1.5*interQuartileRange), (summaryTicVector[5] + 1.5*interQuartileRange))
# Find the outliers based on the fences condition
outliersPosition <- which (ticVector < IQRfences[1] | ticVector > IQRfences[2])
if (length(outliersPosition) > 0) {
	# Remove the correspondent spectra from the dataset
	ticVector <- ticVector [-outliersPosition]
	spectra <- spectra [-outliersPosition]
}	else {spectra <- spectra}
#
return (spectra)
}
################################################################################





################################# REMOVE UNDIGESTED SPECTRA
removeUndigestedSpectra <- function (spectra, massRangeOfInterest=c(0,2500), percentageOfSpectrumInTheRange=75) {
# Make the cluster (one for each core/thread)
cl <- makeCluster (CPUcoreNumber)
clusterEvalQ (cl, {library(MALDIquant)})
################################################# UNDIGESTED FILTERING FUNCTION
undigestedFilteringFunction <- function (spectra, massRangeOfInterest, percentageOfSpectrumInTheRange) {
	## Peak picking
	peaks <- detectPeaks (spectra, method="MAD", SNR=3, halfWindowSize=20)
	# Define the peaks of interest based on the given interval
	peaksOfInterest <- peaks@mass [peaks@mass > massRangeOfInterest[1] & peaks@mass < massRangeOfInterest[2]]
	# Calculate the percentage of the spectrum that is within the given range
	percentageInRange <- length(peaksOfInterest)/length(peaks@mass) * 100
	# The spectrum is to be discarded if not of interest
	if (percentageInRange < percentageOfSpectrumInTheRange) {
		spectra <- NULL
	}
	return (spectra)
}
########################################################################
spectraDigested <- parLapply (cl, spectra, fun=function (spectra) undigestedFilteringFunction (spectra, massRangeOfInterest, percentageOfSpectrumInTheRange))
# Close the processes
stopCluster(cl)
### Keep only the elements that are different from NULL
spectraDigested <- spectraDigested [!sapply(spectraDigested, is.null)]
#
return (spectraDigested)
}





################################################################################





################# GENERATE REPRESENTATIVE SPECTRA, BASED UPON SPECTRA SIMILARITY (ONE IMZML)
generateRepresentativeSpectraHCA <- function (spectra, spectraPerPatient=10, tofMode="linear", algorithm="hierarchicalClustering") {
# Detect and align peaks
if (tofMode=="linear") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=3, halfWindowSize=20)
	peaks <- alignAndFilterPeaks (peaks, tolerancePPM=2000, peaksFiltering=FALSE, frequencyThreshold=0.25)
}
if (tofMode=="reflectron" | tofMode=="reflector") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=3, halfWindowSize=20)
	peaks <- alignAndFilterPeaks (peaks, tolerancePPM=200, peaksFiltering=FALSE, frequencyThreshold=0.25)
}
# Generate the peaklist matrix
peaklist <- intensityMatrix (peaks, spectra)
############################ HCA
if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
	# Compute the distance matrix
	distanceMatrix <- dist (peaklist, method="euclidean")
	hca <- hclust (distanceMatrix)
	# Associate to each row/spectrum the subgroup of the tree to which it belongs
	# k must be between 1 and 56
	if (spectraPerPatient > 56) {
		spectraPerPatient <- 56
	}
	hcaGroups <- cutree (hca, k=spectraPerPatient)
	# Generate the final list of spectra
	spectraGrouped <- list()
	# For each subgroup to be isolated...
	for (s in 1:spectraPerPatient) {
		# Index the spectra under in the selected subgroup of the HCA
		index <- which (hcaGroups == s)
		spectraHCA <- spectra[index]
		# Average the spectra in this HCA subgroup
		spectraAVG <- averageMassSpectra (spectraHCA, method="mean")
		# Add the average to the final list
		spectraGrouped <- append (spectraGrouped, spectraAVG)
	}
}
######################### K-MEANS
if (algorithm == "kMeans" | algorithm == "kmeans" | algorithm == "k-Means") {
	# Compute the k-Means clustering
	kMeans <- kmeans (peaklist, centers=spectraPerPatient)
	# Associate to each row/spectrum the subgroup/cluster to which it belongs
	kMeansGroups <- kMeans$cluster
	# Generate the final list of spectra
	spectraGrouped <- list()
	# For each subgroup to be isolated...
	for (s in 1:spectraPerPatient) {
		# Index the spectra under in the selected subgroup of the HCA
		index <- which (kMeansGroups == s)
		spectraKMeans <- spectra[index]
		# Average the spectra in this HCA subgroup
		spectraAVG <- averageMassSpectra (spectraKMeans, method="mean")
		# Add the average to the final list
		spectraGrouped <- append (spectraGrouped, spectraAVG)
	}
}
return (spectraGrouped)
}





################################################################################





############################################## DEISOTOPING
deIsotopePeaks <- function (peaks) {
##################################################### Deisotoping Function
deisotopingFunction <- function (peaks) {
	# Deisotoped Ions
	deisotopedIonsMass <- numeric()
	deisotopedIonsIntensity <- numeric()
	# Take each ion and evaluate the consecutive ions
	isotopeBinMass <- numeric ()
	isotopeBinIntensity <- numeric ()
	# Scroll the ion list
	for (i in 1:(length(peaks@mass)-1)) {
		# Create the isotope bins, if the consecutive peaks belong to a isotope cluster
		if (abs(peaks@mass[i+1] - peaks@mass[i]) < 1.5) {
			isotopeBinMass <- append (isotopeBinMass, peaks@mass[i])
			isotopeBinIntensity <- append (isotopeBinIntensity, peaks@intensity[i])
		}
		# When a major distance is found... We are out of the cluster
		if (abs(peaks@mass[i+1] - peaks@mass[i]) >= 1.5) {
			# Add this mass to the isotope bins (it will be the last belonging to a this isotope cluster
			isotopeBinMass <- append (isotopeBinMass, peaks@mass[i])
			isotopeBinIntensity <- append (isotopeBinIntensity, peaks@intensity[i])
			# Store the isotope of the isotope Bins with the maximum intensity created so far in a new vector
			mostIntensePeakMass <- isotopeBinMass [which(isotopeBinIntensity == max(isotopeBinIntensity))]
			mostIntensePeakIntensity <- max(isotopeBinIntensity)
			# Prevent adding more masses and only one intensity, kep the balance
			if (length(mostIntensePeakMass) == 1) {
				deisotopedIonsMass <- append (deisotopedIonsMass, mostIntensePeakMass)
				deisotopedIonsIntensity <- append (deisotopedIonsIntensity, mostIntensePeakIntensity)
			}
			if (length(mostIntensePeakMass) > 1) {
				deisotopedIonsMass <- append (deisotopedIonsMass, mostIntensePeakMass[1])
				deisotopedIonsIntensity <- append (deisotopedIonsIntensity, mostIntensePeakIntensity[1])
			}
			# Empty the isotope Bins in order for it to become available for the future cycle
			isotopeBinMass <- numeric()
			isotopeBinIntensity <- numeric()

		}
	}
	# Fill the peaklist back in
	peaks@mass <- deisotopedIonsMass
	peaks@intensity <- deisotopedIonsIntensity
	return (peaks)
}
################################################################################
cl <- makeCluster (CPUcoreNumber)
peaks <- parLapply (cl, peaks, fun = function(peaks) deisotopingFunction(peaks))
stopCluster (cl)
return (peaks)
}





################################################################################





############################################ OUTLIERS REMOVAL
outliersRemoval <- function (vector, replaceWith="") {
summaryVector <- summary (vector)
# Calculate the interquartile range
interQuartileRange <- summaryVector[5] - summaryVector[2]
# Calculate the fences, beyond which the spectrum is an outlier
IQRfences <- c((summaryVector[2] - 1.5*interQuartileRange), (summaryVector[5] + 1.5*interQuartileRange))
# Find the outliers based on the fences condition
outliersPosition <- which (vector < IQRfences[1] | vector > IQRfences[2])
# If the outliers have to be discarded...
if (replaceWith == "") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vector <- vector [-outliersPosition]
	}
}
if (replaceWith == "mean") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vectorNoOutliers <- vector [-outliersPosition]
		# Replace the outliers with the vector mean (without outliers)
		vector [outliersPosition] <- mean (vectorNoOutliers)
	}
}
if (replaceWith == "median") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vectorNoOutliers <- vector [-outliersPosition]
		# Replace the outliers with the vector median
		vector [outliersPosition] <- median (vectorNoOutliers)
	}
}
return (list (vector=vector, outliersPosition=outliersPosition))
}





########################################################################




beep (5)





