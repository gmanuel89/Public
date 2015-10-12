###INSTALL THE REQUIRED PACKAGES
required_packages <- c ("parallel", "MALDIquant", "MALDIquantForeign", "caret", "tcltk", "beepr")
installed_packages <- installed.packages () [,1]
missing_packages <- character()
for (p in 1:length(required_packages)) {
	if ((required_packages[p] %in% installed_packages) == FALSE) {
		missing_packages <- append(missing_packages, required_packages[p])
	}
}
if (length(missing_packages) > 0) {
	install.packages (missing_packages)
}

### LOAD THE REQUIRED PACKAGES
library (parallel)
library (MALDIquant)
library (MALDIquantForeign)
library (tcltk)
library (caret)
#library (e1071)
library (beepr)
#library (gWidgets)
#library (gWidgetstcltk)


##### MEMORY LIMIT
# Set the memory (RAM) limit to 100GB (100000MB)
memory.limit (100000)



######### PARALLELISATION
# Detect the number of cores
cpu_thread_Number <- detectCores (logical=TRUE)
cpu_core_number <- cpu_thread_number/2
##############################




####################################### READ (SOURCE) SCRIPT FROM URL
source_github <- function(url, ...) {
  # Load required package (install it if not present)
  if ("RCurl" %in% installed.packages()[,1] == FALSE) {
	  install.packages ("RCurl")
  }
  require(RCurl)
  # Parse and evaluate each R script (in the global environement)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}





################################################ INSTALL REQUIRED PACKAGES
install_and_load_required_packages <- function (required_packages) {
installed_packages <- installed.packages () [,1]
missing_packages <- character()
for (p in 1:length(required_packages)) {
	if ((required_packages[p] %in% installed_packages) == FALSE) {
		missing_packages <- append(missing_packages, required_packages[p])
	}
}
if (length(missing_packages) > 0) {
	install.packages (missing_packages)
}
for (i in 1:length(required_packages)) {
	library (required_packages[i])
}
}





################################################ SPECTRA QUALITY CONTROL
spectra_quality_control <- function (spectra) {
## Empty spectra
print(paste("Are there any empty spectra?", (any(sapply(spectra, isEmpty)))))
## Same number of data points
data_points_table <- table(sapply(spectra, length))
print("Same number of data points?")
print(data_points_table)
## Same distance between data points?
print(paste("Same distance between data points?", (all(sapply(spectra, isRegular)))))
## Flat spectra
flat_spectra <- findEmptyMassObjects (spectra)
print(paste("Flat spectra in the dataset:", length(flat_spectra)))
}





######################################### REPLACE THE metaData$file WITH THE SAMPLE NAME (SPECTRA or PEAKS)
replace_sample_name <- function (spectra, folder, file_format="imzml") {
setwd(folder)
# Read the files in the folder
folder_files <- read_spectra_files(folder, file_format=file_format, full_path=FALSE)
if (file_format == "imzml" | file_format == "imzML") {
	for (s in 1:length(folder_files)) {
		for (i in 1:length(spectra)) {
			if (length(grep(folder_files[s],spectra[[i]]@metaData$file, fixed=TRUE)) != 0) {
				spectra[[i]]@metaData$file <- folder_files[s]
				spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(".imzML")))[1]
			}
		}
	}
}
#if (file_format == "imzml" | file_format == "imzML") {
#	for (i in 1:length(spectra)) {
#		spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(folder)))[2]
#		spectra[[i]]@metaData$file <- unlist(strsplit(spectra[[i]]@metaData$file, split=(".imzML")))[1]
 #   }
#}
return (spectra)
}





#########################################################################





######################################### REPLACE THE metaData$file WITH THE CLASS NAME (SPECTRA or PEAKS)
replace_class_name <- function (spectra, class_list, file_format="imzml") {
class_list <- sort(class_list)
if (file_format == "imzml" | file_format == "imzML") {
	for (w in 1:length(class_list)) {
		for (i in 1:length(spectra)) {
			if (length(grep(class_list[w],spectra[[i]]@metaData$file)) == 1) {
				spectra[[i]]@metaData$file <- class_list [w]
			}
		}
	}
}
if (file_format == "brukerflex") {
	for (w in 1:length(class_list)) {
		for (i in 1:length(spectra)) {
			if (length(grep(class_list[w],spectra[[i]]@metaData$sampleName)) == 1) {
				spectra[[i]]@metaData$sampleName <- class_list [w]
			}
		}
	}
}
return (spectra)
}





########################################################################





######################################### SPECTRA PRE-PROCESSING
# Preprocess spectra
preprocess_spectra <- function (spectra, tof_mode="linear", smoothing_strength="medium", process_in_packages_of=length(spectra)) {
if (process_in_packages_of <= 0 || preprocess_in_packages_of > length(spectra)) {
	process_in_packages_of <- length(spectra)
}
# Create the list containing the processed spectra
preprocessed_spectra <- list()
# Calculate the k_fold
k_fold <- floor(length(spectra) / process_in_packages_of)
# Determine the spectra to be randomly allocated into processing folds
index <- createFolds (y=seq(1:length(spectra)), k=k_fold)
for (i in 1:length(index)) {
	# Make the cluster (one for each core/thread)
	cl <- makeCluster(cpu_core_number)
	clusterEvalQ (cl, {library(MALDIquant)})
	# Allocate the spectra (in the package) to be processed temporarily
	spectra_temp <- spectra [index[[i]]]
	if (tof_mode == "linear") {
		## Remove flat spectra
		#spectra <- removeEmptyMassObjects (spectra)
		if (smoothing_strength == "small") {
			## Smoothing (Savitzky-Golay filter, with window size 5, 11 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=5))
		}
		if (smoothing_strength == "medium") {
			## Smoothing (Savitzky-Golay filter, with window size 10, 21 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=10))
		}
		if (smoothing_strength == "strong") {
			## Smoothing (Savitzky-Golay filter, with window size 20, 41 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=20))
		}
		if (smoothing_strength == "veryStrong") {
			## Smoothing (Savitzky-Golay filter, with window size 30, 61 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=30))
		}
		## Baseline correction
		spectra_temp <- parLapply(cl, spectra_temp, fun= function (spectra) removeBaseline (spectra, method="TopHat"))
		## Normalisation (TIC)
		spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) calibrateIntensity (spectra, method="TIC"))
	}
	if (tof_mode == "reflectron" || tof_mode == "reflector") {
		## Remove flat spectra
		#spectra <- removeEmptyMassObjects (spectra)
		if (smoothing_strength == "small") {
			## Smoothing (Savitzky-Golay filter, with window size 2, 5 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=1))
		}
		if (smoothing_strength == "medium") {
			## Smoothing (Savitzky-Golay filter, with window size 6, 13 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=3))
		}
		if (smoothing_strength == "strong") {
			## Smoothing (Savitzky-Golay filter, with window size 12, 25 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=6))
		}
		if (smoothing_strength == "veryStrong") {
			## Smoothing (Savitzky-Golay filter, with window size 18, 37 points)
			spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=9))
		}
		## Baseline correction
		spectra_temp <- parLapply(cl, spectra_temp, fun= function (spectra) removeBaseline (spectra, method="TopHat"))
		## Normalisation (TIC)
		spectra_temp <- parLapply(cl, spectra_temp, fun=function (spectra) calibrateIntensity (spectra, method="TIC"))
	}
	# Close the processes
	stopCluster(cl)
	# Append the processed spectra package to the global list
	preprocessed_spectra <- append(preprocessed_spectra, spectra_temp)
}
return (preprocessed_spectra)
}





########################################################################





######################################### AVERAGE SIGNAL NUMBER
# Average number of signals with that S/N in the peaklist set
average_signal_number_with_selected_snr <- function (spectra, snr=5) {
peaks <- detectPeaks(spectra, method="MAD", snr=snr)
average_signal_number <- numeric()
signal_number_vector <- numeric()
# For each peaklist
for (p in 1:length(peaks)) {
	# Determine the number of peaks
	signal_number <- length(peaks[[p]]@mass)
	signal_number_vector <- append(signal_number_vector, signal_number)
}
average_signal_number <- sum(signal_number_vector) / length(signal_number_vector)
return (average_signal_number)
}





#########################################################################





######################################### REPRESENTATIVE SPECTRA (one imzML per patient)
# Number of spectra with more than a defined number of signals with that snr
number_of_representative_spectra <- function (spectra, snr=15, signal_number_threshold=20) {
counter <- 0
## Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=snr)
# For each peaklist
for (p in 1:length(peaks)) {
	# Increase the counter if the number of peaks is more than the defined threshold
	if (length(peaks[[p]]@mass) >= signal_number_threshold) {
		counter <- counter + 1
	}
}
return (counter)
}





########################################################################





######################################### REPRESENTATIVE SPECTRA DATASET (one imzML per patient)
# Average number of spectra with more than a defined number of signals with that snr
average_number_of_representative_spectra <- function (filepath, snr=15, signal_number_threshold=20, tof_mode="linear",
	file_format="imzml") {
counter_vector <- numeric()
for (f in 1:length(filepath)) {
	# Import Spectra
	if (file_format == "imzml" | file_format == "imzML") {
		spectra <- importImzMl(filepath[f])
	}
	if (file_format == "brukerflex") {
		spectra <- importBrukerFlex(filepath[f])
	}
	# Pre-processing
	spectra <- preprocess_spectra(spectra, tof_mode=tof_mode)
	# Counter
	counter <- 0
	## Peak picking
	peaks <- detectPeaks(spectra, method="MAD", snr=snr)
	# For each peaklist
	for (p in 1:length(peaks)) {
		# Increase the counter if the number of peaks is more than the defined threshold
		if (length(peaks[[p]]@mass) >= signal_number_threshold) {
			counter <- counter + 1
		}
	}
	counter_vector <- append(counter_vector, counter)
}
average_counter <- mean(counter_vector)
return (counter_vector)
}





########################################################################





########################################## SPECTRA FILTERING
# Remove all the bad spectra, keep only the best ones (Number of signals with a certain snr)
filter_spectra_signal_number <- function (spectra, snr_filter=15, signal_number_threshold=15) {
	# Make the cluster (one for each core/thread)
	cl <- makeCluster(cpu_core_number)
	clusterEvalQ (cl, {library(MALDIquant)})
##################################################### FILTERING FUNCTION
spectra_filtering_function <- function (spectra) {
	peaks <- detectPeaks (spectra, method="MAD", snr=snr_filter)
	####### Select only the desired peaklists, based upon the threshold characteristics
		if (length(peaks@mass) < signal_number_threshold) {
			# Remove the bad spectrum from the list
			spectra <- NULL
		}
	return (spectra)
}
########################################################################
	purified_spectra <- parLapply(cl, spectra, fun=function(spectra) spectra_filtering_function(spectra))
	# Close the processes
	stopCluster(cl)
	### Keep only the elements that are different from NULL
	purified_spectra <- purified_spectra [!sapply(purified_spectra, is.null)]
	return (purified_spectra)
}






#########################################################################





########################################## SPECTRA GROUPING (PATIENTS) (one imzML file per patient)
# Obtain a certain number of spectra per patient -> Representative spectra
group_spectra <- function (spectra, spectra_per_patient=1, file_format="imzml", tof_mode="linear", seed=0, algorithm="random", discarded_nodes=1, balanced=TRUE, method="mean") {
### Create the file Vector
if (file_format == "imzml" | file_format == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		file_vector <- append(file_vector, spectra[[i]]@metaData$file)
	}
	patient_vector <- unique(file_vector)
}
if (file_format == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		file_vector <- append(file_vector, spectra[[i]]@metaData$sampleName)
	}
	patient_vector <- unique(file_vector)
}
######################################## Balancing
############ The spectra_per_patient get adjusted based on the shortest spectra list
if (balanced == TRUE) {
	########## All the patients should have the same number of spectra, always
	### Determine the shortest patient's spectra number
	#### Set the lowest as the row number of the first patient
	lowest_number_of_observations <- NULL
	for (p in 1:length(patient_vector)) {
		# Isolate the spectra of that patient
		spectra_patient <- list()
		for (s in 1:length(spectra)) {
			if (spectra[[s]]@metaData$file == patient_vector[[p]]) {
				spectra_patient <- append(spectra_patient, spectra[[s]])
			}
		}
		# If its lower than the reference (or not set yet)
		if (lowest_number_of_observations == NULL || length(spectra_patient) < lowest_number_of_observations) {
			lowest_number_of_observations <- length(spectra_patient)
		}
	}
	#### Rows per patient should not be more than the minimum number of rows
	if (spectra_per_patient >= lowest_number_of_observations) {
		spectra_per_patient <- lowest_number_of_observations
	}
}
# If the spectra per patient is equal to zero or one, do the normal averaging
if (spectra_per_patient == 0 || spectra_per_patient == 1) {
	if (method == "mean") {
		patient_spectra_final <- averageMassSpectra (spectra, labels=file_vector, method="mean")
	}
	if (method == "skyline") {
		patient_spectra_final <- group_spectra_skyline (spectra, file_format=file_format)
	}
}
# Run this script if there has to be two or more representative spectra per patient
if (spectra_per_patient > 1) {
	# If there is only one spectrum, use it
	if (length(spectra) == 1) {
		patient_spectra_final <- spectra
	}
	# If there are more spectra...
	if (length(spectra) > 1) {
		# Make the randomness reproducible
		if (seed != 0) {
			set.seed (seed)
		}
		############################# OUTPUTS
		# Generate the final list of spectra
		patient_spectra_final <- list()
		discarded_spectra <- list()
		discarded_spectra_average <- list()
		# Create a new list of spectra for plotting purposes (the intensities will be replaced)
		spectra_for_plotting <- list()
		# List of additional plots
		plots <- list()
		# Generate the final list of MSI images
		msi_plots <- list()
		##############################
		# For each patient (p)
		for (p in 1:length(patient_vector)) {
			# Add the single patient spectra to a list
			patient_spectra <- list()
			for (s in 1:length(spectra)) {
				if (spectra[[s]]@metaData$file == patient_vector[p]) {
					patient_spectra <- append(patient_spectra, spectra[[s]])
				}
			}
			############################################## RANDOMNESS
			if (algorithm == "random") {
				# Do this if the spectra lengthis more than two, otherwise it does not make any sense
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
				index <- createFolds(file_vector[file_vector==patient_vector[p]], k = spectra_per_patient)
				# For each fold
				for (k in 1:length(index)) {
					# Generate a temporary list, where the patient spectra of one fold will be stored (they will be averaged)
					patient_spectra_temp <- list()
					# Create a new list of spectra for plotting purposes (the intensities will be replacd)
					patient_spectra_temp_for_plotting <- list()
					# Take the corresponding spectra subset (of the selected fold)
					for (i in 1:length(index[[k]])) {
						fold_index <- index[[k]] [i]
						patient_spectra_temp <- append(patient_spectra_temp, patient_spectra [[fold_index]])
						patient_spectra_temp_for_plotting <- append(patient_spectra_temp_for_plotting, patient_spectra [[fold_index]])
					}
					# Replace the intensities with the K number for plotting purposes
					for (n in 1:length(patient_spectra_temp_for_plotting)) {
						patient_spectra_temp_for_plotting[[n]]@intensity <- rep (k, length(patient_spectra_temp_for_plotting [[n]]@intensity))
					}
					# Add these modified spectra to the final list of spectra for plotting purposes
					spectra_for_plotting <- append(spectra_for_plotting, patient_spectra_temp_for_plotting)
					# Average them
					if (method == "mean") {
						temporary_patient_average <- averageMassSpectra (patient_spectra_temp, method="mean")
					}
					if (method == "skyline") {
						temporary_patient_average <- generate_skyline_spectrum (patient_spectra_temp)
					}
					# Store them in the final spectra list
					patient_spectra_final <- append(patient_spectra_final, temporary_patient_average)
				}
				# Store the plot into the list
				plotMsiSlice (spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot ()
			}
			######################################### SIMILARITY
			############################ HCA
			if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectra_per_patient <- spectra_per_patient + discarded_nodes
				# Detect and align peaks
				if (tof_mode=="linear") {
					peaks <- detectPeaks (patient_spectra, method="MAD", snr=3, halfWindowSize=20)
					peaks <- align_and_filter_peaks (peaks, tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25)
				}
				if (tof_mode=="reflectron" | tof_mode=="reflector") {
					peaks <- detectPeaks (patient_spectra, method="MAD", snr=3, halfWindowSize=5)
					peaks <- align_and_filter_peaks (peaks, tolerance_ppm=200, peaks_filtering=TRUE, frequency_threshold=0.25)
				}
				# Generate the peaklist matrix
				peaklist <- intensityMatrix (peaks, patient_spectra)
				# Compute the distance matrix
				distance_matrix <- dist (peaklist, method="euclidean")
				hca <- hclust(distance_matrix)
				# Store the plot
				plot(hca, xlab="MSI elements", main="Hierarchical clustering analysis", sub="")
				legend_text <- paste("Cluster method:", hca$method, "\nDistance method: ", hca$dist.method, "\n")
				legend("topright", title="Hierarchical clustering", legend=legend_text, xjust=0.5, border="black")
				plots[[p]] <- recordPlot()
				# Associate to each row/spectrum the subgroup of the tree to which it belongs
				if (spectra_per_patient > 56) {
					spectra_per_patient <- 56
				}
				hca_groups <- cutree(hca, k=spectra_per_patient)
				if (discarded_nodes != 0) {
				# For each subgroup to be isolated...
					for (d in 1:discarded_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep (0, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# They do not get added to the final list of patient spectra
						# Add them to the discarded spectra
						discarded_spectra <- append(discarded_spectra, spectra_hca)
						# Average the spectra
						if (method == "mean") {
							spectra_average <- averageMassSpectra (spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum (spectra_hca)
						}
						# Add the average to the final list of discarded spectra AVG
						discarded_spectra_average <- append(discarded_spectra_average, spectra_average)
					}
					# For each subgroup to be isolated...
					for (d in (discarded_nodes+1):spectra_per_patient) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep (d, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# Average the spectra in this HCA subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra (spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum (spectra_hca)
						}
						# Add the average to the final list
						patient_spectra_final <- append(patient_spectra_final, spectra_average)
					}
				} else {
				# For each subgroup to be isolated...
				for (s in 1:spectra_per_patient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which(hca_groups == s)
					spectra_hca <- patient_spectra[index]
					spectra_hca_for_plotting <- patient_spectra[index]
					# Replace the intensities with the S number for plotting purposes
					for (n in 1:length(spectra_hca_for_plotting)) {
						spectra_hca_for_plotting[[n]]@intensity <- rep (s, length(spectra_hca_for_plotting[[n]]@intensity))
					}
					# Add these modified spectra to the final list of spectra for plotting purposes
					spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
					# Average the spectra in this HCA subgroup
					if (method == "mean") {
						spectra_average <- averageMassSpectra (spectra_hca, method="mean")
					}
					if (method == "skyline") {
						spectra_average <- generate_skyline_spectrum (spectra_hca)
					}
					# Add the average to the final list
					patient_spectra_final <- append(patient_spectra_final, spectra_average)
				}
				}
				# Store the plot into the list
				plotMsiSlice (spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot ()
			}
			######################### K-MEANS
			if (algorithm == "k_means" | algorithm == "kmeans" | algorithm == "k-Means") {
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectra_per_patient <- spectra_per_patient + discarded_nodes
				# Detect and align peaks
				if (tof_mode=="linear") {
					peaks <- detectPeaks (patient_spectra, method="MAD", snr=3, halfWindowSize=20)
					peaks <- align_and_filter_peaks (peaks, tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25)
				}
				if (tof_mode=="reflectron" | tof_mode=="reflector") {
					peaks <- detectPeaks (patient_spectra, method="MAD", snr=3, halfWindowSize=5)
					peaks <- align_and_filter_peaks (peaks, tolerance_ppm=200, peaks_filtering=TRUE, frequency_threshold=0.25)
				}
				# Generate the peaklist matrix
				peaklist <- intensityMatrix (peaks, patient_spectra)
				# Compute the k-Means clustering
				k_means <- kmeans(peaklist, centers=spectra_per_patient)
				# Associate to each row/spectrum the subgroup/cluster to which it belongs
				k_means_groups <- k_means$cluster
				if (discarded_nodes != 0) {
				# For each subgroup to be isolated (to discard)...
					for (d in 1:discarded_nodes) {
						# Index the spectra under in the selected subgroup
						index <- which(k_means_groups == d)
						spectra_k_means <- patient_spectra[index]
						spectra_k_means_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_k_means_for_plotting)) {
							spectra_k_means_for_plotting[[n]]@intensity <- rep (0, length(spectra_k_means_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_k_means_for_plotting)
						# They do not get averaged and added to the final list of patient spectra
						# Add them to the discarded spectra
						discarded_spectra <- append(discarded_spectra, spectra_k_means)
						# Average the spectra
						if (method == "mean") {
							spectra_average <- averageMassSpectra (spectra_k_means, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum (spectra_k_means)
						}
						# Add the average to the final list of discarded spectra AVG
						discarded_spectra_average <- append(discarded_spectra_average, spectra_average)
					}
					# For each subgroup to be isolated (to keep)...
					for (d in (discarded_nodes+1):spectra_per_patient) {
						# Index the spectra under in the selected subgroup
						index <- which(k_means_groups == d)
						spectra_k_means <- patient_spectra[index]
						spectra_k_means_for_plotting <- patient_spectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectra_k_means_for_plotting)) {
							spectra_k_means_for_plotting[[n]]@intensity <- rep (d, length(spectra_k_means_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_k_means_for_plotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra (spectra_k_means, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum (spectra_k_means, method="mean")
						}
						# Add the average to the final list
						patient_spectra_final <- append(patient_spectra_final, spectra_average)
					}
				} else {
					# For each subgroup to be isolated...
					for (s in 1:spectra_per_patient) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(k_means_groups == s)
						spectra_k_means <- patient_spectra[index]
						spectra_k_means_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_k_means_for_plotting)) {
							spectra_k_means_for_plotting[[n]]@intensity <- rep (s, length(spectra_k_means_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectra_for_plotting <- append(spectra_for_plotting, spectra_k_means_for_plotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra (spectra_k_means, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum (spectra_k_means, method="mean")
						}
						# Add the average to the final list
						patient_spectra_final <- append(patient_spectra_final, spectra_average)
					}
					}
				# Store the plot into the list
				plotMsiSlice (spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot ()
			}
			#####################
		}
	}
}
# Processing before returning
patient_spectra_final <- removeBaseline (patient_spectra_final, method="TopHat")
patient_spectra_final <- calibrateIntensity (patient_spectra_final, method="TIC")
#
return (list(spectra=patient_spectra_final, ms_images=msi_plots, discarded_spectra=discarded_spectra, discarded_average_spectra=discarded_spectra_average, plots=plots))
}
################################################################################





########################################## SPECTRA GROUPING (PATIENTS) - SKYLINE
# Obtain one spectrum per patient (average)
group_spectra_skyline <- function (spectra, file_format="imzml") {
if (file_format == "imzml" | file_format == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		file_vector <- append(file_vector, spectra[[i]]@metaData$file)
	}
	# Put the class names into a vector
	path_spectra <- unique(file_vector)
	# Create the empty average_spectrumList
	spectra_grouped <- list ()
	# For each filepath...
	for (p in 1:length(path_spectra)) {
		# Create the empty list in which spectra will be allocated
		spectra_list <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length(spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$file == path_spectra [p]) {
				spectra_list <- append(spectra_list, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same path)
		skyline_spectrum  <- generate_skyline_spectrum (spectra_list)
		# Add the average spectrum to another final list of average spectra
		spectra_grouped <- append(spectra_grouped, skyline_spectrum)
	}
}
if (file_format == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		file_vector <- append(file_vector, spectra[[i]]@metaData$sampleName)
	}
	# Put the file names into a vector
	path_spectra <- unique(file_vector)
	# Create the empty average_spectrumList
	spectra_grouped <- list ()
	# For each filepath...
	for (p in 1:length(path_spectra)) {
		# Create the empty list in which spectra will be allocated
		spectra_list <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length(spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$sampleName == path_spectra [p]) {
				spectra_list <- append(spectra_list, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same path)
		skyline_spectrum  <- generate_skyline_spectrum (spectra_list)
		# Add the average spectrum to another final list of average spectra
		spectra_grouped <- append(spectra_grouped, skyline_spectrum )
	}
}
return (spectra_grouped)
}





########################################################################





########################################## SPECTRA GROUPING (CLASSES)
# Obtain one spectrum per class (average)
group_spectra_class <- function (spectra, class_list, file_format="imzml") {
class_list <- sort(class_list)
if (file_format == "imzml" | file_format == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name (spectra, class_list=class_list, file_format="imzml")
	# Put the filenames/classes in a vector
	# Create the empty vector
	class_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		class_vector <- append(class_vector, spectra[[i]]@metaData$file)
	}
	class_spectra_grouped <- averageMassSpectra (spectra, labels=class_vector, method="mean")
}
if (file_format == "brukerflex") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name (spectra, class_list=class_list, file_format="brukerflex")
	# Put the filenames/classes in a vector
	# Create the empty vector
	class_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		class_vector <- append(class_vector, spectra[[i]]@metaData$sampleName)
	}
	class_spectra_grouped <- averageMassSpectra (spectra, labels=class_vector, method="mean")
}
return (class_spectra_grouped)
}





########################################################################





########################################## SPECTRA GROUPING (CLASSES) - SKYLINE
# Obtain one spectrum per class (average)
group_spectra_class_skyline <- function (spectra, class_list, file_format="imzml") {
class_list <- sort(class_list)
if (file_format == "imzml" | file_format == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name (spectra, class_list=class_list, file_format="imzml")
	# Put the filenames/classes in a vector
	# Create the empty vector
	class_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		class_vector <- append(class_vector, spectra[[i]]@metaData$file)
	}
	# Put the class names into a vector
	class_spectra <- unique(class_vector)
	# Create the empty average_spectrumList
	class_spectra_grouped <- list ()
	# For each class...
	for (p in 1:length(class_spectra)) {
		# Create the empty list in which spectra will be allocated
		spectra_list <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length(spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$file == class_spectra [p]) {
				spectra_list <- append(spectra_list, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same class)
		skyline_spectrum  <- generate_skyline_spectrum (spectra_list)
		# Add the average spectrum to another final list of average spectra
		class_spectra_grouped <- append(class_spectra_grouped, skyline_spectrum )
	}
}
if (file_format == "brukerflex") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name (spectra, class_list=class_list, file_format="brukerflex")
	# Put the filenames/classes in a vector
	# Create the empty vector
	class_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		class_vector <- append(class_vector, spectra[[i]]@metaData$sampleName)
	}
	# Put the class names into a vector
	class_spectra <- unique(class_vector)
	# Create the empty average_spectrumList
	class_spectra_grouped <- list ()
	# For each class...
	for (p in 1:length(class_spectra)) {
		# Create the empty list in which spectra will be allocated
		spectra_list <- list ()
		# Scroll the entire spectral dataset
		for (i in 1:length(spectra)) {
			# Add the spectra with the same filepath to a list of spectra
			if (spectra[[i]]@metaData$sampleName == class_spectra [p]) {
				spectra_list <- append(spectra_list, spectra[[i]])
			}
		}
		# Create the average spectrum of the spectra in the list (with the same class)
		skyline_spectrum  <- generate_skyline_spectrum (spectra_list)
		# Add the average spectrum to another final list of average spectra
		class_spectra_grouped <- append(class_spectra_grouped, skyline_spectrum )
	}
}
return (class_spectra_grouped)
}





################################################################################





######################################### ADD THE CLASS AND THE SAMPLE NAME TO THE MATRIX
matrix_add_class_and_sample <- function (signal_matrix, peaks=list(), class_list=list(), file_format="imzml", sample_output=TRUE, class_output=TRUE) {
signal_matrix <- as.matrix(signal_matrix)
number_of_spectra <- length(peaks)
### The name of the rows will be either the sample name or the class name (depending on the function parameter)
# If the rows are named according to the sample name, an additional column for the class is added
# Otherwise, use the class names as the row names
if ((class_output==FALSE && sample_output==TRUE) || (class_output==TRUE && length(class_list)==0 && sample_output==TRUE)) {
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length(if averaged, it contains the name of all the single spectra)
		if (file_format == "imzml" | file_format == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
			}
		}
		if (file_format == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	# Create the sample matrix column and appendit to the global matrix
	sample_column <- matrix (0, ncol=1, nrow=number_of_spectra)
	colnames(sample_column) <- "Sample"
	sample_column [,1] <- cbind(file_vector)
	signal_matrix <- cbind(signal_matrix, sample_column)
}
if (class_output==TRUE && length(class_list)>=1 && sample_output==TRUE) {
	### Put the sample names in the first column of the matrix
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length(if averaged, it contains the name of all the single spectra)
		if (file_format == "imzml" | file_format == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
			}
		}
		if (file_format == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	# Create the sample matrix column and appendit to the global matrix
	sample_column <- matrix (0, ncol=1, nrow=number_of_spectra)
	colnames(sample_column) <- "Sample"
	sample_column [,1] <- cbind(file_vector)
	signal_matrix <- cbind(signal_matrix, sample_column)
	### Add the class column
	class_list <- sort(class_list)
	class_column <- matrix (0, ncol=1, nrow=number_of_spectra)
	colnames(class_column) <- "Class"
	# Rename the classes according to the class_list vector
	class_vector <- file_vector
	for (p in 1:length(class_vector)) {
		for (w in 1:length(class_list)) {
			if (length(grep(class_list[w],class_vector[p], ignore.case=TRUE)) == 1) {
				class_vector[p] <- class_list [w]
			}
		}
	}
	# Fill in the matrix column with the file_vector classes and samples
	class_column [,1] <- cbind(class_vector)
	signal_matrix <- cbind(signal_matrix, class_column)
}
if (class_output==TRUE && length(class_list)>=1 && sample_output==FALSE) {
	class_list <- sort(class_list)
	### Put the sample names in the first column of the matrix
	# Create the empty vector
	file_vector <- vector(length=0)
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		# Check the length(if averaged, it contains the name of all the single spectra)
		if (file_format == "imzml" | file_format == "imzML") {
			if (length(peaks[[i]]@metaData$file) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file)
			}
			if (length(peaks[[i]]@metaData$file) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
			}
		}
		if (file_format == "brukerflex") {
			if (length(peaks[[i]]@metaData$sampleName) == 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
			}
			if (length(peaks[[i]]@metaData$sampleName) > 1) {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	}
	### Add the class column
	class_column <- matrix (0, ncol=1, nrow=number_of_spectra)
	colnames(class_column) <- "Class"
	# Rename the classes according to the class_list vector
	class_vector <- vector(length=0)
	for (w in 1:length(class_list)) {
		for (p in 1:length(file_vector)) {
			if (length(grep(class_list[w],file_vector[p], ignore.case=TRUE)) !=0) {
				class_vector[p] <- class_list [w]
			}
		}
	}
	# Fill in the matrix column with the file_vector classes and samples
	class_column [,1] <- cbind(class_vector)
	signal_matrix <- cbind(signal_matrix, class_column)
}
### Add these matrix columns to the peaklist matrix
return (signal_matrix)
}





########################################################################





###################################### NO HEMOGLOBIN
no_hemoglobin <- function (spectra, hemoglobin_mass=15400, tolerance_ppm=2000, average_spectrum_evaluation=FALSE) {
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
#################################################### FILTERING FUNCTION
hemoglobin_removal_function <- function (spectra, hemoglobin_mass, tolerance_ppm) {
## Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=3)
hemoglobin_bin <- numeric()
# Do this if there are peaks
if (length(peaks@mass) > 0) {
for (s in 1:length(peaks@mass)) {
	# Select the Hemoglobin peaks (based on the difference from the exact mass)
	if ((abs(peaks@mass[s] - hemoglobin_mass)*10^6/hemoglobin_mass) <= tolerance_ppm) {
		# That is the peak of interest, add its intensity to a bin
		hemoglobin_bin <- append(hemoglobin_bin, peaks@intensity[s])
	}
}
# Evaluate the peaks in the hemoglobin Bin (if there are any)...
if (length(hemoglobin_bin) > 0) {
	# Start assuming that the spectrum is good
	is_spectrum_to_be_removed <- FALSE
	# If one of this is the base peak...
	for (h in 1:length(hemoglobin_bin)) {
		if (hemoglobin_bin[h] == max(peaks@intensity)) {
			# The spectrum is bad and to be removed
			is_spectrum_to_be_removed <- TRUE
		}
	}
	# Remove the spectrum if it is to be removed
	if (is_spectrum_to_be_removed == TRUE) {
		spectra <- NULL
	}
}
}
return (spectra)
}
########################################################################
spectra_no_hemoglobin <- parLapply(cl, spectra, fun=function (spectra) hemoglobin_removal_function (spectra, hemoglobin_mass, tolerance_ppm))
# Close the processes
stopCluster(cl)
### Keep only the elements that are different from NULL
spectra_no_hemoglobin <- spectra_no_hemoglobin [!sapply(spectra_no_hemoglobin, is.null)]
##### Now we have discarded all the spectra with Hemoglobin, but the average can still be bad
# Evaluate the average spectrum
if (average_spectrum_evaluation == TRUE) {
	# Compute the average spectrum of this filtered dataset
	average_spectrum <- averageMassSpectra (spectra_no_hemoglobin, method="mean")
	## Peak picking on the AVG
	peaks_average <- detectPeaks (average_spectrum, method="MAD", snr=3)
	hemoglobin_bin <- numeric()
	# Do this if there are peaks
	if (length(peaks_average@mass) > 0) {
	for (s in 1:length(peaks_average@mass)) {
		# Select the Hemoglobin peaks (based on the difference from the exact mass)
		if ((abs(peaks_average@mass[s] - hemoglobin_mass)*10^6/hemoglobin_mass) <= tolerance_ppm) {
			# That is the peak of interest, add its intensity to a bin
			hemoglobin_bin <- append(hemoglobin_bin, peaks_average@intensity[s])
		}
	}
	# Evaluate the peaks in the hemoglobin Bin (if there are any)...
	if (length(hemoglobin_bin) > 0) {
	# Start assuming that the spectrum is good
	is_spectrum_to_be_removed <- FALSE
	# If one of this is the base peak...
	for (h in 1:length(hemoglobin_bin)) {
		if (hemoglobin_bin[h] == max(peaks_average@intensity)) {
			# The spectrum is bad and to be removed
			is_spectrum_to_be_removed <- TRUE
		}
	}
	# Remove the spectra if they are to be removed
	if (is_spectrum_to_be_removed == TRUE) {
		spectra_no_hemoglobin <- list()
	}
	}
	}
}
return (spectra_no_hemoglobin)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE
biotyper_performance <- function (sample_vector, class_list, disease_name="Diseased", healthy_name="Controls") {
class_list <- sort(class_list)
# Declare the parameters
true_positive <- 0
true_negative <- 0
false_positive <- 0
false_negative <- 0
# For each tested sample
for (s in 1:number_of_samples) {
	# For each class taken into account
	for (w in 1:length(class_list)) {
		# If the sample and the class match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) & (length(grep(class_list[w],disease_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true positive number increases
			true_positive <- true_positive + 1
		}
		# If the sample and the class match, the class is HP and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) & (length(grep(class_list[w],healthy_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true negative number increases
			true_negative <- true_negative + 1
		}
		# If the sample and the class do not match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) & (length(grep(class_list[w],disease_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false positive number increases
			false_positive <- false_positive + 1
		}
		# If the sample and the class do not match, the class is HP and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) & (length(grep(class_list[w],healthy_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false negative number increases
			false_negative <- false_negative + 1
		}
    }
}
# Calculate the performance values
sensitivity <- (true_positive / (true_positive + false_negative))*100
specificity <- (true_negative / (false_positive + true_negative))*100
ppv <- true_positive / (true_positive + false_positive)
npv <- true_negative / (true_negative + false_negative)
fdr <- 1 - ppv
for_value <- 1 - npv
# Create the output matrix
test_performance <- matrix (0, nrow=6, ncol=1)
test_performance_rows <- c("Sensitivity", "Specificity", "PPV", "NPV", "FDR", "FOR")
rownames(test_performance) <- test_performance_rows
colnames(test_performance) <- "Biotyper"
test_performance [1,1] <- paste(sensitivity,"%")
test_performance [2,1] <- paste(specificity, "%")
test_performance [3,1] <- ppv
test_performance [4,1] <- npv
test_performance [5,1] <- fdr
test_performance [6,1] <- for_value
return (test_performance)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE 2
biotyper_performance_2 <- function (sample_vector, class_list, disease_name="Diseased", healthy_name="Controls") {
class_list <- sort(class_list)
# Declare the parameters
true_positive <- 0
true_negative <- 0
false_positive <- 0
false_negative <- 0
# For each tested sample
for (s in 1:number_of_samples) {
	# For each class taken into account
	for (w in 1:length(class_list)) {
		# If the sample and the class match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) & (length(grep(class_list[w],disease_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true positive number increases
			true_positive <- true_positive + 1
		}
		# If the sample and the class match, the class is HP and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) & (length(grep(class_list[w],healthy_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The true negative number increases
			true_negative <- true_negative + 1
		}
		# If the sample and the class do not match, the class is PTC and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) & (length(grep(class_list[w],disease_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false positive number increases
			false_positive <- false_positive + 1
		}
		# If the sample and the class do not match, the class is HP and it is a YES (score >=2)
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) & (length(grep(class_list[w],healthy_name, ignore.case=TRUE)) == 1) & (score[s,w]>=2)) {
			# The false negative number increases
			false_negative <- false_negative + 1
		}
	}
}
# Calculate the performance values
sensitivity <- (true_positive / (true_positive + false_negative))*100
specificity <- (true_negative / (false_positive + true_negative))*100
ppv <- true_positive / (true_positive + false_positive)
npv <- true_negative / (true_negative + false_negative)
fdr <- 1 - ppv
for_value <- 1 - npv
# Create the output matrix
test_performance <- matrix (0, nrow=3, ncol=3)
rownames(test_performance) <- c("Test outcome positive", "Test outcome negative", "")
colnames(test_performance) <- c("Condition positive", "Condition negative", "")
test_performance [1,1] <- paste("True positive =", true_positive)
test_performance [1,2] <- paste("False positive =", false_positive)
test_performance [2,1] <- paste("False negative =", false_negative)
test_performance [2,2] <- paste("True negative =", true_negative)
test_performance [3,1] <- paste(sensitivity,"%")
test_performance [3,2] <- paste(specificity, "%")
test_performance [1,3] <- paste("Positive predictive value =", ppv)
test_performance [2,3] <- paste("Negative predictive value =", npv)

return (test_performance)
}





########################################################################





################# PLOT THE MEAN SPECTRUM WITH THE SD BARS ON THE AVERAGE
average_spectrum_bars <- function (spectra, snr=5, tolerance_ppm=2000, mass_range_plot=c(4000,15000), graph_title="Spectrum", average_spectrum_colour="black", peak_points="yes", points_colour="red", bar_width=40, bar_colour="blue") {
# Generate the average spectrum
average_spectrum <- averageMassSpectra (spectra, method="mean")
average_spectrum <- removeBaseline (average_spectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaks_average <- detectPeaks (average_spectrum, method="MAD", snr=snr)
# Peak picking on the dataset
peaks <- detectPeaks (spectra, method="MAD", snr=snr)
peaks <- removeEmptyMassObjects (peaks)
# Alignment: merge the two lists and align them all
peaks_global <- list()
peaks_global <- peaks
peaks_global <- append(peaks_global, peaks_average)
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
# Empty the lists
peaks_average <- list()
peaks <- list()
# Re-fill the lists with the aligned peaklists
for (i in 1:(length(peaks_global) - 1)) {
	peaks <- append(peaks, peaks_global[[i]])
}
peaks_average <- peaks_global[[length(peaks_global)]]
######## Generate a vector with the standard deviations of the peaks in the average peaklist
st_dev_intensityensity <- vector(length=0)
# For each peak in the average peaklist
for (a in 1:length(peaks_average@mass)) {
	intensity_vector <- vector(length=0)
	# Search for it in the dataset peaklist
	for (p in 1:length(peaks)) {
		for (i in 1:length(peaks[[p]]@mass)) {
			if (abs(peaks[[p]]@mass[i] == peaks_average@mass[a])) {
				intensity_vector <- append(intensity_vector, peaks[[p]]@intensity[i])
			}
		}
	}
	st_dev_intensityensity <- append(st_dev_intensityensity, sd(intensity_vector))
}
####### Plot
# Spectrum
plot(average_spectrum, main=graph_title, col.main=average_spectrum_colour, xlab="m/z", ylab="Intensity (a.i.)", xlim=mass_range_plot, col=average_spectrum_colour)
# Peaks
if (peak_points == "yes") {
	points(peaks_average, pch=4, col=points_colour)
}
# Bars
# epsilon: lengthof the horizontal segment
epsilon = bar_width
for (i in 1:length(peaks_average@mass)) {
	# Define the upper and the lower limit of the vertical segment
	up <- peaks_average@intensity[i] + st_dev_intensityensity [i]
	low <- peaks_average@intensity[i] - st_dev_intensityensity [i]
	# Vertical bar (x,y x,y)
	segments(peaks_average@mass[i], low, peaks_average@mass[i], up, col=bar_colour)
	# Horizontal segments(x,y , x,y)
	segments(peaks_average@mass[i]-epsilon, low, peaks_average@mass[i]+epsilon, low, col=bar_colour)
	segments(peaks_average@mass[i]-epsilon, up, peaks_average@mass[i]+epsilon, up, col=bar_colour)
}
}





########################################################################





################# PLOT THE SIGNALS OF INTEREST WITH THE SD BARS ON THE AVERAGE
average_spectrum_bars_signals_of_interest <- function (spectra, snr=5, signals_of_interest=peaks_average@mass, tolerance_ppm=2000, half_window_plot=1000, graph_title="Spectrum", average_spectrum_colour="black", peak_points=TRUE, points_colour="red", bar_width=40, bar_colour="blue") {
# Generate the average spectrum
average_spectrum <- averageMassSpectra (spectra, method="mean")
average_spectrum <- removeBaseline (average_spectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaks_average <- detectPeaks (average_spectrum, method="MAD", snr=snr)
# Peak picking on the dataset
peaks <- detectPeaks (spectra, method="MAD", snr=snr)
peaks <- removeEmptyMassObjects (peaks)
######## Generate a vector with the standard deviations of the peaks of interest
# For each peak in the average peaklist
for (a in 1:length(signals_of_interest)) {
	# Search for the custom mass value in the average spectrum, in order to assign the real value
	for (j in 1:length(peaks_average@mass)) {
		if ((abs(signals_of_interest[a] - peaks_average@mass[j])*10^6/peaks_average@mass[j]) <= tolerance_ppm) {
			signals_of_interest[a] <- peaks_average@mass[j]
			# Position of the peak of interest in the average_spectrum peak list
			peak_index <- j
		}
	}
	# Generate the intensity vector
	intensity_vector <- vector(length=0)
	# Search for it in the dataset peaklist
	for (p in 1:length(peaks)) {
		for (i in 1:length(peaks[[p]]@mass)) {
			if ((abs((signals_of_interest[a] - peaks[[p]]@mass[i])/peaks[[p]]@mass[i])*10^6) <= tolerance_ppm) {
				intensity_vector <- append(intensity_vector, peaks[[p]]@intensity[i])
			}
		}
	}
	# Calculate the standard deviation
	st_dev_intensityensity <- sd(intensity_vector)
	####### Plot
	# Define the limits on the x-axis
	low_mass_plot <- signals_of_interest[a] - half_window_plot
	up_mass_plot <- signals_of_interest[a] + half_window_plot
	# Spectrum
	#X11() #Opens a new graphic instance
	plot(average_spectrum, main=graph_title, col.main=average_spectrum_colour, xlab="m/z", ylab="Intensity (a.i.)", xlim=c(low_mass_plot,up_mass_plot), col=average_spectrum_colour)
	# Peaks
	if (peak_points == TRUE) {
		points(peaks_average[peak_index], pch=4, col=points_colour)
	}
	# Bars
	# epsilon: lengthof the horizontal segment
	epsilon = bar_width
	# Define the upper and the lower limit of the vertical segment
	up <- peaks_average@intensity[peak_index] + st_dev_intensityensity
	low <- peaks_average@intensity[peak_index] - st_dev_intensityensity
	# Vertical bar (x,y x,y)
	segments(signals_of_interest[a], low, signals_of_interest[a], up, col=bar_colour)
	# Horizontal segments(x,y , x,y)
	segments(signals_of_interest[a]-epsilon, low, signals_of_interest[a]+epsilon, low, col=bar_colour)
	segments(signals_of_interest[a]-epsilon, up, signals_of_interest[a]+epsilon, up, col=bar_colour)
}
}





########################################################################





################################### MEMORY EFFICIENT IMPORTING
################## FOR EACH PATIENT: IMPORT, PREPROCESS, FILTER, AVERAGE, REPRESENTATIVE SPECTRA
memory_efficient_import <- function (folder, tof_mode="linear", tic_purification=FALSE, absolute_tic_threshold=0, smoothing_strength="medium", mass_range=c(0,0), preprocess_spectra=TRUE, process_in_packages_of=length(spectra), generate_representative_spectra=FALSE, spectra_per_patient=5, algorithm_for_representative_spectra="hca", discarded_nodes=1, average_patient=FALSE, skyline=FALSE, align_spectra=FALSE, alignment_tolerance_ppm=2000, file_format="imzml") {
setwd(folder)
folder_files <- read_spectra_files (folder, file_format=file_format, full_path=FALSE)
spectra_dataset <- list ()
# For each imzML file (patient)...
for (i in 1:length(folder_files)) {
	### Load the spectra
	if (file_format == "imzml" | file_format == "imzML") {
		spectra <- importImzMl(folder_files[i])
	}
	if (file_format == "brukerflex") {
		spectra <- importBrukerFlex(folder_files[i])
	}
	if (tic_purification == TRUE) {
		spectra <- spectra_tic_purification (spectra, tof_mode=tof_mode, absolute_tic_threshold=absolute_tic_threshold)
	}
	if (average_patient == FALSE) {
	### Preprocessing
		if (preprocess_spectra == TRUE) {
			spectra <- preprocess_spectra(spectra, tof_mode=tof_mode, smoothing_strength=smoothing_strength, process_in_packages_of=process_in_packages_of)
		}
	}
	### Average the spectra
	if (average_patient == TRUE) {
		# At least two spectra for the averaging
		if (length(spectra) > 1) {
			if (tof_mode == "linear") {
				spectra <- averageMassSpectra (spectra, method="mean")
				spectra <- smoothIntensity (spectra, method="SavitzkyGolay", halfWindowSize=10)
				spectra <- removeBaseline (spectra, method="TopHat")
				spectra <- calibrateIntensity (spectra, method="TIC")
			}
			if (tof_mode == "reflector" || tof_mode == "reflectron") {
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
			spectra <- generate_skyline_spectrum (spectra)
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
		spectra_dataset <- append(spectra_dataset, spectra)
	}
	### Free the memory
	rm(spectra)
	gc()
}
### Trimming the entire dataset
	if (mass_range[1] == 0 && mass_range[2] == 0) {
		spectra_dataset <- trim(spectra_dataset)
	}
	if ((mass_range[1] != 0 || mass_range[2] == 0) && mass_range[2] != 0) {
		spectra_dataset <- trim(spectra_dataset, range=mass_range)
	}
	if (mass_range[1] != 0 && mass_range[2] == 0) {
		upperCroppingBoundary <- Inf
		spectra_dataset <- trim(spectra_dataset, range=mass_range)
	}
	spectra_dataset <- calibrateIntensity (spectra_dataset, method="TIC")
### Spectral alignment
if (align_spectra==TRUE) {
	# Determine the warping function
	peaks <- detectPeaks (spectra_dataset, snr=5)
	reference_for_alignment <- referencePeaks (peaks, method="relaxed", tolerance=(alignment_tolerance_ppm/10^6), minFrequency=0.9)
	if (length(reference_for_alignment@mass) > 0) {
		if (tof_mode == "linear") {
			spectra_dataset <- align_spectra (spectra_dataset, halfWindowSize=20, noiseMethod="MAD", snr=3,
					reference=reference_for_alignment, tolerance=(alignment_tolerance_ppm/10^6), warpingMethod="lowess")
		}
		if (tof_mode == "reflector" | tof_mode == "reflectron") {
			spectra_dataset <- align_spectra (spectra_dataset, halfWindowSize=5, noiseMethod="MAD", snr=3,
					reference=reference_for_alignment, tolerance=(alignment_tolerance_ppm/10^6), warpingMethod="lowess")
		}
	}
}
### Generate representative spectra
if (generate_representative_spectra == TRUE) {
	if (length(spectra) > 1) {
		spectra <- group_spectra (spectra, spectra_per_patient=spectra_per_patient, file_format=file_format, tof_mode=tof_mode, seed=0, algorithm=algorithm_for_representative_spectra, discarded_nodes=discarded_nodes)
	}
	# If there is only one spectrum left, use it
	if (length(spectra) == 1) {
		spectra <- spectra
	}
}
return (spectra_dataset)
}





########################################################################





#################### SIGNAL COUNTER IN THE AVERAGE SPECTRUM
##### Count signals in a mass range
countSignalsRangeAverageSpectrum <- function (spectra, snr=5, lower=4000, interval_width=1000, max=20000) {
average_spectrum <- averageMassSpectra (spectra, method="mean")
average_spectrum <- removeBaseline (average_spectrum, method="TopHat")
# Peak picking on the average_spectrum
peaks_average <- detectPeaks (average_spectrum, method="MAD", snr=snr)
# Vector with signals
peak_vector <- peaks_average@mass
#Create the result_matrix
result_matrix <- matrix (ncol=1, nrow=1)
result_matrix [1,1] <- "Signal counter"
colnames(result_matrix) <- "Signals"
# Count the signals in the specified mass range, with a counter
counter <- 0
upper <- lower + interval_width
for (n in 1:((max-lower)/interval_width)) {
	for (i in 1:length(peak_vector)) {
		if (peak_vector [i] >= lower & peak_vector [i] <= upper & upper <= max) {
			counter <- counter + 1
		}
	}
	# Put the counter into a matrix column
	signal_counter <- matrix (ncol=1, nrow=1)
	colnames(signal_counter) <- paste(lower, "-", upper)
	signal_counter [1,1] <- counter
	# Add the signal_counter column to the result_matrix, each time
	result_matrix <- cbind(result_matrix, signal_counter)
	# Move to the next interval
	lower <- lower + interval_width
	upper <- upper + interval_width
	counter <- 0
}
return (result_matrix)
}





########################################################################





#################### SIGNAL COUNTER IN THE AVERAGE SPECTRUM
##### Count signals in a mass range
count_signals_range_spectra <- function (spectra, snr=5, lower=4000, interval_width=1000, max=20000) {
peaks <- detectPeaks (spectra, method="MAD", snr=snr)
#Create the result_matrix
result_matrix <- matrix ("Number of signals", ncol=(((max-lower)/interval_width)+1), nrow=1)
# Count the signals in the specified mass range, with a counter
counter <- 0
# The first upper value would be
upper <- lower + interval_width
# Then it will increase during the loop
# Scroll the "peaks" list and create a vector each time containing the peaklist of that spectrum
for (s in 1:length(peaks)) {
	peak_vector <- peaks[[s]]@mass
	#Create the result_matrix_row (this would be referred to each spectrum)
	result_matrix_row <- matrix (ncol=1, nrow=1)
	# Reiterate for every interval (from 1 to the number of intervals)
	for (n in 1:((max-lower)/interval_width)) {
		# Check every value of the peak_vector and adjust the counter based on the condition
		for (i in 1:length(peak_vector)) {
			if (peak_vector [i] >= lower & peak_vector [i] <= upper & upper <= max) {
				counter <- counter + 1
			}
		}
		# Put the counter into a matrix column
		signal_counter <- matrix (ncol=1, nrow=1)
		colnames(signal_counter) <- paste(lower, "-", upper)
		signal_counter [1,1] <- counter
		# Add the signal_counter column to the result_matrix, each time
		result_matrix_row [1,1] <- s
		result_matrix_row <- cbind(result_matrix_row, signal_counter)
		# Move to the next interval
		lower <- lower + interval_width
		upper <- upper + interval_width
		counter <- 0
	}
	# Add this row to the result_matrix
	result_matrix <- rbind (result_matrix, result_matrix_row)
	# Restore the initial values, in order to reiterate for the next spectrum
	lower <- 4000
	counter <- 0
	upper <- lower + interval_width
}
result_matrix [1,1] <- "Spectrum number"
return (result_matrix)
}





########################################################################





################################## GENERATE A SKYLINE SPECTRUM
generate_skyline_spectrum <- function (spectra) {
# Generate the matrix containing all the intensities of the data points for each spectrum
spectra_matrix <- matrix (0, nrow=(length(spectra)+1), ncol=length(spectra[[1]]@mass))
colnames(spectra_matrix) <- spectra[[1]]@mass
# Fill in the matrix rows with the intensity values of the data points for each spectrum
for (s in 1:length(spectra)) {
	spectra_matrix [s,] <- spectra[[s]]@intensity
}
# In the last row, calculate the maximum intensity value of all the data points
for (i in 1:length(spectra_matrix [(length(spectra)+1),])) {
	spectra_matrix [(length(spectra)+1),i] <- max (spectra_matrix[,i])
}
# Generate the mean spectrum with the R function (it will be replaced by the skyline)
skyline_spectrum  <- averageMassSpectra (spectra, method="mean")
# Replace the mean intensities with the maximum intensities
skyline_spectrum @intensity <- spectra_matrix [(length(spectra)+1),]
# Return the function
return (skyline_spectrum )
}





########################################################################





####################### LOWEST NUMBER OF DATAPOINTS
estimate_lowest_data_points <- function (spectra) {
lowest_data_points <- NULL
# Compare it with all the others
for (s in 1:length(spectra)) {
	if (lowest_data_points == NULL || length(spectra[[s]]@mass) < lowest_data_points) {
		lowest_data_points <- length(spectra[[s]]@mass)
	}
}
return (lowest_data_points)
}





########################################################################





###################### SPECTRA BINNING
resample_spectra <- function (spectra, final_data_points=lowest_data_points, binning_method="sum") {
####################################################### BINNING FUNCTION
binning_function <- function (spectra, final_data_points, binning_method) {
	# Create the new spectra_binned list
	spectra_binned <- spectra
	# Empty the mass and intensity values
	for (s in 1:length(spectra_binned)) {
		spectra_binned@mass <- numeric()
		spectra_binned@intensity <- numeric()
	}
	# Calculate the number of datapoints per bin
	data_points_per_bin <- length(spectra@mass) / final_data_points
	data_points_per_bin <- floor (data_points_per_bin)
	# Define the indexes
	index1 <- 1
	index2 <- data_points_per_bin
	# For each bin...
	for (d in 1:final_data_points) {
		# Create the (temporary) bin vectors (mass and intensity), where the data points will be stored for the binning
		bin_mass <- numeric()
		bin_intensity <- numeric()
		# Scroll the data points, grouping them by bins
		for (i in index1:index2) {
			bin_mass <- append(bin_mass, spectra@mass[i])
			bin_intensity <- append(bin_intensity, spectra@intensity[i])
		}
		# Calculate the value of the new data point
		data_point_mass <- mean(bin_mass)
		if (binning_method == "sum") {
			data_point_intensity <- sum(bin_intensity)
		}
		if (binning_method == "mean") {
			data_point_intensity <- mean(bin_intensity)
		}
		if (binning_method == "median") {
			data_point_intensity <- median(bin_intensity)
		}
		if (method == "RMS" | method == "rms") {
			data_point_intensity <- sqrt(sum(bin_intensity^2))
		}
		# Store it in the new spectra_binned list
		spectra_binned@mass <- append(spectra_binned@mass, data_point_mass)
		spectra_binned@intensity <- append(spectra_binned@intensity, data_point_intensity)
		# Increase the indexes
		index1 <- index1 + data_points_per_bin
		index2 <- index2 + data_points_per_bin
	}
	return (spectra_binned)
}
#######################################################################
# Calculate the lowest amount of data points, that corresponds to the maximum
# number of data points that can be used for the binning
# Use this value as default if final_data_points is not specified (Function)
lowest_data_points <- estimate_lowest_data_points (spectra)
# Check the qualities of the spectral dataset
data_points_table <- table(sapply(spectra, length))
datapoints_dataset <- as.numeric(names(data_points_table))
equality_data_points <- length(data_points_table)
######## Do not bin if all the spectra are of the same lengthand have the same number of datapoints as defined
cl <- makeCluster(cpu_core_number)
if ((equality_data_points == 1) && (datapoints_dataset == final_data_points)) {
	spectra_binned <- spectra
}
######## Perform the binning if the spectra are not of the same lengthor
# they are of the same lengthbut with a different number of datapoints than
# the desired one
if (!(equality_data_points == 1) || ((equality_data_points == 1) && !(datapoints_dataset == final_data_points))) {
	# Do the binning only if the number of final data points is lower than the
	# lowest number of original data points
	if (final_data_points > lowest_data_points) {
		finalDatapoints <- lowest_data_points
		print("Binning at this sample rate is not possible, the highest number of data points possible will be used")
		spectra_binned <- parLapply(cl, spectra, fun=function (spectra) binning_function (spectra, final_data_points, binning_method))
	}
	if (final_data_points <= lowest_data_points) {
		spectra_binned <- parLapply(cl, spectra, fun=function (spectra) binning_function (spectra, final_data_points, binning_method))
	}
}
# Close the processes
stopCluster(cl)
return (spectra_binned)
print(table(sapply(spectra_binned, length)))
print(paste("Equal distance between datapoints", (all(sapply(spectra_binned, isRegular)))))
}





########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N MOST INTENSE
most_intense_signals <- function (spectra, signals_to_take=20) {
####################################################### PICKING FUNCTION
picking_function <- function (peaks, signals_to_take) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Sort the dataframe according to the snr
	peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$snr),]
	# Select only the first most intense signals
	selected_signals <- peaks_data_frame [1:signals_to_take,]
	# Sort the dataframe back according to mass
	selected_signals <- selected_signals [order(selected_signals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selected_signals$mass
	peaks@intensity <- selected_signals$intensity
	peaks@snr <- selected_signals$snr
	return (peaks)
}
########################################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=3)
most_intense_peaks <- mclapply(peaks, FUN= function (peaks) picking_function (peaks, signals_to_take=signals_to_take))
return (most_intense_peaks)
}






########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N LESS INTENSE
less_intense_signals <- function (spectra, base_snr=2, signals_to_take=20) {
############################################################
picking_function <- function (peaks, signals_to_take) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Sort the dataframe according to the snr
	peaks_data_frame <- peaks_data_frame [order(peaks_data_frame$snr),]
	# Select only the first most intense signals
	selected_signals <- peaks_data_frame [1:signals_to_take,]
	# Sort the dataframe back according to mass
	selected_signals <- selected_signals [order(selected_signals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selected_signals$mass
	peaks@intensity <- selected_signals$intensity
	peaks@snr <- selected_signals$snr
	return (peaks)
}
###########################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=base_snr)
less_intense_peaks <- mclapply(peaks, FUN= function (peaks) picking_function (peaks, signals_to_take=signals_to_take))
return (less_intense_peaks)
}





########################################################################
####### LESS AND MOST INTENSE SIGNALS
highest_and_lowest_intensity_signals <- function (spectra, base_snr=3, signals_to_take) {
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=base_snr)
# For each peaklist...
for (p in 1:length(peaks)) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks[[p]]@mass, intensity = peaks[[p]]@intensity, snr = peaks[[p]]@snr)
	# Sort the dataframe according to the snr
	most_intense_peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$snr),]
	less_intense_peaks_data_frame <- peaks_data_frame [order(peaks_data_frame$snr),]
	# Select only the first most intense signals
	selected_signals_most_intense <- most_intense_peaks_data_frame [1:(signals_to_take/2),]
	selected_signals_less_intense <- less_intense_peaks_data_frame [1:(signals_to_take/2),]
	# Sort the dataframe back according to mass
	selected_signals_most_intense <- selected_signals_most_intense [order(selected_signals_most_intense$mass),]
	selected_signals_less_intense <- selected_signals_less_intense [order(selected_signals_less_intense$mass),]
	# Put these signals back into the peaklist
	peaks[[p]]@mass <- selected_signals_less_intense$mass
	peaks[[p]]@intensity <- selected_signals_less_intense$intensity
	peaks[[p]]@snr <- selected_signals_less_intense$snr
	peaks[[p]]@mass <- append(peaks[[p]]@mass, selected_signals_most_intense$mass)
	peaks[[p]]@intensity <- append(peaks[[p]]@intensity, selected_signals_most_intense$intensity)
	peaks[[p]]@snr <- append(peaks[[p]]@snr, selected_signals_most_intense$snr)
}
return (peaks)
}





########################################################################





################################### DIAGNOSTIC TEST PERFORMANCE
biotyper_performance_wide <- function (biotyper_score) {
class_list <- colnames(biotyper_score)
sample_vector <- rownames(biotyper_score)
number_of_samples <- length(sample_vector)
# Declare the parameters
correctly_identified <- 0
falsely_identified <- 0
not_identified <- 0
# For each tested sample
for (s in 1:number_of_samples) {
	# For each class taken into account
	for (w in 1:length(class_list)) {
		# If the sample and the class match and it is a YES
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) && (length(grep("YES", biotyper_score[s,w], ignore.case=TRUE)) == 1)) {
			# The correctly identified number increases
			correctly_identified <- correctly_identified + 1
		}
		# If the sample and the class match and it is a NI
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) && (length(grep("NI", biotyper_score[s,w], ignore.case=TRUE)) == 1)) {
			# The correctly identified number increases
			correctly_identified <- correctly_identified + 0.5
		}
		# If the sample and the class do not match and it is a YES
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) && (length(grep("YES", biotyper_score[s,w], ignore.case=TRUE)) == 1)) {
			# The falsely identified number increases
			falsely_identified <- falsely_identified + 1
		}
		# If the sample and the class do not match and it is a NI
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 0) && (length(grep("NI", biotyper_score[s,w], ignore.case=TRUE)) == 1)) {
			# The falsely identified number increases
			falsely_identified <- falsely_identified + 0.5
		}
		# If the sample and the class match and it is a NO
		if ((length(grep(class_list[w],sample_vector[s], ignore.case=TRUE)) == 1) & (length(grep("NO", biotyper_score[s,w], ignore.case=TRUE)) == 1)) {
			# The not identified number increases
			not_identified <- not_identified + 1
		}
	}
}
# Calculate the performance values
sensitivity <- (correctly_identified / number_of_samples)*100
false_positive_rate <- (falsely_identified / number_of_samples)*100
missed_identifications <- (not_identified / number_of_samples)*100
# Create the output matrix
test_performance <- matrix (0, nrow=3, ncol=1)
test_performance_rows <- c("Sensitivity", "False positive rate", "Missed Identifications")
rownames(test_performance) <- test_performance_rows
colnames(test_performance) <- "Biotyper"
test_performance [1,1] <- paste(sensitivity,"%")
test_performance [2,1] <- paste(false_positive_rate, "%")
test_performance [3,1] <- paste(missed_identifications, "%")
return (test_performance)
}





########################################################################





####################### SPECTRA LIST CLASS
spectra_list_class <- function (spectra, class_name, file_format="imzml") {
spectra_list <- list()
if (file_format == "imzml" | file_format == "imzML") {
	for (s in 1:length(spectra)) {
		if (length(grep(class_name,spectra[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
			spectra_list <- append(spectra_list, spectra[[s]])
		}
	}
}
if (file_format == "brukerflex") {
	for (s in 1:length(spectra)) {
		if (length(grep(class_name,spectra[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
			spectra_list <- append(spectra_list, spectra[[s]])
		}
	}
}
return (spectra_list)
}




########################################################################





########################### ALIGN THE DATAPOINTS AFTER TRUNCATION OR BINNING
align_data_points_after_truncation <- function (spectra, tof_mode="linear") {
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
########################################################################
data_points_alignment_function <- function (spectra, reference_spectrum, tof_mode) {
data_points_difference <- spectra@mass[1] - reference_spectrum@mass[1]
if (tof_mode == "linear") {
	if (data_points_difference <= 4) {
		for (m in 1:length(spectra@mass)) {
			spectra@mass[m] <- spectra@mass[m] - data_points_difference
		}
	}
}
if (tof_mode == "reflector") {
	if (data_points_difference <= 1) {
		for (m in 1:length(spectra@mass)) {
			spectra@mass[m] <- spectra@mass[m] - data_points_difference
		}
	}
}
return (spectra)
}
########################################################################
#### Align the data points after the truncation
reference_spectrum <- spectra[[1]]
spectra <- parLapply(cl, spectra, fun=function(spectra) data_points_alignment_function (spectra, reference_spectrum, tof_mode=tof_mode))
# Close the processes
stopCluster(cl)
return (spectra)
}





##################################### BIOTYPER SCORE
biotyper_score <- function (filepath_samples, peaks_test, peaks_library, class_list, tolerance_ppm=2000, intensity_tolerance_percent=50, file_format="brukerflex", spectra_path_output=TRUE, replicates_average=FALSE, score_only=TRUE, standard_deviation_evaluation=FALSE, number_of_st_dev=1) {
number_of_samples <- length(peaks_test)
library_size <- length(peaks_library)
# Create the sample vector
sample_vector <- vector(length=number_of_samples)
if (file_format == "brukerflex") {
	if (replicates_average == FALSE) {
		for (s in 1:number_of_samples) {
			sample_vector[s] <- peaks_test[[s]]@metaData$sampleName
		}
	}
	# If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	if (replicates_average == TRUE) {
		for (s in 1:number_of_samples) {
			sample_vector[s] <- peaks_test[[s]]@metaData$sampleName[1]
		}
	}
}
if (file_format == "imzml" | file_format == "imzML") {
	for (s in 1:number_of_samples) {
		sample_vector[s] <- peaks_test[[s]]@metaData$file
	}
}
# Read the list of spectra folders in the sample mother folder
spectra_path_vector <- read_spectra_files (filepath_samples, file_format="brukerflex")
################# PEAK ALIGNMENT
# Create a global list for the alignment (with library and tests)
peaks_global <- peaks_library
peaks_global <- append(peaks_global, peaks_test)
# Alignment (R function)
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_global <- binPeaks(peaks_global, method="relaxed", tolerance=(tolerance_ppm/10^6))
# Empty the two original lists
peaks_library <- list()
peaks_test <- list()
# Re-fill the peaklists with the aligned peaks
for (i in 1:library_size) {
	peaks_library <- append(peaks_library, peaks_global[[i]])
}
for (i in (library_size+1):length(peaks_global)) {
	peaks_test <- append(peaks_test, peaks_global[[i]])
}
############################################################ SCORE (FRI)
if (standard_deviation_evaluation == FALSE) {
	# Compare the peaks in the single sample peaklists with the peaks in each peaklist in the library
	# COUNTER 1 - FIT
	# Create a counter, symmetrical to the database Peaklist
	counter1 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each sample
	for (s in 1:number_of_samples) {
		# For each peaklist in the Library
		for (l in 1:library_size) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_test[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_library[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
					counter1 [s,l] <- counter1 [s,l] + 1
					}
				}
			}
		}
	}
	for (s in 1:number_of_samples) {
		for (l in 1:library_size) {
			counter1[s,l] <- counter1[s,l] / length(peaks_test[[s]]@mass)
		}
	}
	# Compare the peaks in the library peaklists with the peaks in each sample
	# COUNTER 2 - RETRO FIT
	# Create a counter, symmetrical to the database Peaklist
	counter2 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each peaklist in the library
	for (l in 1:library_size) {
		# For each sample
		for (s in 1:number_of_samples) {
			# Scroll the peaks in the library
			for (j in 1:length(peaks_library[[l]]@mass)) {
				# Scroll the peaklist in the sample
				for (i in 1:length(peaks_test[[s]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
						counter2 [s,l] <- counter2 [s,l] + 1
					}
				}
			}
		}
	}
	for (l in 1:library_size) {
		for (s in 1:number_of_samples) {
			counter2[s,l] <- counter2[s,l] / length(peaks_library[[l]]@mass)
		}
	}
	# COUNTER 3
	# Symmetry -> comparison between intensities
	# Create a counter, symmetrical to the database Peaklist
	counter3 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each sample
	for (s in 1:number_of_samples) {
		# For each peaklist in the Library
		for (l in 1:library_size) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_test[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_library[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
						# Evaluate the difference in intensity
						if ((abs(peaks_test[[s]]@intensity[i] - peaks_library[[l]]@intensity[j])*100/peaks_library[[l]]@intensity[j]) < intensity_tolerance_percent) {
							counter3[s,l] <- counter3[s,l] + 1
						}
					}
				}
			}
		}
	}
	for (s in 1:number_of_samples) {
		for (l in 1:library_size) {
			counter3[s,l] <- counter3[s,l] / length(peaks_test[[s]]@mass)
		}
	}
	### Score calculation
	score <- matrix (0, nrow=number_of_samples, ncol=library_size)
	for (i in 1:number_of_samples) {
		for (j in 1:length(peaks_library)) {
			score[i,j] <- log10(counter1[i,j]*counter2[i,j]*counter3[i,j]*1000)
		}
	}
	#### Output the classification
	output <- matrix ("NO", nrow=number_of_samples, ncol=library_size)
	colnames(output) <- sort(class_list)
	if ((spectra_path_output == FALSE) || (replicates_average == TRUE)) {
		rownames(output) <- sample_vector
	}
	if ((spectra_path_output == TRUE) && (replicates_average == FALSE)) {
		rownames(output) <- spectra_path_vector
	}
	if (score_only == TRUE) {
		for (r in 1:number_of_samples) {
			for (w in 1:length(peaks_library)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste("YES","(", round(score[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste("NO", "(", round(score[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste("NI","(", round(score[r,w], digits=3), ")")
				}
			}
		}
	}
	if (score_only == FALSE) {
		for (r in 1:number_of_samples) {
			for (w in 1:length(peaks_library)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste("YES","(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste("NO", "(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste("NI","(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
				}
			}
		}
	}
}
####################################################### SCORE (FR-STDEV)
############################ peaks_library SHULD HAVE snr REPLACED BY SD!
if (standard_deviation_evaluation == TRUE) {
	# Compare the peaks in the single sample peaklists with the peaks in each peaklist in the library
	# COUNTER 1 - FIT
	# Create a counter, symmetrical to the database Peaklist
	counter1 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each sample
	for (s in 1:number_of_samples) {
		# For each peaklist in the Library
		for (l in 1:library_size) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_test[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_library[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
					counter1 [s,l] <- counter1 [s,l] + 1
					}
				}
			}
		}
	}
	for (s in 1:number_of_samples) {
		for (l in 1:library_size) {
			counter1[s,l] <- counter1[s,l] / length(peaks_test[[s]]@mass)
		}
	}
	# Compare the peaks in the library peaklists with the peaks in each sample
	# COUNTER 2 - RETRO FIT
	# Create a counter, symmetrical to the database Peaklist
	counter2 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each peaklist in the library
	for (l in 1:library_size) {
		# For each sample
		for (s in 1:number_of_samples) {
			# Scroll the peaks in the library
			for (j in 1:length(peaks_library[[l]]@mass)) {
				# Scroll the peaklist in the sample
				for (i in 1:length(peaks_test[[s]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
						counter2 [s,l] <- counter2 [s,l] + 1
					}
				}
			}
		}
	}
	for (l in 1:library_size) {
		for (s in 1:number_of_samples) {
			counter2[s,l] <- counter2[s,l] / length(peaks_library[[l]]@mass)
		}
	}
	# COUNTER 3
	# Symmetry -> comparison between intensities
	# Create a counter, symmetrical to the database Peaklist
	counter3 <- matrix (0, nrow=number_of_samples, ncol=library_size)
	# For each sample
	for (s in 1:number_of_samples) {
		# For each peaklist in the Library
		for (l in 1:library_size) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_test[[s]]@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_library[[l]]@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_test[[s]]@mass[i] == peaks_library[[l]]@mass[j]) {
						# Evaluate the difference in intensity
						if ((peaks_test[[s]]@intensity[i] <= (peaks_library[[l]]@intensity[j] + (number_of_st_dev * peaks_library[[l]]@snr[j]))) && (peaks_test[[s]]@intensity[i] >= (peaks_library[[l]]@intensity[j] - (number_of_st_dev * peaks_library[[l]]@snr[j])))) {
							counter3[s,l] <- counter3[s,l] + 1
						}
					}
				}
			}
		}
	}
	for (s in 1:number_of_samples) {
		for (l in 1:library_size) {
			counter3[s,l] <- counter3[s,l] / length(peaks_test[[s]]@mass)
		}
	}
	### Score calculation
	score <- matrix (0, nrow=number_of_samples, ncol=library_size)
	for (i in 1:number_of_samples) {
		for (j in 1:length(peaks_library)) {
			score[i,j] <- log10(counter1[i,j]*counter2[i,j]*counter3[i,j]*1000)
		}
	}
	#### Output the classification
	output <- matrix ("NO", nrow=number_of_samples, ncol=library_size)
	colnames(output) <- sort(class_list)
	if ((spectra_path_output == FALSE) || (replicates_average == TRUE)) {
		rownames(output) <- sample_vector
	}
	if ((spectra_path_output == TRUE) && (replicates_average == FALSE)) {
		rownames(output) <- spectra_path_vector
	}
	if (score_only == TRUE) {
		for (r in 1:number_of_samples) {
			for (w in 1:length(peaks_library)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste("YES","(", round(score[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste("NO", "(", round(score[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste("NI","(", round(score[r,w], digits=3), ")")
				}
			}
		}
	}
	if (score_only == FALSE) {
		for (r in 1:number_of_samples) {
			for (w in 1:length(peaks_library)) {
				if (score[r,w]>=2) {
					output[r,w] <- paste("YES","(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
				}
				if (score[r,w]<1.5) {
					output[r,w] <- paste("NO", "(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
				}
				if (score[r,w]>=1.5 && score[r,w]<2) {
					output[r,w] <- paste("NI","(Score:", round(score[r,w], digits=3), "), ","(Fit:", round(counter1[r,w], digits=3), "), ","(Retrofit:", round(counter2[r,w], digits=3), "), ","(Intensity:", round(counter3[r,w], digits=3), ")")
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
signal_follower_statistics <- function (filepath, signal_list, mass_labels=list(), snr=5, file_format="imzml") {
# Define the column and row headers (if mass labels are provided or not)
if (length(mass_labels) != 0) {
	column_names_st_dev <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_st_dev[i] <- paste(signal_list[i], "-", mass_labels[i], "- StDev")
	}
	column_names_coeff_var <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_coeff_var[i] <- paste(signal_list[i], "-", mass_labels[i], "- CV")
	}
	column_names_mean_abs_dev <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_mean_abs_dev[i] <- paste(signal_list[i], "-", mass_labels[i], "- MeanAbsDev")
	}
}
if (length(mass_labels) == 0) {
	column_names_st_dev <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_st_dev[i] <- paste(signal_list[i], "- StDev")
	}
	column_names_coeff_var <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_coeff_var[i] <- paste(signal_list[i], "- CV")
	}
	column_names_mean_abs_dev <- vector(length=length(signal_list))
	for (i in 1:length(signal_list)) {
		column_names_mean_abs_dev[i] <- paste(signal_list[i], "- MeanAbsDev")
	}
}
sample_names <- filepath
### Create the result matrices
st_dev_result_matrix <- matrix (ncol=length(signal_list), nrow=length(filepath))
colnames(st_dev_result_matrix) <- column_names_st_dev
rownames(st_dev_result_matrix) <- sample_names
coeff_var_result_matrix <- matrix (ncol=length(signal_list), nrow=length(filepath))
colnames(coeff_var_result_matrix) <- column_names_coeff_var
rownames(coeff_var_result_matrix) <- sample_names
mean_abs_dev_result_matrix <- matrix (ncol=length(signal_list), nrow=length(filepath))
colnames(mean_abs_dev_result_matrix) <- column_names_mean_abs_dev
rownames(mean_abs_dev_result_matrix) <- sample_names
# The script is run for every imzML file
for (lib in 1:length(filepath)) {
	# Spectra import
	spectra <- importImzMl(filepath[lib])
	### Smoothing (Savitzky-Golay filter, with window size 10, 21 points)
	spectra <- smoothIntensity (spectra, method="SavitzkyGolay")
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
	### Normalisation (TIC)
	spectra <- calibrateIntensity (spectra, method="TIC")
	### Peak Picking and alignment
	peaks <- detectPeaks (spectra, method="MAD", snr=snr)
	#### Create the partial result matrices
	st_dev_matrix <- matrix (ncol=length(signal_list), nrow=1)
	colnames(st_dev_matrix) <- column_names_st_dev
	rownames(st_dev_matrix) <- spectra[[1]]@metaData$file
	coeff_var_matrix <- matrix (0, ncol=length(signal_list), nrow=1)
	colnames(coeff_var_matrix) <- column_names_coeff_var
	rownames(coeff_var_matrix) <- spectra[[1]]@metaData$file
	mean_abs_dev_matrix <- matrix (0, ncol=length(signal_list), nrow=1)
	colnames(mean_abs_dev_matrix) <- column_names_mean_abs_dev
	rownames(mean_abs_dev_matrix) <- spectra[[1]]@metaData$file
	# For each signal in the signal_list
	for (s in 1:length(signal_list)) {
		intensity_vector <- vector(length=0)
		# For each peaklist
		for (p in 1:length(peaks)) {
			# Scroll the peaklist
			for (m in 1:length(peaks[[p]]@mass)) {
				# If there is a match
				if (abs(signal_list[s] - peaks[[p]]@mass[m]) <= 4) {
					# Put the intensity of that peak in the peaklist in a vector
					intensity_vector <- append(intensity_vector, peaks[[p]]@intensity[m])
				}
			}
		}
		# Calculate the standard deviation and the coefficient of variation
		st_dev_intensity <- sd(intensity_vector, na.rm=TRUE)
		mean_abs_dev_int <- mad(intensity_vector, na.rm=TRUE)
		mean_intensity <- mean(intensity_vector, na.rm=TRUE)
		median_intensity <- median(intensity_vector, na.rm=FALSE)
		coeff_var_intensity <- (st_dev_intensity / mean_intensity)*100
		# Fill in the partial result matrices with the values
		st_dev_matrix [1,s] <- paste(mean_intensity, "+/-", st_dev_intensity)
		coeff_var_matrix [1,s] <- paste(coeff_var_intensity, "%", "( mean =", mean_intensity, ")")
		mean_abs_dev_matrix [1,s] <- paste(median_intensity, "+/-", mean_abs_dev_int)
	}
	# Put the partial result matrices together in one matrix
	st_dev_result_matrix [lib,] <- st_dev_matrix
	coeff_var_result_matrix [lib,] <- coeff_var_matrix
	mean_abs_dev_result_matrix [lib,] <- mean_abs_dev_matrix
}
result_matrix <- st_dev_result_matrix
result_matrix <- cbind(result_matrix, coeff_var_result_matrix)
result_matrix <- cbind(result_matrix, mean_abs_dev_result_matrix)
return (result_matrix)
}





########################################################################





################# READ THE LIST OF imzML FILES FROM A FOLDER
read_spectra_files <- function (folder, file_format="imzml", full_path=FALSE) {
if (file_format=="imzml") {
	# Read all the files
	folder_all_files <- list.files (folder, full.names=full_path)
	# Create the empty vector in which only the imzML files will be listed
	spectra_files <- character()
	# Put the imzML files in the new vector discarding the ibd files
	for (l in 1:length(folder_all_files)) {
		if (length(grep(".imzML", folder_all_files [l], fixed=TRUE)) == 1) {
			spectra_files <- append(spectra_files, folder_all_files [l])
		}
	}
}
if (file_format=="brukerflex") {
	# Read all the folders
	sample_directory_list_all_folders <- dir (folder, ignore.case=TRUE, full.names=full_path, recursive=TRUE, include.dirs=TRUE)
	# Take only the fid paths (spectra)
	spectra_files <- character()
	for (i in 1:length(sample_directory_list_all_folders)) {
		if ((length(grep("fid", sample_directory_list_all_folders[i], ignore.case=TRUE)) == 1)) {
			spectra_files <- append(spectra_files, sample_directory_list_all_folders[i])
		}
	}
}
return (spectra_files)
}





########################################################################





################### AVERAGE THE REPLICATES (BY NAME)
average_replicates_by_name <- function (spectra, pattern_for_replicates="name_of_the_file", number_of_characters_for_similarity=10, file_format="brukerflex") {
## Put the filenames in a vector
# Create the empty vector
sample_vector <- vector(length=0)
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(spectra)) {
	sample_vector <- append(sample_vector, spectra[[i]]@metaData$sampleName)
}
if (pattern_for_replicates == "name_of_the_file") {
	for (s in 1: length(sample_vector)) {
		reference_sample <- sample_vector [s]
		for (j in s: length(sample_vector)) {
			# Determine the matching with each other sample
			matching <- pmatch(unlist(strsplit(reference_sample,"")), unlist(strsplit(sample_vector[j],"")))
			# If the FIRST characters match but the LAST character does not match, it is a replicate
			true_matching <- seq (from=1, to=number_of_characters_for_similarity)
			if ((length(intersect(true_matching, matching)) == length(true_matching)) && is.na (matching[length(matching)])) {
				sample_vector [j] <- reference_sample
			}
		}
	}
}
if ((pattern_for_replicates != "name_of_the_file") || (length(pattern_for_replicates) != 1)) {
	for (p in 1:length(pattern_for_replicates)) {
		for (s in 1: length(sample_vector)) {
			reference_sample <- sample_vector [s]
			for (j in s: length(sample_vector)) {
				# Determine the matching with each other sample
				matching <- pmatch(unlist(strsplit(reference_sample,"")), unlist(strsplit(sample_vector[j],"")))
				# If the FIRST characters match but the LAST character does not match, it is a replicate
				true_matching <- seq (from=1, to=number_of_characters_for_similarity)
				if ((length(intersect(true_matching, matching)) == length(true_matching)) && is.na (matching[length(matching)]) && length(grep(pattern_for_replicates[p], reference_sample, ignore.case=TRUE)) == 1 && length(grep(pattern_for_replicates[p], sample_vector [j], ignore.case=TRUE))) {
					sample_vector [j] <- reference_sample
				}
			}
		}
	}
}
spectra_replicates_averaged <- averageMassSpectra (spectra, labels=sample_vector)
return (spectra_replicates_averaged)
}





###########################################################################





########################################## PEAKLIST FILTERING
# Remove all the bad peaklists, keep only the best ones
filter_peaklist <- function (peaks, snr_filter=15, signal_number_threshold=15) {
if (snr_filter != 0 && signal_number_threshold != 0) {
	# Create the clusters
	cl <- makeCluster(cpu_core_number)
	############################################# FILTERING FUNCTION
	peak_filtering_function <- function (peaks, snr_filter, signal_number_threshold) {
		### Create a dataframe with the snr
		peaks_data_frame <- data.frame (snr = peaks@snr)
		### elect only the snr more than the threshold
		peaks_data_frame <- peaks_data_frame [peaks_data_frame$snr >= snr_filter,]
		### Select only the desired peaklists, based upon the threshold characteristics
			if (length(peaks_data_frame) < signal_number_threshold) {
				# Remove the bad spectrum from the list
				peaks <- NULL
			}
		return (peaks)
	}
	################################################################
	peaks_purified <- parLapply(cl, peaks, fun=function(peaks) peak_filtering_function(peaks, snr_filter=snr_filter, signal_number_threshold=signal_number_threshold))
	### Keep only the elements that are different from NULL
	peaks_purified <- peaks_purified [!sapply(peaks_purified, is.null)]
	# Stop the clusters
	stopCluster(cl)
	return (peaks_purified)
}
}





########################################################################





#################################### SPECTRA GROUPING (PEAKLIST MERGING)
merge_peaklist <- function (peaks, file_format="imzml", tof_mode="linear") {
### Alignment
if (tof_mode == "linear") {
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
}
if (tof_mode == "reflectron" || tof_mode == "reflector") {
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
}
if (file_format == "imzml" | file_format == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		file_vector <- append(file_vector, peaks[[i]]@metaData$file)
	}
	peaks_merged <- mergeMassPeaks (peaks, labels=file_vector, method="mean", ignore.na=FALSE)
}
if (file_format == "brukerflex") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
	}
	peaks_merged <- mergeMassPeaks (peaks, labels=file_vector, method="mean", ignore.na=FALSE)
}
return (peaks_merged)
}





########################################################################





############## SPECTRA GROUPING CLASS (PEAKLIST MERGING)
merge_peaklist_class <- function (peaks, class_list, file_format="imzml", tof_mode="linear") {
### Alignment
if (tof_mode == "linear") {
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.002)
}
if (tof_mode == "reflectron" || tof_mode == "reflector") {
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
peaks <- binPeaks(peaks, method="relaxed", tolerance=0.0002)
}
if (file_format == "imzml" | file_format == "imzML") {
	# Replace the name with the class name
	peaks <- replace_class_name (peaks, class_list=class_list, file_format="imzml")
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		file_vector <- append(file_vector, peaks[[i]]@metaData$file)
	}
	peaks_merged <- mergeMassPeaks (peaks, labels=file_vector, method="mean", ignore.na=FALSE)
}
if (file_format == "brukerflex") {
	# Replace the name with the class name
	peaks <- replace_class_name (peaks, class_list=class_list, file_format="brukerflex")
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(peaks)) {
		file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
	}
	peaks_merged <- mergeMassPeaks (peaks, labels=file_vector, method="mean", ignore.na=FALSE)
}
return (peaks_merged)
}






########################################################################





######################### PEAK INTENSITY PLOT
peak_intensity_plot <- function (peaks, peak_of_interest, class_list=list(), tolerance_ppm=2000) {
if (length(class_list) > 2) {
	print("This algorithm can be computed only with no classes or one class")
}
# If there are one or no classes, compute the plot on the entire dataset
if (length(class_list) == 0 || length(class_list) == 1) {
	intensity_vector <- numeric()
	# Scroll the peaklists and add the corresponding intensities (of the peak of interese) to a vector
	for (s in 1:length(peaks)) {
		for (p in 1:length(peaks[[s]]@mass)) {
			if (abs(peaks[[s]]@mass[p] - peak_of_interest)*10^6/peak_of_interest <= tolerance_ppm) {
				intensity_vector <- append(intensity_vector, peaks[[s]]@intensity[p])
			}
		}
	}
	# Plot the results
	peak_intensity_plot <- plot(intensity_vector)
}
# If there are one or no classes, compute the plot on the entire dataset
if (length(class_list) == 2) {
	class_names <- class_list
	intensity_vector <- list(class1=numeric(), class2=numeric())
	for (w in 1:length(class_list)) {
		# Scroll the peaklists and add the corresponding intensities (of the peak of interese) to a vector
		for (s in 1:length(peaks)) {
			# Check if the peaklist belongs to the class
			if (length(grep(class_list[w], peaks[[s]]@metaData$file)) == 1) {
				for (p in 1:length(peaks[[s]]@mass)) {
					if (abs(peaks[[s]]@mass[p] - peak_of_interest)*10^6/peak_of_interest <= tolerance_ppm) {
						intensity_vector[[w]] <- append(intensity_vector[[w]], peaks[[s]]@intensity[p])
					}
				}
			}
		}
	}
	# Plot the results
	peak_intensity_plot <- plot(x=intensity_vector[[1]], xlab=class_list[1], y=intensity_vector[[2]], ylab=class_list[2])
}
return (peak_intensity_plot)
}





########################################################################





##################### AVERAGE THE REPLICATES (BY FOLDER)
average_replicates_by_folder <- function (spectra, folder, file_format="brukerflex") {
# List the spectra files
folder_files <- read_spectra_files (folder, file_format=file_format, full_path=TRUE)
# Split the path into individual folders
folder_files_splitted <- list()
# Split the paths into folders
for (f in 1:length(folder_files)) {
	folder_files_splitted[f] <- strsplit(folder_files[f], "/")
}
# Check for replicates in the same folder
for (f in 1:(length(folder_files_splitted)-1)) {
	# If the spectrum folder name is different between two samples but the folder above is the same, they are replicates (based on the brukerflex folder structure)
	if ((folder_files_splitted[[f]] [length(folder_files_splitted[[f]]) - 4] != folder_files_splitted[[f+1]] [length(folder_files_splitted[[f+1]]) - 4]) && (folder_files_splitted[[f]] [length(folder_files_splitted[[f]]) - 5] == folder_files_splitted[[f+1]] [length(folder_files_splitted[[f+1]]) - 5])) {
		folder_files_splitted[[f+1]] [length(folder_files_splitted[[f+1]]) - 4] <- folder_files_splitted[[f]] [length(folder_files_splitted[[f]]) - 4]
	}
}
# Recreate the sample Vector
sample_vector <- character(length=length(folder_files_splitted))
for (f in 1:length(folder_files_splitted)) {
	sample_vector[f] <- folder_files_splitted[[f]] [length(folder_files_splitted[[f]]) -4]
}
spectra_replicates_averaged <- averageMassSpectra (spectra, labels=sample_vector)
return (spectra_replicates_averaged)
}





########################################################################





######################################## LIBRARY CREATION (for Biotyper)
library_creation <- function (filepath_library, class_list=list(), class_grouping=TRUE, mass_range=c(3000,15000), replicates_average=FALSE, patients_average=FALSE, snr=5, most_intense_peaks=TRUE, signals_to_take=20, standard_deviation=FALSE, coeff_of_var=FALSE, low_intensity_peak_removal=FALSE, intensity_threshold_percent=0.1, tof_mode="linear", file_format="brukerflex") {
if (tof_mode == "linear") {
	tolerance_ppm <- 2000
}
if (tof_mode == "reflectron" || tof_mode == "reflector") {
	tolerance_ppm <- 200
}
if (file_format == "brukerflex") {
	### Load the spectra
	spectra <- importBrukerFlex(filepath_library)
}
if (file_format == "imzml" | file_format == "imzML") {
	### Load the spectra
	spectra <- importImzMl(filepath_library)
}
### Truncation
spectra <- trim(spectra, range = mass_range)
### Preprocessing
spectra <- preprocess_spectra(spectra, tof_mode=tof_mode)
### Average the replicates
if (replicates_average == TRUE) {
	spectra <- average_replicates_by_folder (spectra, filepath_library, file_format=file_format)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
if (patients_average == TRUE) {
	spectra <- group_spectra (spectra, spectra_per_patient=1, file_format=file_format, tof_mode=tof_mode)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
### Peak picking on the individual spectra
if (class_grouping == TRUE) {
	### Replace the sample name with the class name
	spectra <- replace_class_name (spectra, class_list=class_list, file_format=file_format)
	### Spectra grouping (class)
	spectra <- group_spectra_class (spectra, class_list, file_format=file_format)
	### Baseline correction
	spectra <- removeBaseline (spectra, method="TopHat")
}
if (class_grouping == FALSE) {
	spectra <- spectra
}
##########################
### Peak picking
if (most_intense_peaks == TRUE) {
	peaks_library <- most_intense_signals (spectra, signals_to_take=signals_to_take)
}
if (most_intense_peaks == FALSE) {
	peaks_library <- detectPeaks (spectra, method="MAD", snr=snr)
}
###########################
if (standard_deviation == TRUE && coeff_of_var == FALSE) {
	# For each peaklist in the library, compute the SD
	for (j in 1:length(peaks_library)) {
		peaks_library[[j]] <- replace_snr_with_st_dev_in_peaklist (peaks_library[[j]], peaks, tolerance_ppm=tolerance_ppm, file_format=file_format)
	}
}
if (standard_deviation == FALSE && coeff_of_var == TRUE) {
	# For each peaklist in the library, compute the CV
	for (j in 1:length(peaks_library)) {
		peaks_library[[j]] <- replace_snr_with_cv_in_peaklist (peaks_library[[j]], peaks, tolerance_ppm=tolerance_ppm, file_format=file_format)
	}
}
if (standard_deviation == TRUE && coeff_of_var == TRUE) {
	# For each peaklist in the library, compute the SD
	for (j in 1:length(peaks_library)) {
		peaks_library[[j]] <- replace_snr_with_st_dev_in_peaklist (peaks_library[[j]], peaks, tolerance_ppm=tolerance_ppm, file_format=file_format)
	}
}
############
if (low_intensity_peak_removal == TRUE) {
	peaks_library <- remove_low_intensity_peaks (peaks_library, intensity_threshold_percent=intensity_threshold_percent)
}
####
library_list <- list (spectra = spectra, peaks = peaks_library)
return (library_list)
}





########################################################################





######################## PEAK ALIGNMENT
align_and_filter_peaks <- function (peaks, tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1) {
peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
if (peaks_filtering == TRUE) {
	peaks_aligned <- filterPeaks(peaks_aligned, minFrequency=frequency_threshold)
}
if (low_intensity_peaks_removal == TRUE) {
	peaks_aligned <- remove_low_intensity_peaks (peaks_aligned, intensity_threshold_percent=intensity_threshold_percent)
}
return (peaks_aligned)
}





########################################################################





############################ PEAK FILTERING
discard_peaks <- function (peaks, list_of_peaks=list(), mass_type="monoisotopic", cluster_size=4, tof_mode="reflectron") {
if (mass_type == "monoisotopic") {
	# For each peak in the peaklist of interest
	for (p in 1:length(list_of_peaks)) {
		# Add the peaks of the cluster
		for (i in 1:(cluster_size-1)) {
			list_of_peaks <- append(list_of_peaks, list_of_peaks[p]+i)
		}
	}
}
if (mass_type == "average") {
	# For each peak in the peaklist of interest
	for (p in 1:length(list_of_peaks)) {
		# Add the peaks of the cluster
		for (i in 1:(cluster_size/2)) {
			list_of_peaks <- append(list_of_peaks, list_of_peaks[p]-i)
			list_of_peaks <- append(list_of_peaks, list_of_peaks[p]+i)
		}
	}
}
if (tof_mode == "linear") {
	tolerance_ppm <- 2000
}
if (tof_mode == "reflectron" || tof_mode == "reflector") {
	tolerance_ppm <- 1000
}
###################################################### ERASING FUNCTION
peaks_erasing_function <- function (peaks, list_of_peaks) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, snr = peaks@snr)
	# Select the peaks of interest and remove them from the peaklist
	for (p in 1: length(list_of_peaks)) {
		for (s in 1:length(peaks_data_frame$mass)) {
			if (abs(peaks_data_frame$mass[s] - list_of_peaks[p])*10^6/list_of_peaks[p] <= tolerance_ppm) {
				peaks_data_frame$mass[s] <- 0
			}
		}
	}
	# Keep only the signals that are not zero
	peaks_data_frame <- peaks_data_frame [peaks_data_frame$mass != 0, ]
	# Put these signals back into the peaklist
	peaks@mass <- peaks_data_frame$mass
	peaks@intensity <- peaks_data_frame$intensity
	peaks@snr <- peaks_data_frame$snr
	return (peaks)
}
#######################################################################
if (length(list_of_peaks) == 0) {
	peaks <- peaks
}
if (length(list_of_peaks) > 0) {
	cl <- makeCluster(cpu_core_number)
	peaks <- parLapply(cl, peaks, fun=function(peaks) peaks_erasing_function (peaks, list_of_peaks))
	stopCluster(cl)
}
return (peaks)
}





########################################################################





########################## FILTER SPECTRA
filter_spectra_noise <- function (spectra, snr_filter=20, percentage_above_threshold=10) {
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
################################################## NOISY SPECTRA REMOVAL
noisy_spectra_removal_function <- function (spectra, snr_filter, percentage_above_threshold) {
	# Peak picking
	peaks <- detectPeaks (spectra, method="MAD", snr=3)
	# Keep only the spectra with a certain percentage of signals above the threshold
	peaks_new <- detectPeaks (spectra, method="MAD", snr=snr_filter)
	if ((length(peaks_new@mass) / length(peaks@mass))*100 < percentage_above_threshold) {
		spectra <- NULL
	}
	return (spectra)
}
########################################################################
spectra <- parLapply(cl, spectra, fun=function (spectra) noisy_spectra_removal_function (spectra, snr_filter, percentage_above_threshold))
stopCluster(cl)
### Keep only the elements that are different from NULL
spectra <- spectra [!sapply(spectra, is.null)]
return (spectra)
}





########################################################################





######################### REPLACE THE snr WITH THE STDEV IN THE PEAKS
### Compute the standard deviation of each peak of a average spectrum peaklist,
# by replacing the existing snr slot with the SD (each peak of the average peaklist
# is searched across the dataset
replace_snr_with_st_dev_in_peaklist <- function (peaks_average, peaks, tolerance_ppm=2000, file_format="imzml") {
	# Scroll the peaks of the average peaklist
	for (p in 1:length(peaks_average@mass)) {
		# Create an empty vector where to allocate the intensity of this peak into the dataset
		intensity_vector <- numeric()
		# Pass every peaklist of the dataset (one for each spectrum)
		for (s in 1:length(peaks)) {
			# Do the operation only for the peaklists that matches the AVG
			if (file_format == "brukerflex") {
				if (length(grep(peaks_average@metaData$sampleName, peaks[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaks_average@mass[p])*10^6/peaks_average@mass[p]) <= tolerance_ppm) {
							# Record its intensity into a vector
							intensity_vector <- append(intensity_vector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
			if (file_format == "imzml" | file_format == "imzML") {
				if (length(grep(peaks_average@metaData$file, peaks[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaks_average@mass[p])*10^6/peaks_average@mass[p]) <= tolerance_ppm) {
							# Record its intensity into a vector
							intensity_vector <- append(intensity_vector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
		}
		# Calculate the standard deviation of the peak (whose intensities are in the vector)
		st_dev <- sd(intensity_vector)
		# Replace the snr slot with the st_dev (for the peak)
		peaks_average@snr[p] <- st_dev
	}
return (peaks_average)
}




################################################################################





######################### REPLACE THE snr WITH THE CV IN THE PEAKS
### Compute the coefficient of variation of each peak of a average spectrum peaklist,
# by replacing the existing snr slot with the CV (each peak of the average peaklist
# is searched across the dataset
replace_snr_with_cv_in_peaklist <- function (peaks_average, peaks, tolerance_ppm=2000, file_format="imzml") {
	# Scroll the peaks of the average peaklist
	for (p in 1:length(peaks_average@mass)) {
		# Create an empty vector where to allocate the intensity of this peak into the dataset
		intensity_vector <- numeric()
		# Pass every peaklist of the dataset (one for each spectrum)
		for (s in 1:length(peaks)) {
			# Do the operation only for the peaklists that matches the AVG
			if (file_format == "brukerflex") {
				if (length(grep(peaks_average@metaData$sampleName, peaks[[s]]@metaData$sampleName, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaks_average@mass[p])*10^6/peaks_average@mass[p]) <= tolerance_ppm) {
							# Record its intensity into a vector
							intensity_vector <- append(intensity_vector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
			if (file_format == "imzml" | file_format == "imzML") {
				if (length(grep(peaks_average@metaData$file, peaks[[s]]@metaData$file, ignore.case=TRUE)) == 1) {
					# Scroll the peaklist
					for (m in 1:length(peaks[[s]]@mass)) {
						# Identify the peak
						if ((abs(peaks[[s]]@mass[m] - peaks_average@mass[p])*10^6/peaks_average@mass[p]) <= tolerance_ppm) {
							# Record its intensity into a vector
							intensity_vector <- append(intensity_vector, peaks[[s]]@intensity[m])
						}
					}
				}
			}
		}
		# Calculate the standard deviation and the mean of the peak (whose intensities are in the vector)
		st_dev <- sd(intensity_vector)
		mean_intensity <- mean(intensity_vector)
		# Calculate the coefficient of variation
		coeff_var_intensity <- (st_dev / mean_intensity) * 100
		# Replace the snr slot with the CV (for the peak)
		peaks_average@snr[p] <- coeff_var_intensity
	}
return (peaks_average)
}





################################################################################





######## RICHEST PEAKLIST AS REFERENCE FOR THE SPECTRA ALIGNMENT/WARPING
richest_peaklist_for_alignment <- function (spectra, snr=5, tolerance_ppm=2000) {
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", snr=snr)
# Define the indexes
index <- 0
max_length <- 0
# Scroll the peaklist list
for (i in 1:length(peaks)) {
# Determine the longest peaklist and record its position
	if (length(peaks[[i]]@mass) >= max_length) {
	max_length <- length(peaks[[i]]@mass)
	index <- i
	}
}
# Create the vector containing the richest peaklist
reference_peaklist <- peaks[[index]]
# SPECTRA WARPING / ALIGNMENT
spectra <- align_spectra (spectra, reference=reference_peaklist, tolerance=(tolerance_ppm/10^6), noiseMethod="MAD", snr=snr, warpingMethod="quadratic")
return (spectra)
}





################################################################################





############################# MATRIX ZERO FILLING
# Replace the NA values with a zero
matrix_zero_filling <- function (input_matrix) {
cl <- makeCluster(cpu_core_number)
######################################################## NA REPLACEMENT FUNCTION
na_replacement_function <- function (x) {
	if (is.na(x) == TRUE) {
		x <- 0
	}
return (x)
}
################################################################################
input_matrix <- parApply(cl, input_matrix, c(1,2), FUN=function(input_matrix) na_replacement_function (input_matrix))
stopCluster(cl)
return (input_matrix)
}





###############################################################################





######################################### ADD THE THY CLASS TO THE MATRIX (THYROID)
matrix_add_thy <- function (signal_matrix, peaks, file_format="imzml") {
number_of_spectra <- length(peaks)
# Create the empty vector
file_vector <- character()
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(peaks)) {
	if (file_format == "imzml" | file_format == "imzML") {
		file_vector <- append(file_vector, peaks[[i]]@metaData$file)
	}
	if (file_format == "brukerflex") {
		file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
	}
}
# Create the thy vector, symmetrical to the file vector
thy_vector <- file_vector
# Find the "THY" in the file name and store its value
for (t in 1:length(thy_vector)) {
	############## THY
	# THY1
	if (length(grep("THY1", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 1
	}
	# THY2
	if (length(grep("THY2", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 2
	}
	# THY3
	if (length(grep("THY3", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 3
	}
	# THY4
	if (length(grep("THY4", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 4
	}
	# THY5
	if (length(grep("THY5", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 5
	}
	############## TIR
	# TIR1
	if (length(grep("TIR1", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 1
	}
	# TIR2
	if (length(grep("TIR2", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 2
	}
	# TIR3
	if (length(grep("TIR3", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 3
	}
	# TIR4
	if (length(grep("TIR4", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 4
	}
	# TIR5
	if (length(grep("TIR5", thy_vector[t], ignore.case=TRUE)) == 1) {
		thy_vector[t] <- 5
	}
}
# Create the sample matrix column and appendit to the global matrix
thy_column <- matrix (0, ncol=1, nrow=number_of_spectra)
colnames(thy_column) <- "THY"
# Fill in the matrix thy column with the thy_vector and attach it to the matrix
thy_column [,1] <- cbind(thy_vector)
signal_matrix <- cbind(signal_matrix, thy_column)
#
return (signal_matrix)
}





################################################################################





############################## PEAK STATISTICS (on processed Spectra)
peak_statistics <- function (spectra, snr=3, class_list=list(), tof_mode="linear", file_format="imzml", tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25, remove_outliers=TRUE) {
# Determine the number of classes
if (length(class_list) == 0 | length(class_list) == 1) {
number_of_classes <- 1
}
if (length(class_list) > 1) {
	number_of_classes <- length(class_list)
}
# Detect and Align Peaks
if (tof_mode == "linear") {
	peaks <- detectPeaks (spectra, method="MAD", snr=snr, halfWindowSize=20)
}
if (tof_mode == "reflector" | tof_mode == "reflectron") {
	peaks <- detectPeaks (spectra, method="MAD", snr=snr, halfWindowSize=20)
}
peaks <- align_and_filter_peaks (peaks, tolerance_ppm=tolerance_ppm, peaks_filtering=peaks_filtering, frequency_threshold=frequency_threshold)
# Generate the matrix (and convert it into a data frame)
signal_matrix <- intensityMatrix (peaks, spectra)
# Peak vector
#peak_vector <- as.numeric(names(signal_matrix))
############################################################## ONE CLASS
if (number_of_classes == 1) {
	# Fix the signal_matrix (Add the sample column)
	signal_matrix <- matrix_add_class_and_sample (signal_matrix, peaks=peaks,
		file_format=file_format, sample_output=TRUE, class_output=FALSE)
	# Output matrix
	peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-1), ncol=9)
	rownames(peak_stat_matrix) <- as.numeric(colnames(signal_matrix)[1:(ncol(signal_matrix)-1)])
	colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter", "Sample")
	# For each peak
	for (p in 1:(ncol(signal_matrix)-1)) {
		intensity_vector <- as.numeric(signal_matrix[,p])
		if (remove_outliers == TRUE) {
			intensity_vector <- outliers_removal (intensity_vector)
			intensity_vector <- intensity_vector$vector
		}
		# Calculate the statistical parameters on the intensity values in the vector
		# Normality
		if (length(intensity_vector) >= 3 & length(intensity_vector) <= 5000) {
		shapiro_test <- shapiro.test (intensity_vector)
			if (shapiro_test$p.value < 0.05) {
			distribution_type <- "Non-normal"
			}
			if (shapiro_test$p.value >= 0.05) {
			distribution_type <- "Normal"
			}
		}
		if (length(intensity_vector) < 3) {
		distribution_type <- "Not determinable, number of samples too low"
		}
		if (length(intensity_vector) > 5000) {
		distribution_type <- "Number of samples too high, assume it is normal"
		}
		# Other parameters
		spectra_names <- levels (as.factor(signal_matrix[,ncol(signal_matrix)]))
		st_dev_intensityensity <- sd(intensity_vector)
		summary_intensity_vector <- summary (intensity_vector)
		mean_intensity <- summary_intensity_vector [4]
		coeff_variation <- (st_dev_intensityensity / mean_intensity) *100
		median_intensityensity <- summary_intensity_vector [3]
		first_quartile <- summary_intensity_vector [2]
		third_quartile <- summary_intensity_vector [5]
		inter_quartile_range <- third_quartile - first_quartile
		spectra_counter <- length(intensity_vector)
		# Fill the matrix with the values
		peak_stat_matrix [p,1] <- distribution_type
		peak_stat_matrix [p,2] <- mean_intensity
		peak_stat_matrix [p,3] <- st_dev_intensityensity
		peak_stat_matrix [p,4] <- coeff_variation
		peak_stat_matrix [p,5] <- median_intensityensity
		peak_stat_matrix [p,6] <- inter_quartile_range
		peak_stat_matrix [p,7] <- paste("1st quartile", first_quartile, "; 3rd quartile", third_quartile)
		peak_stat_matrix [p,8] <- spectra_counter
		peak_stat_matrix [p,9] <- spectra_names
	}
}
############################################################ TWO OR MORE CLASSES
# Every variable now is a list, each element of which corresponds to a certain value from a class
# So every variable is a list with the same lengthof the class list (each element of the list
# is referred to a class
if (number_of_classes > 1) {
	# Fix the signal_matrix (Add the sample column)
	signal_matrix <- matrix_add_class_and_sample (signal_matrix, peaks=peaks, class_list=class_list,
		file_format=file_format, sample_output=TRUE, class_output=TRUE)
	# Output matrix
	peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-2), ncol=14)
	rownames(peak_stat_matrix) <- as.numeric(names(signal_matrix)[1:(ncol(signal_matrix)-2)])
	colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %",
		"Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)",
		"t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
	# For each peak
	for (p in 1:(ncol(signal_matrix)-2)) {
		# Put the intensity of that peak into one vector per class (in a global list)
		intensity_vector <- list()
		# Scroll the peaklists and Add the peak intensity to a vector(one for each class)
		for (l in 1:length(class_list)) {
			# Allocate in the intensity vector the rows for that peak belonging to the certain class
			intensity_vector [[l]] <- as.numeric(signal_matrix [signal_matrix[,ncol(signal_matrix)] == class_list[l],p])
		}
		if (remove_outliers == TRUE) {
			for (i in 1:length(intensity_vector)) {
				intensity_vector[[l]] <- outliers_removal (intensity_vector[[l]])
				intensity_vector[[l]] <- intensity_vector[[l]]$vector
			}
		}
		######################## STATISTICAL PARAMETERS
		############################################### Normality for each class
		shapiro_test <- list()
		distribution_type <- list()
		for (l in 1:length(class_list)) {
			if (length(intensity_vector[[l]]) >= 3 && length(intensity_vector[[l]]) <=5000) {
				shapiro_test[[l]] <- shapiro.test (intensity_vector[[l]])
				if (shapiro_test[[l]]$p.value < 0.05) {
				distribution_type[[l]] <- "Non-normal"
				}
				if (shapiro_test[[l]]$p.value >= 0.05) {
				distribution_type[[l]] <- "Normal"
				}
			}
			if (length(intensity_vector[[l]]) < 3) {
			distribution_type[[l]] <- "Not determinable, number of samples too low"
			}
			if (length(intensity_vector) > 5000) {
			distribution_type[[l]] <- "Number of samples too high, assume it is normal"
			}
		}
		##################################################### Homoscedasticity
		if (length(class_list) == 2) {
			variance_test_parametric <- var.test (intensity_vector[[1]], intensity_vector[[2]])
		}
		if (length(class_list) >= 2) {
			variance_test_non_parametric <- bartlett.test (as.numeric(signal_matrix[,p]), g=as.factor(signal_matrix[,ncol(signal_matrix)]))
		}
		########################################### Other parameters (per class)
		st_dev_intensityensity <- list()
		summary_intensity_vector <- list()
		mean_intensity <- list()
		coeff_variation <- list()
		median_intensityensity <- list()
		first_quartile <- list()
		third_quartile <- list()
		inter_quartile_range <- list()
		spectra_counter <- list()
		variance <- list()
		for (l in 1:length(class_list)) {
			st_dev_intensityensity[[l]] <- sd(intensity_vector[[l]])
			summary_intensity_vector [[l]] <- summary (intensity_vector[[l]])
			mean_intensity[[l]] <- summary_intensity_vector[[l]] [4]
			coeff_variation[[l]] <- (st_dev_intensityensity[[l]] / mean_intensity[[l]]) *100
			median_intensityensity[[l]] <- summary_intensity_vector[[l]] [3]
			first_quartile[[l]] <- summary_intensity_vector[[l]] [2]
			third_quartile[[l]] <- summary_intensity_vector[[l]] [5]
			inter_quartile_range[[l]] <- third_quartile[[l]] - first_quartile[[l]]
			spectra_counter[[l]] <- length(intensity_vector[[l]])
			variance[[l]] <- var (intensity_vector[[l]])
		}
		############################################# Parameters between classes
		# T-test
		if (length(class_list) == 2) {
			t_test <- t.test (intensity_vector[[1]], intensity_vector[[2]])
		}
		# ANOVA TEST
		if (length(class_list) >= 2) {
		anova_test <- aov (signal_matrix[,p] ~ signal_matrix[,ncol(signal_matrix)])
		}
		# WILCOXON - MANN-WHITNEY TEST
		if (length(class_list) == 2) {
			wilcoxon_test <- wilcox.test (intensity_vector[[1]], intensity_vector[[2]])
		}
		# KRUSKAL-WALLIS TEST
		if (length(class_list) >= 2) {
			kruskal_wallis_test <- kruskal.test (signal_matrix[,p], g=as.factor(signal_matrix[,ncol(signal_matrix)]))
		}
		######################################## Fill the matrix with the values
		# Distribution Type
		distribution_type_name <- character()
		for (l in length(class_list)) {
			distribution_type_name <- paste(distribution_type_name, " ", distribution_type[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,1] <- paste(distribution_type_name)
		# Mean
		mean_intensity_name <- character()
		for (l in length(class_list)) {
			mean_intensity_name <- paste(mean_intensity_name, " ", mean_intensity[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,2] <- mean_intensity_name
		# Standard Deviation
		st_dev_name <- character()
		for (l in length(class_list)) {
			st_dev_name <- paste(st_dev_name, " ", st_dev_intensityensity[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,3] <- st_dev_name
		# Coefficient of Variation
		coeff_variation_name <- character()
		for (l in length(class_list)) {
			coeff_variation_name <- paste(coeff_variation_name, " ", coeff_variation[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,4] <- coeff_variation_name
		# Median
		median_intensityensity_name <- character()
		for (l in length(class_list)) {
			median_intensityensity_name <- paste(median_intensityensity_name, " ", median_intensityensity[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,5] <- median_intensityensity_name
		# Interquartile Range (IQR)
		inter_quartile_range_name <- character()
		for (l in length(class_list)) {
			inter_quartile_range_name <- paste(inter_quartile_range_name, " ", inter_quartile_range[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,6] <- inter_quartile_range_name
		# Spectra counter
		spectra_counter_name <- character()
		for (l in length(class_list)) {
			spectra_counter_name <- paste(inter_quartile_range_name, " ", spectra_counter[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,7] <- spectra_counter_name
		# Class
		class_name <- character()
		for (l in length(class_list)) {
			class_name <- paste(class_name, " ", class_list[[l]], " - ", class_list[l], sep="")
		}
		peak_stat_matrix [p,8] <- class_name
		# Homoscedasticity (Parametric)
		if (variance_test_parametric$p.value < 0.05) {
		homoscedasticity_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
		}
		if (variance_test_parametric$p.value >= 0.05) {
		homoscedasticity_parametric <- paste("Homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
		}
		if (variance_test_non_parametric$p.value < 0.05) {
		homoscedasticity_non_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_non_parametric$p.value, ")")
		}
		if (variance_test_non_parametric$p.value >= 0.05) {
		homoscedasticity_non_parametric <- paste("Homoscedastic data", "(p-value:", variance_test_non_parametric$p.value, ")")
		}
		peak_stat_matrix [p,9] <- homoscedasticity_parametric
		peak_stat_matrix [p,10] <- homoscedasticity_non_parametric
		# t-Test
		peak_stat_matrix [p,11] <- t_test$p.value
		# ANOVA
		peak_stat_matrix [p,12] <- summary(anova_test)[[1]]$"Pr(>F)"[1]
		# Wilcoxon / Mann-Whitney test
		peak_stat_matrix [p,13] <- wilcoxon_test$p.value
		# Kruskal-Wallis test
		peak_stat_matrix [p,14] <- kruskal_wallis_test$p.value
	}
}
return (peak_stat_matrix)
}



################################################################################





########################## REMOVE LOW INTENSITY PEAKS
remove_low_intensity_peaks <- function (peaks, intensity_threshold_percent=0.1) {
cl <- makeCluster(cpu_core_number)
########################################### INTENSITY FILTERING FUNCTION
intensity_filtering_function <- function (peaks, intensity_threshold_percent) {
	# Filter out the peaks whose intensity is below a certain threshold
	for (p in 1:length(peaks@mass)) {
		if ((abs(peaks@intensity[p] - max(peaks@intensity,na.rm=TRUE))*100/max(peaks@intensity,na.rm=TRUE)) < intensity_threshold_percent) {
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
peaks_filtered <- parLapply(cl, peaks, fun = function (peaks) intensity_filtering_function (peaks, intensity_threshold_percent=intensity_threshold_percent))
return (peaks_filtered)
stopCluster(cl)
}





################################################################################





###################### SPECTRA PURIFICATION BASED ON THE TIC (Before processing)
spectra_tic_purification <- function (spectra, tof_mode="reflectron", absolute_tic_threshold=0) {
# If not specified, set it based on the TOF mode
if (absolute_tic_threshold == 0) {
	if (tof_mode == "linear") {
		absolute_tic_threshold <- 2000
	}
	if (tof_mode == "reflectron" | tof_mode == "reflector") {
		absolute_tic_threshold <- 200
	}
}	else {absolute_tic_threshold <- absolute_tic_threshold}
# Before preprocessing (and thus normalisation), evaluate the TIC
tic_vector <- numeric()
for (s in 1:length(spectra)) {
	# Calculate the spectrum TIC
	tic_spectrum <- sum(spectra[[s]]@intensity)
	# Add this value to the global TIC vector
	tic_vector <- append(tic_vector, tic_spectrum)
}
# Remove the spectra with an absolute TIC value of less than ...
low_tic_position <- which(tic_vector <= absolute_tic_threshold)
if (length(low_tic_position) > 0) {
	# Remove the bad spectra both from the spectra list and from the tic_vector
	spectra <- spectra [-low_tic_position]
	tic_vector <- tic_vector [-low_tic_position]
}	else {spectra <- spectra}
###### Outliers detection
# Output the summary of the vector(quartiles, mean, median, max and min)
summary_tic_vector <- summary (tic_vector)
# Calculate the interquartile range
inter_quartile_range <- summary_tic_vector[5] - summary_tic_vector[2]
# Calculate the fences, beyond which the spectrum is an outlier
iqr_fences <- c((summary_tic_vector[2] - 1.5*inter_quartile_range), (summary_tic_vector[5] + 1.5*inter_quartile_range))
# Find the outliers based on the fences condition
outliers_position <- which(tic_vector < iqr_fences[1] | tic_vector > iqr_fences[2])
if (length(outliers_position) > 0) {
	# Remove the correspondent spectra from the dataset
	tic_vector <- tic_vector [-outliers_position]
	spectra <- spectra [-outliers_position]
}	else {spectra <- spectra}
#
return (spectra)
}
################################################################################





################################# REMOVE UNDIGESTED SPECTRA
remove_undigested_spectra <- function (spectra, mass_range_of_interest=c(0,2500), percentage_of_spectrum_in_the_range=75) {
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
################################################# UNDIGESTED FILTERING FUNCTION
undigested_filtering_function <- function (spectra, mass_range_of_interest, percentage_of_spectrum_in_the_range) {
	## Peak picking
	peaks <- detectPeaks (spectra, method="MAD", snr=3, halfWindowSize=20)
	# Define the peaks of interest based on the given interval
	peaks_of_interest <- peaks@mass [peaks@mass > mass_range_of_interest[1] & peaks@mass < mass_range_of_interest[2]]
	# Calculate the percentage of the spectrum that is within the given range
	percentage_in_range <- length(peaks_of_interest)/length(peaks@mass) * 100
	# The spectrum is to be discarded if not of interest
	if (percentage_in_range < percentage_of_spectrum_in_the_range) {
		spectra <- NULL
	}
	return (spectra)
}
########################################################################
spectra_digested <- parLapply(cl, spectra, fun=function (spectra) undigested_filtering_function (spectra, mass_range_of_interest, percentage_of_spectrum_in_the_range))
# Close the processes
stopCluster(cl)
### Keep only the elements that are different from NULL
spectra_digested <- spectra_digested [!sapply(spectra_digested, is.null)]
#
return (spectra_digested)
}





################################################################################





################# GENERATE REPRESENTATIVE SPECTRA, BASED UPON SPECTRA SIMILARITY (ONE IMZML)
generate_representative_spectra_hca <- function (spectra, spectra_per_patient=10, tof_mode="linear", algorithm="hierarchicalClustering") {
# Detect and align peaks
if (tof_mode=="linear") {
	peaks <- detectPeaks (spectra, method="MAD", snr=3, halfWindowSize=20)
	peaks <- align_and_filter_peaks (peaks, tolerance_ppm=2000, peaks_filtering=FALSE, frequency_threshold=0.25)
}
if (tof_mode=="reflectron" | tof_mode=="reflector") {
	peaks <- detectPeaks (spectra, method="MAD", snr=3, halfWindowSize=20)
	peaks <- align_and_filter_peaks (peaks, tolerance_ppm=200, peaks_filtering=FALSE, frequency_threshold=0.25)
}
# Generate the peaklist matrix
peaklist <- intensityMatrix (peaks, spectra)
############################ HCA
if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
	# Compute the distance matrix
	distance_matrix <- dist (peaklist, method="euclidean")
	hca <- hclust(distance_matrix)
	# Associate to each row/spectrum the subgroup of the tree to which it belongs
	# k must be between 1 and 56
	if (spectra_per_patient > 56) {
		spectra_per_patient <- 56
	}
	hca_groups <- cutree(hca, k=spectra_per_patient)
	# Generate the final list of spectra
	spectra_grouped <- list()
	# For each subgroup to be isolated...
	for (s in 1:spectra_per_patient) {
		# Index the spectra under in the selected subgroup of the HCA
		index <- which(hca_groups == s)
		spectra_hca <- spectra[index]
		# Average the spectra in this HCA subgroup
		spectra_average <- averageMassSpectra (spectra_hca, method="mean")
		# Add the average to the final list
		spectra_grouped <- append(spectra_grouped, spectra_average)
	}
}
######################### K-MEANS
if (algorithm == "k_means" | algorithm == "kmeans" | algorithm == "k-Means") {
	# Compute the k-Means clustering
	k_means <- kmeans(peaklist, centers=spectra_per_patient)
	# Associate to each row/spectrum the subgroup/cluster to which it belongs
	k_means_groups <- k_means$cluster
	# Generate the final list of spectra
	spectra_grouped <- list()
	# For each subgroup to be isolated...
	for (s in 1:spectra_per_patient) {
		# Index the spectra under in the selected subgroup of the HCA
		index <- which(k_means_groups == s)
		spectra_k_means <- spectra[index]
		# Average the spectra in this HCA subgroup
		spectra_average <- averageMassSpectra (spectra_k_means, method="mean")
		# Add the average to the final list
		spectra_grouped <- append(spectra_grouped, spectra_average)
	}
}
return (spectra_grouped)
}





################################################################################





############################################## DEISOTOPING
de_isotope_peaks <- function (peaks) {
##################################################### Deisotoping Function
deisotoping_function <- function (peaks) {
	# Deisotoped Ions
	deisotoped_ions_mass <- numeric()
	deisotoped_ions_intensity <- numeric()
	# Take each ion and evaluate the consecutive ions
	isotope_bin_mass <- numeric ()
	isotope_bin_intensity <- numeric ()
	# Scroll the ion list
	for (i in 1:(length(peaks@mass)-1)) {
		# Create the isotope bins, if the consecutive peaks belong to a isotope cluster
		if (abs(peaks@mass[i+1] - peaks@mass[i]) < 1.5) {
			isotope_bin_mass <- append(isotope_bin_mass, peaks@mass[i])
			isotope_bin_intensity <- append(isotope_bin_intensity, peaks@intensity[i])
		}
		# When a major distance is found... We are out of the cluster
		if (abs(peaks@mass[i+1] - peaks@mass[i]) >= 1.5) {
			# Add this mass to the isotope bins (it will be the last belonging to a this isotope cluster
			isotope_bin_mass <- append(isotope_bin_mass, peaks@mass[i])
			isotope_bin_intensity <- append(isotope_bin_intensity, peaks@intensity[i])
			# Store the isotope of the isotope Bins with the maximum intensity created so far in a new vector
			most_intense_peak_mass <- isotope_bin_mass [which(isotope_bin_intensity == max(isotope_bin_intensity))]
			most_intense_peak_intensity <- max(isotope_bin_intensity)
			# Prevent adding more masses and only one intensity, kep the balance
			if (length(most_intense_peak_mass) == 1) {
				deisotoped_ions_mass <- append(deisotoped_ions_mass, most_intense_peak_mass)
				deisotoped_ions_intensity <- append(deisotoped_ions_intensity, most_intense_peak_intensity)
			}
			if (length(most_intense_peak_mass) > 1) {
				deisotoped_ions_mass <- append(deisotoped_ions_mass, most_intense_peak_mass[1])
				deisotoped_ions_intensity <- append(deisotoped_ions_intensity, most_intense_peak_intensity[1])
			}
			# Empty the isotope Bins in order for it to become available for the future cycle
			isotope_bin_mass <- numeric()
			isotope_bin_intensity <- numeric()

		}
	}
	# Fill the peaklist back in
	peaks@mass <- deisotoped_ions_mass
	peaks@intensity <- deisotoped_ions_intensity
	return (peaks)
}
################################################################################
cl <- makeCluster(cpu_core_number)
peaks <- parLapply(cl, peaks, fun = function(peaks) deisotoping_function(peaks))
stopCluster(cl)
return (peaks)
}





################################################################################





############################################ OUTLIERS REMOVAL
outliers_removal <- function (vector, replace_with="") {
summary_vector <- summary (vector)
# Calculate the interquartile range
inter_quartile_range <- summary_vector[5] - summary_vector[2]
# Calculate the fences, beyond which the spectrum is an outlier
iqr_fences <- c((summary_vector[2] - 1.5*inter_quartile_range), (summary_vector[5] + 1.5*inter_quartile_range))
# Find the outliers based on the fences condition
outliers_position <- which(vector < iqr_fences[1] | vector > iqr_fences[2])
# If the outliers have to be discarded...
if (replace_with == "") {
	if (length(outliers_position) > 0) {
		# Remove the correspondent spectra from the dataset
		vector <- vector [-outliers_position]
	}
}
if (replace_with == "mean") {
	if (length(outliers_position) > 0) {
		# Remove the correspondent spectra from the dataset
		vector_no_outliers <- vector [-outliers_position]
		# Replace the outliers with the vector mean(without outliers)
		vector [outliers_position] <- mean(vector_no_outliers)
	}
}
if (replace_with == "median") {
	if (length(outliers_position) > 0) {
		# Remove the correspondent spectra from the dataset
		vector_no_outliers <- vector [-outliers_position]
		# Replace the outliers with the vector median
		vector [outliers_position] <- median(vector_no_outliers)
	}
}
return (list (vector=vector, outliers_position=outliers_position))
}





########################################################################





beep(5)
