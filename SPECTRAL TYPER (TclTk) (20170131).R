###################### FUNCTIONS - MASS SPECTROMETRY 2017.01.31

# Update the packages
update.packages(repos="http://cran.mirror.garr.it/mirrors/CRAN/", ask=FALSE)


########################################################################## MISC

##################################################### INSTALL REQUIRED PACKAGES
# This function installs and loads the selected packages
install_and_load_required_packages <- function(required_packages, repository="http://cran.mirror.garr.it/mirrors/CRAN/") {
# Retrieve the installed packages
installed_packages <- installed.packages() [,1]
# Determine the missing packages
missing_packages <- character()
for (p in 1:length(required_packages)) {
	if ((required_packages[p] %in% installed_packages) == FALSE) {
		missing_packages <- append(missing_packages, required_packages[p])
	}
}
# If a repository is specified
if (repository != "" || !is.null(repository)) {
	if (length(missing_packages) > 0) {
		install.packages(missing_packages, repos=repository)
	}
} else {
	if (length(missing_packages) > 0) {
		install.packages(missing_packages)
	}
}
# Load the packages
for (i in 1:length(required_packages)) {
	library (required_packages[i], character.only=TRUE)
}
}





###############################################################################





############################### ADD CUSTOM FEATURES TO THE PEAKLIST INTENSITY MATRIX
# The function takes a list of spectra and a vector of custom features to be included in the generation of the final peaklist intensity matrix. If a matrix is specified as input, the columns corresponding to the features to be searched for are appended to the matrix itself. The input matrix must have the number of spectra as the number of rows.
custom_peaklist_intensity_matrix <- function (spectra, features_to_add=numeric(), final_sample_matrix=NULL, multicore_processing=TRUE, tolerance_ppm=2000) {
	## Install the required packages
	install_and_load_required_packages("MALDIquant")
	#################### Multiple spectra
	if (isMassSpectrumList(spectra)) {
		### If there are features to add...
		if (length(features_to_add > 0)) {
			## Generate a fake spectrum and a fake peaklist with the features to add
			fake_spectrum <- createMassSpectrum(mass=spectra[[1]]@mass, intensity=spectra[[1]]@intensity, metaData=list(name="Fake spectrum"))
			fake_peaks <- createMassPeaks(mass=as.numeric(features_to_add), intensity=rep(1, length(features_to_add)), snr=rep(3, length(features_to_add)), metaData=list(name="Fake peaklist"))
			## Detect the peaks in the spectra
			peaks <- detectPeaks(spectra, method="SuperSmoother", SNR=3)
			## Append the fake spectrum and the fake peaklist to he original lists
			spectra_all <- append(spectra, fake_spectrum)
			peaks_all <- append(peaks, fake_peaks)
			## Generate the intensity matrix
			intensity_matrix_all <- intensityMatrix(peaks_all, spectra_all)
			## Remove the last row (corresponding to the fake spectrum)
			intensity_matrix_all <- intensity_matrix_all[1:(nrow(intensity_matrix_all)-1),]
			## Keep only the columns that are corresponding to the desired features
			final_intensity_matrix <- intensity_matrix_all[,features_to_add]
			# If the final matrix does not exist yet and it is null, the final matrix becomes the feature column
			if (is.null(final_sample_matrix)) {
				final_sample_matrix <- final_intensity_matrix
			} else {
				# If the final matrix exists, append the feature column to the matrix
				final_sample_matrix <- cbind(final_sample_matrix, final_intensity_matrix)
			}
		}
	} else if (isMassSpectrum(spectra)) {
		### If there are features to add...
		if (length(features_to_add > 0)) {
			# Scroll the features to add (in the model but not in the sample)
			for (f in 1:length(features_to_add)) {
				# Initialise the output
				x_intensity <- NA
				# Scroll the mass list of each spectrum
				for (m in 1:length(spectra@mass)) {
					# If there is a match
					if (abs(spectra@mass[m]-as.numeric(features_to_add[f]))*10^6/as.numeric(features_to_add[f]) <= tolerance_ppm) {
						# Store the corresponding intensity and generate the matrix column
						x_intensity <- as.matrix(cbind(spectra@intensity[m]))
						colnames(x_intensity) <- features_to_add[f]
						# If the final matrix does not exist yet and it is null, the final matrix becomes the feature column
						if (is.null(final_sample_matrix)) {
							final_sample_matrix <- cbind(x_intensity)
						} else {
							# If the final matrix exists, append the feature column to the matrix
							final_sample_matrix <- cbind(final_sample_matrix, x_intensity)
						}
						# Do not keep searching
						break
					}
				}
			}
		}
	}
	##### Return the final matrix with the custom features (plus the original ones if a matrix is specified as input)
	return(final_sample_matrix)
}





###############################################################################






############################### ADD THE CLASS AND THE SAMPLE NAME TO THE MATRIX
# This function adds two column to the peaklist matrix (rows: spectra/patients, columns: aligned peaks): Sample and Class, according to the file name.
### The name of the rows will be either the sample name or the class name (depending on the function parameter).
# If the rows are named according to the sample name, an additional column for the class is added
matrix_add_class_and_sample <- function (signal_matrix, peaks=list(), class_list=list(), spectra_format="imzml", sample_output=TRUE, class_output=TRUE, row_labels="Sample") {
	# Convert the input matrix/dataframe into a matrix
	signal_matrix <- as.matrix(signal_matrix)
	# Determine the number of spectra/peaklists
	if (isMassPeaksList(peaks)) {
		number_of_spectra <- length(peaks)
	} else if (isMassPeaks(peaks)) {
		number_of_spectra <- 1
	}
	####################################### FILE VECTOR
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	if (isMassPeaksList(peaks)) {
		for (i in 1:length(peaks)) {
			if (spectra_format == "imzml" || spectra_format == "imzML") {
				file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
			} else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
				file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName[1])
			}
		}
	} else if (isMassPeaks(peaks)) {
		if (spectra_format == "imzml" || spectra_format == "imzML") {
			file_vector <- append(file_vector, peaks@metaData$file[1])
		} else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
			file_vector <- append(file_vector, peaks@metaData$sampleName[1])
		}
	}
	################################## Only the sample
	if ((class_output == FALSE && sample_output == TRUE) || (class_output == TRUE && length(class_list) == 0 && sample_output == TRUE)) {
		# Create the sample matrix column and append it to the global matrix
		# Sample as rownames
		if ("sample" %in% row_labels) {
			rownames(signal_matrix) <- file_vector
		} else {
			sample_column <- matrix("", ncol = 1, nrow = number_of_spectra)
			colnames(sample_column) <- "Sample"
			sample_column[,1] <- cbind(file_vector)
			signal_matrix <- cbind(signal_matrix, sample_column)
		}
	}
	################################## Both the class and the sample
	if (class_output == TRUE && length(class_list) >= 1 && sample_output == TRUE) {
		# Create the sample matrix column and append it to the global matrix
		# Sample as rownames
		if (row_labels == "sample") {
		rownames(signal_matrix) <- file_vector
		} else {
			sample_column <- matrix("", ncol = 1, nrow = number_of_spectra)
			colnames(sample_column) <- "Sample"
			sample_column[,1] <- cbind(file_vector)
			signal_matrix <- cbind(signal_matrix, sample_column)
		}
		### Add the class column
		class_list <- sort(class_list)
		# Rename the classes according to the class_list vector (the match should be /class/ to avoid catching the name of the class in previous folders)
		class_vector <- file_vector
		for (p in 1:length(class_vector)) {
			for (w in 1:length(class_list)) {
				if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
					if (length(grep(paste("/", class_list[w], "/", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
						class_vector[p] <- class_list[w]
					}
				} else if (Sys.info()[1] == "Windows") {
					if (length(grep(paste("\\\\", class_list[w], "\\\\", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
						class_vector[p] <- class_list[w]
					}
				}
			}
		}
		# Class as rownames
		if (row_labels == "class") {
			rownames(signal_matrix) <- class_vector
		} else {
			class_column <- matrix("", ncol = 1, nrow = number_of_spectra)
			colnames(class_column) <- "Class"
			# Fill in the matrix column with the file_vector classes and samples
			class_column[,1] <- cbind(class_vector)
			signal_matrix <- cbind(signal_matrix, class_column)
		}
	}
	################################## Only the class
	if (class_output == TRUE && length(class_list) >= 1 && sample_output == FALSE) {
		class_list <- sort(class_list)
		# Rename the classes according to the class_list vector (the match should be /class/ to avoid catching the name of the class in previous folders)
		class_vector <- file_vector
		for (p in 1:length(class_vector)) {
			for (w in 1:length(class_list)) {
				if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
					if (length(grep(paste("/", class_list[w], "/", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
						class_vector[p] <- class_list[w]
					}
				} else if (Sys.info()[1] == "Windows") {
					if (length(grep(paste("\\\\", class_list[w], "\\\\", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
						class_vector[p] <- class_list[w]
					}
				}
			}
		}
		if ("class" %in% row_labels) {
			rownames(signal_matrix) <- class_vector
		} else {
			### Add the class column
			class_column <- matrix("", ncol = 1, nrow = number_of_spectra)
			colnames(class_column) <- "Class"
			# Fill in the matrix column with the file_vector classes and samples
			class_column[,1] <- cbind(class_vector)
			signal_matrix <- cbind(signal_matrix, class_column)
		}
	}
	### Add these matrix columns to the peaklist matrix
	return(signal_matrix)
}





################################################################################





###################################### ADD THE THY CLASS TO THE MATRIX (THYROID)
# This function adds the THY column to the peaklist matrix, by reading the THY value from the sample name (imzML)
matrix_add_thy <- function (signal_matrix, peaks, spectra_format="imzml") {
number_of_spectra <- length(peaks)
# Create the empty vector
file_vector <- character()
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(peaks)) {
	if (spectra_format == "imzml" || spectra_format == "imzML") {
		file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
	}
	if (spectra_format == "brukerflex" || spectra_format == "xmass") {
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





################################################################ SIGNAL FOLLOWER
# This function takes a folder containing imzML files in it and a list of signals of interest (along with their possible names). It imports the spectra and returns a matrix listing the signal of interests along with their statistics.
signal_follower_statistics <- function (filepath, signal_list, mass_labels=list(), SNR=5, spectra_format="imzml", tof_mode="linear", smoothing_strength_preprocessing="medium", process_in_packages_of=length(spectra), tolerance_ppm=2000, peak_picking_algorithm="SuperSmoother") {
# Load the required libraries
install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
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
} else {
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
	# Spectra import and preprocessing
	if (!is.null(mass_range)) {
		spectra <- importImzMl(filepath[lib], massRange=mass_range)
	} else {
		spectra <- importImzMl(filepath[lib])
	}
	spectra <- preprocess_spectra(spectra, tof_mode=tof_mode, smoothing_strength=smoothing_strength_preprocessing, process_in_packages_of=process_in_packages_of)
	### Peak Picking and alignment
	peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	peaks <- align_and_filter_peaks(peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
	# Generate the intensity matrix
	intensity_matrix <- intensityMatrix(peaks, spectra)
	#### Create the partial result matrices
	st_dev_matrix <- matrix (ncol=length(signal_list), nrow=1)
	colnames(st_dev_matrix) <- column_names_st_dev
	rownames(st_dev_matrix) <- spectra[[1]]@metaData$file[1]
	coeff_var_matrix <- matrix (0, ncol=length(signal_list), nrow=1)
	colnames(coeff_var_matrix) <- column_names_coeff_var
	rownames(coeff_var_matrix) <- spectra[[1]]@metaData$file[1]
	mean_abs_dev_matrix <- matrix (0, ncol=length(signal_list), nrow=1)
	colnames(mean_abs_dev_matrix) <- column_names_mean_abs_dev
	rownames(mean_abs_dev_matrix) <- spectra[[1]]@metaData$file[1]
	# For each signal in the signal_list
	for (s in 1:length(signal_list)) {
		# Identify the column of interest
		signal_column <- which(abs((as.numeric(colnames(intensity_matrix))-signal_list[s])*10^6/signal_list[s]) <= tolerance_ppm)
		# Isolate the column of interest (there might not be one!!)
		if (length(signal_column) != 0) {
			intensity_vector <- intensity_matrix[,signal_column]
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
	}
	# Put the partial result matrices together in one matrix
	st_dev_result_matrix [lib,] <- st_dev_matrix
	coeff_var_result_matrix [lib,] <- coeff_var_matrix
	mean_abs_dev_result_matrix [lib,] <- mean_abs_dev_matrix
}
result_matrix <- cbind(st_dev_result_matrix, coeff_var_result_matrix, mean_abs_dev_result_matrix)
return (result_matrix)
}





################################################################################





##################################################### REMOVE LOW INTENSITY PEAKS
# This function removes low-intensity peaks (in terms of level of intensity compared with the most intense peak in the peaklist) from the list of provided peaks (MALDIquant).
# If the method is selected to be "element-wise", each element of the peaklist is evaluated, and the intensity threshold is calculated over the peaks of only that element. Otherwise, if "whole" is selected, the threshold is calculated on all the peaks in the dataset.
remove_low_intensity_peaks <- function (peaks, intensity_threshold_percent=0.1, intensity_threshold_method="element-wise", multicore_processing=TRUE) {
# Load the required libraries
install_and_load_required_packages("parallel")
# If there is only one peaklist, there is no point in doing the whole method, but only the element-wise.
if (!isMassPeaksList(peaks)) {
	intensity_threshold_method <- "element-wise"
}
################################################################## ELEMENT-WISE
if (intensity_threshold_method == "element-wise") {
	############### INTENSITY FILTERING FUNCTION (ELEMENT-WISE)
	intensity_filtering_subfunction_element <- function (peaks, intensity_threshold_percent) {
		# Filter out the peaks whose intensity is below a certain threshold
		# Store mass and intensity into vectors
		intensity_values <- peaks@intensity
		mass_values <- peaks@mass
		snr_values <- peaks@snr
		# Identify the positions of the values to be discarded
		values_to_be_discarded <- intensity_values[((intensity_values*100/max(intensity_values,na.rm=TRUE)) < intensity_threshold_percent)]
		# If there are values to be discarded...
		if (length(values_to_be_discarded) > 0) {
			# Identify the positions
			positions_to_be_discarded <- numeric()
			for (i in 1:length(values_to_be_discarded)) {
				value_position <- which(intensity_values == values_to_be_discarded[i])
				positions_to_be_discarded <- append(positions_to_be_discarded, value_position)
			}
			# Discard the values from the vectors
			intensity_values <- intensity_values [-positions_to_be_discarded]
			mass_values <- mass_values [-positions_to_be_discarded]
			snr_values <- snr_values [-positions_to_be_discarded]
			# Put the values back into the MALDIquant list
			peaks@mass <- mass_values
			peaks@intensity <- intensity_values
			peaks@snr <- snr_values
		} else {
			# If there aren't any values to be discarded...
			peaks@mass <- mass_values
			peaks@intensity <- intensity_values
			peaks@snr <- snr_values
		}
		return (peaks)
	}
	######################################### Multiple peaks elements
	if (isMassPeaksList(peaks)) {
		########## MULTICORE
		if (multicore_processing == TRUE) {
			# Detect the number of cores
			cpu_thread_number <- detectCores(logical=TRUE) - 1
			if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
				peaks_filtered <- mclapply(peaks, FUN = function (peaks) intensity_filtering_subfunction_element(peaks, intensity_threshold_percent), mc.cores=cpu_thread_number)
			} else if (Sys.info()[1] == "Windows") {
				# Make the CPU cluster for parallelisation
				cl <- makeCluster(cpu_thread_number)
				# Make the cluster use the custom functions and the package functions along with their parameters
				clusterEvalQ(cl, {library(MALDIquant)})
				# Pass the variables to the cluster for running the function
				clusterExport(cl=cl, varlist=c("peaks", "intensity_threshold_percent", "intensity_filtering_subfunction_element"), envir=environment())
				# Apply the multicore function
				peaks_filtered <- parLapply(cl, peaks, fun = function (peaks) intensity_filtering_subfunction_element(peaks, intensity_threshold_percent))
				stopCluster(cl)
			}
		} else {
			########## SINGLE CORE
			peaks_filtered <- lapply(peaks, FUN = function (peaks) intensity_filtering_subfunction_element(peaks, intensity_threshold_percent))
		}
	} else {
		######################################### Single peaks element
		peaks_filtered <- intensity_filtering_subfunction_element(peaks, intensity_threshold_percent)
	}
}
################################################################# WHOLE DATASET
if (intensity_threshold_method == "whole") {
	############### INTENSITY FILTERING FUNCTION (ELEMENT-WISE)
	intensity_filtering_subfunction_whole <- function (peaks, intensity_threshold_percent, highest_intensity) {
		# Filter out the peaks whose intensity is below a certain threshold
		# Store mass and intensity into vectors
		intensity_values <- peaks@intensity
		mass_values <- peaks@mass
		snr_values <- peaks@snr
		# Identify the positions of the values to be discarded
		values_to_be_discarded <- intensity_values[((intensity_values*100/highest_intensity) < intensity_threshold_percent)]
		# If there are values to be discarded...
		if (length(values_to_be_discarded) > 0) {
			# Identify the positions
			positions_to_be_discarded <- numeric()
			for (i in 1:length(values_to_be_discarded)) {
				value_position <- which(intensity_values == values_to_be_discarded[i])
				positions_to_be_discarded <- append(positions_to_be_discarded, value_position)
			}
			# Discard the values from the vectors
			intensity_values <- intensity_values [-positions_to_be_discarded]
			mass_values <- mass_values [-positions_to_be_discarded]
			snr_values <- snr_values [-positions_to_be_discarded]
			# Put the values back into the MALDIquant list
			peaks@mass <- mass_values
			peaks@intensity <- intensity_values
			peaks@snr <- snr_values
		} else {
			# If there aren't any values to be discarded...
			peaks@mass <- mass_values
			peaks@intensity <- intensity_values
			peaks@snr <- snr_values
		}
		return (peaks)
	}
	############### Determine the highest peak in the dataset
	highest_peak <- NULL
	highest_intensity <- 0
	for (p in 1:length(peaks)) {
		if (length(peaks[[p]]@mass) > 0) {
			for (m in 1:length(peaks[[p]]@mass)) {
				if (highest_intensity == 0 || peaks[[p]]@intensity[m]>highest_intensity) {
					highest_intensity <- peaks[[p]]@intensity[m]
					highest_peak <- peaks[[p]]@mass[m]
				}
			}
		}
	}
	# Filter the peaks
	######################################### Multiple peaks elements
	if (isMassPeaksList(peaks)) {
		########## MULTICORE
		if (multicore_processing == TRUE) {
			# Detect the number of cores
			cpu_thread_number <- detectCores(logical=TRUE) - 1
			if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
				peaks_filtered <- mclapply(peaks, FUN = function(peaks) intensity_filtering_subfunction_whole(peaks, intensity_threshold_percent, highest_intensity), mc.cores=cpu_thread_number)
			} else if (Sys.info()[1] == "Windows") {
				# Make the CPU cluster for parallelisation
				cl <- makeCluster(cpu_thread_number)
				# Make the cluster use the custom functions and the package functions along with their parameters
				clusterEvalQ(cl, {library(MALDIquant)})
				# Pass the variables to the cluster for running the function
				clusterExport(cl=cl, varlist=c("peaks", "intensity_threshold_percent", "intensity_filtering_subfunction_whole", "highest_intensity"), envir=environment())
				# Apply the multicore function
				peaks_filtered <- parLapply(cl, peaks, fun = function(peaks) intensity_filtering_subfunction_whole(peaks, intensity_threshold_percent, highest_intensity))
				stopCluster(cl)
			}
		} else {
			########## SINGLE CORE
			peaks_filtered <- lapply(peaks, FUN = function(peaks) intensity_filtering_subfunction_whole(peaks, intensity_threshold_percent, highest_intensity))
		}
	} else {
		######################################### Single peaks element
		peaks_filtered <- intensity_filtering_subfunction_whole(peaks, intensity_threshold_percent, highest_intensity)
	}
}
return (peaks_filtered)
}





################################################################################





########################################################### SPECTRA FILES READER
# This function reads all the files from a folder and returns only the imzML files or the fid files.
read_spectra_files <- function (folder, spectra_format="imzml", full_path=TRUE) {
if (spectra_format == "imzml" || spectra_format == "imzML") {
	# Read all the files
	folder_all_files <- list.files(folder, full.names=full_path, recursive=TRUE)
	# Create the empty vector in which only the imzML files will be listed
	spectra_files <- character()
	# Put the imzML files in the new vector discarding the ibd files
	for (l in 1:length(folder_all_files)) {
		if (length(grep(".imzML", folder_all_files[l], fixed=TRUE)) == 1) {
			spectra_files <- append(spectra_files, folder_all_files[l])
		}
	}
}
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	# Read all the folders
	sample_directory_list_all_folders <- dir(folder, ignore.case=TRUE, full.names=full_path, recursive=TRUE, include.dirs=TRUE)
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





################################################################################





###################### SPECTRA PURIFICATION BASED ON THE TIC (Before processing)
# This function discards spectra with a total ion current (TIC) below the selected threshold.
# A specific threshold can be specified; if just the TOF MS mode is specified, a default TIC threshold value is used.
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





######################### REPLACE THE SNR WITH THE STDEV IN THE PEAKS
# This function computes the standard deviation of each peak of an average spectrum peaklist, by replacing the existing SNR slot with the SD or CV: all the peaks (average and dataset) are aligned and each peak of the average peaklist is searched across the dataset thanks to the intensity matrix.
replace_SNR_in_avg_peaklist <- function (spectra, SNR=5, tof_mode="linear", tolerance_ppm=2000, spectra_format="imzml", replace_snr_with="std") {
install_and_load_required_packages(c("MALDIquant", "stats"), peak_picking_algorithm="SuperSmoother")
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Check if it is not a single spectrum
# If the spectra are many...
if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
	# Average the spectra
	avg_spectrum <- averageMassSpectra(spectra, method="mean")
	avg_spectrum <- removeBaseline(avg_spectrum, method="TopHat")
	avg_spectrum <- calibrateIntensity(avg_spectrum, method="TIC")
	# Peak picking
	peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	peaks_avg <- peak_picking(avg_spectrum, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	# Merge for the alignment
	global_peaks <- append(peaks_avg, peaks)
	global_peaks <- align_and_filter_peaks(global_peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
	peaks_avg <- global_peaks[[1]]
	peaks <- global_peaks [2:length(global_peaks)]
	# Compute the intensity matrix
	intensity_matrix <- intensityMatrix(peaks, spectra)
	# Scroll the peaks of the average peaklist...
	for (p in 1:length(peaks_avg@mass)) {
		# Search for the column in the intensity_matrix
		for (z in 1:ncol(intensity_matrix)) {
			# Match...
			if (peaks_avg@mass[p] == colnames(intensity_matrix)[z]) {
				# Create an empty vector where to allocate the intensity of this peak into the dataset
				intensity_vector <- intensity_matrix[,z]
				# Compute the standard deviation
				if (replace_snr_with == "std" || replace_snr_with == "sd") {
					std_int <- sd(intensity_vector)
					# Replacement
					peaks_avg@snr[p] <- std_int
				}
				if (replace_snr_with == "cv" || replace_snr_with == "CV") {
					std_int <- sd(intensity_vector)
					mean_int <- mean(intensity_vector)
					cv_int <- std_int / mean_int
					# Replacement
					peaks_avg@snr[p] <- cv_int
				}
				# Save time avoiding to scroll to the end when found
				break
			} else {peaks_avg@snr[p] <- NaN}
		}
	}
	return (peaks_avg)
} else {
	return(print("The spectra list is empty or the spectra list contains only one spectrum"))
}
}





################################################################################





############################################################### OUTLIERS REMOVAL
# This function takes a vector on values as input and returns the same vector without the outliers, calculated based upon the interquartile range (fence rule).
# The outliers can be replaced with nothing (so they are removed from the vector) or with some value
outliers_removal <- function (v, replace_with="") {
summary_vector <- summary(v)
# Calculate the interquartile range
inter_quartile_range <- summary_vector[5] - summary_vector[2]
# Calculate the fences, beyond which the spectrum is an outlier
iqr_fences <- c((summary_vector[2] - 1.5*inter_quartile_range), (summary_vector[5] + 1.5*inter_quartile_range))
# Find the outliers based on the fences condition
outliers_position <- which(v < iqr_fences[1] | v > iqr_fences[2])
# If the outliers have to be discarded...
if (replace_with == "") {
	if (length(outliers_position) > 0) {
		# Remove the correspondent elements from the dataset
		v <- v [-outliers_position]
	}
}
# If the outliers have to be replaced...
if (is.numeric(replace_with)) {
	if (length(outliers_position) > 0) {
        # Replace the outliers with the replacement
		v [outliers_position] <- replace_with
	}
}
if (replace_with == 0 || replace_with == "zero") {
	if (length(outliers_position) > 0) {
        # Replace the outliers with the replacement
		v [outliers_position] <- 0
	}
}
if (replace_with == "NA" || is.na(replace_with)) {
	if (length(outliers_position) > 0) {
        # Replace the outliers with the replacement
		v [outliers_position] <- NA
	}
}
if (replace_with == "mean") {
	if (length(outliers_position) > 0) {
        # Remove the correspondent elements from the dataset
		vector_no_outliers <- v [-outliers_position]
		# Replace the outliers with the vector mean (no outliers)
		v [outliers_position] <- mean(vector_no_outliers)
	}
}
if (replace_with == "median") {
	if (length(outliers_position) > 0) {
        # Remove the correspondent elements from the dataset
		vector_no_outliers <- v [-outliers_position]
		# Replace the outliers with the vector median (no outliers)
		v [outliers_position] <- median(vector_no_outliers)
	}
}
return (list (vector=v, outliers_position=outliers_position))
}





################################################################################





######################################### PEAK STATISTICS (on processed Spectra)
# This function computes the peak statistics onto a selected spectra dataset (or to the provided peaks), both when the spectra belong to no (or one) class and more classes.
# It returns a NULL value if the peak statistics cannot be performed.
peak_statistics <- function (spectra, peaks=NULL, SNR=3, peak_picking_algorithm="SuperSmoother", class_list=NULL, class_in_file_name=TRUE, tof_mode="linear", spectra_format="imzml", exclude_spectra_without_peak=FALSE, alignment_iterations=5, peaks_filtering=TRUE, frequency_threshold_percent=25, remove_outliers=TRUE, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1, intensity_threshold_method="element_wise") {
	########## Load the required libraries
	install_and_load_required_packages(c("MALDIquant", "stats"))
	########## Rename the trim function
	trim_spectra <- get(x="trim", pos="package:MALDIquant")
	########## Define the tolerance in PPM
	if (tof_mode =="linear" || tof_mode=="Linear" || tof_mode=="L") {
		tolerance_ppm <- 2000
	} else if (tof_mode =="reflector" || tof_mode=="reflectron" || tof_mode=="R") {
		tolerance_ppm <- 200
	}
	########## Determine the number of classes
	if (length(class_list) == 0 || length(class_list) == 1 || is.null(class_list)) {
	number_of_classes <- 1
	} else if (length(class_list) > 1) {
		number_of_classes <- length(class_list)
	}
	########## Detect (if not already provided) and Align Peaks
	if (is.null(peaks)) {
		peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	}
	peaks <- align_and_filter_peaks(peaks, tof_mode=tof_mode, alignment_iterations=alignment_iterations, peaks_filtering=peaks_filtering, frequency_threshold_percent=frequency_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=intensity_threshold_percent, intensity_threshold_method=intensity_threshold_method, reference_peaklist=NULL, spectra=spectra)
	# Generate the matrix (and convert it into a data frame)
	if (exclude_spectra_without_peak == FALSE) {
		signal_matrix <- intensityMatrix(peaks, spectra)
	} else if (exclude_spectra_without_peak == TRUE) {
		signal_matrix <- intensityMatrix(peaks)
	}
	# Peak vector
	if (is.matrix(signal_matrix)) {
		peak_vector <- as.numeric(colnames(signal_matrix))
	} else if (is.data.frame(signal_matrix)) {
		peak_vector <- as.numeric(names(signal_matrix))
	}
	############################################################## ONE CLASS
	if (number_of_classes == 1) {
		################################# FUNCTION for matrix APPLY (it will applied for each matrix column, for each peak)
		peak_statistcs_function <- function (signal_matrix_column, signal_matrix, remove_outliers) {
			# Generate the output matrix row
			peak_stat_matrix_row <- matrix (0, nrow=1, ncol=7)
			rownames(peak_stat_matrix_row) <- as.numeric(colnames(signal_matrix_column))
			colnames(peak_stat_matrix_row) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles")
			# Start the calculation
			intensity_vector <- as.numeric(signal_matrix_column)
			if (remove_outliers == TRUE) {
				intensity_vector <- outliers_removal(intensity_vector)$vector
			}
			# Calculate the statistical parameters on the intensity values in the vector
			# Normality
			if (length(intensity_vector) >= 3 & length(intensity_vector) <= 5000) {
				shapiro_test <- shapiro.test(intensity_vector)
				if (shapiro_test$p.value < 0.05) {
					distribution_type <- paste("Non-normal", "(Shapiro p-value:", round(shapiro_test$p.value,3), ")")
				}
				if (shapiro_test$p.value >= 0.05) {
					distribution_type <- paste("Normal", "(Shapiro p-value:", round(shapiro_test$p.value,3), ")")
				}
			} else if (length(intensity_vector) < 3) {
				distribution_type <- "Not determinable, number of samples too low"
			} else if (length(intensity_vector) > 5000) {
				distribution_type <- "Number of samples too high, assume it is normal"
			}
			# Other parameters
			st_dev_intensity <- sd(intensity_vector, na.rm=TRUE)
			summary_intensity_vector <- summary(intensity_vector)
			mean_intensity <- summary_intensity_vector [4]
			coeff_variation <- (st_dev_intensity / mean_intensity) *100
			median_intensityensity <- summary_intensity_vector [3]
			first_quartile <- summary_intensity_vector [2]
			third_quartile <- summary_intensity_vector [5]
			inter_quartile_range <- third_quartile - first_quartile
			# Fill the matrix with the values
			peak_stat_matrix_row [,1] <- distribution_type
			peak_stat_matrix_row [,2] <- as.numeric(mean_intensity)
			peak_stat_matrix_row [,3] <- as.numeric(st_dev_intensity)
			peak_stat_matrix_row [,4] <- as.numeric(coeff_variation)
			peak_stat_matrix_row [,5] <- as.numeric(median_intensityensity)
			peak_stat_matrix_row [,6] <- as.numeric(inter_quartile_range)
			peak_stat_matrix_row [,7] <- paste("1st quartile", first_quartile, "; 3rd quartile", third_quartile)
			return (peak_stat_matrix_row)
		}
		###############
		# Fix the signal_matrix (Add the sample column)
		signal_matrix <- matrix_add_class_and_sample(signal_matrix, peaks=peaks, spectra_format=spectra_format, sample_output=TRUE, class_output=FALSE)
		# Output matrix
		peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-1), ncol=8)
		rownames(peak_stat_matrix) <- as.numeric(colnames(signal_matrix)[1:(ncol(signal_matrix)-1)])
		colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Sample")
		# Only peaks
		signal_matrix_peaks <- signal_matrix [,1:(ncol(signal_matrix)-1)]
		# Apply the function (transpose the result matrix)
		peak_stat_matrix <- t(apply(signal_matrix_peaks, MARGIN=2, FUN=function(x) peak_statistcs_function(x, signal_matrix, remove_outliers=remove_outliers)))
		# Generate the intensity matrix with NA if the peak is not present in the spectra
		intensity_matrix_with_na <- intensityMatrix(peaks)
		spectra_counter_vector <- numeric()
		for (pk in 1:ncol(intensity_matrix_with_na)) {
			intensity_vector <- intensity_matrix_with_na[,pk]
			spectra_counter_vector <- append(spectra_counter_vector, length(intensity_vector[!is.na(intensity_vector)]))
		}
		peak_stat_matrix <- cbind(peak_stat_matrix, spectra_counter_vector)
		# Fix the column names
		colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter")
		## Return
		return(peak_stat_matrix)
	} else if (number_of_classes > 1) {
		############################################################ TWO OR MORE CLASSES
		# Every variable now is a list, each element of which corresponds to a certain value from a class
		# So every variable is a list with the same length of the class list (each element of the list
		# is referred to a class
		# Fix the signal_matrix (Add the sample column)
		signal_matrix <- matrix_add_class_and_sample(signal_matrix, peaks=peaks, class_list=class_list, spectra_format=spectra_format, sample_output=TRUE, class_output=TRUE)
		# Check if there is a sufficient number of observations per class
		observations_per_class <- numeric()
		for (i in 1:length(class_list)) {
			observations_per_class <- append(observations_per_class, length(which(signal_matrix[,ncol(signal_matrix)] == class_list[i])))
		}
		sufficient_number_of_observations_per_class <- TRUE
		for (i in 1:length(observations_per_class)) {
			if (observations_per_class[i] < 3) {
				sufficient_number_of_observations_per_class <- FALSE
			}
		}
		##### Run only if there is a sufficient number of samples
		if (sufficient_number_of_observations_per_class == TRUE) {
			# Output matrix
			peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-2), ncol=14)
			rownames(peak_stat_matrix) <- as.numeric(colnames(signal_matrix)[1:(ncol(signal_matrix)-2)])
			colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)", "t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
			# For each peak
			for (p in 1:(ncol(signal_matrix)-2)) {
				# Put the intensity of that peak into one vector per class (in a global list)
				intensity_vector <- list()
				# Scroll the peaklists and Add the peak intensity to a vector(one for each class)
				for (l in 1:length(class_list)) {
					# Allocate in the intensity vector the rows for that peak belonging to the certain class
					intensity_vector[[l]] <- as.numeric(signal_matrix [signal_matrix[,ncol(signal_matrix)] == class_list[l],p])
				}
				if (remove_outliers == TRUE) {
					for (i in 1:length(intensity_vector)) {
						intensity_vector[[i]] <- outliers_removal(intensity_vector[[i]])
						intensity_vector[[i]] <- intensity_vector[[i]]$vector
					}
				}
				######################## STATISTICAL PARAMETERS
				############################################### Normality for each class
				shapiro_test <- list()
				distribution_type <- list()
				for (l in 1:length(class_list)) {
					if (length(intensity_vector[[l]]) >= 3 && length(intensity_vector[[l]]) <= 5000) {
						shapiro_test[[l]] <- shapiro.test(intensity_vector[[l]])
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
					variance_test_parametric <- var.test(intensity_vector[[1]], intensity_vector[[2]])
				}
				if (length(class_list) >= 2) {
					variance_test_non_parametric <- bartlett.test(as.numeric(signal_matrix[,p]), g=as.factor(signal_matrix[,ncol(signal_matrix)]))
				}
				########################################### Other parameters (per class)
				st_dev_intensity <- list()
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
					st_dev_intensity[[l]] <- sd(intensity_vector[[l]])
					summary_intensity_vector [[l]] <- summary(intensity_vector[[l]])
					mean_intensity[[l]] <- summary_intensity_vector[[l]] [4]
					coeff_variation[[l]] <- (st_dev_intensity[[l]] / mean_intensity[[l]]) *100
					median_intensityensity[[l]] <- summary_intensity_vector[[l]] [3]
					first_quartile[[l]] <- summary_intensity_vector[[l]] [2]
					third_quartile[[l]] <- summary_intensity_vector[[l]] [5]
					inter_quartile_range[[l]] <- third_quartile[[l]] - first_quartile[[l]]
					spectra_counter[[l]] <- length(intensity_vector[[l]])
					variance[[l]] <- var(intensity_vector[[l]])
				}
				############################################# Parameters between classes
				# T-test
				if (length(class_list) == 2) {
					t_test <- t.test(intensity_vector[[1]], intensity_vector[[2]])
				}
				# ANOVA TEST
				if (length(class_list) >= 2) {
				anova_test <- aov(signal_matrix[,p] ~ signal_matrix[,ncol(signal_matrix)])
				}
				# WILCOXON - MANN-WHITNEY TEST
				if (length(class_list) == 2) {
					wilcoxon_test <- wilcox.test(intensity_vector[[1]], intensity_vector[[2]])
				}
				# KRUSKAL-WALLIS TEST
				if (length(class_list) >= 2) {
					kruskal_wallis_test <- kruskal.test(signal_matrix[,p], g=as.factor(signal_matrix[,ncol(signal_matrix)]))
				}
				######################################## Fill the matrix with the values
				# Distribution Type
				distribution_type_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(distribution_type_name)) {
						distribution_type_name <- paste(distribution_type[[l]], " - ", class_list[l], sep="")
					} else {
						distribution_type_name <- paste(distribution_type_name, " , ", distribution_type[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,1] <- paste(distribution_type_name)
				# Mean
				mean_intensity_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(mean_intensity_name)) {
						mean_intensity_name <- paste(mean_intensity[[l]], " - ", class_list[l], sep="")
					} else {
						mean_intensity_name <- paste(mean_intensity_name, " , ", mean_intensity[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,2] <- mean_intensity_name
				# Standard Deviation
				st_dev_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(st_dev_name)) {
						st_dev_name <- paste(st_dev_intensity[[l]], " - ", class_list[l], sep="")
					} else {
						st_dev_name <- paste(st_dev_name, " , ", st_dev_intensity[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,3] <- st_dev_name
				# Coefficient of Variation
				coeff_variation_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(coeff_variation_name)) {
						coeff_variation_name <- paste(coeff_variation[[l]], " - ", class_list[l], sep="")
					} else {
						coeff_variation_name <- paste(coeff_variation_name, " , ", coeff_variation[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,4] <- coeff_variation_name
				# Median
				median_intensityensity_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(median_intensityensity_name)) {
						median_intensityensity_name <- paste(median_intensityensity[[l]], " - ", class_list[l], sep="")
					} else {
						median_intensityensity_name <- paste(median_intensityensity_name, " , ", median_intensityensity[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,5] <- median_intensityensity_name
				# Interquartile Range (IQR)
				inter_quartile_range_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(inter_quartile_range_name)) {
						inter_quartile_range_name <- paste(inter_quartile_range[[l]], " - ", class_list[l], sep="")
					} else {
						inter_quartile_range_name <- paste(inter_quartile_range_name, " , ", inter_quartile_range[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,6] <- inter_quartile_range_name
				# Spectra counter
				spectra_counter_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(spectra_counter_name)) {
						spectra_counter_name <- paste(spectra_counter[[l]], " - ", class_list[l], sep="")
					} else {
						spectra_counter_name <- paste(spectra_counter_name, " , ", spectra_counter[[l]], " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,7] <- spectra_counter_name
				# Class
				class_name <- NULL
				for (l in 1:length(class_list)) {
					if (is.null(class_name)) {
						class_name <- class_list[l]
					} else {
						class_name <- paste(class_name, " - ", class_list[l], sep="")
					}
				}
				peak_stat_matrix [p,8] <- class_name
				# Homoscedasticity (Parametric)
				if (variance_test_parametric$p.value < 0.05) {
				homoscedasticity_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
				} else if (variance_test_parametric$p.value >= 0.05) {
				homoscedasticity_parametric <- paste("Homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
				}
				if (variance_test_non_parametric$p.value < 0.05) {
				homoscedasticity_non_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_non_parametric$p.value, ")")
				} else if (variance_test_non_parametric$p.value >= 0.05) {
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
			## Return
			return(peak_stat_matrix)
		} else {
			## Return NULL
			return(NULL)
		}
	}
}





################################################################################


























































































############################################# SPECTRA

################## PLOT THE SIGNALS OF INTEREST WITH THE SD BARS ON THE AVERAGE
average_spectrum_bars_signals_of_interest <- function (spectra, SNR=5, signals_of_interest=peaks_average@mass, tolerance_ppm=2000, tof_mode="linear", half_window_plot=1000, graph_title="Spectrum", average_spectrum_colour="black", peak_points=TRUE, points_colour="red", bar_width=40, bar_colour="blue", peak_picking_algorithm="SuperSmoother") {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Outputs
spectrum_images <- list()
# Generate the average spectrum
average_spectrum <- averageMassSpectra(spectra, method="mean")
average_spectrum <- removeBaseline(average_spectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaks_average <- peak_picking(average_spectrum, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
# Peak picking on the dataset
peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
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
	st_dev_intensity <- sd(intensity_vector)
	####### Plot
	# Define the limits on the x-axis
	low_mass_plot <- signals_of_interest[a] - half_window_plot
	up_mass_plot <- signals_of_interest[a] + half_window_plot
	# Spectrum
	plot(average_spectrum, main=graph_title, col.main=average_spectrum_colour, xlab="m/z", ylab="Intensity (a.i.)", xlim=c(low_mass_plot,up_mass_plot), col=average_spectrum_colour)
	# Peaks
	if (peak_points == TRUE) {
		points(peaks_average[peak_index], pch=4, col=points_colour)
	}
	# Bars
	# epsilon: lengthof the horizontal segment
	epsilon = bar_width
	# Define the upper and the lower limit of the vertical segment
	up <- peaks_average@intensity[peak_index] + st_dev_intensity
	low <- peaks_average@intensity[peak_index] - st_dev_intensity
	# Vertical bar (x,y x,y)
	segments(signals_of_interest[a], low, signals_of_interest[a], up, col=bar_colour)
	# Horizontal segments(x,y , x,y)
	segments(signals_of_interest[a]-epsilon, low, signals_of_interest[a]+epsilon, low, col=bar_colour)
	segments(signals_of_interest[a]-epsilon, up, signals_of_interest[a]+epsilon, up, col=bar_colour)
    # Store the zoomed part of the spectrum with the signal
    avg_spectrum_zoom <- recordPlot()
    # Add this to a final list of images
    spectrum_images <- append(spectrum_images, avg_spectrum_zoom)
}
return (spectrum_images)
}





################################################################################





##################################################### MEMORY EFFICIENT IMPORTING
# This function imports the spectra in a memory efficient way: it reads spectra from one imzML file at a time, it can discard spectra according to their TIC, it runs the preprocessing of the spectra from the imzML file into packages of spectra, it can generate a set of representative average spectra (by grouping spectra randomly or according to a clustering algorithm). After this, it stores all the spectra from all the imzML files into a variable and from here it can align the spectra with the peaklist of the average spectrum of the dataset and it can crop all the spectra to a selected mass range.
# It relies upon other functions.
# The functions returns (the user can select what to compute) a list of elements: all the spectra, the representative spectra, the MS images after clustering.
memory_efficient_import <- function (folder, tof_mode="linear", tic_purification=FALSE, absolute_tic_threshold=0, smoothing_strength="medium", crop_spectra=FALSE, mass_range=NULL, spectra_preprocessing=TRUE, multicore_processing=TRUE, data_transformation=FALSE, transformation_algorithm="sqrt", peak_picking_algorithm="SuperSmoother", process_in_packages_of=length(spectra), generate_representative_spectra=FALSE, spectra_per_patient=1, algorithm_for_representative_spectra="hca", clustering_method_for_hca="agglomerative", discarded_nodes=1, skyline=FALSE, spectra_alignment=FALSE, spectra_alignment_method="cubic", spectra_format="imzml", seed=NULL, output_list=c("spectra","average","representative")) {
##### Load the required libraries
install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign"))
##### If it is not an imzML file...
if (length(grep(".imzML",folder, fixed=TRUE)) == 0) {
	##### Read the spectra files in the folder and in the subfolders
	folder_files <- read_spectra_files(folder, spectra_format=spectra_format, full_path=TRUE)
} else {
	folder_files <- folder
}
##### Tolerance - TOF-mode
if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
	tolerance_ppm <- 2000
} else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
	tolerance_ppm <- 200
}
##### Output initialisation
spectra_dataset <- list ()
spectra_dataset_grouped <- list()
spectra_dataset_average <- list()
spectra_dataset_skyline <- list()
spectra_dataset_representative <- list()
########## For each imzML file (patient)...
for (ff in 1:length(folder_files)) {
	##### Load the spectra
	if (spectra_format == "imzml" || spectra_format == "imzML") {
		# Crop the spectra while importing
		if (!is.null(mass_range)) {
			spectra <- importImzMl(folder_files[ff], massRange=mass_range)
		} else {
			spectra <- importImzMl(folder_files[ff])
		}
	}
	if (spectra_format == "brukerflex" || spectra_format == "xmass") {
		if (!is.null(mass_range)) {
			spectra <- importBrukerFlex(folder_files[ff], massRange=mass_range)
		} else {
			spectra <- importBrukerFlex(folder_files[ff])
		}
	}
	##### TIC purification
	if (tic_purification == TRUE) {
		spectra <- spectra_tic_purification(spectra, tof_mode=tof_mode, absolute_tic_threshold=absolute_tic_threshold)
	}
	##### Preprocessing (only if there is still some spectra left and only if the spectral dataset has to be returned) and add this purified spectra list to a global list (if there is still something after the TIC purification)
	if ("spectra" %in% output_list) {
		if (spectra_preprocessing == TRUE) {
			if (length(spectra) > 0) {
				spectra <- preprocess_spectra(spectra, tof_mode=tof_mode, smoothing_strength=smoothing_strength, process_in_packages_of=process_in_packages_of, align_spectra=spectra_alignment, spectra_alignment_method=spectra_alignment_method, multicore_processing=multicore_processing, data_transformation=data_transformation, transformation_algorithm=transformation_algorithm, crop_spectra=crop_spectra, mass_range=mass_range)
			}
		}
		spectra_dataset <- append(spectra_dataset, spectra)
	}
	##### Average the spectra (preprocessing not needed)
	if (generate_representative_spectra == TRUE && spectra_per_patient == 1 && "average" %in% output_list) {
		# At least two spectra for the averaging (the averaging will fail if there is only one spectrum, so use it if there is only one; after TIC purification no spectra can be left!) Check if there is any spectra and if there are one or many spectra.
		if (tof_mode == "linear") {
			if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
				spectra_avg <- averageMassSpectra(spectra, method="mean")
				spectra_avg <- smoothIntensity(spectra_avg, method="SavitzkyGolay", halfWindowSize=10)
				spectra_avg <- removeBaseline(spectra_avg, method="TopHat")
				spectra_avg <- calibrateIntensity(spectra_avg, method="TIC")
			}
		}
		if (tof_mode == "reflector" || tof_mode == "reflectron") {
			if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
				spectra_avg <- averageMassSpectra(spectra, method="mean")
				spectra_avg <- smoothIntensity(spectra_avg, method="SavitzkyGolay", halfWindowSize=2.5)
				spectra_avg <- removeBaseline(spectra_avg, method="TopHat")
				spectra_avg <- calibrateIntensity(spectra_avg, method="TIC")
			}
		}
		spectra_dataset_average <- append(spectra_dataset_average, spectra_avg)
	}
	##### Skyline spectrum (preprocessing not needed)
	if (generate_representative_spectra == TRUE && skyline == TRUE && "skyline" %in% output_list) {
		# At least two spectra for the skyline (the process will fail if there is only one spectrum, so use it if there is only one; after TIC purification no spectra can be left!)
		if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
			spectra_skyline <- generate_skyline_spectrum(spectra)
			spectra_skyline <- removeBaseline(spectra_skyline, method="TopHat")
			spectra_skyline <- calibrateIntensity(spectra_skyline, method="TIC")
		}
		if (length(spectra_syline) > 0) {
			spectra_dataset_syline <- append(spectra_dataset_skyline, spectra_skyline)
		}
	}
	##### Generate representative spectra (after TIC purification no spectra can be left!)
	if (generate_representative_spectra == TRUE && "representative" %in% output_list) {
		if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
			spectra_dataset_representative[[i]] <- group_spectra(spectra, spectra_per_patient=spectra_per_patient, spectra_format=spectra_format, tof_mode=tof_mode, seed=ifelse(is.null(seed), 0, seed), algorithm=algorithm_for_representative_spectra, clustering_method=clustering_method_for_hca, discarded_nodes=discarded_nodes)
			spectra_dataset_grouped <- append(spectra_dataset_grouped,spectra_dataset_representative[[i]]$spectra)
		}
	}
	### Free the memory
	rm(spectra)
	gc()
}
############################################################# SPECTRA ALIGNMENT
if (spectra_alignment == TRUE && "spectra" %in% output_list) {
	if (length(spectra_dataset) > 0 && isMassSpectrumList(spectra_dataset)) {
		# Average the spectra: the avg spectrum peaks will be used as reference
	    spectra_avg_ref <- averageMassSpectra(spectra_dataset, method="mean")
		peaks_avg_ref <- peak_picking(spectra_avg_ref, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
		reference_for_alignment <- peaks_avg_ref@mass
		if (length(reference_for_alignment) > 0) {
			if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
				spectra_dataset <- alignSpectra(spectra_dataset, halfWindowSize=20, noiseMethod=peak_picking_algorithm, SNR=3, reference=peaks_avg_ref, tolerance=(tolerance_ppm/10^6), warpingMethod=spectra_alignment_method)
			}
			if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
				spectra_dataset <- alignSpectra(spectra_dataset, halfWindowSize=5, noiseMethod=peak_picking_algorithm, SNR=3, reference=peaks_avg_ref, tolerance=(tolerance_ppm/10^6), warpingMethod=spectra_alignment_method)
			}
		}
	}
}
########## Return the values
return (list(spectra=spectra_dataset, spectra_dataset_grouped=spectra_dataset_grouped, spectra_average=spectra_dataset_average, spectra_skyline=spectra_dataset_skyline, spectra_dataset_representative=spectra_dataset_representative))
}





################################################################################





################################################################ SPECTRA BINNING
# The function performs the binning onto a selected spectra dataset (list of MALDIquant spectra objects)
resample_spectra <- function (spectra, final_data_points=lowest_data_points, binning_method="sum", multicore_processing=TRUE) {
####################################################### BINNING FUNCTION
binning_subfunction <- function (spectra, final_data_points, binning_method) {
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
		if (binning_method == "RMS" | binning_method == "rms") {
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
############# More spectra
if (isMassSpectrumList(spectra)) {
	# Load the required libraries
	install_and_load_required_packages("parallel")
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	########################
	# Calculate the lowest amount of data points, that corresponds to the maximum
	# number of data points that can be used for the binning
	# Use this value as default if final_data_points is not specified (Function)
	lowest_data_points <- NULL
	# Compare it with all the others
	for (s in 1:length(spectra)) {
		if (lowest_data_points == NULL || length(spectra[[s]]@mass) < lowest_data_points) {
			lowest_data_points <- length(spectra[[s]]@mass)
		}
	}
	# Check the qualities of the spectral dataset
	data_points_table <- table(sapply(spectra, length))
	datapoints_dataset <- as.numeric(names(data_points_table))
	equality_data_points <- length(data_points_table)
	######## Do not bin if all the spectra are of the same lengthand have the same number of datapoints as defined
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
			if (multicore_processing == TRUE) {
				if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
					spectra_binned <- mclapply(spectra, FUN=function (spectra) binning_subfunction(spectra, final_data_points, binning_method), mc.cores=cpu_thread_number)
				} else if (Sys.info()[1] == "Windows") {
					cl <- makeCluster(cpu_thread_number)
					# Pass the variables to the cluster for running the function
					clusterExport(cl=cl, varlist=c("final_data_points", "binning_method"), envir=environment())
					spectra_binned <- parLapply(cl, spectra, fun=function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
					stopCluster(cl)
				}
			} else {
				spectra_binned <- lapply(spectra, FUN=function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
			}
		}
		if (final_data_points <= lowest_data_points) {
			if (multicore_processing == TRUE) {
				if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
					spectra_binned <- mclapply(spectra, fun=function (spectra) binning_subfunction(spectra, final_data_points, binning_method), mc.cores=cpu_thread_number)
				} else if (Sys.info()[1] == "Windows") {
					cl <- makeCluster(cpu_thread_number)
					# Pass the variables to the cluster for running the function
					clusterExport(cl=cl, varlist=c("final_data_points", "binning_method"), envir=environment())
					spectra_binned <- parLapply(cl, spectra, fun=function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
					stopCluster(cl)
				}
			} else {
				spectra_binned <- lapply(spectra, fun=function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
			}
		}
	}
	print(table(sapply(spectra_binned, length)))
	print(paste("Equal distance between datapoints", (all(sapply(spectra_binned, isRegular)))))
} else {
	# Retrieve the number of datapoints
	lowest_data_points <- length(spectra@mass)
	# Perform the binning only if necessary
	if (final_data_points < lowest_data_points) {
		spectra_binned <- binning_subfunction(spectra, final_data_points, binning_method)
	} else {spectra_binned <- spectra}
}
return (spectra_binned)
}





################################################################################





######################################################### SAMPLE NAME REPLACING
# This function replaces the sample name field in the spectrum with the actual sample name (keeping only the last part of the file path and discarding the folder tree)
# The input can be both spectra or peaks (MALDIquant)
replace_sample_name <- function (spectra, spectra_format="imzml") {
#### imzML
if (spectra_format == "imzml" || spectra_format == "imzML") {
	# Scroll the spectra
	if (isMassSpectrumList(spectra) || isMassPeaksList(spectra)) {
		for (i in 1:length(spectra)) {
			if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
				# Split the filepath at /
				sample_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"/"))
				# The sample name is the last part of the path
				sample_name <- sample_name[length(sample_name)]
				# Detach the file extension
				sample_name <- unlist(strsplit(sample_name, ".imzML"))
				sample_name <- sample_name[1]
				# Put the name back into the spectra
				spectra[[i]]@metaData$file <- sample_name
			} else if (Sys.info()[1] == "Windows") {
				# Split the filepath at \
				sample_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"\\\\"))
				# The sample name is the last part of the path
				sample_name <- sample_name[length(sample_name)]
				# Detach the file extension
				sample_name <- unlist(strsplit(sample_name, ".imzML"))
				sample_name <- sample_name[1]
				# Put the name back into the spectra
				spectra[[i]]@metaData$file <- sample_name
			} else {
				# Keep the name unaltered
				sample_name <- spectra[[i]]@metaData$file[1]
			}
		}
	} else {
		if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
			# Split the filepath at /
			sample_name <- unlist(strsplit(spectra@metaData$file[1],"/"))
			# The sample name is the last part of the path
			sample_name <- sample_name[length(sample_name)]
			# Detach the file extension
			sample_name <- unlist(strsplit(sample_name, ".imzML"))
			sample_name <- sample_name[1]
			# Put the name back into the spectra
			spectra@metaData$file <- sample_name
		} else if (Sys.info()[1] == "Windows") {
			# Split the filepath at \
			sample_name <- unlist(strsplit(spectra@metaData$file[1],"\\\\"))
			# The sample name is the last part of the path
			sample_name <- sample_name[length(sample_name)]
			# Detach the file extension
			sample_name <- unlist(strsplit(sample_name, ".imzML"))
			sample_name <- sample_name[1]
			# Put the name back into the spectra
			spectra@metaData$file <- sample_name
		} else {
			# Keep the name unaltered
			sample_name <- spectra@metaData$file[1]
		}
	}
}
#### Xmass
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	# Scroll the spectra
	if (isMassSpectrumList(spectra) || isMassPeaksList(spectra)) {
		for (i in 1:length(spectra)) {
			sample_name <- spectra[[i]]@metaData$sampleName[1]
			# Put the name back into the spectra
			spectra[[i]]@metaData$file <- sample_name
		}
	} else {
		sample_name <- spectra@metaData$sampleName[1]
		# Put the name back into the spectra
		spectra@metaData$file <- sample_name
	}
}
return (spectra)
}





#########################################################################





########################################################## CLASS NAME REPLACING
# This function replaces the sample name field in the spectrum with the class the spectrum belongs to. If the filename contains the class, it replaces the filename with the class name, otherwise a class list symmetrical to the spectra must be provided.
# The input can be both spectra or peaks (MALDIquant)
replace_class_name <- function (spectra, class_list=NULL, class_in_file_name=TRUE, spectra_format="imzml") {
# If a class list is provided...
if (!is.null(class_list)) {
	# If the class name is in the name
	if (class_in_file_name == TRUE) {
		# Extract the file name from the path
		#spectra <- replace_sample_name(spectra, spectra_format=spectra_format)
		# Scroll the class list...
		for (w in 1:length(class_list)) {
    		# Scroll the spectra...
			for (i in 1:length(spectra)) {
        		# If there is a match between the file name and the class
				if (length(grep(class_list[w],spectra[[i]]@metaData$file[1], fixed=TRUE)) != 0) {
            		# Replace the file name with the class name
					spectra[[i]]@metaData$file <- class_list [w]
				}
			}
		}
	} else {
		# If the class name is not in the name
		# The class list and the spectra are symmetrical
		for (w in 1:length(class_list)) {
			spectra[[w]]@metaData$file <- class_list [w]
		}
	}
}
# If a class list is not provided...
if (is.null(class_list)) {
	# It's guessed from the filepath
	if (class_in_file_name == TRUE) {
		#### imzML
		if (spectra_format == "imzml" || spectra_format=="imzML") {
		    # Scroll the spectra...
		    for (i in 1:length(spectra)) {
		        # Split the filepath at / (universal for different OSs)
		        class_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"/"))
		        # The sample name is the last part of the path
		        class_name <- class_name[(length(class_name)-1)]
		        # Put the name back into the spectra
		        spectra[[i]]@metaData$file <- class_name
		    }
		}
		#### Xmass
		if (spectra_format == "brukerflex" || spectra_format == "xmass") {
		    # Scroll the spectra...
		    for (i in 1:length(spectra)) {
		        # Split the filepath at / (universal for different OSs)
		        class_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"/"))
		        # The sample name is the last part of the path
		        class_name <- class_name[length(class_name)-6]
		        # Put the name back into the spectra
		        spectra[[i]]@metaData$file <- class_name
		    }
		}
	} else {
		# The worst case returns the list of spectra as it is.
		spectra <- spectra
	}
}
return (spectra)
}





###############################################################################





######################################################## SPECTRA PRE-PROCESSING
# The function runs the preprocessing on the selected spectra (smoothing, baseline subtraction and normalisation). The function can be applied both to a spectra list or a single spectrum, allowing parallel computation.
# The function allows to select some additional parameters of the preprocessing.
# This version of the function whould be faster because each element of the spectral list is subjected to all the preprocessing step.
preprocess_spectra <- function (spectra, tof_mode="linear", preprocessing_parameters=list(crop_spectra=FALSE, mass_range=NULL, data_transformation=FALSE, transformation_algorithm="sqrt", smoothing_algorithm="SavitzkyGolay", smoothing_strength="medium", baseline_subtraction_algorithm="SNIP", baseline_subtraction_iterations=100, normalisation_algorithm="TIC", normalisation_mass_range=NULL), process_in_packages_of=length(spectra), align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=TRUE) {
##### Load the required libraries
install_and_load_required_packages(c("MALDIquant", "parallel"))
##### Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
##### Extract the parameters from the input list
crop_spectra <- preprocessing_parameters$crop_spectra
mass_range <- preprocessing_parameters$mass_range
data_transformation <- preprocessing_parameters$data_transformation
transformation_algorithm <- preprocessing_parameters$transformation_algorithm
smoothing_algorithm <- preprocessing_parameters$smoothing_algorithm
smoothing_strength <- preprocessing_parameters$smoothing_strength
baseline_subtraction_algorithm <- preprocessing_parameters$baseline_subtraction_algorithm
baseline_subtraction_iterations <- preprocessing_parameters$baseline_subtraction_iterations
normalisation_algorithm <- preprocessing_parameters$normalisation_algorithm
normalisation_mass_range <- preprocessing_parameters$normalisation_mass_range
##### Define the smoothing half wondow size
# Initialisation
smoothing_half_window_size <- NULL
if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
	#if (!is.null(smoothing_strength) && smoothing_strength == "small") {
		#if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			#smoothing_half_window_size <- 5
		#} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			#smoothing_half_window_size <- 1
		#}
	if (!is.null(smoothing_strength) && smoothing_strength == "medium") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 10
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 2
		}
	} else if (!is.null(smoothing_strength) && smoothing_strength == "strong") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 20
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 4
		}
	} else if (!is.null(smoothing_strength) && smoothing_strength == "stronger") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 30
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 6
		}
	}
} else if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
	#if (!is.null(smoothing_strength) && smoothing_strength == "small") {
		#if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			#smoothing_half_window_size <- 1
		#} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			#smoothing_half_window_size <- 0.2
		#}
	if (!is.null(smoothing_strength) && smoothing_strength == "medium") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 3
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 0.6
		}
	} else if (!is.null(smoothing_strength) && smoothing_strength == "strong") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 6
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 1.2
		}
	} else if (!is.null(smoothing_strength) && smoothing_strength == "stronger") {
		if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
			smoothing_half_window_size <- 9
		} else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
			smoothing_half_window_size <- 1.8
		}
	}
}
##### Generate the preprocessing function to be applied to every element of the spectra_temp list (x = spectrum)
preprocessing_subfunction <- function(x, crop_spectra, mass_range, data_transformation, transformation_algorithm, smoothing_algorithm, smoothing_half_window_size, baseline_subtraction_algorithm, baseline_subtraction_iterations, normalisation_algorithm, normalisation_mass_range) {
	### Remove flat spectra
	# x <- removeEmptyMassObjects (x)
	### Trimming
	if (crop_spectra == TRUE) {
		# Mass range specified
		if (!is.null(mass_range)) {
			x <- trim_spectra(x, range=mass_range)
		}
	}
	### Transformation
	if (data_transformation == TRUE) {
		x <- transformIntensity(x, method=transformation_algorithm)
	}
	### Smoothing
	if (!is.null(smoothing_algorithm)) {
		x <- smoothIntensity(x, method=smoothing_algorithm, halfWindowSize=smoothing_half_window_size)
	}
	### Baseline removal
	if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "TopHat") {
		x <- removeBaseline(x, method=baseline_subtraction_algorithm)
	} else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "SNIP") {
		# Default value for the number of iterations
		if (baseline_subtraction_iterations <= 0) {
			baseline_subtraction_iterations <- 100
		}
		x <- removeBaseline(x, method=baseline_subtraction_algorithm, iterations=baseline_subtraction_iterations)
	}
	### Normalisation
	if (normalisation_algorithm == "TIC") {
		if(!is.null(normalisation_mass_range)) {
			x <- calibrateIntensity(x, method=normalisation_algorithm, range=normalisation_mass_range)
		} else {
			x <- calibrateIntensity(x, method=normalisation_algorithm)
		}
	} else {
		x <- calibrateIntensity(x, method=normalisation_algorithm)
	}
	### Return the preprocessed spectrum (x)
	return(x)
}
######################################### Multiple spectra
if (isMassSpectrumList(spectra)) {
	### Trimming (same mass range for all the dataset)
	if (crop_spectra == TRUE && is.null(mass_range)) {
			spectra <- trim_spectra(spectra)
	}
	##### Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	##### Packages of preprocessing
	if (process_in_packages_of <= 0 || process_in_packages_of > length(spectra)) {
		process_in_packages_of <- length(spectra)
	}
	##### Create the list containing the processed spectra
	preprocessed_spectra <- list()
	index1 <- 1
	index2 <- process_in_packages_of
	spectra_packages <- ceiling(length(spectra) / process_in_packages_of)
	for (p in 1:spectra_packages) {
		## If the index 2 is more than the length of the spectra list, it has to be equal to the length of the list, it is not possible to go beyond the last element of the list
		if (index2 < length(spectra)) {
			spectra_temp <- spectra [index1:index2]
		} else {spectra_temp <- spectra [index1:length(spectra)]}
		## Fix the indexes at every cycle
		index1 <- index2 + 1
		index2 <- index2 + process_in_packages_of
		##################### Process the selected spectra (spectra_temp)
		##### Apply the function to the list of spectra_temp
		if (multicore_processing == TRUE) {
			if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
				spectra_temp <- mclapply(spectra_temp, FUN=function(spectra_temp) preprocessing_subfunction(spectra_temp, crop_spectra=crop_spectra, mass_range=mass_range, data_transformation=data_transformation, transformation_algorithm=transformation_algorithm, smoothing_algorithm=smoothing_algorithm, smoothing_half_window_size=smoothing_half_window_size, baseline_subtraction_algorithm=baseline_subtraction_algorithm, baseline_subtraction_iterations=baseline_subtraction_iterations, normalisation_algorithm=normalisation_algorithm, normalisation_mass_range=normalisation_mass_range), mc.cores=cpu_thread_number)
			} else if (Sys.info()[1] == "Windows") {
				cl <- makeCluster(cpu_thread_number)
				clusterEvalQ(cl, {library(MALDIquant)})
				clusterExport(cl=cl, varlist=c("crop_spectra", "mass_range", "data_transformation", "transformation_algorithm", "smoothing_algorithm", "smoothing_half_window_size", "baseline_subtraction_algorithm", "baseline_subtraction_iterations", "normalisation_algorithm", "normalisation_mass_range", "preprocessing_subfunction"), envir=environment())
				spectra_temp <- parLapply(cl, spectra_temp, fun=function(spectra_temp) preprocessing_subfunction(spectra_temp, crop_spectra=crop_spectra, mass_range=mass_range, data_transformation=data_transformation, transformation_algorithm=transformation_algorithm, smoothing_algorithm=smoothing_algorithm, smoothing_half_window_size=smoothing_half_window_size, baseline_subtraction_algorithm=baseline_subtraction_algorithm, baseline_subtraction_iterations=baseline_subtraction_iterations, normalisation_algorithm=normalisation_algorithm, normalisation_mass_range=normalisation_mass_range))
				stopCluster(cl)
			}
		} else {
				spectra_temp <- preprocessing_subfunction(spectra_temp, crop_spectra=crop_spectra, mass_range=mass_range, data_transformation=data_transformation, transformation_algorithm=transformation_algorithm, smoothing_algorithm=smoothing_algorithm, smoothing_half_window_size=smoothing_half_window_size, baseline_subtraction_algorithm=baseline_subtraction_algorithm, baseline_subtraction_iterations=baseline_subtraction_iterations, normalisation_algorithm=normalisation_algorithm, normalisation_mass_range=normalisation_mass_range)
		}
		########## Add to the final preprocessed spectral dataset
		preprocessed_spectra <- append(preprocessed_spectra, spectra_temp)
	}
}
######################################### Single spectra
if (!isMassSpectrumList(spectra) && isMassSpectrum(spectra)) {
	spectra <- preprocessing_subfunction(spectra, crop_spectra, mass_range, data_transformation, transformation_algorithm, smoothing_algorithm, smoothing_half_window_size, baseline_subtraction_algorithm, baseline_subtraction_iterations, normalisation_algorithm, normalisation_mass_range)
	########## Add to the final preprocessed spectral dataset
	preprocessed_spectra <- spectra
}
######################################### SPECTRAL ALIGNMENT
if (align_spectra == TRUE) {
	if (isMassSpectrumList(preprocessed_spectra)) {
		if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
			half_window_alignment <- 20
			tolerance_ppm <- 2000
		} else if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
			half_window_alignment <- 5
			tolerance_ppm <- 200
		}
		preprocessed_spectra <- alignSpectra(preprocessed_spectra, halfWindowSize=half_window_alignment, SNR=3, tolerance=(tolerance_ppm/10^6), warpingMethod=spectra_alignment_method)
	}
}
return (preprocessed_spectra)
}





###############################################################################





################################################### SPECTRA GROUPING (PATIENTS)
# The functions takes a list of already preprocessed spectra (MALDIquant) and generates a list of representative spectra, averaging spectra randomly or according to a cluster algorithm.
# It is advisable to use a list of spectra coming from one patient (one imzML)
# Spectra_per_patient = 1 returns the average spectrum of the spectra dataset.
# The function returns a list of elements: representative spectra, MS images of pixels under the same nodes coloured the same way, list of spectra under the discarded nodes (with their average spectrum), list of spectra grouped according to the node they belong to.
group_spectra <- function(spectra, spectra_per_patient=1, spectra_format="imzml", tof_mode="linear", seed=NULL, algorithm="random", clustering_method="agglomerative", discarded_nodes=1, balanced=TRUE, method="mean", peak_picking_algorithm="SuperSmoother") {
# Load the required libraries
install_and_load_required_packages (c("MALDIquant", "caret", "stats"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Output variables
patient_spectra_final <- list()
msi_plots <- list()
discarded_spectra <- list()
discarded_spectra_average <- list()
plots <- list()
spectra_hca_grouped <- list()
# Check if it is not a single spectrum
### Create the file Vector
if (spectra_format == "imzml" | spectra_format == "imzML") {
	# Put the filenames in a vector
	# Create the empty vector
	file_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	for (i in 1:length(spectra)) {
		file_vector <- append(file_vector, spectra[[i]]@metaData$file[1])
	}
	patient_vector <- unique(file_vector)
}
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
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
			if (spectra[[s]]@metaData$file[1] == patient_vector[[p]]) {
				spectra_patient <- append(spectra_patient, spectra[[s]])
			}
		}
		# If its lower than the reference (or not set yet)
		if (is.null(lowest_number_of_observations) || length(spectra_patient) < lowest_number_of_observations) {
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
		patient_spectra_final <- averageMassSpectra(spectra, labels=file_vector, method="mean")
	}
	if (method == "skyline") {
		patient_spectra_final <- group_spectra_skyline(spectra, spectra_format=spectra_format)
	}
}
# Run this script if there has to be two or more representative spectra per patient
if (spectra_per_patient > 1) {
	# If there is only one spectrum, use it
	if (length(spectra) > 0 && !isMassSpectrumList(spectra)) {
		patient_spectra_final <- spectra
	}
	# If there are more spectra...
	if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
		# Make the randomness reproducible
		if (!is.null(seed)) {
			set.seed (seed)
		}
		########################################################### OUTPUTS
		# Generate the final list of spectra
		patient_spectra_final <- list()
		discarded_spectra <- list()
		discarded_spectra_average <- list()
        spectra_hca_grouped <- list()
        spectra_kmeans_grouped <- list()
		# Create a new list of spectra for plotting purposes (the intensities will be replaced)
		spectra_for_plotting <- list()
		# List of additional plots
		plots <- list()
		# Generate the final list of MSI images
		msi_plots <- list()
		################################################# For each patient (p)
		for (p in 1:length(patient_vector)) {
			# Add the single patient spectra to a list
			patient_spectra <- list()
			for (s in 1:length(spectra)) {
				if (spectra[[s]]@metaData$file[1] == patient_vector[p]) {
					patient_spectra <- append(patient_spectra, spectra[[s]])
				}
			}
			############################################## RANDOMNESS
			if ("random" %in% algorithm) {
				# Do this if the spectra lengthis more than two, otherwise it does not make any sense
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
                # Make the randomness reproducible
        		if (!is.null(seed)) {
        			set.seed (seed)
        		}
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
				plotMsiSlice(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot()
			}
			######################################### SIMILARITY
			############################ HCA
			# Agglomerative
			if ((algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") && (clustering_method == "agglomerative" || is.null(clustering_method))) {
                spectra_hca_grouped_patient <- list()
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectra_per_patient <- spectra_per_patient + discarded_nodes
				# Detect and align peaks
				peaks <- peak_picking(patient_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
				peaks <- align_and_filter_peaks(peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
				# Generate the peaklist matrix
				peaklist <- intensityMatrix(peaks, patient_spectra)
				# Compute the distance matrix
				distance_matrix <- dist(peaklist, method="euclidean")
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
                ########## With discarded nodes
				if (discarded_nodes != 0) {
				# For each subgroup to be isolated...
					for (d in 1:discarded_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep(0, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# They do not get added to the final list of patient spectra
						# Add them to the discarded spectra
						discarded_spectra <- append(discarded_spectra, spectra_hca)
						# Average the spectra
						if (method == "mean") {
							spectra_average <- averageMassSpectra(spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum(spectra_hca)
						}
						# Add the average to the final list of discarded spectra AVG
						discarded_spectra_average <- append(discarded_spectra_average, spectra_average)
					}
					# For each subgroup to be isolated...
					for (d in (discarded_nodes+1):spectra_per_patient) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
                        spectra_hca_grouped_patient[[d]] <- spectra_hca
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep(d, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# Average the spectra in this HCA subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra(spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum(spectra_hca)
						}
						# Add the average to the final list
						patient_spectra_final <- append(patient_spectra_final, spectra_average)
					}
				} else {
                ########## Without discarded nodes
				# For each subgroup to be isolated...
				for (s in 1:spectra_per_patient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which(hca_groups == s)
					spectra_hca <- patient_spectra[index]
                    spectra_hca_grouped_patient[[s]] <- spectra_hca
					spectra_hca_for_plotting <- patient_spectra[index]
					# Replace the intensities with the S number for plotting purposes
					for (n in 1:length(spectra_hca_for_plotting)) {
						spectra_hca_for_plotting[[n]]@intensity <- rep (s, length(spectra_hca_for_plotting[[n]]@intensity))
					}
					# Add these modified spectra to the final list of spectra for plotting purposes
					spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
					# Average the spectra in this HCA subgroup
					if (method == "mean") {
						spectra_average <- averageMassSpectra(spectra_hca, method="mean")
					}
					if (method == "skyline") {
						spectra_average <- generate_skyline_spectrum(spectra_hca)
					}
					# Add the average to the final list
					patient_spectra_final <- append(patient_spectra_final, spectra_average)
				}
				}
				# Store the plot into the list
				plotMsiSlice(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot()
                # Spectra clustered
                spectra_hca_grouped[[p]] <- spectra_hca_grouped_patient
			}
			# K-MEANS
			if ((algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") && (clustering_method == "k_means" | clustering_method == "kmeans"| clustering_method == "k-Means")) {
                spectra_hca_grouped_patient <- list()
				# Adjust the spectra per patient value accordingo to how many nodes should be discarded
				spectra_per_patient <- spectra_per_patient + discarded_nodes
				# Detect and align peaks
				peaks <- peak_picking(patient_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
				peaks <- align_and_filter_peaks(peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
				# Generate the peaklist matrix
				peaklist <- intensityMatrix(peaks, patient_spectra)
				# Compute the k-Means clustering
				hca <- kmeans(peaklist, centers=spectra_per_patient)
				# Associate to each row/spectrum the subgroup/cluster to which it belongs
				hca_groups <- hca$cluster
				if (discarded_nodes != 0) {
				# For each subgroup to be isolated (to discard)...
					for (d in 1:discarded_nodes) {
						# Index the spectra under in the selected subgroup
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep(0, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# They do not get averaged and added to the final list of patient spectra
						# Add them to the discarded spectra
						discarded_spectra <- append(discarded_spectra, spectra_hca)
						# Average the spectra
						if (method == "mean") {
							spectra_average <- averageMassSpectra(spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum(spectra_hca)
						}
						# Add the average to the final list of discarded spectra AVG
						discarded_spectra_average <- append(discarded_spectra_average, spectra_average)
					}
					# For each subgroup to be isolated (to keep)...
					for (d in (discarded_nodes+1):spectra_per_patient) {
						# Index the spectra under in the selected subgroup
						index <- which(hca_groups == d)
						spectra_hca <- patient_spectra[index]
                        spectra_hca_grouped_patient[[d]] <- spectra_hca
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the D number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep(d, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra(spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum(spectra_hca, method="mean")
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
                        spectra_hca_grouped_patient[[s]] <- spectra_hca
						spectra_hca_for_plotting <- patient_spectra[index]
						# Replace the intensities with the S number for plotting purposes
						for (n in 1:length(spectra_hca_for_plotting)) {
							spectra_hca_for_plotting[[n]]@intensity <- rep(s, length(spectra_hca_for_plotting[[n]]@intensity))
						}
						# Add these modified spectra to the final list of spectra for plotting purposes
						spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_for_plotting)
						# Average the spectra in this subgroup
						if (method == "mean") {
							spectra_average <- averageMassSpectra(spectra_hca, method="mean")
						}
						if (method == "skyline") {
							spectra_average <- generate_skyline_spectrum(spectra_hca, method="mean")
						}
						# Add the average to the final list
						patient_spectra_final <- append(patient_spectra_final, spectra_average)
					}
					}
				# Store the plot into the list
				plotMsiSlice(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
				msi_plots [[p]] <- recordPlot()
                # Clustered spectra
                spectra_hca_grouped <- spectra_hca_grouped_patient
			}
			#####################
		}
	}
}
# Processing before returning
patient_spectra_final <- removeBaseline(patient_spectra_final, method="TopHat")
patient_spectra_final <- calibrateIntensity(patient_spectra_final, method="TIC")
#
return (list(spectra=patient_spectra_final, ms_images=msi_plots, discarded_spectra=discarded_spectra, discarded_average_spectra=discarded_spectra_average, plots=plots, spectra_hca_grouped=spectra_hca_grouped))
}





################################################################################





#################################################### SPECTRA GROUPING (CLASSES)
# The functions takes a list of already preprocessed spectra (MALDIquant) and generates a list of representative spectra, averaging spectra according to the class they belong to, generating one average spectrum per class.
group_spectra_class <- function (spectra, class_list, method="mean", spectra_format="imzml", class_in_file_name=TRUE) {
class_list <- sort(class_list)
####### IMZML
if (spectra_format == "imzml" || spectra_format == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name(spectra, class_list=class_list, spectra_format=spectra_format, class_in_file_name=class_in_file_name)
	# Put the filenames/classes in a vector
	# Create the empty vector
	class_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	if (isMassSpectrumList(spectra)) {
		for (i in 1:length(spectra)) {
			class_vector <- append(class_vector, spectra[[i]]@metaData$file[1])
		}
	} else {
		class_vector <- spectra@metaData$file[1]
	}
}
####### XMASS
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	# Replace the names in the spectra with the class
	for (cl in 1:length(class_list)) {
		if (isMassSpectrumList(spectra)) {
			for (s in 1:length(spectra)) {
				if (length(grep(class_list[cl], spectra[[s]]@metaData$file[1], fixed=TRUE)) > 0) {
					spectra[[s]]@metaData$file <- class_list[cl]
				}
			}
		} else {
			if (length(grep(class_list[cl], spectra@metaData$file[1], fixed=TRUE)) > 0) {
				spectra@metaData$file <- class_list[cl]
			}
		}
	}
	# Create the empty vector
	class_vector <- character()
	# Add the file names recursively, scrolling the whole spectral dataset
	if (isMassSpectrumList(spectra)) {
		for (i in 1:length(spectra)) {
			class_vector <- append(class_vector, spectra[[i]]@metaData$file[1])
		}
	} else {
		class_vector <- spectra@metaData$file[1]
	}
}
# Average
if (method == "mean") {
	class_spectra_grouped <- averageMassSpectra(spectra, labels=class_vector, method="mean")
}
# Skyline
if (method == "skyline") {
	class_spectra_grouped <- list()
	# For each class
	for (class in class_list) {
		# Add the spectra for each class to a list
		for (spectrum in spectra) {
			if (spectrum@metaData$file[1] == class) {
				class_spectra <- append(class_spectra, spetrum)
			}
		}
		# Skyline spectrum
		class_spectra_skyline <- generate_skyline_spectrum(class_spectra)
		# Add this to the final list
		class_spectra_grouped <- append(class_spectra_grouped, class_spectra_skyline)
	}
}
return (class_spectra_grouped)
}





###############################################################################





######################## PLOT THE MEAN SPECTRUM WITH THE SD BARS ON THE AVERAGE
# This function takes a list of spectra (MALDIquant) as input and returns the average spectrum of the provided dataset with bars onto the peaks, after calculating the standard deviation.
average_spectrum_bars <- function (spectra, SNR=5, peak_picking_algorithm="SuperSmoother", tolerance_ppm=2000, mass_range_plot=c(4000,15000), graph_title="Spectrum", average_spectrum_colour="black", peak_points="yes", points_colour="red", bar_width=40, bar_colour="blue") {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Generate the average spectrum
average_spectrum <- averageMassSpectra(spectra, method="mean")
average_spectrum <- removeBaseline(average_spectrum, method="TopHat")
# Peak picking on the average spectrum (for plotting)
peaks_average <- peak_picking(average_spectrum, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
# Peak picking on the dataset
peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
# Alignment: merge the two lists and align them all
peaks_all <- append(peaks, peaks_average)
peaks_all <- align_and_filter_peaks(peaks_all, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1)
# Empty the lists
peaks_average <- list()
peaks <- list()
# Re-fill the lists with the aligned peaklists
for (i in 1:(length(peaks_all) - 1)) {
	peaks <- append(peaks, peaks_all[[i]])
}
peaks_average <- peaks_all[[length(peaks_all)]]
######## Generate a vector with the standard deviations of the peaks in the average peaklist
st_dev_intensity <- vector(length=0)
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
	st_dev_intensity <- append(st_dev_intensity, sd(intensity_vector))
}
####### Plot
# Spectrum
plot(average_spectrum, main=graph_title, col.main=average_spectrum_colour, xlab="m/z", ylab="Intensity (a.i.)", xlim=mass_range_plot, col=average_spectrum_colour)
# Peaks
if (peak_points == "yes") {
	points(peaks_average, pch=4, col=points_colour)
}
# Bars
# epsilon: length of the horizontal segment
epsilon = bar_width
for (i in 1:length(peaks_average@mass)) {
	# Define the upper and the lower limit of the vertical segment
	up <- peaks_average@intensity[i] + st_dev_intensity [i]
	low <- peaks_average@intensity[i] - st_dev_intensity [i]
	# Vertical bar (x,y x,y)
	segments(peaks_average@mass[i], low, peaks_average@mass[i], up, col=bar_colour)
	# Horizontal segments(x,y , x,y)
	segments(peaks_average@mass[i]-epsilon, low, peaks_average@mass[i]+epsilon, low, col=bar_colour)
	segments(peaks_average@mass[i]-epsilon, up, peaks_average@mass[i]+epsilon, up, col=bar_colour)
}
avg_spectrum_with_bars <- recordPlot()
return(avg_spectrum_with_bars)
}





################################################################################





############################################# MOST INTENSE PEAKS IN PEAK PICKING
# This function returns a peak list containing only the most intense peaks per spectrum. If the input is a list of spectra, the function computes the peak picking and keeps only the most intense ones, if it's a list of peaklists, it applies the filtering function directly on the peaks.
most_intense_signals <- function (spectra, signals_to_take=20, tof_mode="linear", peak_picking_algorithm="SuperSmoother", multicore_processing=TRUE) {
# Load the required libraries
install_and_load_required_packages(c("parallel", "MALDIquant"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
####################################################### PICKING FUNCTION
picking_subfunction <- function (peaks, signals_to_take) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, SNR = peaks@snr)
	# Check if the provided number does not exceed the number of available signals
	if (signals_to_take > nrow(peaks_data_frame) || signals_to_take <= 0) {
		signals_to_take <- nrow(peaks_data_frame)
	}
	# Sort the dataframe according to the SNR
	peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$SNR),]
	# Select only the first most intense signals
	selected_signals <- peaks_data_frame [1:signals_to_take,]
	# Sort the dataframe back according to mass
	selected_signals <- selected_signals [order(selected_signals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selected_signals$mass
	peaks@intensity <- selected_signals$intensity
	peaks@snr <- selected_signals$SNR
	return (peaks)
}
########################################################################
# Peak picking
if (isMassSpectrumList(spectra)) {
	peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3, multicore_processing=multicore_processing)
} else if (isMassPeaksList(spectra)) {
peaks <- spectra
}
# Most intense signals
if (isMassPeaksList(peaks)) {
	if (multicore_processing == TRUE) {
		# Detect the number of cores
		cpu_thread_number <- detectCores(logical=TRUE) - 1
		if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
			most_intense_peaks <- mclapply(peaks, FUN = function(peaks) picking_subfunction(peaks, signals_to_take=signals_to_take), mc.cores=cpu_thread_number)
		} else if (Sys.info()[1] == "Windows") {
			# Make the CPU cluster for parallelisation
			cl <- makeCluster(cpu_thread_number)
			# Pass the variables to the cluster for running the function
			clusterExport(cl=cl, varlist="signals_to_take", envir=environment())
			most_intense_peaks <- parLapply(cl, peaks, fun = function(peaks) picking_subfunction(peaks, signals_to_take=signals_to_take))
			stopCluster(cl)
		}
	} else {
		most_intense_peaks <- lapply(peaks, FUN = function(peaks) picking_subfunction(peaks, signals_to_take=signals_to_take))
	}
} else if (isMassPeaks(peaks)) {
	most_intense_peaks <- picking_subfunction(peaks, signals_to_take)
}
return (most_intense_peaks)
}





################################################################################





############################################# AVERAGE THE REPLICATES (BY FOLDER)
# This function averages the spectra contained in the same folder (more suitable for brukerflex format)
average_replicates_by_folder <- function (spectra, folder, spectra_format="brukerflex") {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Count the number of spectra
if (isMassSpectrumList(spectra)) {
	number_of_spectra <- length(spectra)
} else {
	number_of_spectra <- 1
}
# List the spectra files (list the imzML files or the FID files)
folder_files <- read_spectra_files(folder, spectra_format=spectra_format, full_path=TRUE)
if (spectra_format == "brukerflex") {
	# Split the path into individual folders (list, each element is a vector with the path splitted for that spectrum)
	folder_files_splitted <- list()
	# Split the paths into folders
	for (f in 1:number_of_spectra) {
		folder_files_splitted[f] <- strsplit(folder_files[f], "/")
	}
	#### SITUATION #1
	# The n-5 folder (where n is the length of the folder_splitted) is not unique, since it can represent the treatment folder; the n-4 folder is not unique, since one can call the replicates all in the same way (e.g. rep1, rep2), because they are placed in the sample folder; the combination n-5/n-4 should be almost unique, but to guarantee the uniqueness it is better to use the n-6/n-5/n-4 combination, so that it can be folder/sample/replicate or sample/treatment/replicate.
	#### SITUATION #2
	# The second case is when the sample folder is the same, but the coordinates change, so the replicates are the spectra under the same coordinate folder, while different coordinates mean different samples.
	# The folder in which the samples are should be the same, there shouldn't be a more complicate folder structure, in that case this method should be revised.
	# n-6 = sample name
	sample_name <- character()
	for (f in 1:number_of_spectra) {
		sample_name <- append(sample_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-6])
	}
	# n-5 = treatment / class
	treatment_name <- character()
	for (f in 1:number_of_spectra) {
		treatment_name <- append(treatment_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-5])
	}
	# n-4 = replicate_name
	#replicate_name <- character()
	#for (f in 1:number_of_spectra) {
		#replicate_name <- append(replicate_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-4])
	#}
	# n-3 = coordinates
	#coordinate_name <- character()
	#for (f in 1:number_of_spectra) {
		#coordinate_name <- append(coordinate_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-3])
	#}
	# n-2 = replicate_coordinates
	#replicate_coordinates <- character()
	#for (f in 1:number_of_spectra) {
		#replicate_coordinates <- append(replicate_coordinates, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-2])
	#}
	# Generate a file_vector with only the folder to identify the replicates
	unique_sample_name <- character(length=number_of_spectra)
	# Check for replicates in the same folder
	for (f in 1:(number_of_spectra-1)) {
		### Situation #1 or Situation #2
		if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] == treatment_name[f])) {
			unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], sep="/")
			unique_sample_name[f+1] <- paste(sample_name[f], treatment_name[f], sep="/")
		} else {
			unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], sep="/")
			unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], sep="/")
		}
		#if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] == treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
			### Situation #1
			#replicate_name[f+1] <- replicate_name[f]
			#unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], replicate_name[f], sep="/")
			#unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep="/")
		#} else if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] != treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
			#unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], replicate_name[f], sep="/")
			#unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep="/")
		#} else if ((treatment_name[f+1] == treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
			# Situation #1 + Situation #2
			#unique_sample_name[f] <- replicate_name[f]
			#unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep="/")
		#}
	}
} else if (spectra_format == "imzML" || spectra_format == "imzml") {
	unique_sample_name <- character()
	for (f in 1:number_of_spectra) {
		unique_sample_name <- append(unique_sample_name, spectra[[f]]@metaData$file[1])
	}
}
# Average the mass spectra, grouping them according to the sample_vector
spectra_replicates_averaged <- averageMassSpectra(spectra, labels=unique_sample_name)
return (spectra_replicates_averaged)
}





################################################################################





################################################################### PEAK PICKING
# This function takes a list of spectra (MALDIquant) and computes the peak picking.
peak_picking <- function(spectra, peak_picking_algorithm="SuperSmoother", tof_mode="linear", SNR=3, multicore_processing=TRUE) {
	########## Load the required libraries
	install_and_load_required_packages(c("MALDIquant", "parallel"))
	########## Multi-core
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	##### TOF-MODE
	if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
		half_window_size <- 20
	} else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
		half_window_size <- 5
	}
	########## SINGLE SPECTRUM
	if (isMassSpectrum(spectra)) {
		peaks <- detectPeaks(spectra, method=peak_picking_algorithm, halfWindowSize=half_window_size, SNR=SNR)
	}
	########## MULTIPLE SPECTRA
	if (isMassSpectrumList(spectra)) {
		# Peak detection
		peaks <- list()
		if (multicore_processing == TRUE) {
			if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
				peaks <- mclapply(spectra, FUN=function (spectra) detectPeaks(spectra, method=peak_picking_algorithm, halfWindowSize=half_window_size, SNR=SNR), mc.cores=cpu_thread_number)
			} else if (Sys.info()[1] == "Windows") {
				# Make the cluster (one for each core/thread)
				cl <- makeCluster(cpu_thread_number)
				clusterEvalQ(cl, {library(MALDIquant)})
				clusterExport(cl=cl, varlist=c("peak_picking_algorithm", "half_window_size", "SNR"), envir=environment())
				peaks <- parLapply(cl, spectra, fun=function (spectra) detectPeaks(spectra, method=peak_picking_algorithm, halfWindowSize=half_window_size, SNR=SNR))
				stopCluster(cl)
			}
		} else {
			peaks <- detectPeaks(spectra, method=peak_picking_algorithm, halfWindowSize=half_window_size, SNR=SNR)
		}

	}
	return (peaks)
}




################################################################# PEAK ALIGNMENT
# This function takes a list of peaks (MALDIquant) and computes the peak alignment, along with the false positive removal and the removal of low-intensity peaks.
align_and_filter_peaks <- function (peaks, peak_picking_algorithm="SuperSmoother", tof_mode="linear", peaks_filtering=TRUE, frequency_threshold_percent=25, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1, intensity_threshold_method="element-wise", reference_peaklist=NULL, spectra=NULL, alignment_iterations=5, multicore_processing=TRUE) {
########## Determine the tolerance in PPM
if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
	tolerance_ppm <- 2000
} else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
	tolerance_ppm <- 200
}
########## Align only if there are many peaklists
if (isMassPeaksList(peaks)) {
	##### Load the required libraries
	install_and_load_required_packages(c("MALDIquant", "parallel"))
	##### Rename the trim function
	trim_spectra <- get(x="trim", pos="package:MALDIquant")
	###### Peak alignment
	# Fix the iteration number
	if (alignment_iterations <= 0) {
		alignment_iterations <- 5
	}
	# Initialize the variable
	peaks_aligned <- NULL
	# For each iteration...
	for (iter in 1:alignment_iterations) {
		# If it is the first time that the alignment is run...
		if (is.null(peaks_aligned)) {
			# Run the alignment
			peaks_aligned <- binPeaks(peaks, method="relaxed", tolerance=(tolerance_ppm/10^6))
		} else {
			# Otherwise run another alignment on the already aligned peaks
			peaks_aligned <- binPeaks(peaks_aligned, method="relaxed", tolerance=(tolerance_ppm/10^6))
		}
	}
	##### False positive removal
	if (peaks_filtering == TRUE) {
		peaks_aligned <- filterPeaks(peaks_aligned, minFrequency=(frequency_threshold_percent/100))
	}
	##### Low-intensity peaks removal
	if (low_intensity_peaks_removal == TRUE) {
		peaks_aligned <- remove_low_intensity_peaks(peaks_aligned, intensity_threshold_percent=intensity_threshold_percent, intensity_threshold_method=intensity_threshold_method, multicore_processing=multicore_processing)
	}
	##### Align to a reference peaklist: AVERAGE SPECTRUM (if a spectra list is provided)
	if (is.character(reference_peaklist) && reference_peaklist == "average" && !is.null(spectra)) {
		# Average the spectra
		average_spectrum <- averageMassSpectra(spectra, method="mean")
		# Peak picking
		average_spectrum_peaks <- peak_picking(average_spectrum, peak_picking_algorithm=peak_picking_algorithm, SNR=5, multicore_processing=multicore_processing)
		reference_peaklist <- average_spectrum_peaks@mass
	} else if (is.character(reference_peaklist) && reference_peaklist == "average" && is.null(spectra)) {
		reference_peaklist <- NULL
	}
	##### Align to the reference peaklist
	if (!is.null(reference_peaklist)) {
		############ Function for lapply
		align_peaks_subfunction <- function (peaks, reference_peaklist, tolerance_ppm) {
			mass_vector <- peaks@mass
			# For each reference peak
			for (ref in reference_peaklist) {
				# Replace the value in the vector with the reference value
				mass_vector[which(abs(mass_vector - ref)*10^6/ref <= tolerance_ppm)] <- ref
			}
			# Put the fixed mass vector back into the MALDIquant peaklist
			peaks@mass <- mass_vector
			return (peaks)
		}
		############# If there are many peaklists or one peaklist (use multicore)
		if (isMassPeaksList(peaks_aligned)) {
			if (multicore_processing == TRUE) {
				# Detect the number of cores
				cpu_thread_number <- detectCores(logical=TRUE) - 1
				if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
					peaks_aligned <- mclapply(peaks_aligned, FUN=function(peaks_aligned) align_peaks_subfunction(peaks_aligned, reference_peaklist, tolerance_ppm), mc.cores=cpu_thread_number)
				} else if (Sys.info()[1] == "Windows") {
					# Make the CPU cluster for parallelisation
					cl <- makeCluster(cpu_thread_number)
					# Apply the multicore function
					# Pass the variables to the cluster for running the function
					clusterExport(cl=cl, varlist=c("reference_peaklist", "tolerance_ppm"), envir=environment())
					peaks_aligned <- parLapply(cl, peaks_aligned, fun=function(peaks_aligned) align_peaks_subfunction(peaks_aligned, reference_peaklist, tolerance_ppm))
					stopCluster(cl)
				}
			} else {
				peaks_aligned <- lapply(peaks_aligned, FUN=function(peaks_aligned) align_peaks_subfunction(peaks_aligned, reference_peaklist, tolerance_ppm))
			}
		} else {
			peaks_aligned <- align_peaks_subfunction(peaks_aligned, reference_peaklist, tolerance_ppm)
		}
	}
	return (peaks_aligned)
} else {
	# Low-intensity peaks removal
	if (low_intensity_peaks_removal == TRUE) {
		peaks <- remove_low_intensity_peaks(peaks, intensity_threshold_percent=intensity_threshold_percent, intensity_threshold_method=intensity_threshold_method, multicore_processing=multicore_processing)
	}
	return (peaks)
}
}





################################################################################





################################################ LIBRARY CREATION
# This function reads the files contained in a provided folder (no memory efficient importing), it can average them according to the class they belong to or according to the folder they are into (average replicate/patient). It also computes the peak picking and can replace the SNR field in the peaklist with the standard deviation of the peaks or the coefficient of variation.
# It returns a list containing: the spectra and the peaks (NOT aligned).
# It is used to create the library/database.
library_creation <- function (filepath_database, peak_picking_algorithm="SuperSmoother", class_grouping=TRUE, mass_range=c(3000,15000), spectra_preprocessing=list(smoothing_strength="medium", preprocess_in_packages_of=length(spectra)), average_replicates=FALSE, average_patients=FALSE, SNR=5, most_intense_peaks=FALSE, signals_to_take=20, reference_peaklist_for_alignment=NULL, tof_mode="linear", spectra_format="brukerflex") {
# Load the required libraries
install_and_load_required_packages(c("MALDIquantForeign", "MALDIquant"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
### Define the classes (the classes are the folders in the library directory)
class_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
# TOF mode
if (tof_mode == "linear") {
	tolerance_ppm <- 2000
}
if (tof_mode == "reflectron" || tof_mode == "reflector") {
	tolerance_ppm <- 200
}
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	### Load the spectra
	if (!is.null(mass_range)) {
		spectra <- importBrukerFlex(filepath_database, massRange=mass_range)
	} else {
		spectra <- importBrukerFlex(filepath_database)
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	### Load the spectra
	if (!is.null(mass_range)) {
		spectra <- importImzMl(filepath_database, massRange=mass_range)
	} else {
		spectra <- importImzMl(filepath_database)
	}
}
### Average the replicates
if (average_replicates == TRUE) {
	spectra <- average_replicates_by_folder(spectra, filepath_database, spectra_format=spectra_format)
}
# Average the patients
if (average_patients == TRUE) {
	spectra <- group_spectra(spectra, spectra_per_patient=1, spectra_format=spectra_format, tof_mode=tof_mode)
}
### Class grouping
if (class_grouping == TRUE) {
	### Spectra grouping (class)
	spectra <- group_spectra_class(spectra, class_list=class_list, spectra_format=spectra_format, class_in_file_name=TRUE)
} else {spectra <- spectra}
### Preprocessing
spectra <- preprocess_spectra(spectra, tof_mode=tof_mode, smoothing_strength=spectra_preprocessing$smoothing_strength, process_in_packages_of=spectra_preprocessing$preprocess_in_packages_of, sampling_method="sequential")
##########################
### Peak picking on the individual spectra
if (most_intense_peaks == TRUE) {
	peaks_database <- most_intense_signals(spectra, signals_to_take=signals_to_take)
} else {
	peaks_database <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	#peaks_database <- align_and_filter_peaks(peaks_database, tof_mode=tof_mode, peaks_filtering=FALSE, low_intensity_peaks_removal=FALSE, reference_peaklist=reference_peaklist_for_alignment)
}
####
library_list <- list(spectra = spectra, peaks = peaks_database)
return(library_list)
}





################################################################################





































































































################################################################ CLASSIFICATION

#################################### CLASSIFICATION VIA HIERARCHICAL CLUSTERING
# The input is a S4 list of spectra (MALDIquant), both for the classification and for the database, the class names (if in the spectrum file name, the filenames will be replaced with the class name and averaged to create representative spectra for the classes; if not in the files, one representative spectrum per class should be provided in the same order as the class list), the TOF-MS mode, the algorithm for the clustering and the nodes to discard.
# The script returns a list containing: the hierarchical analysis dendrograms, a matrix with the classification (of pixels/spectra and samples) and the MS image with pixel or area classification, along with the classification of the patients.
hierarchical_clustering_classification <- function (spectra_to_be_classified, spectra_database, peak_picking_algorithm="SuperSmoother", class_list=c("HP","PTC"), class_in_file_name=TRUE, tof_mode="linear", clustering_method="agglomerative", nodes=4, discarded_nodes=0, seed=NULL, classification_of=c("average","subareas","pixels"), true_class_in_file_name=FALSE, true_class_list=NULL) {
# Do everything if the user wants some outputs
if (!is.null(classification_of) && length(classification_of) > 0) {
################################################## Load the required libraries
install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign", "stats", "caret", "pROC"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
########################################################### Outputs
classification_hca_results <- NULL
classification_hca_results_pixel_by_pixel <- NULL
classification_hca_results_avg <- NULL
msi_plots <- list()
msi_plots_pixel_by_pixel <- list()
hca_classification_performances <- list()
hca_roc <- NULL
hca_classification_performances <- NULL
################### DATABASE
# Replace the filename with the class name in the spectra for the database
if (class_in_file_name == TRUE) {
    class_list <- sort(class_list)
    # Replace the filename with the class name in the spectra database
    spectra_database <- replace_class_name(spectra_database, class_list)
}
# Generate representative spectra for each class
spectra_database <- group_spectra_class(spectra_database, class_list=class_list, spectra_format="imzML")
########################## PATIENT SPECTRA
# Fix the filenames
spectra_to_be_classified <- replace_sample_name(spectra_to_be_classified, spectra_format="imzML")
# Create the empty vector
file_vector <- character()
# Add the file names recursively, scrolling the whole spectral dataset
for (i in 1:length(spectra_to_be_classified)) {
	file_vector <- append(file_vector, spectra_to_be_classified[[i]]@metaData$file[1])
}
# Compute the patient (imzML file) vector
patient_vector <- unique(file_vector)
######################## For each patient...
for (p in 1:length(patient_vector)) {
    patient_spectra <- NULL
    patient_spectra_number <- 0
    # Isolate the patient spectra
    for (s in 1:length(spectra_to_be_classified)) {
        # If the spectrum filename corresponds to the patient
        if (spectra_to_be_classified[[s]]@metaData$file == patient_vector[p]) {
            # Add the selected spectra to the list of the patient spectra
			if (is.null(patient_spectra)) {
				patient_spectra <- spectra_to_be_classified[[s]]
				patient_spectra_number <- patient_spectra_number + 1
			} else {
				patient_spectra <- append(patient_spectra, spectra_to_be_classified[[s]])
				patient_spectra_number <- patient_spectra_number + 1
			}
        }
    }
    ################## CLASSIFICATION OF ENTIRE SUB-AREAS
	if ("subareas" %in% classification_of) {
    #### Group spectra
    patient_spectra_grouped <- group_spectra(patient_spectra, spectra_per_patient=nodes, spectra_format="imzml", tof_mode=tof_mode, seed=seed, algorithm="hca", clustering_method=clustering_method, discarded_nodes=discarded_nodes, balanced=TRUE, method="mean")
    patient_spectra_clustered_average <- patient_spectra_grouped$spectra
    patient_spectra_clustered <- patient_spectra_grouped$spectra_hca_grouped[[1]]
    # Put the database and the spectra together
    global_spectra <- append(spectra_database, patient_spectra_clustered_average)
	global_peaks <- peak_picking(global_spectra, SNR=SNR, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode)
	global_peaks <- align_and_filter_peaks(global_peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
    # Compute the intensity matrix
    intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
    intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks=global_peaks, class_list=class_list, spectra_format="imzml", sample_output=TRUE, class_output=FALSE)
    rownames(intensity_matrix) <- intensity_matrix[,"Sample"]
	# Compute the agglomerative hierarchical clustering
	if (clustering_method == "agglomerative") {
	    distance_matrix <- dist(intensity_matrix[,1:(ncol(intensity_matrix)-1)])
	    hierarchical_clustering <- hclust(distance_matrix)
	    plot(hierarchical_clustering)
	    hca_graph <- recordPlot()
	    distance_matrix <- as.matrix(distance_matrix)
	    # Remove the first rows (database spectra) and keep only the first columns (database spectra)
	    distance_matrix <- distance_matrix [(length(spectra_database)+1):nrow(distance_matrix),(1:length(spectra_database))]
	    ################### Classification output
	    ## Function for matrix apply
	    classification_hca_function <- function (dist_row) {
	        # Set the closest value to NULL
	        closest <- NULL
	        # Scroll the row...
	        for (l in 1:length(dist_row)) {
	            # If closest has not been set yet or the value in the cell is the closest...
	            if (is.null(closest) || dist_row[l] < closest){
	                # Set the closest value to be the one in the cell
	                closest <- dist_row[l]
	            }
	        }
	        # Find the position of the row where the closest is
	        closest_position <- which(dist_row == closest)
	        # The belonging class is corresponding to the element of the database which the spectrum is closer to
	        belonging_class <- names(dist_row)[closest_position]
	        return (belonging_class)
	    }
	    # Patient classification
	    classification_hca  <- cbind(apply(distance_matrix, MARGIN=1, FUN=function(x) classification_hca_function(x)))
	    ### Generate the final output results
	    if (is.null(classification_hca_results)) {
	        classification_hca_results <- classification_hca
	    } else {classification_hca_results <- rbind(classification_hca_results, classification_hca)}
		colnames(classification_hca_results) <- "Classification"
	    ####################################### MS IMAGES
	    spectra_for_plotting <- list()
	    # Correspondence between the classificatio_hca matrix classification and the spectra under the hca nodes
	    for (n in 1:nrow(classification_hca)) {
	        # Isolate the spectra under the hca node...
	        spectra_hca_node <- patient_spectra_clustered[[n]]
	        # Replace the intensity with the classification
	        for (s in 1:length(spectra_hca_node)) {
	            for (z in 1:length(class_list)) {
	                if (classification_hca[n,1] == class_list[z]) {
	                    spectra_hca_node[[s]]@intensity <- rep(z, length(spectra_hca_node[[s]]@intensity))
	                }
	            }
	        }
	        spectra_for_plotting <- append(spectra_for_plotting, spectra_hca_node)
	    }
	    plotMsiSlice(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
	    legend(x="bottomright", legend=class_list, fill=c("green","red"))
	    legend(x="topright", legend=spectra_for_plotting[[1]]@metaData$file[1])
	    msi_plots [[p]] <- recordPlot()
		}
	}
	######################### PIXEL-BY-PIXEL CLASSIFICATION
	if ("pixels" %in% classification_of) {
    # Put the database and the spectra together
    global_spectra <- append(spectra_database, patient_spectra)
	global_peaks <- peak_picking(global_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
        global_peaks <- align_and_filter_peaks(global_peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
    # Compute the intensity matrix
    intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
    intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks=global_peaks, class_list=class_list, spectra_format="imzml", sample_output=TRUE, class_output=FALSE)
    rownames(intensity_matrix) <- intensity_matrix[,"Sample"]
    # Compute the agglomerative hierarchical clustering
	if (clustering_method == "agglomerative") {
		distance_matrix <- dist(intensity_matrix[,1:(ncol(intensity_matrix)-1)])
	    hierarchical_clustering <- hclust(distance_matrix)
	    plot(hierarchical_clustering)
	    hca_graph <- recordPlot()
	    distance_matrix <- as.matrix(distance_matrix)
	    # Remove the first rows (database spectra) and keep only the first columns (database spectra)
	    distance_matrix <- distance_matrix [(length(spectra_database)+1):nrow(distance_matrix),(1:length(spectra_database))]
	    ################### Classification output
		## Function for matrix apply
	    classification_hca_function <- function (dist_row) {
	        # Set the closest value to NULL
	        closest <- NULL
	        # Scroll the row...
	        for (l in 1:length(dist_row)) {
	            # If closest has not been set yet or the value in the cell is the closest...
	            if (is.null(closest) || dist_row[l] < closest){
	                # Set the closest value to be the one in the cell
	                closest <- dist_row[l]
	            }
	        }
	        # Find the position of the row where the closest is
	        closest_position <- which(dist_row == closest)
	        # The belonging class is corresponding to the element of the database which the spectrum is closer to
	        belonging_class <- names(dist_row)[closest_position]
	        return (belonging_class)
	    }
	    classification_hca  <- cbind(apply(distance_matrix, MARGIN=1, FUN=function(x) classification_hca_function(x)))
	    ### Generate the final output results
	    if (is.null(classification_hca_results_pixel_by_pixel)) {
	        classification_hca_results_pixel_by_pixel <- classification_hca
	    } else {classification_hca_results_pixel_by_pixel <- rbind(classification_hca_results_pixel_by_pixel, classification_hca)}
		colnames(classification_hca_results_pixel_by_pixel) <- "Pixel-by-pixel Classification"
	    ####################################### MS IMAGES
	    spectra_for_plotting <- patient_spectra
	    # Correspondence between the classificatio_hca matrix classification and the spectra under the hca nodes
	    for (s in 1:length(spectra_for_plotting)) {
	        for (z in 1:length(class_list)) {
	            if (classification_hca[s,1] == class_list[z]) {
	                spectra_for_plotting[[s]]@intensity <- rep(z, length(spectra_for_plotting[[s]]@intensity))
	            }
	        }
	    }
	    plotMsiSlice(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=1, legend=FALSE)
	    legend(x="bottomright", legend=class_list, fill=c("green","red"))
	    legend(x="topright", legend=spectra_for_plotting[[1]]@metaData$file[1])
	    msi_plots_pixel_by_pixel [[p]] <- recordPlot()
		}
	}
	############################# AVERAGE SPECTRUM
	if ("average" %in% classification_of) {
		# If the spectra per patient are more, average them for patient classification...
		if (isMassSpectrumList(patient_spectra)) {
			patient_spectra <- averageMassSpectra(patient_spectra, method="mean")
		}
		# Put the database and the spectra together
	    global_spectra <- append(spectra_database, patient_spectra)
	global_peaks <- peak_picking(global_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=SNR)
	global_peaks <- align_and_filter_peaks(global_peaks, tof_mode=tof_mode, peaks_filtering=TRUE, frequency_threshold_percent=25)
	    # Compute the intensity matrix
	    intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
	    intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks=global_peaks, class_list=class_list, spectra_format="imzml", sample_output=TRUE, class_output=FALSE)
	    rownames(intensity_matrix) <- intensity_matrix[,"Sample"]
		# Compute the agglomerative hierarchical clustering
		if (clustering_method == "agglomerative") {
		    distance_matrix <- dist(intensity_matrix[,1:(ncol(intensity_matrix)-1)])
		    hierarchical_clustering <- hclust(distance_matrix)
		    plot(hierarchical_clustering)
		    hca_graph <- recordPlot()
		    distance_matrix <- as.matrix(distance_matrix)
		    # Remove the first rows (database spectra) and keep only the first columns (database spectra)
		    distance_matrix <- distance_matrix [(length(spectra_database)+1):nrow(distance_matrix),(1:length(spectra_database))]
		    ################### Classification output
		    # Set the closest value to NULL
		    closest <- NULL
		    # Scroll the row...
		    for (l in 1:length(distance_matrix)) {
		        # If closest has not been set yet or the value in the cell is the closest...
		        if (is.null(closest) || distance_matrix[l] < closest){
		            # Set the closest value to be the one in the cell
		            closest <- distance_matrix[l]
		        }
		    }
		    # Find the position of the row where the closest is
		    closest_position <- which(distance_matrix == closest)
		    # The belonging class is corresponding to the element of the database which the spectrum is closer to
		    classification_hca <- as.matrix(names(distance_matrix)[closest_position])
			rownames(classification_hca) <- patient_vector[p]
		}
	    ### Generate the final output results
	    if (is.null(classification_hca_results_avg)) {
	        classification_hca_results_avg <- classification_hca
	    } else {classification_hca_results_avg <- rbind(classification_hca_results_avg, classification_hca)}
		colnames(classification_hca_results_avg) <- "Patient classification (average spectrum)"
	}
}
############################## PERFORMANCES
if ("average" %in% classification_of && ((is.null(true_class_list) || length(true_class_list) == 0) && true_class_in_file_name == TRUE)) {
	# Add the true class to the patient classification matrix
	true_class_list <- character()
	spectra_to_be_classified <- replace_class_name(spectra_to_be_classified, class_list=class_list)
	for (sp in spectra_to_be_classified) {
		true_class_list <- append(true_class_list, sp@metaData$file[1])
	}
	classification_hca_results_avg <- cbind(classification_hca_results_avg, true_class_list)
	hca_classification_performances <- confusionMatrix(classification_hca_results_avg[,1], classification_hca_results_avg[,2], positive="HP")
	#### ROC analysis
	#hca_roc <- roc(response=as.numeric(classification_hca_results_avg[,2]), predictor=as.numeric(classification_hca_results_avg[,1]))
}
# Output returning
return (list(classification_hca_results_avg=classification_hca_results_avg, classification_hca_results=classification_hca_results, classification_hca_results_pixel_by_pixel=classification_hca_results_pixel_by_pixel, hca_graph=hca_graph, ms_images=msi_plots, ms_images_pixel_by_pixel=msi_plots_pixel_by_pixel, hca_classification_performances=hca_classification_performances, hca_roc=hca_roc))
}
}





################################################################################





################################################ PROFILE CLASSIFICATION (ENSEMBLE, MULTICORE)
# The function takes a folder in which there are imzML files (one for each patient) or an imzML file or a list of MALDIquant spectra files, the R workspace containing the models with the name of the model objects in the workspace, and allows the user to specify something regarding the preprocessing of the spectra to be classified.
# The features in the model must be aligned to the features in the dataset.
# The function outputs a list containing: a matrix with the classification (patient's average spectrum), the model list and the average spectrum of the patients with red bars on the signals used by the models to classify it, a matrix with the ensemble classification (patient's average spectrum).
spectral_classification_profile <- function(spectra_path, filepath_R, list_of_models = c("svm", "pls", "nbc", "knn"), model_names = list(svm="RSVM_model",pls="pls_model",nbc="nbc_model",knn="knn_model"), spectra_preprocessing=TRUE, preprocessing_parameters=list(crop_spectra=TRUE, mass_range=c(4000,15000), data_transformation=FALSE, transformation_algorithm="sqrt", smoothing_algorithm="SavitzkyGolay", smoothing_strength="medium", baseline_subtraction_algorithm="SNIP", baseline_subtraction_iterations=100, normalisation_algorithm="TIC", normalisation_mass_range=NULL), spectral_alignment=FALSE, tof_mode="linear", peak_picking_algorithm="SuperSmoother", preprocess_spectra_in_packages_of=length(sample_spectra), multicore_processing=TRUE, decision_method_ensemble="majority", vote_weights_ensemble="equal") {
########## Load the required packages
install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign","stats", "parallel", "kernlab", "MASS", "klaR", "pls", "randomForest","nnet"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
#### TOF-MODE
if (tof_mode == "linear") {
	tolerance_ppm <- 2000
} else if (tof_mode == "reflectron" || tof_mode == "reflector") {
	tolerance_ppm <- 200
}
#### Mass range
mass_range <- preprocessing_parameters$mass_range
########### List the imzML files in the selected folder (if the path provided is a folder): check if its is folder, imzML file or spectra list
## Multiple imzML file path provided (folder)
if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed=TRUE)) == 0) {
	filepath_test_imzml <- read_spectra_files(spectra_path, spectra_format="imzml", full_path=TRUE)
}
## imzML file path provided (path)
if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed=TRUE)) != 0) {
	filepath_test_imzml <- spectra_path
}
## List of spectra
if (is.list(spectra_path)) {
	# Set the value of the filepath to a string (for the future IF cycles)
	filepath_test_imzml <- "List of spectra"
	sample_spectra <- spectra_path
}
######################################## Global OUTPUT Initialisation
final_result_matrix_all <- NULL
model_list <- list()
classification_ensemble_matrix_all <- NULL
average_spectra_with_bars_list <- list()
######################################## SPECTRA
########## Process the sample to be classified
##### For each imzML file...
for (p in 1:length(filepath_test_imzml)) {
	## Output initialisation (patient)
	final_result_matrix <- NULL
	average_spectra_with_bars <- list()
	## Import the spectra
	if(!isMassSpectrumList(spectra_path)) {
		# Import the spectra (one imzML at a time)
		if (!is.null(mass_range)) {
			sample_spectra <- importImzMl(filepath_test_imzml[p], massRange=mass_range)
		} else {
			sample_spectra <- importImzMl(filepath_test_imzml[p])
		}
	}
	## Replace the sample name (path) with the actual sample name
	sample_spectra <- replace_sample_name(sample_spectra)
	## Preprocess spectra
	if (spectra_preprocessing == TRUE) {
		sample_spectra <- preprocess_spectra(sample_spectra, tof_mode=tof_mode, preprocessing_parameters=preprocessing_parameters, process_in_packages_of=preprocess_spectra_in_packages_of, align_spectra=spectral_alignment, spectra_alignment_method="cubic", multicore_processing=multicore_processing)
	}
	## Generate the average spectrum
	sample_spectra <- averageMassSpectra(sample_spectra, method="mean")
	##### Sample name
	sample_name <- sample_spectra@metaData$file[[1]]
	## Preprocess average spectrum
	if (spectra_preprocessing == TRUE) {
		sample_spectra <- preprocess_spectra(sample_spectra, tof_mode=tof_mode, preprocessing_parameters=preprocessing_parameters, process_in_packages_of=preprocess_spectra_in_packages_of, align_spectra=spectral_alignment, spectra_alignment_method="cubic", multicore_processing=multicore_processing)
	}
	######################################## SUPPORT VECTOR MACHINE
	##### Detect if a SVM model is in the model list
	if ("svm" %in% list_of_models) {
		svm_is_present <- TRUE
	} else {
		svm_is_present <- FALSE
	}
	##### SVM
	if (svm_is_present == TRUE) {
		# Identify the SVM models
		svm_models <- model_names$svm
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each SVM model...
		for (md in 1:length(svm_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the svm model from the workspace
			svm_model <- get(svm_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, svm_model)
			# Class list
			class_list <- levels(factor(svm_model@fitted))
			##### Isolate the peaks used to create the model
			features_model <- colnames(svm_model@xmatrix[[1]])
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(svm_model, newdata = final_sample_matrix)
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class SVM"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
							break
						}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
		##### Return an error message in case of incompatibility
		#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## PARTIAL LEAST SQUARES
	##### Detect if a PLS model is in the model list
	if ("pls" %in% list_of_models) {
		pls_is_present <- TRUE
	} else {
		pls_is_present <- FALSE
	}
	##### pls
	if (pls_is_present == TRUE) {
		# Identify the PLS models
		pls_models <- model_names$pls
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each pls model...
		for (md in 1:length(pls_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the pls model from the workspace
			pls_model <- get(pls_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, pls_model)
			# Class list
			class_list <- levels(factor(colnames(pls_model$finalModel$model$y)))
			##### Isolate the peaks used to create the model
			features_model <- pls_model$finalModel$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(pls_model, newdata = final_sample_matrix)
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class PLS"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
			}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## NAIVE BAYES CLASSIFIER
	##### Detect if a NBC model is in the model list
	if ("nbc" %in% list_of_models) {
		nbc_is_present <- TRUE
	} else {
		nbc_is_present <- FALSE
	}
	##### NBC
	if (nbc_is_present == TRUE) {
		# Identify the NBC models
		nbc_models <- model_names$nbc
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each nbc model...
		for (md in 1:length(nbc_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the pls model from the workspace
			nbc_model <- get(nbc_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, nbc_model)
			# Class list
			class_list <- levels(factor(nbc_model$finalModel$levels))
			##### Isolate the peaks used to create the model
			features_model <- nbc_model$finalModel$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(nbc_model, newdata = final_sample_matrix)
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class NBC"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
					}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## K-NEAREST NEIGHBOUR (KNN)
	##### Detect if a KNN model is in the model list
	if ("knn" %in% list_of_models) {
		knn_is_present <- TRUE
	} else {
		knn_is_present <- FALSE
	}
	##### KNN
	if (knn_is_present == TRUE) {
		# Identify the KNN models
		knn_models <- model_names$knn
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each KNN model...
		for (md in 1:length(knn_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the knn model from the workspace
			knn_model <- get(knn_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, knn_model)
			# Class list
			class_list <- levels(factor(knn_model$obsLevels))
			##### Isolate the peaks used to create the model
			features_model <- knn_model$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(knn_model, newdata = final_sample_matrix, type="class")
			# Generate a matrix with the results
			result_matrix <- matrix(nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix[,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class KNN"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
					}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## LINEAR DISCRIMINANT ANALYSIS
	##### Detect if a LDA model is in the model list
	if ("lda" %in% list_of_models) {
		lda_is_present <- TRUE
	} else {
		lda_is_present <- FALSE
	}
	##### LDA
	if (lda_is_present == TRUE) {
		# Identify the LDA models
		lda_models <- model_names$lda
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each LDA model...
		for (md in 1:length(lda_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the svm model from the workspace
			lda_model <- get(lda_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, lda_model)
			# Class list
			class_list <- lda_model$lev
			##### Isolate the peaks used to create the model
			features_model <- lda_model$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(lda_model, newdata = final_sample_matrix, type="class")$class
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class LDA"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
					}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message incase of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## RANDOM FOREST
	##### Detect if a RF model is in the model list
	if ("rf" %in% list_of_models) {
		rf_is_present <- TRUE
	} else {
		rf_is_present <- FALSE
	}
	##### RF
	if (rf_is_present == TRUE) {
		# Identify the RF models
		rf_models <- model_names$rf
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each RF model...
		for (md in 1:length(rf_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the svm model from the workspace
			rf_model <- get(rf_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, rf_model)
			# Class list
			class_list <- rf_model$classes
			##### Isolate the peaks used to create the model
			features_model <- rf_model$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(rf_model, newdata = final_sample_matrix, type="class")
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class RF"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
					}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## NEURAL NETWORK
	##### Detect if a NN model is in the model list
	if ("nn" %in% list_of_models) {
		nn_is_present <- TRUE
	} else {
		nn_is_present <- FALSE
	}
	##### RF
	if (nn_is_present == TRUE) {
		# Identify the NN models
		nn_models <- model_names$nn
		########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
		##### For each NN model...
		for (md in 1:length(nn_models)) {
			##### LOAD THE R WORKSPACE
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_R, envir=temporary_environment)
			# Get the svm model from the workspace
			nn_model <- get(nn_models[md], pos=temporary_environment)
			# Add it to the output list
			model_list <- append(model_list, nn_model)
			# Class list
			class_list <- levels(factor(nn_model$obsLevels))
			##### Isolate the peaks used to create the model
			features_model <- nn_model$xNames
			# Remove the X
			for (f in 1:length(features_model)) {
				name_splitted <- unlist(strsplit(features_model[f],""))
				feature_def <- name_splitted [2]
				for (i in 3:length(name_splitted)) {
					feature_def <- paste(feature_def, name_splitted[i], sep="")
				}
				features_model[f] <- feature_def
			}
			##### Check if the model features and the dataset features are compatible
			# IMPLEMENT THIS AFTER USING THE FUNCTION GENERATE_CUSTOM_INTENSITY_MATRIX
			### Generate the intensity matrix with the features from the model
			final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
			# Put the X at the beginning of the peak names
			for (n in 1:length(colnames(final_sample_matrix))) {
				name <- paste("X", colnames(final_sample_matrix)[n], sep="")
				colnames(final_sample_matrix)[n] <- name
			}
			# Predictions
			predicted_classes <- predict(nn_model, newdata = final_sample_matrix, type="class")
			# Generate a matrix with the results
			result_matrix <- matrix (nrow=1, ncol=1)
			rownames(result_matrix) <- sample_name
			result_matrix [,1] <- as.character(predicted_classes)
			colnames(result_matrix) <- "Predicted Class NN"
			#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
			if (is.null(final_result_matrix)) {
				final_result_matrix <- result_matrix
			} else {
				final_result_matrix <- cbind(final_result_matrix, result_matrix)
			}
			################# Average spectrum with bars onto the signals used by the model
			# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
			# Detect peaks in the avg (SNR=1)
			sample_peaks_avg_for_bars <- peak_picking(sample_spectra, peak_picking_algorithm=peak_picking_algorithm, tof_mode=tof_mode, SNR=3)
			# Determine the coordinates of the bars
			coordinates_of_bars <- list(x=numeric(), y=numeric())
			# Check if the features used for the model are in the spectrum
			for (f in features_model) {
				presence_in_the_avg <- FALSE
				# Scroll the peaks
				for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
					# If there is a match...
					if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
						# Add the intensity of this peak to the y coordinates of the bars
						coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
						coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
						# Set the presence in the peaklist to true
						presence_in_the_avg <- TRUE
						# Break the for cycle to avoid duplicates and to continue
						break
					}
				}
				if (presence_in_the_avg == FALSE) {
					# If the feature is not in the peaklist, scroll the datapoints in the spectrum
					for (j in 1:length(sample_spectra@mass)) {
						# If there is a match...
						if ((abs(sample_spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= tolerance_ppm) {
							# Add the intensity of this peak to the y coordinates of the bars
							coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra@mass[j])
							coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra@intensity[j])
							# Break the for cycle to avoid duplicates and to continue
							break
						}
					}
				}
			}
			plot(sample_spectra, xlab="m/z", ylab="Intensity (a.i.)")
			legend(x="topright", legend=sample_name)
			# Draw the bars
			for (s in 1:length(coordinates_of_bars$x)) {
				# Vertical bars (x,y x,y)
				segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
				# Horizontal segments(x,y , x,y)
				segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
				segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
			}
			average_spectra_with_bars[[md]] <- recordPlot()
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
		}
		##### Append the list of average spectra with bars to the final global list
		average_spectra_with_bars_list <- append(average_spectra_with_bars_list, average_spectra_with_bars)
	}
	######################################## Add the result matrix for the patient to the global matrix
	if (is.null(final_result_matrix_all)) {
		final_result_matrix_all <- final_result_matrix
	} else {
		final_result_matrix_all <- rbind(final_result_matrix_all, final_result_matrix)
	}
	######################################## ENSEMBLE VOTE
	if (length(model_names) > 1 && !is.null(final_result_matrix)) {
		########## Ensemble results
		classification_ensemble_matrix <- final_result_matrix
		########## Vote
		##### Majority vote
		if (decision_method_ensemble == "majority" && vote_weights_ensemble == "equal") {
			# Function for matrix apply (x=row)
			majority_vote_function <- function (x, class_list) {
				# Sort the class list for reproducibility
				class_list <- sort(class_list)
				# Generate the vote vector (same length as the class list, with the number of the votes for each class, labeled)
				votes <- integer(length=length(class_list))
				names(votes) <- class_list
				# Record the vote numbers
				for (class in class_list) {
					votes [which(class_list==class)] <- length(which(x == class))
				}
				# Establish the final vote
				final_vote <- names(votes)[which(votes==max(votes))]
				# Even vote
				if (length(final_vote) != 1) {
					final_vote <- class_list[1]
				}
				# Return
				return(final_vote)
			}
			# For each spectrum (matrix row), establish the final majority vote
			classification_ensemble_matrix <- cbind(apply(X=final_result_matrix, MARGIN=1, FUN=function(x) majority_vote_function(x, class_list)))
			colnames(classification_ensemble_matrix) <- "Ensemble classification"
			# Store the ensemble classification matrix in the final output list
			if (is.null(classification_ensemble_matrix_all)) {
				classification_ensemble_matrix_all <- classification_ensemble_matrix
			} else {
				classification_ensemble_matrix_all <- rbind(classification_ensemble_matrix_all, classification_ensemble_matrix)
			}
		}
	}
}
return (list(final_result_matrix=final_result_matrix_all, average_spectra_with_bars_list=average_spectra_with_bars_list, model_list=model_list, classification_ensemble_matrix=classification_ensemble_matrix_all))
}





################################################################################





######################################### CLASSIFICATION: PIXEL-BY-PIXEL (MULTICORE, ENSEMBLE)
# The function takes a folder in which there are imzML files (one for each patient) or an imzML file or a list of MALDIquant spectra files, the R workspace containing the models with the name of the model objects in the workspace, and allows the user to specify something regarding the preprocessing of the spectra to be classified.
# The function outputs a list containing: a matrix with the classification (pixel-by-pixel), MS images with the pixel-by-pixel classification, the model list, a matrix with the ensemble classification (pixel-by-pixel) and MS images with the pixel-by-pixel ensemble classification.
# Parallel computation implemented
spectral_classification_pixelbypixel <- function (spectra_path, filepath_R, list_of_models = c("svm", "pls", "nbc", "knn"), model_names = list(svm="RSVM_model",pls="pls_model",nbc="nbc_model",knn="knn_model"), peak_picking_algorithm="SuperSmoother", preprocessing_parameters=list(crop_spectra=TRUE, mass_range=c(4000,15000), data_transformation=FALSE, transformation_algorithm="sqrt", smoothing_algorithm="SavitzkyGolay", smoothing_strength="medium", baseline_subtraction_algorithm="SNIP", baseline_subtraction_iterations=100, normalisation_algorithm="TIC", normalisation_mass_range=NULL), spectral_alignment=FALSE, tof_mode="linear", spectra_preprocessing=TRUE, preprocess_spectra_in_packages_of=length(sample_spectra), multicore_processing=TRUE, decision_method_ensemble="majority", vote_weights_ensemble="equal", pixel_grouping=c("single","moving window average","graph","hca"), moving_window_size = 10, number_of_hca_nodes = 10) {
	# Install and load the required packages
	install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign","stats", "parallel", "kernlab", "MASS", "klaR", "pls", "randomForest", "lda"))
	# Rename the trim function
	trim_spectra <- get(x="trim", pos="package:MALDIquant")
	##### TOF-MODE
	if (tof_mode == "linear") {
		tolerance_ppm <- 2000
	} else if (tof_mode == "reflector" || tof_mode == "reflectron") {
		tolerance_ppm <- 200
	}
	##### Mass range
	mass_range <- preprocessing_parameters$mass_range
	########## List the imzML files in the selected folder (if the path provided is a folder): check if its is folder, imzML file or spectra list
	## Multiple imzML filepath (path)
	if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed=TRUE)) == 0) {
		filepath_test_imzml <- read_spectra_files(spectra_path, spectra_format="imzml", full_path=TRUE)
	}
	## Single imzML filepath
	if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed=TRUE)) != 0) {
		filepath_test_imzml <- spectra_path
	}
	## Spectra list
	if (is.list(spectra_path)) {
		# Assign a string to the filepath for the future IF cycles
		filepath_test_imzml <- "List of spectra"
		sample_spectra <- spectra_path
	}
	######################################## Global OUTPUT Initialisation
	final_result_matrix_list <- list()
	classification_msi_list <- list()
	model_list <- list()
	classification_ensemble_matrix_list <- list()
	classification_ensemble_msi_list <- list()
	######################################## SPECTRA
	########## Process the sample to be classified
	##### For each imzML file...
	for (p in 1:length(filepath_test_imzml)) {
		## Output initialisation (patient)
		final_result_matrix <- NULL
		classification_msi_patient <- list()
		## Import the spectra
		if(!isMassSpectrumList(spectra_path)) {
			# Import the spectra (one imzML at a time)
			if (!is.null(mass_range)) {
				sample_spectra <- importImzMl(filepath_test_imzml[p], massRange=mass_range)
			} else {
				sample_spectra <- importImzMl(filepath_test_imzml[p])
			}
		}
		## Replace the sample name (path) with the actual sample name
		sample_spectra <- replace_sample_name(sample_spectra)
		## Preprocess spectra
		if (spectra_preprocessing == TRUE) {
			sample_spectra <- preprocess_spectra(sample_spectra, tof_mode=tof_mode, preprocessing_parameters=preprocessing_parameters, process_in_packages_of=preprocess_spectra_in_packages_of, align_spectra=spectral_alignment, spectra_alignment_method="cubic", multicore_processing=multicore_processing)
		}
		######################################## SUPPORT VECTOR MACHINE
		##### Detect if a SVM model is in the model list
		if ("svm" %in% list_of_models) {
			svm_is_present <- TRUE
		} else {
			svm_is_present <- FALSE
		}
		##### SVM
		if (svm_is_present == TRUE) {
			# Identify the SVM models
			svm_models <- model_names$svm
			########## Outputs
			classification_msi_svm <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each SVM model...
			for (md in 1:length(svm_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the svm model from the workspace
				svm_model <- get(svm_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, svm_model)
				# Class list
				class_list <- levels(factor(svm_model@fitted))
				##### Isolate the peaks used to create the model
				features_model <- colnames(svm_model@xmatrix[[1]])
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						colnames(final_sample_matrix)[n] <- paste("X", colnames(final_sample_matrix)[n], sep = "")
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(svm_model, newdata = final_sample_matrix)
					# Generate a matrix with the results
					result_matrix_svm <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_svm) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_svm [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_svm) <- "Predicted Class SVM"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_svm
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_svm)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_svm[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
				########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_svm <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(svm_model, newdata = sample_bin_matrix)
						# Generate a matrix with the results
						result_matrix_svm_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_svm_bin) <- sample_name
						result_matrix_svm_bin [1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_svm_bin) <- "Predicted Class SVM"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_svm)) {
							result_matrix_svm <- result_matrix_svm_bin
						} else {
							result_matrix_svm <- rbind(result_matrix_svm, result_matrix_svm_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_svm
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_svm)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_svm[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_svm[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
				########## HCA
					### Initialise the model result matrix
					result_matrix_svm <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_svm) <- "Predicted Class SVM"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_svm) <- rep(sample_name, nrow(result_matrix_svm))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(svm_model, newdata = sample_hca_matrix)
						# Fill the model matrix with the results
						result_matrix_svm[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_svm[[md]] <- recordPlot()
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_svm
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_svm)
					}
				}
				##### Return an error message
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (SVM)
			classification_msi_patient <- append (classification_msi_patient, classification_msi_svm)
		}
		######################################## PARTIAL LEAST SQUARES
		##### Detect if a PLS model is in the model list
		if ("pls" %in% list_of_models) {
			pls_is_present <- TRUE
		} else {
			pls_is_present <- FALSE
		}
		##### PLS
		if (pls_is_present == TRUE) {
			# Identify the PLS models
			pls_models <- model_names$pls
			########## Outputs
			classification_msi_pls <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			#### For each PLS model...
			for (md in 1:length(pls_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the pls model from the workspace
				pls_model <- get(pls_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, pls_model)
				# Class list
				class_list <- levels(factor(colnames(pls_model$finalModel$model$y)))
				##### Isolate the peaks used to create the model
				features_model <- pls_model$finalModel$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					### Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(pls_model, newdata = final_sample_matrix)
					# Generate a matrix with the results
					result_matrix_pls <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_pls) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_pls [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_pls) <- "Predicted Class PLS"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_pls
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_pls)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images
					classification_msi_pls[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_pls <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(pls_model, newdata = sample_bin_matrix)
						# Generate a matrix with the results
						result_matrix_pls_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_pls_bin) <- sample_name
						result_matrix_pls_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_pls_bin) <- "Predicted Class PLS"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_pls)) {
							result_matrix_pls <- result_matrix_pls_bin
						} else {
							result_matrix_pls <- rbind(result_matrix_pls, result_matrix_pls_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_pls
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_pls)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_pls[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_pls[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_pls <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_pls) <- "Predicted Class PLS"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_pls) <- rep(sample_name, nrow(result_matrix_pls))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(pls_model, newdata = sample_hca_matrix)
						# Fill the model matrix with the results
						result_matrix_pls[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_pls[[md]] <- recordPlot()
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_pls
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_pls)
					}
				}
			##### Return an error message in case of incompatibility
			#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (PLS)
			classification_msi_patient <- append(classification_msi_patient, classification_msi_pls)
		}
		######################################## NAIVE BAYES CLASSIFIER
		##### Detect if a NBC model is in the model list
		if ("nbc" %in% list_of_models) {
			nbc_is_present <- TRUE
		} else {
			nbc_is_present <- FALSE
		}
		##### NBC
		if (nbc_is_present == TRUE) {
			# Identify the NBC models
			nbc_models <- model_names$nbc
			########## Outputs
			classification_msi_nbc <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each NBC model...
			for (md in 1:length(nbc_models)) {
				### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the nbc model from the workspace
				nbc_model <- get(nbc_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, nbc_model)
				# Class list
				class_list <- levels(factor(nbc_model$finalModel$levels))
				##### Isolate the peaks used to create the model
				features_model <- nbc_model$finalModel$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(nbc_model, newdata = final_sample_matrix)
					# Generate a matrix with the results
					result_matrix_nbc <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_nbc) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_nbc [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_nbc) <- "Predicted Class Bayes"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nbc
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nbc)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green", "red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images
					classification_msi_nbc[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_nbc <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(nbc_model, newdata = sample_bin_matrix)
						# Generate a matrix with the results
						result_matrix_nbc_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_nbc_bin) <- sample_name
						result_matrix_nbc_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_nbc_bin) <- "Predicted Class NBC"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_nbc)) {
							result_matrix_nbc <- result_matrix_nbc_bin
						} else {
							result_matrix_nbc <- rbind(result_matrix_nbc, result_matrix_nbc_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nbc
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nbc)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_nbc[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_nbc[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_nbc <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_nbc) <- "Predicted Class NBC"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_nbc) <- rep(sample_name, nrow(result_matrix_nbc))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(nbc_model, newdata = sample_hca_matrix)
						# Fill the model matrix with the results
						result_matrix_nbc[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_nbc[[md]] <- recordPlot()
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nbc
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nbc)
					}
				}
				##### Return an error message in case of incompatibility
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (NBC)
			classification_msi_patient <- append(classification_msi_patient, classification_msi_nbc)
		}
		######################################## K-NEAREST NEIGHBOUR (KNN)
		##### Detect if a KNN model is in the model list
		if ("knn" %in% list_of_models) {
			knn_is_present <- TRUE
		} else {
			knn_is_present <- FALSE
		}
		##### KNN
		if (knn_is_present == TRUE) {
			# Identify the KNN models
			knn_models <- model_names$knn
			########## Outputs
			classification_msi_knn <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each KNN model...
			for (md in 1:length(knn_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the knn model from the workspace
				knn_model <- get(knn_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, knn_model)
				# Class list
				class_list <- levels(factor(knn_model$obsLevels))
				##### Isolate the peaks used to create the model
				features_model <- knn_model$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(knn_model, newdata = final_sample_matrix, type = "class")
					# Generate a matrix with the results
					result_matrix_knn <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_knn) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_knn[,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_knn) <- "Predicted Class KNN"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_knn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_knn)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for KNN)
					classification_msi_knn[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_knn <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(knn_model, newdata = sample_bin_matrix, type = "class")
						# Generate a matrix with the results
						result_matrix_knn_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_knn_bin) <- sample_name
						result_matrix_knn_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_knn_bin) <- "Predicted Class KNN"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_knn)) {
							result_matrix_knn <- result_matrix_knn_bin
						} else {
							result_matrix_knn <- rbind(result_matrix_knn, result_matrix_knn_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_knn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_knn)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_knn[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_knn[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_knn <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_knn) <- "Predicted Class KNN"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_knn) <- rep(sample_name, nrow(result_matrix_knn))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(knn_model, newdata = sample_hca_matrix, type = "class")
						# Fill the model matrix with the results
						result_matrix_knn[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_knn[[md]] <- recordPlot()
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_knn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_knn)
					}
				}
				##### Return an error message in case of incompatibility
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (KNN)
			classification_msi_patient <- append (classification_msi_patient, classification_msi_knn)
		}
		######################################## LINEAR DISCRIMINANT ANALYSIS
		##### Detect if a LDA model is in the model list
		if ("lda" %in% list_of_models) {
			lda_is_present <- TRUE
		} else {
			lda_is_present <- FALSE
		}
		##### LDA
		if (lda_is_present == TRUE) {
			# Identify the LDA models
			lda_models <- model_names$lda
			########## Outputs
			classification_msi_lda <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each LDA model...
			for (md in 1:length(lda_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the svm model from the workspace
				lda_model <- get(lda_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, lda_model)
				# Class list
				class_list <- lda_model$lev
				##### Isolate the peaks used to create the model
				features_model <- lda_model$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(lda_model, newdata = final_sample_matrix, type="class")$class
					# Generate a matrix with the results
					result_matrix_lda <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_lda) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_lda [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_lda) <- "Predicted Class LDA"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_lda
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_lda)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for LDA)
					classification_msi_lda[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_lda <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(lda_model, newdata = sample_bin_matrix, type="class")$class
						# Generate a matrix with the results
						result_matrix_lda_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_lda_bin) <- sample_name
						result_matrix_lda_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_lda_bin) <- "Predicted Class LDA"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_lda)) {
							result_matrix_lda <- result_matrix_lda_bin
						} else {
							result_matrix_lda <- rbind(result_matrix_lda, result_matrix_lda_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_lda
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_lda)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_lda[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_lda[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_lda <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_lda) <- "Predicted Class LDA"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_lda) <- rep(sample_name, nrow(result_matrix_lda))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(lda_model, newdata = as.data.frame(sample_hca_matrix), type="class")$class
						# Fill the model matrix with the results
						result_matrix_lda[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_lda[[md]] <- recordPlot()
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_lda
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_lda)
					}
				}
				##### Return an error message in case of incompatibility
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (LDA)
			classification_msi_patient <- append (classification_msi_patient, classification_msi_lda)
		}
		######################################## RANDOM FOREST
		##### Detect if a RF model is in the model list
		if ("rf" %in% list_of_models) {
			rf_is_present <- TRUE
		} else {
			rf_is_present <- FALSE
		}
		##### RF
		if (rf_is_present == TRUE) {
			# Identify the RF models
			rf_models <- model_names$rf
			########## Outputs
			classification_msi_rf <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each RF model...
			for (md in 1:length(rf_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the svm model from the workspace
				rf_model <- get(rf_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, rf_model)
				# Class list
				class_list <- rf_model$classes
				##### Isolate the peaks used to create the model
				features_model <- rf_model$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(rf_model, newdata = final_sample_matrix, type="class")
					# Generate a matrix with the results
					result_matrix_rf <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_rf) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_rf [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_rf) <- "Predicted Class RF"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_rf
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_rf)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(predicted_classes)
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for LDA)
					classification_msi_rf[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_rf <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(rf_model, newdata = sample_bin_matrix, type = "class")
						# Generate a matrix with the results
						result_matrix_rf_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_rf_bin) <- sample_name
						result_matrix_rf_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_rf_bin) <- "Predicted Class RF"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_rf)) {
							result_matrix_rf <- result_matrix_rf_bin
						} else {
							result_matrix_rf <- rbind(result_matrix_rf, result_matrix_rf_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_rf
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_rf)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_rf[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_rf[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_rf <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_rf) <- "Predicted Class RF"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_rf) <- rep(sample_name, nrow(result_matrix_rf))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(rf_model, newdata = sample_hca_matrix, type="class")
						# Fill the model matrix with the results
						result_matrix_rf[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_rf[[md]] <- recordPlot()
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_rf
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_rf)
					}
				}
				##### Return an error message in case of incompatibility
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (RF)
			classification_msi_patient <- append (classification_msi_patient, classification_msi_rf)
		}
		######################################## NEURAL NETWORK
		##### Detect if a NN model is in the model list
		if ("nn" %in% list_of_models) {
			nn_is_present <- TRUE
		} else {
			nn_is_present <- FALSE
		}
		##### NN
		if (nn_is_present == TRUE) {
			# Identify the NN models
			nn_models <- model_names$nn
			########## Outputs
			classification_msi_nn <- list()
			########## The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
			##### For each NN model...
			for (md in 1:length(nn_models)) {
				##### LOAD THE R WORKSPACE
				# Create a temporary environment
				temporary_environment <- new.env()
				# Load the workspace
				load(filepath_R, envir=temporary_environment)
				# Get the svm model from the workspace
				nn_model <- get(nn_models[md], pos=temporary_environment)
				# Add it to the output list
				model_list <- append(model_list, nn_model)
				# Class list
				class_list <- levels(factor(nn_model$obsLevels))
				##### Isolate the peaks used to create the model
				features_model <- nn_model$xNames
				# Remove the X
				for (f in 1:length(features_model)) {
					name_splitted <- unlist(strsplit(features_model[f],""))
					feature_def <- name_splitted [2]
					for (i in 3:length(name_splitted)) {
						feature_def <- paste(feature_def, name_splitted[i], sep="")
					}
					features_model[f] <- feature_def
				}
				########## SINGLE PIXEL CLASSIFICATION
				if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(sample_spectra)) {
					### Generate the intensity matrix with the features from the model
					final_sample_matrix <- generate_custom_intensity_matrix(sample_spectra, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
					# Put the X at the beginning of the peak names
					for (n in 1:length(colnames(final_sample_matrix))) {
						name <- paste("X", colnames(final_sample_matrix)[n], sep="")
						colnames(final_sample_matrix)[n] <- name
					}
					##### Predictions (spectra by spectra)
					predicted_classes <- predict(nn_model, newdata = final_sample_matrix, type="class")
					# Generate a matrix with the results
					result_matrix_nn <- matrix (nrow=length(predicted_classes), ncol=1)
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_nn) <- cbind(rep(sample_name, length(sample_spectra)))
					result_matrix_nn [,1] <- cbind(as.character(predicted_classes))
					colnames(result_matrix_nn) <- "Predicted Class NN"
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nn)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					class_as_number <- as.numeric(factor(predicted_classes))
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for NN)
					classification_msi_nn[[md]] <- recordPlot()
				} else if (pixel_grouping == "moving window average") {
					########## MOVING WINDOW AVERAGE
					##### Rearrange the spectra accoring to the space coordinates
					sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
					##### Initialise the model result matrix
					result_matrix_nn <- NULL
					##### For each spectrum...
					for (s in 1:length(sample_spectra)) {
						### Define the indices
						index1 <- s - floor(moving_window_size/2)
						index2 <- s + floor(moving_window_size/2)
						### Check the indices
						if (index1 <= 0) {
							index1 <- 1
						} else if (index1 > length(sample_spectra)) {
							index1 <- length(sample_spectra)
						}
						if (index2 <= 0) {
							index2 <- 1
						} else if (index2 > length(sample_spectra)) {
							index2 <- length(sample_spectra)
						}
						### Isolate the spectra from the bin
						spectra_bin <- sample_spectra[index1:index2]
						### Generate the average spectrum for the bin
						average_spectrum_bin <- averageMassSpectra(spectra_bin)
						### Preprocessing the AVG spectrum
						average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (n in 1:length(colnames(sample_bin_matrix))) {
							colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(nn_model, newdata = sample_bin_matrix)
						# Generate a matrix with the results
						result_matrix_nn_bin <- matrix(nrow=1, ncol=1)
						sample_name <- average_spectrum_bin@metaData$file[1]
						rownames(result_matrix_nn_bin) <- sample_name
						result_matrix_nn_bin[1,1] <- as.character(predicted_class_avg)
						colnames(result_matrix_nn_bin) <- "Predicted Class NN"
						### Add the result matrix to the result matrix of the model
						if (is.null(result_matrix_nn)) {
							result_matrix_nn <- result_matrix_nn_bin
						} else {
							result_matrix_nn <- rbind(result_matrix_nn, result_matrix_nn_bin)
						}
					}
					### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nn)
					}
					########## Generate a molecular image of the classification
					# Replace the spectra intensities with the class number for plotting purposes
					predicted_classes <- result_matrix_nn[,1]
					spectra_for_plotting <- sample_spectra
					for (s in 1:length(spectra_for_plotting)) {
						spectra_for_plotting[[s]]@intensity <- rep(ifelse(predicted_classes[s] == class_list[1],1,2)/2, length(spectra_for_plotting[[s]]@intensity))
					}
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_nn[[md]] <- recordPlot()
				} else if (pixel_grouping == "hca") {
					########## HCA
					### Initialise the model result matrix
					result_matrix_nn <- matrix("class", nrow=length(sample_spectra), ncol = 1)
					colnames(result_matrix_nn) <- "Predicted Class RF"
					sample_name <- sample_spectra[[1]]@metaData$file[1]
					rownames(result_matrix_nn) <- rep(sample_name, nrow(result_matrix_nn))
					### Detect and align peaks
					peaks <- peak_picking(sample_spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
					peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
					### Generate the peaklist matrix
					peaklist <- intensityMatrix(peaks, sample_spectra)
					### Compute the distance matrix
					distance_matrix <- dist(peaklist, method = "euclidean")
					# Generate the dendrogram
					hca <- hclust(distance_matrix)
					### Cut the tree to generate K number of sub-clusters
					hca_groups <- cutree(hca, k = number_of_hca_nodes)
					## Generate the spectra for plotting list
					spectra_for_plotting <- list()
					## For each HCA node...
					for (n in 1:number_of_hca_nodes) {
						# Index the spectra under in the selected subgroup of the HCA
						index <- which(hca_groups == n)
						spectra_hca <- sample_spectra[index]
						# Generate the average spectrum for these spectra under the node
						average_spectrum_hca <- averageMassSpectra(spectra_hca)
						### Preprocessing the AVG spectrum
						average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=FALSE)
						### Peak picking on the AVG spectrum
						sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, spectra_preprocessing = FALSE, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 5, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 0.1, intensity_threshold_method = "element-wise", process_in_packages_of = preprocess_spectra_in_packages_of, multicore_processing = multicore_processing)
						### Put the X at the beginning of the peak names
						for (x in 1:length(colnames(sample_hca_matrix))) {
							colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
						}
						### Predictions (AVG spectrum)
						predicted_class_avg <- predict(nn_model, newdata = sample_hca_matrix)
						# Fix the character conversion (the class is a number, conver it back to the character) (the class number is one less than the converted number, numbers start from 0, the conversion from 1)
						for (cls in 1:length(class_list)) {
							try({
								if ((as.numeric(predicted_class_avg) + 1) == cls) {
									predicted_class_avg <- as.character(class_list[cls])
								}}, silent = TRUE)
						}
						# Fill the model matrix with the results
						result_matrix_nn[index,1] <- as.character(predicted_class_avg)
						### Fix the spectra for plotting list (edit the intensity values)
						spectra_for_plotting_hca <- spectra_hca
						for (s in 1:length(spectra_for_plotting_hca)) {
							spectra_for_plotting_hca[[s]]@intensity <- rep(ifelse(predicted_class_avg == class_list[1],1,2)/2, length(spectra_for_plotting_hca[[s]]@intensity))
						}
						# Append to the final list of spectra for plotting
						spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_hca)
					}
					# MSI
					slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
					plotMsiSlice(slices, legend=FALSE, scale=F)
					legend(x="bottomright", legend=class_list, fill=c("green","red"))
					legend(x="topright", legend=sample_name)
					# Store the plot into the list of images (for SVM)
					classification_msi_nn[[md]] <- recordPlot()
					##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
					if (is.null(final_result_matrix)) {
						final_result_matrix <- result_matrix_nn
					} else {
						final_result_matrix <- cbind(final_result_matrix, result_matrix_nn)
					}
				}
				##### Return an error message in case of incompatibility
				#tkmessageBox (title = "Classification not possible", message = "The classification cannot be performed due to incompatibility between the model's features and the dataset's features: the sample's mass range must be wider than the model's", icon = "warning")
			}
			# Append the classification image list to the final list (NN)
			classification_msi_patient <- append (classification_msi_patient, classification_msi_nn)
		}
		######################################## Store the matrix with the classification from all the models in the final list of matrices
		final_result_matrix_list[[p]] <- final_result_matrix
		######################################## Store all the MS images in the element of the final list of MS images
		classification_msi_list[[p]] <- classification_msi_patient
		######################################## ENSEMBLE VOTE
		if (length(model_names) > 1 && !is.null(final_result_matrix)) {
			########## Ensemble results
			classification_ensemble_matrix <- final_result_matrix
			########## Vote
			##### Majority vote
			if (decision_method_ensemble == "majority" && vote_weights_ensemble == "equal") {
				# Function for matrix apply (x=row)
				majority_vote_function <- function (x, class_list) {
					# Sort the class list for reproducibility
					class_list <- sort(class_list)
					# Generate the vote vector (same length as the class list, with the number of the votes for each class, labeled)
					votes <- integer(length=length(class_list))
					names(votes) <- class_list
					# Count the votes for each class
					for (class in class_list) {
						votes [which(class_list==class)] <- length(which(x == class))
					}
					# Determine the final majority vote
					final_vote <- names(votes)[which(votes==max(votes))]
					# Even vote
					if (length(final_vote) != 1) {
						final_vote <- class_list[1]
					}
					# Return the vote
					return(final_vote)
				}
				# For each spectrum (matrix row), establish the final majority vote
				classification_ensemble_matrix <- cbind(apply(X=final_result_matrix, MARGIN=1, FUN=function(x) majority_vote_function(x, class_list)))
				colnames(classification_ensemble_matrix) <- "Ensemble classification"
				# Store the ensemble classification matrix in the final output list
				classification_ensemble_matrix_list[[p]] <- classification_ensemble_matrix
				########## Generate a molecular image of the classification
				# Generate the "predicted classes" vector from the ensemble classification matrix
				predicted_classes <- as.factor(classification_ensemble_matrix)
				# Replace the spectra intensities with the class number for plotting purposes
				class_as_number <- as.numeric(predicted_classes)
				spectra_for_plotting <- sample_spectra
				for (s in 1:length(spectra_for_plotting)) {
					spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s]/2, length(spectra_for_plotting[[s]]@intensity))
				}
				slices <- msiSlices(spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance=1, adjust=TRUE, method="median")
				plotMsiSlice(slices, legend=FALSE, scale=F)
				legend(x="bottomright", legend=class_list, fill=c("green", "red"))
				legend(x="topright", legend=sample_name)
				# Store the plot into the list of images
				classification_ensemble_msi_list[[p]] <- recordPlot()
			}
		}
	}
	# Return the results
	return (list(final_result_matrix_list=final_result_matrix_list, classification_msi_list=classification_msi_list, model_list=model_list, classification_ensemble_matrix_list=classification_ensemble_matrix_list, classification_ensemble_msi_list=classification_ensemble_msi_list))
}











################################################################################











































##################################################################### STATISTICS

#################################################### PEARSON CORRELATION p-value
# This function returns a p-value for the significance of a correlation coefficient.
correlation_pvalue <- function (correlation_coefficient, number_of_samples, correlation_type="pearson", tails=2) {
    degrees_of_freedom <- number_of_samples - 2
    if (correlation_type == "pearson") {
        t_value = correlation_coefficient * sqrt((number_of_samples-2)/(1-(correlation_coefficient)^2))
        p_value <- pt(t_value, df=degrees_of_freedom, lower.tail = FALSE, log.p = FALSE)
    }
    if (tails == 1) {
        return (p_value)
    }
    if (tails == 2) {
        return (p_value*2)
    }
}





################################################################################





##################################### DEFINED SPECTRA GROUPING (PEAKLIST MATRIX)
# This function takes a peaklist matrix as input and generates a peaklist dataframe with a certain number of rows per patient, by averaging spectra randomly or according to a clustering algorithm.
group_peaklist <- function (peaklist, rows_per_patient=100, seed=NULL, balanced=TRUE, discard_poor_samples=FALSE, discard_if_lower_than=100, non_features=c("Sample","Class"), algorithm="random", clustering_method="agglomerative", grouping_variable="Sample") {
# Load the required libraries
install_and_load_required_packages(c("caret", "stats"))
# Results
result_data_frame <- data.frame()
# Create the file vector and the patient vector
file_vector <- peaklist[,grouping_variable]
patient_vector <- levels(file_vector)
#### If the dataset has to be balanced...
if (balanced == TRUE) {
	########## All the patients should have the same number of rows, always
	### Determine the shortest patient's spectra number
	lowest_observation_number <- NULL
	# Scroll the other patients
	for (p in 1:length(patient_vector)) {
		# Determine the number of observations for this patient
		observation_number <- nrow(peaklist [peaklist[,grouping_variable] == patient_vector[p],])
		# Set the lowest as the row number of the first patient and then check
		if (lowest_observation_number == NULL || observation_number < lowest_observation_number) {
			lowest_observation_number <- observation_number
		}
	}
	#### Rows per patient should not be more than the minimum number of rows
	if (rows_per_patient >= lowest_observation_number) {
		rows_per_patient <- lowest_observation_number
	} else {rows_per_patient <- rows_per_patient}
}
# If the spectra per patient is equal to one, do the normal averaging
if (rows_per_patient == 0 | rows_per_patient == 1) {
	# For each patient
	for (p in 1:length(patient_vector)) {
		# Subset the entire dataframe, isolating only the part corresponding to that patient
		patient_data_frame <- peaklist [peaklist[,grouping_variable] == patient_vector[p],]
		patient_data_frame_rows <- nrow(patient_data_frame)
		# Continue only if nothing has to be discarded or the dimension of the patient dataset is higher than the threshold
		if (discard_poor_samples == FALSE || (discard_poor_samples == TRUE && patient_data_frame_rows >= discard_if_lower_than)) {
		# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
		patient_data_frame_matrix <- as.matrix(patient_data_frame [,!(names(patient_data_frame) %in% non_features)])
		# Generate a single matrix row with the means of each column (feature)
		patient_data_frame_average <- apply(patient_data_frame_matrix, 2, mean)
		# Put this row together again with the Sample and the Class back into a dataframe
		patient_data_frame_final <- data.frame (rbind(patient_data_frame_average))
		for (i in 1:length(non_features)) {
			patient_data_frame_final <- data.frame (patient_data_frame_final, patient_data_frame[,non_features[i]][1])
		}
		# Fix the column names
		names(patient_data_frame_final) <- names(patient_data_frame)
		# Attach this row to the final dataframe
		result_data_frame <- rbind(result_data_frame, patient_data_frame_final)
		}
	}
}
if (rows_per_patient > 1) {
	# For each patient
	for (p in 1:length(patient_vector)) {
		# Subset the entire dataframe, isolating only the part corresponding to that patient
		patient_data_frame <- peaklist [peaklist[,grouping_variable] == patient_vector[p],]
		patient_data_frame_rows <- nrow(patient_data_frame)
		# Continue only if nothing has to be discarded or the dimension of the patient dataset is higher than the threshold
		if (discard_poor_samples == FALSE || (discard_poor_samples == TRUE && patient_data_frame_rows >= discard_if_lower_than)) {
			########################### RANDOMNESS
			if (algorithm == "random") {
				# Plant the seed only if a specified value is entered
				if (!is.null(seed)) {
					# Make the randomness reproducible
					set.seed(seed)
				}
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
				index <- createFolds(patient_data_frame[,grouping_variable], k = rows_per_patient)
				# For each fold
				for (k in 1:length(index)) {
					# Take the corresponding spectra subset (of the selected fold)
					fold_index <- index[[k]]
					patient_data_frame_temp <- patient_data_frame [foldIndex,]
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patient_data_frame_matrix <- as.matrix(patient_data_frame_temp [,!(names(patient_data_frame_temp) %in% non_features)])
					# Generate a single matrix row with the means of each column (feature)
					patient_data_frame_average <- apply(patient_data_frame_matrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patient_data_frame_final <- data.frame (rbind(patient_data_frame_average))
					for (i in 1:length(non_features)) {
						patient_data_frame_final <- data.frame (patient_data_frame_final, patient_data_frame[,non_features[i]][1])
					}
					# Fix the column names
					names(patient_data_frame_final) <- names(patient_data_frame_temp)
					# Attach this row to the final dataframe
					result_data_frame <- rbind(result_data_frame, patient_data_frame_final)
				}
			}
			############################ HCA
			if ((algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") && (clustering_method == "agglomerative" || is.null(clustering_method))) {
				# Compute the distance matrix
				distance_matrix <- dist(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], method="euclidean")
				# Plant the seed only if a specified value is entered
				if (!is.null(seed)) {
					# Make the randomness reproducible
					set.seed(seed)
				}
				hca <- hclust(distance_matrix)
				# Associate to each row/spectrum the subgroup of the tree to which it belongs
				# k must be between 1 and 56
				if (rows_per_patient > 56) {
					rows_per_patient <- 56
				}
				hca_groups <- cutree(hca, k=rows_per_patient)
				# For each subgroup to be isolated...
				for (s in 1:rows_per_patient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which(hca_groups == s)
					patient_data_frame_temp <- patient_data_frame [index,]
					# Average the spectra in this HCA subgroup
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patient_data_frame_matrix <- as.matrix(patient_data_frame_temp [,!(names(patient_data_frame_temp) %in% non_features)])
					# Generate a single matrix row with the means of each column (feature)
					patient_data_frame_average <- apply(patient_data_frame_matrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patient_data_frame_final <- data.frame (rbind(patient_data_frame_average))
					for (i in 1:length(non_features)) {
						patient_data_frame_final <- data.frame (patient_data_frame_final, patient_data_frame[,non_features[i]][1])
					}
					# Fix the column names
					names(patient_data_frame_final) <- names(patient_data_frame_temp)
					# Attach this row to the final dataframe
					result_data_frame <- rbind(result_data_frame, patient_data_frame_final)
				}
			}
			######################### K-MEANS
			if ((algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") && (clustering_method == "kMeans" | clustering_method == "kmeans" | clustering_method == "k-Means")) {
				# Plant the seed only if a specified value is entered
				if (!is.null(seed)) {
					# Make the randomness reproducible
					set.seed(seed)
				}
				# Compute the k-Means clustering
				k_means <- kmeans(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], centers=rows_per_patient)
				# Associate to each row/spectrum the subgroup/cluster to which it belongs
				k_means_groups <- k_means$cluster
				# For each subgroup to be isolated...
				for (s in 1:rows_per_patient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which(k_means_groups == s)
					patient_data_frame_temp <- patient_data_frame [index,]
					# Average the spectra in this kMeans subgroup
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patient_data_frame_matrix <- as.matrix(patient_data_frame_temp [,!(names(patient_data_frame_temp) %in% non_features)])
					# Generate a single matrix row with the means of each column (feature)
					patient_data_frame_average <- apply(patient_data_frame_matrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patient_data_frame_final <- data.frame (rbind(patient_data_frame_average))
					for (i in 1:length(non_features)) {
						patient_data_frame_final <- data.frame (patient_data_frame_final, patient_data_frame[,non_features[i]][1])
					}
					# Fix the column names
					names(patient_data_frame_final) <- names(patient_data_frame_temp)
					# Attach this row to the final dataframe
					result_data_frame <- rbind(result_data_frame, patient_data_frame_final)
				}
			}
		}
	}
}
return (result_data_frame)
}





################################################################################





##################### MATRIX/DATAFRAME SPLITTING FUNCTION (TRAINING AND TESTING)
# This function takes a peaklist matrix as input and outputs a list containing the part of the dataset used for the training and the part to be used as test dataset. The split happens onto the discriminant (class) variable.
# The split is made according to the user desire and according to a selected column (usually the class column).
# All the entries (rows) are considered independent observations.
matrix_splitting_training_test <- function(peaklist, discriminant_feature="Class", seed=123456, percentage_of_observation_for_training=66) {
	# Install the required packages
	install_and_load_required_packages("caret")
	if (!is.null(discriminant_feature)) {
		##### Determine the class list
		class_list <- levels(as.factor(peaklist[,discriminant_feature]))
	} else {
		class_list <- character()
	}
	### If there are no classes or one class, just randomly select rows on the entire dataset
	if (length(class_list) <= 1) {
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		index_training <- sample(nrow(peaklist), size=round(nrow(peaklist)*percentage_of_observation_for_training/100))
	} else if (length(class_list) > 1) {
		### If there are more classes, randomly select the rows for each class
		# Create a list of dataframes, one for each class
		class_data_frame_list <- list()
		for (j in 1:length(class_list)) {
			class_data_frame_list[[j]] <- peaklist[peaklist[,discriminant_feature] == class_list[j],]
		}
		### Create a list containing the training and testing indexes of the class dataframes
		index_training <- list()
		for (i in 1:length(class_data_frame_list)) {
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			index_training[[i]] <- createDataPartition(y = class_data_frame_list[[i]][,discriminant_feature], p = (percentage_of_observation_for_training/100), list = FALSE)
		}
		names(index_training) <- class_list
	}
	### Now we have randomly selected some patients for the training and some others for the testing, for each class
	### The indexTraining is integer if there are no classes or one class, and rows are selected randomly
	if (is.integer(index_training)) {
		training_dataset <- peaklist[index_training,]
		test_dataset <- peaklist[-index_training,]
	} else if (is.list(index_training)) {
		# The indexTraining is a list if there are more than one classes (one element per class), and rows are selected randomly for each class
		### Create the training and the testing patient dataframes
		training_data_frame_list <- list()
		test_data_frame_list <- list()
		# For each class (and so for each element of the list with the indexes per class)
		for (i in 1:length(index_training)) {
			# Create the two complementary dataframes: randomly select the patients for each class
			# for training and testing
			training_data_frame_list[[i]] <- class_data_frame_list[[i]][index_training[[i]],]
			test_data_frame_list[[i]] <- class_data_frame_list[[i]][-index_training[[i]],]
		}
		### Merge the two traininf and test sets
		training_dataset <- NULL
		test_dataset <- NULL
		for (p in 1:length(training_data_frame_list)) {
			if (is.null(training_dataset)) {
				training_dataset <- training_data_frame_list[[p]]
			} else {
				training_dataset <- rbind(training_dataset, training_data_frame_list[[p]])
			}
		}
		for (p in 1:length(test_data_frame_list)) {
			if (is.null(test_dataset)) {
				test_dataset <- test_data_frame_list[[p]]
			} else {
				test_dataset <- rbind(test_dataset, test_data_frame_list[[p]])
			}
		}
	}
	return(list(training_dataset = training_dataset, test_dataset = test_dataset))
}





################################################################################





############################################################## FEATURE SELECTION
# This function runs the feature selection algorithm onto the peaklist matrix, returning the peaklist without the redundant/non-informative features, the original peaklist and the list of selected features.
# The function allows for the use of several feature selection algorithms.
feature_selection <- function (peaklist, feature_selection_method="rfe", features_to_select=20, selection_method="pls", selection_metric="Kappa", correlation_method="pearson", correlation_threshold=0.75, auc_threshold=0.7, cv_repeats_control=5, k_fold_cv_control=10, discriminant_attribute="Class", non_features=c("Sample", "Class", "THY"), seed=NULL, automatically_select_features=FALSE, generate_plots=TRUE, preprocessing=c("center","scale"), multicore_processing=TRUE, feature_reranking=FALSE) {
# Load the required libraries
install_and_load_required_packages(c("caret", "pls", "stats", "randomForest", "pROC"))
### PARALLEL BACKEND
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE) - 1
if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
	install_and_load_required_packages("doMC")
	# Register the foreach backend
	registerDoMC(cores = cpu_thread_number)
} else if (Sys.info()[1] == "Windows") {
	install_and_load_required_packages("doParallel")
	# Register the foreach backend
	cl <- makeCluster(cpu_thread_number, type='PSOCK')
	registerDoParallel(cl)
}
# Initialization
feature_weights <- NULL
variable_importance <- NULL
fs_model_performance <- NULL
########################################################### RFE MODEL
if (feature_selection_method == "rfe" || feature_selection_method == "recursive feature elimination") {
	# The simulation will fit models with subset sizes: (the subset size is the number of predictors to use)
	if (automatically_select_features == TRUE) {
		if (selection_method != "") {
			# Define the feature subset sizes
			subset_sizes <- seq(2, features_to_select, by = 1)
			# Set the seeds (for each resampling iteration in the cv)
			number_of_iterations <- (k_fold_cv_control * cv_repeats_control) + 1
			seeds <- vector(mode = "list", length = number_of_iterations)
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			for (s in 1:(number_of_iterations - 1)) {
				seeds[[s]] <- sample.int(nrow(peaklist), length(subset_sizes) + 1, replace = TRUE)
			}
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			seeds[number_of_iterations] <- sample.int(nrow(peaklist), 1)
			# Define the control function of the RFE
			rfe_ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, seeds=seeds, saveDetails=TRUE, allowParallel = multicore_processing, rerank = feature_reranking)
			# Run the RFE
			rfe_model <- rfe(x = peaklist[,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, method = selection_method, metric = selection_metric, preProcess = preprocessing)
			# Variable importance
			variable_importance <- varImp(rfe_model$fit)
			# Feature weights
			feature_weights <- rfe_model$fit
			# Model performances
			if (selection_metric == "kappa" || selection_metric == "Kappa") {
				fs_model_performance <- max(rfe_model$fit$results$Kappa)
			} else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
				fs_model_performance <- max(rfe_model$fit$results$Accuracy)
			}
		} else {
			# Define the feature subset sizes
			subset_sizes <- seq(2, features_to_select, by = 1)
			# Set the seeds (for each resampling iteration in the cv)
			number_of_iterations <- (k_fold_cv_control * cv_repeats_control) +1
			seeds <- vector(mode = "list", length = number_of_iterations)
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			for (s in 1:(number_of_iterations - 1)) {
				seeds[[s]] <- sample.int(nrow(peaklist), length(subset_sizes) + 1, replace = TRUE)
			}
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			seeds[number_of_iterations] <- sample.int(nrow(peaklist), 1)
			# Define the control function of the RFE
			rfe_ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, seeds=seeds, saveDetails = TRUE, allowParallel = multicore_processing, rerank = feature_reranking)
			# Run the RFE
			rfe_model <- rfe(x = peaklist[,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, metric = selection_metric, preProcess = preprocessing)
			# Variable importance
			variable_importance <- varImp(rfe_model$fit)
			# Feature weights
			feature_weights <- rfe_model$fit
			# Model performances
			if (selection_metric == "kappa" || selection_metric == "Kappa") {
				fs_model_performance <- max(rfe_model$fit$results$Kappa)
			} else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
				fs_model_performance <- max(rfe_model$fit$results$Accuracy)
			}
		}
		# Output the best predictors after the RFE
		predictors_rfe <- predictors(rfe_model)
	} else {
		############# NO AUTOMATIC DETERMINATION
		if (selection_method != "") {
			# Define the feature subset sizes
			subset_sizes <- features_to_select
			# Set the seeds (for each resampling iteration in the cv)
			number_of_iterations <- (k_fold_cv_control * cv_repeats_control) + 1
			seeds <- vector(mode = "list", length = number_of_iterations)
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			for (s in 1:(number_of_iterations - 1)) {
				seeds[[s]] <- sample.int(nrow(peaklist), length(subset_sizes) + 1, replace = TRUE)
			}
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			seeds[number_of_iterations] <- sample.int(nrow(peaklist), 1)
			# Define the control function of the RFE
			rfe_ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, seeds = seeds, saveDetails = TRUE, allowParallel = multicore_processing, rerank = feature_reranking)
			# Run the RFE
			rfe_model <- rfe(x = peaklist[,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, method = selection_method, metric = selection_metric, preProcess = preprocessing)
			# Variable importance
			variable_importance <- varImp(rfe_model$fit)
			# Feature weights
			feature_weights <- rfe_model$fit
			# Model performances
			if (selection_metric == "kappa" || selection_metric == "Kappa") {
				fs_model_performance <- max(rfe_model$fit$results$Kappa)
			} else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
				fs_model_performance <- max(rfe_model$fit$results$Accuracy)
			}
		} else {
			# Define the feature subset sizes
			subset_sizes <- features_to_select
			# Set the seeds (for each resampling iteration in the cv)
			number_of_iterations <- (k_fold_cv_control * cv_repeats_control) + 1
			seeds <- vector(mode = "list", length = number_of_iterations)
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			for (s in 1:(number_of_iterations - 1)) {
				seeds[[s]] <- sample.int(nrow(peaklist), length(subset_sizes) + 1, replace = TRUE)
			}
			# Plant the seed only if a specified value is entered
			if (!is.null(seed)) {
				# Make the randomness reproducible
				set.seed(seed)
			}
			seeds[number_of_iterations] <- sample.int(nrow(peaklist), 1)
			# Define the control function of the RFE
			rfe_ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, seeds = seeds, saveDetails = TRUE, allowParallel = multicore_processing, rerank = feature_reranking)
			# Run the RFE
			rfe_model <- rfe(x = peaklist[,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, metric = selection_metric, preProcess = preprocessing)
			# Variable importance
			variable_importance <- varImp(rfe_model$fit)
			# Feature weights
			feature_weights <- rfe_model$fit
			# Model performances
			if (selection_metric == "kappa" || selection_metric == "Kappa") {
				fs_model_performance <- max(rfe_model$fit$results$Kappa)
			} else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
				fs_model_performance <- max(rfe_model$fit$results$Accuracy)
			}
		}
		# Output the best predictors after the RFE
		predictors_rfe <- predictors(rfe_model) [1:features_to_select]
	}
	# Predictors
	predictors_feature_selection <- predictors_rfe
}
######################################################################### ANOVA
if (feature_selection_method == "ANOVA") {
	feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
	features_ANOVA <- character()
	# For each feature, calculate the impact on the classification capability by fitting an ANOVA
	for (feature in feature_list){
		# Isolate the peak intensities
		intensity_vector <- peaklist [,feature]
		# Compute an ANOVA based upon the discriminant attribute
		feature_ANOVA <- aov(intensity_vector ~ peaklist[,discriminant_attribute])
		# Extract the p-value
		feature_ANOVA_pvalue <- summary(feature_ANOVA)[[1]]$"Pr(>F)"[1]
		# Keep the feature if with a great impact
		if (feature_ANOVA_pvalue <= 0.05) {
			features_ANOVA <- append(features_ANOVA, feature)
		}
	}
	# Predictors
	predictors_feature_selection <- features_ANOVA
}
################################################################ KRUSKAL-WALLIS
if (feature_selection_method == "kruskal" || feature_selection_method == "kruskal-wallis" || feature_selection_method == "Kruskal-Wallis") {
	feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
	features_kruskal <- character()
	# For each feature, calculate the impact on the classification capability by fitting an ANOVA
	for (feature in feature_list){
		# Isolate the peak intensities
		intensity_vector <- peaklist [,feature]
		# Compute an ANOVA based upon the discriminant attribute
		feature_kruskal <- kruskal.test(intensity_vector ~ peaklist[,discriminant_attribute])
		# Extract the p-value
		feature_kruskal_pvalue <- feature_kruskal$p.value
		# Keep the feature if with a great impact
		if (feature_kruskal_pvalue <= 0.05) {
			features_kruskal <- append(features_kruskal, feature)
		}
	}
	# Predictors
	predictors_feature_selection <- features_kruskal
}
################################################# CORRELATION FEATURE SELECTION
if (feature_selection_method == "correlation") {
	# Take only the part of the matrix without Class and Sample
	peaklist_features <- peaklist [,!(names(peaklist) %in% non_features)]
	# Compute the correlation
	feature_correlation <- cor(peaklist_features, method=correlation_method)
	# Output the highly correlated features
	highly_correlated <- findCorrelation(feature_correlation, correlation_threshold)
	# List the highly correlated features
	highly_correlated_features <- names(peaklist_features[,highly_correlated])
	# Features to keep
	low_correlation_features <- names(peaklist_features[,-highly_correlated])
	print(paste("The number of selected features is", length(low_correlation_features), "out of", length(peaklist_features)))
	# Predictors
	predictors_feature_selection <- low_correlation_features
}
#################################################################### IMPORTANCE
if (feature_selection_method == "importance") {
	feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
	# For each feature, calculate the impact on the classification capability
	model_control <- trainControl(method="repeatedcv", number=k_fold_cv_control, repeats=cv_repeats_control)
	# Compute a model based upon the discriminant attribute
	feature_model <- train(x=peaklist [,!(names(peaklist) %in% non_features)], y=peaklist[,discriminant_attribute], method="pls", preProcess="scale", trControl=model_control)
	# Estimate variable importance
	feature_importance <- varImp(feature_model, scale=FALSE)
	# Isolate the most important features (rank them according to their importance first!)
	feature_importance_df <- feature_importance$importance
	feature_importance_df$Features <- rownames(feature_importance_df)
	feature_importance_df <- feature_importance_df [order(-feature_importance_df$Overall),]
	# Predictors
	predictors_feature_selection <- feature_importance_df$Features [1:features_to_select]
}
########################################################################### ROC
if (feature_selection_method == "ROC" || feature_selection_method == "roc") {
	feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
	# List of important features
	features_ROC <- character()
	##### Automatically select features
	if (automatically_select_features == TRUE) {
		# For each feature, calculate the impact on the classification capability by computing a ROC
		for (feature in feature_list){
			# Compute the ROC of the feature
			feature_ROC <- roc(response=peaklist[,discriminant_attribute], predictor=peaklist[,feature])
			# Extract the AUC
			feature_ROC_AUC <- feature_ROC$auc
			# Keep the feature if with a great impact
			if (feature_ROC_AUC >= auc_threshold) {
				features_ROC <- append(features_ROC, feature)
			}
		}
	} else {
		##### Select the most N important features
		# For each feature, calculate the impact on the classification capability by computing a ROC
		feature_ROC_vector <- numeric(length=length(feature_list))
		names(feature_ROC_vector) <- feature_list
		for (f in 1:length(feature_list)) {
			# Compute the ROC of the feature
			feature_ROC <- roc(response=peaklist[,discriminant_attribute], predictor=peaklist[,feature_list[f]])
			# Append the AUC to a vector
			feature_ROC_vector[f] <- feature_ROC$auc
		}
		# Sort the vector
		feature_ROC_vector_sorted <- sort(feature_ROC_vector, decreasing=TRUE)
		# Take the first N features
		features_ROC <- names(feature_ROC_vector_sorted)[1:features_to_select]
	}
	# Predictors
	predictors_feature_selection <- features_ROC
}
########################################################################## Plot
if (generate_plots == TRUE) {
	# Initialize
	feature_selection_graphics <- NULL
	if (feature_selection_method == "rfe" || feature_selection_method == "recursive feature elimination") {
		feature_selection_graphics <- plot(rfe_model, type=c("g","o"))
	}
	if (feature_selection_method == "correlation") {
	}
	if (feature_selection_method == "importance") {
		feature_selection_graphics <- plot(feature_importance)
	}
} else {feature_selection_graphics <- NULL}
#################################################### Take the selected features
peaklist_feature_selection <- peaklist [,predictors_feature_selection]
# Add the non features back
for (i in 1:length(non_features)) {
	peaklist_feature_selection <- cbind(peaklist_feature_selection, peaklist[,non_features[i]])
}
names(peaklist_feature_selection) <- c(as.character(predictors_feature_selection), non_features)
# Return the values
return (list(peaklist_feature_selection=peaklist_feature_selection, predictors_feature_selection=predictors_feature_selection, feature_selection_graphics=feature_selection_graphics, feature_weights=feature_weights, variable_importance=variable_importance, fs_model_performance=fs_model_performance))
}





################################################################################





####################################################### TRUNCATE PEAKLIST MATRIX
# This function operates the truncation at the peaklist level, by removing peaks (columns) out of a certain range.
truncate_peaklist <- function (peaklist, range=c(4000, 15000), non_features=c("Sample", "Class")){
features_to_keep <- numeric()
# Do not do anything if the values are set to zero
if (range[1] == 0 & range[2] == 0) {
	# Keep all the features
	features_to_keep <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
}
#
if (range[1] != 0 & range[2] == 0) {
	feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
	# Keep only the features that are above the lower threshold
	for (f in 1:length(feature_list)) {
		if (feature_list[f] >= range[1]) {
			features_to_keep <- append(features_to_keep, feature_list[f])
		}
	}
}
#
if (range[1] != 0 & range[2] != 0) {
	feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
	# Keep only the features that are above the lower threshold and below the upper threshold
	for (f in 1:length(feature_list)) {
		if (feature_list[f] >= range[1] && feature_list[f] <= range[2]) {
			features_to_keep <- append(features_to_keep, feature_list[f])
		}
	}
}
# Keep only the features of interest
peaklist_truncated <- data.frame (peaklist [,as.character(features_to_keep)])
# Add the non features back
for (i in 1:length(non_features)) {
	peaklist_truncated <- data.frame (peaklist_truncated, peaklist[,non_features[i]])
}
names(peaklist_truncated) <- c(features_to_keep, non_features)
#
return (peaklist_truncated)
}





################################################################################





############################ CROSS-VALIDATION FOR A SUPPORT VECTOR MACHINE MODEL
# This function operates a cross-validation to evaluate the performances of the SVM model, trained (and tuned) onto the dataset. The cross-validation is made by training the SVM model (with the same parameters) onto the training subset and evaluating its performances onto the testing subset.
# It returns a matrix with the cross-validation results (the average performance of all the validations performed) along with the results from each validation (for each k of the k-fold). It returns also the seed used for making the randomness reproducible and the list of features.
# It is advisable to have one row per patient in the peaklist matrix.
cross_validation_svm <- function (training_dataset, seed=NULL, k_fold=10, repeats=5, svm_model, non_features=c("Sample","Class","THY"), positive_class=levels(training_dataset$Class)[1]) {
# Plant the seed only if a specified value is entered
if (!is.null(seed)) {
	# Make the randomness reproducible
	set.seed(seed)
}
### Extract the svm model parameters
# Kernel
svm_kernel <- svm_model$kernel
if (svm_kernel == 0) {
	svm_kernel <- "linear"
}
if (svm_kernel == 1) {
	svm_kernel <- "polynomial"
}
if (svm_kernel == 2) {
	svm_kernel <- "radial"
}
if (svm_kernel == 3) {
	svm_kernel <- "sigmoid"
}
svm_cost <- svm_model$cost
svm_degree <- svm_model$degree
svm_gamma <- svm_model$gamma
### Output matrix
confusion_matrix_list <- list()
result_matrix <- matrix (NA, ncol=14, nrow=0)
colnames(result_matrix) <- c("Accuracy", "Kappa", "No information rate", "Accuracy p-value", "McNiemar p-value", "Sensitivity", "Specificity", "PPV", "NPV", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy", "ROC AUC")
#### For each repetition...
for (i in 1:repeats) {
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed*i)
	}
	# Index of the spectra to be allocated in the folds, randomly selected
	# k equal-size partitions are created
	index <- createFolds(training_dataset$Class, k = k_fold)
	# Each single fold has to be used as a training and the rest as testing
	for (k in 1:length(index)) {
		# Generate the training and the testing subsets (based upon the index of the kFold)
		# The kFold indexed with k is used each time as testing
		test_subset <- training_dataset [index[[k]],]
		test_predictors <- test_subset [,!(names(test_subset) %in% non_features)]
		test_outcomes <- test_subset$Class
		# All the other ones are used as training
		training_subset <- training_dataset [-index[[k]],]
		training_predictors <- training_subset [,!(names(training_subset) %in% non_features)]
		training_outcomes <- training_subset$Class
		# Now the model has to be built using the training and tested onto the testing dataset
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		model <- svm (x=training_predictors, y=training_outcomes, scale=TRUE, kernel=svm_kernel, degree=svm_degree, gamma=svm_gamma, cost=svm_cost)
		# Use the model to predict the testing subset
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes <- predict(model, newdata=test_predictors)
		# Compute the performances
		performances <- confusionMatrix(data=predicted_classes, test_outcomes, positive=positive_class)
		confusion_matrix_list [[k*i]] <- performances
		#### ROC analysis
		svm_roc <- roc(response=test_outcomes, predictor=as.numeric(predicted_classes))
		# Create a matrix row with the partial results
		result_matrix_row <- matrix (0, ncol=ncol(result_matrix), nrow=1)
		result_matrix_row [1,1] <- as.numeric(performances$overall[1])
		result_matrix_row [1,2] <- as.numeric(performances$overall[2])
		result_matrix_row [1,3] <- as.numeric(performances$overall[5])
		result_matrix_row [1,4] <- as.numeric(performances$overall[6])
		result_matrix_row [1,5] <- as.numeric(performances$overall[7])
		result_matrix_row [1,6] <- as.numeric(performances$byClass[1])
		result_matrix_row [1,7] <- as.numeric(performances$byClass[2])
		result_matrix_row [1,8] <- as.numeric(performances$byClass[3])
		result_matrix_row [1,9] <- as.numeric(performances$byClass[4])
		result_matrix_row [1,10] <- as.numeric(performances$byClass[5])
		result_matrix_row [1,11] <- as.numeric(performances$byClass[6])
		result_matrix_row [1,12] <- as.numeric(performances$byClass[7])
		result_matrix_row [1,13] <- as.numeric(performances$byClass[8])
		result_matrix_row [1,14] <- as.numeric(svm_roc$auc)
		# Attach the row to the result matrix
		result_matrix <- rbind(result_matrix, result_matrix_row)
	}
}
# Compute the average of each column
result_matrix_average <- apply(result_matrix, 2, mean, na.rm=TRUE)
#return (list(confusion_matrix=confusion_matrix_list, result_matrix=result_matrix, features=names(training_dataset)[!(names(training_dataset) %in% non_features)]))
return (list(result_matrix_average=result_matrix_average))
}





################################################################################





###################################################### SVM TUNING AND VALIDATION
# This function operates the tuning of the Support Vector Machine (SVM), by testing all the parameters of the SVM (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
svm_tuning_and_validation2 <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(gamma=10^(-5:5), cost=10^(-5:5), epsilon=seq(1,2,by=1), degree=1:5, kernel="radial"), k_fold_cv=10, repeats_cv=2, parameters=list(gamma=0.1, cost=10, epsilon=0.1, degree=3, kernel="radial"), positive_class_cv="HP", seed=NULL, pca=FALSE, number_of_components=3) {
# Load the required libraries
install_and_load_required_packages(c("caret", "kernlab", "e1071", "pROC"))
### PARALLEL BACKEND
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE) - 1
if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
	install_and_load_required_packages("doMC")
	# Register the foreach backend
	registerDoMC(cores = cpu_thread_number)
} else if (Sys.info()[1] == "Windows") {
	install_and_load_required_packages("doParallel")
	# Register the foreach backend
	cl <- makeCluster(cpu_thread_number, type='PSOCK')
	registerDoParallel(cl)
}
######################################## PCA
if (pca == TRUE) {
	# Compute the PCs
	pca_training <- prcomp (as.matrix(peaklist_training [,!(names(peaklist_training) %in% non_features)]))
	pca_test <- prcomp (as.matrix(peaklist_test [,!(names(peaklist_test) %in% non_features)]))
	#################### SVM tuning
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	if (autotuning == TRUE) {
		svm_tuning <- tune.svm(as.data.frame(pca_training$x[,1:number_of_components]), factor(peaklist_training$Class), gamma=tuning_parameters$gamma, cost=tuning_parameters$cost, kernel=tuning_parameters$kernel, epsilon=tuning_parameters$epsilon)
		# Select automatically the best model after the tuning
		svm_model <- svm_tuning$best.model
		# Parameters output
		svm_kernel <- svm_tuning$best.model$kernel
		if (svm_kernel == 0) {
			svm_kernel <- "linear"
		}
		if (svm_kernel == 1) {
			svm_kernel <- "polynomial"
		}
		if (svm_kernel == 2) {
			svm_kernel <- "radial"
		}
		if (svm_kernel == 3) {
			svm_kernel <- "sigmoid"
		}
		parameters_output <- data.frame (kernel=svm_kernel, cost=svm_tuning$best.model$cost, degree=svm_tuning$best.model$degree, epsilon=svm_tuning$best.model$epsilon, gamma=svm_tuning$best.model$gamma)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		svm_model <- svm(pca_training, factor(peaklist_training$Class), kernel=parameters$kernel, cost=parameters$cost, epsilon=parameters$epsilon, gamma=parameters$gamma)
		parameters_output <- list (kernel=parameters$kernel, cost=parameters$cost, degree=parameters$degree, epsilon=parameters$epsilon, gamma=parameters$gamma)
	}
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predictthe outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_svm <- predict(svm_model, newdata = as.data.frame(pca_test$x[,1:number_of_components]))
		# Create the outcomes dataframe
		classification_results_svm <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		svm_roc <- list()
		roc_curve <- roc(response=classification_results_svm$true, predictor=as.numeric(classification_results_svm$predicted))
		svm_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste ("ROC area under the curve:", roc_curve$auc)
		legend ("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		svm_roc[[2]] <- recordPlot()
		# Output the results
		return (list (model=svm_model, svm_features=colnames(svm_model$SV), svm_parameters=parameters_output, classification_results=classification_results_svm, performances=test_performances_svm, roc=svm_roc))
	}	else {return (list (model=svm_model, svm_features=colnames(svm_model$SV), svm_parameters=parameters_output, cross_validation=cv_svm_model))}
} else {
	######################################## NO PCA
	################# FEATURES
	if (autotuning == TRUE) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma=tuning_parameters$gamma, cost=tuning_parameters$cost, kernel=tuning_parameters$kernel, epsilon=tuning_parameters$epsilon, degree=tuning_parameters$degree)
		# Select automatically the best model after the tuning
		svm_model <- svm_tuning$best.model
		# Parameters output
		svm_kernel <- svm_tuning$best.model$kernel
		if (svm_kernel == 0) {
			svm_kernel <- "linear"
		}
		if (svm_kernel == 1) {
			svm_kernel <- "polynomial"
		}
		if (svm_kernel == 2) {
			svm_kernel <- "radial"
		}
		if (svm_kernel == 3) {
			svm_kernel <- "sigmoid"
		}
		parameters_output <- data.frame (kernel=svm_kernel, cost=svm_tuning$best.model$cost, degree=svm_tuning$best.model$degree, epsilon=svm_tuning$best.model$epsilon, gamma=svm_tuning$best.model$gamma)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		svm_model <- svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), kernel=parameters$kernel, cost=parameters$cost, epsilon=parameters$epsilon, gamma=parameters$gamma)
		parameters_output <- data.frame (kernel=parameters$kernel, cost=parameters$cost, degree=parameters$degree, epsilon=parameters$epsilon, gamma=parameters$gamma)
	}
	#################### CROSS-VALIDATION
	cv_svm_model <- cross_validation_svm(peaklist_training, k_fold=k_fold_cv, repeats=repeats_cv, svm_model=svm_model, non_features=non_features, positive_class=positive_class_cv, seed=seed)
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_svm <- predict(svm_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
		# Create the outcomes dataframe
		classification_results_svm <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		svm_roc <- list()
		roc_curve <- roc(response=as.numeric(classification_results_svm$true), predictor=as.numeric(classification_results_svm$predicted))
		svm_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
		legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		svm_roc[[2]] <- recordPlot()
		###################### Pie chart classification
		correctly_classified <- 0
		misclassified <- 0
		for (i in 1:nrow(classification_results_svm)) {
			if (classification_results_svm$predicted[i] == classification_results_svm$true[i]) {
				correctly_classified <- correctly_classified + 1
			} else {
				misclassified <- misclassified + 1
			}
		}
		classification_pie <- c(correctly_classified, misclassified)
		pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
		pie_chart_classification <- recordPlot()
		# Output the results
		return (list(model=svm_model, svm_features=colnames(svm_model$SV), cross_validation=cv_svm_model, classification_results=classification_results_svm, parameters=parameters_output, performances=test_performances_svm, roc=svm_roc, pie_chart_classification=pie_chart_classification))
	} else {return (list(model=svm_model, svm_features=colnames(svm_model$SV), parameters=parameters_output, cross_validation=cv_svm_model))}
	}
}





################################################################################





###################################################### SVM TUNING AND VALIDATION
# This function operates the tuning of the Support Vector Machine (SVM), by testing all the parameters of the SVM (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
svm_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(sigma=10^(-5:5), cost=10^(-5:5), epsilon=seq(1,2,by=1), degree=1:5, scale=1, kernel="radial"), k_fold_cv=10, repeats_cv=2, preprocessing=c("center","scale"), parameters=list(sigma=0.001, scale=1, gamma=0.1, cost=10, epsilon=0.1, degree=3, kernel="radial"), positive_class_cv="HP", seed=NULL, evaluation_method="Accuracy") {
# Load the required libraries
install_and_load_required_packages(c("caret", "pROC", "kernlab"))
### PARALLEL BACKEND
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE) - 1
if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
	install_and_load_required_packages("doMC")
	# Register the foreach backend
	registerDoMC(cores = cpu_thread_number)
} else if (Sys.info()[1] == "Windows") {
	install_and_load_required_packages("doParallel")
	# Register the foreach backend
	cl <- makeCluster(cpu_thread_number, type='PSOCK')
	registerDoParallel(cl)
}
################# FEATURES
if (autotuning == TRUE) {
	# Find the best tuning parameters for the SVM
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	# Define the control function for training
	training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
	# Train and tune the model
	if (tuning_parameters$kernel == "radial") {
		if (!is.null(tuning_parameters$sigma)) {
			# Define the tune grid for the model
			svm_tune_grid <- expand.grid(sigma=tuning_parameters$sigma, C=tuning_parameters$cost)
			# Training and tuning
			svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmRadial", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
		} else {
			# Define the tune grid for the model
			svm_tune_grid <- expand.grid(C=tuning_parameters$cost)
			# Training and tuning
			svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmRadialCost", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
		}
	} else if (tuning_parameters$kernel == "polynomial") {
		# Define the tune grid for the model
		svm_tune_grid <- expand.grid(scale=tuning_parameters$scale, C=tuning_parameters$cost, degree=tuning_parameters$degree)
		# Training and tuning
		svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmPoly", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
	} else if (tuning_parameters$kernel == "linear") {
		# Define the tune grid for the model
		svm_tune_grid <- expand.grid(cost=tuning_parameters$cost, gamma=tuning_parameters$gamma)
		# Training and tuning
		svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmLinear2", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
	}
	# Extract the cross-validation information
	cv_svm_model <- confusionMatrix(svm_tuning)
	# Select automatically the best model after the tuning
	svm_model <- svm_tuning$finalModel
	# Parameters output
	if (tuning_parameters$kernel == "radial") {
		parameters_output <- data.frame(kernel=tuning_parameters$kernel, cost=svm_model@param$C, sigma=svm_model@kernelf@kpar[[1]])
	} else if (tuning_parameters$kernel == "polynomial") {
		parameters_output <- data.frame(kernel=tuning_parameters$kernel, cost=svm_model@param$C, degree=svm_model@kernelf@kpar$degree, scale=svm_model@kernelf@kpar$scale)
	} else if (tuning_parameters$kernel == "linear") {
		parameters_output <- data.frame(kernel=tuning_parameters$kernel, cost=svm_model$cost, gamma=svm_model$gamma)
	}
	# Plots
	tuning_plot <- plot(svm_tuning)
}
#################### SVM with defined parameters
if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
	# Find the best tuning parameters for the SVM
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	# Define the control function for training
	training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
	# Train the model
	if (tuning_parameters$kernel == "radial") {
		if (!is.null(tuning_parameters$sigma)) {
			# Define the tune grid for the model
			svm_tune_grid <- expand.grid(sigma=parameters$sigma, C=parameters$cost)
			# Training and tuning
			svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmRadial", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
		} else {
			# Define the tune grid for the model
			svm_tune_grid <- expand.grid(C=parameters$cost)
			# Training and tuning
			svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmRadialCost", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
		}
	} else if (tuning_parameters$kernel == "polynomial") {
		# Define the tune grid for the model
		svm_tune_grid <- expand.grid(scale=parameters$scale, C=parameters$cost, degree=parameters$degree)
		# Training and tuning
		svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmPoly", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
	} else if (tuning_parameters$kernel == "linear") {
		# Define the tune grid for the model
		svm_tune_grid <- expand.grid(cost=parameters$cost, gamma=parameters$gamma)
		# Training and tuning
		svm_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="svmLinear2", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=svm_tune_grid)
	}
	# Extract the cross-validation information
	cv_svm_model <- confusionMatrix(svm_tuning)
	# Select automatically the best model after the tuning
	svm_model <- svm_tuning$finalModel
	# Parameters output
	if (tuning_parameters$kernel == "radial") {
		parameters_output <- data.frame(kernel=parameters$kernel, cost=svm_model@param$C, sigma=svm_model@kernelf@kpar[[1]])
	} else if (tuning_parameters$kernel == "polynomial") {
		parameters_output <- data.frame(kernel=parameters$kernel, cost=svm_model@param$C, degree=svm_model@kernelf@kpar$degree, scale=svm_model@kernelf@kpar$scale)
	} else if (tuning_parameters$kernel == "linear") {
		parameters_output <- data.frame(kernel=parameters$kernel, cost=svm_model$cost, gamma=svm_model$gamma)
	}
}
######. Close the cluster for parallel computation
if (Sys.info()[1] == "Windows") {
	stopCluster(cl)
}
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	predicted_classes_svm <- predict(svm_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
	# Create the outcomes dataframe
	classification_results_svm <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	svm_roc <- list()
	roc_curve <- roc(response=as.numeric(classification_results_svm$true), predictor=as.numeric(classification_results_svm$predicted))
	svm_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	svm_roc[[2]] <- recordPlot()
	###################### Pie chart classification
	correctly_classified <- 0
	misclassified <- 0
	for (i in 1:nrow(classification_results_svm)) {
		if (as.character(classification_results_svm$predicted[i]) == as.character(classification_results_svm$true[i])) {
			correctly_classified <- correctly_classified + 1
		} else {
			misclassified <- misclassified + 1
		}
	}
	classification_pie <- c(correctly_classified, misclassified)
	pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
	pie_chart_classification <- recordPlot()
	# Output the results
	return (list(model=svm_model, svm_features=colnames(svm_model@xmatrix[[1]]), cross_validation=cv_svm_model, classification_results=classification_results_svm, parameters=parameters_output, tuning_plot=tuning_plot, performances=test_performances_svm, roc=svm_roc, pie_chart_classification=pie_chart_classification))
} else {return (list(model=svm_model, svm_features=colnames(svm_model@xmatrix[[1]]), parameters=parameters_output, tuning_plot=tuning_plot, cross_validation=cv_svm_model))}
}





################################################################################





###################################################### PLS TUNING AND VALIDATION
# This function operates the training of a Partial Least Squares (PLS) model, by testing all the parameters of the PLS and selecting the best (in terms of accuracy). The tuning is performed via cross-validation onto the training dataset.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
pls_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), tuning_parameters=data.frame(ncomp=1:5), k_fold_cv=10, repeats_cv=2, positive_class_cv="HP", seed=NULL, preprocessing=c("center","scale"), selection_criteria="Accuracy", maximise_selection_criteria_values=TRUE) {
# Load the required libraries
install_and_load_required_packages(c("caret", "e1071"))
### PARALLEL BACKEND
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE) - 1
if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
	install_and_load_required_packages("doMC")
	# Register the foreach backend
	registerDoMC(cores = cpu_thread_number)
} else if (Sys.info()[1] == "Windows") {
	install_and_load_required_packages("doParallel")
	# Register the foreach backend
	cl <- makeCluster(cpu_thread_number, type='PSOCK')
	registerDoParallel(cl)
}
################ Tuning
# A tune grid has to be generated and passed to the tuning algorithm
if (is.null(tuning_parameters) || is.null(tuning_parameters$ncomp)) {
	tune_grid <- NULL
} else {
	tune_grid <- expand.grid(tuning_parameters)
}
############ Fit the model (and tune it)
# Fit the model with the features and tune it with cross-validation
train_control_pls <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv)
if (!is.null(seed)) {
	# Make the randomness reproducible
	set.seed(seed)
}
pls_model <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="pls", trControl=train_control_pls, preProcess=preprocessing, metric=selection_criteria, maximize=maximise_selection_criteria_values, tuneGrid=tune_grid)
# Output the parameters
pls_performances <- pls_model$results
# Plots
plots <- list()
plot(pls_model)
#plot_pls_accuracy <- recordPlot()
#plots <- append(plots, plot_pls_accuracy)
######. Close the cluster for parallel computation
if (Sys.info()[1] == "Windows") {
	stopCluster(cl)
}
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	predicted_classes_pls <- predict(pls_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
	# Create the outcomes dataframe
	classification_results_pls <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_pls, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_pls <- confusionMatrix(data = predicted_classes_pls, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	pls_roc <- list()
	roc_curve <- roc(response=classification_results_pls$true, predictor=as.numeric(classification_results_pls$predicted))
	pls_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	pls_roc[[2]] <- recordPlot()
	###################### Pie chart classification
	correctly_classified <- 0
	misclassified <- 0
	for (i in 1:nrow(classification_results_pls)) {
		if (as.character(classification_results_pls$predicted[i]) == as.character(classification_results_pls$true[i])) {
			correctly_classified <- correctly_classified + 1
		} else {
			misclassified <- misclassified + 1
		}
	}
	classification_pie <- c(correctly_classified, misclassified)
	pie(x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
	pie_chart_classification <- recordPlot()
	plots <- append(plots, pie_chart_classification)
}
# Output the results
return (list(model=pls_model, classification_results=classification_results_pls, performances_pls_cv=pls_performances, pls_performances=test_performances_pls, roc=pls_roc, plots=plots, pls_features=pls_model$finalModel$xNames))
}





################################################################################





################################### NAIVE BAYES CLASSIFIER TUNING AND VALIDATION
# This function operates the training of a Naive Bayes Classifier (NBC) model, by testing all the parameters of the NBC and selecting the best (in terms of accuracy). The tuning is performed via cross-validation onto the training dataset.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
nbc_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), tuning_parameters=data.frame(fL=NULL,usekernel=NULL), k_fold_cv=10, repeats_cv=2, positive_class_cv="HP", seed=NULL, preprocessing=c("center","scale"), selection_criteria="Accuracy", maximise_selection_criteria_values=TRUE) {
# Load the required libraries
install_and_load_required_packages(c("caret", "e1071", "klaR", "MASS"))
### PARALLEL BACKEND
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE) - 1
if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
	install_and_load_required_packages("doMC")
	# Register the foreach backend
	registerDoMC(cores = cpu_thread_number)
} else if (Sys.info()[1] == "Windows") {
	install_and_load_required_packages("doParallel")
	# Register the foreach backend
	cl <- makeCluster(cpu_thread_number, type='PSOCK')
	registerDoParallel(cl)
}
# Fit the model with the features and tune it with cross-validation
train_control_nbc <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv)
################ Tuning
# A tune grid has to be generated and passed to the tuning algorithm
if (is.null(tuning_parameters) || is.null(tuning_parameters$fL) && is.null(tuning_parameters$usekernel)) {
	tune_grid <- NULL
} else {
	tune_grid <- expand.grid(tuning_parameters)
}
############ Fit the model (and tune it)
# Make the randomness reproducible
if (!is.null(seed)) {
	set.seed(seed)
}
nbc_model <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="nb", trControl=train_control_nbc, preProcess=preprocessing, metric=selection_criteria, maximize=maximise_selection_criteria_values, tuneGrid=tune_grid)
# Output the parameters
nbc_performances <- nbc_model$results
# Plots
plots <- list()
#plot(nbc_model)
#plot_nbc_accuracy <- recordPlot()
#plots <- append(plots, plot_nbc_accuracy)
######. Close the cluster for parallel computation
if (Sys.info()[1] == "Windows") {
	stopCluster(cl)
}
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	# Plant the seed only if a specified value is entered
	if (!is.null(seed)) {
		# Make the randomness reproducible
		set.seed(seed)
	}
	predicted_classes_nbc <- predict(nbc_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
	# Create the outcomes dataframe
	classification_results_nbc <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_nbc, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_nbc <- confusionMatrix(data = predicted_classes_nbc, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	nbc_roc <- list()
	roc_curve <- roc(response=classification_results_nbc$true, predictor=as.numeric(classification_results_nbc$predicted))
	nbc_roc[[1]] <- roc_curve$auc
	plot(roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	nbc_roc[[2]] <- recordPlot()
	###################### Pie chart classification
	correctly_classified <- 0
	misclassified <- 0
	for (i in 1:nrow(classification_results_nbc)) {
		if (as.character(classification_results_nbc$predicted[i]) == as.character(classification_results_nbc$true[i])) {
			correctly_classified <- correctly_classified + 1
		} else {
			misclassified <- misclassified + 1
		}
	}
	classification_pie <- c(correctly_classified, misclassified)
	pie(x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
	pie_chart_classification <- recordPlot()
	plots <- append(plots, pie_chart_classification)
}
# Output the results
return (list(model=nbc_model, classification_results=classification_results_nbc, performances_nbc_cv=nbc_performances, nbc_performances=test_performances_nbc, roc=nbc_roc, plots=plots, nbc_features=nbc_model$finalModel$xNames))
}





################################################################################





###################################################### KNN TUNING AND VALIDATION
# This function operates the tuning of the k-Nearest Neighbour (KNN), by testing all the parameters of the KNN (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
knn_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(k=seq(1,15, by=1)), k_fold_cv=10, repeats_cv=2, preprocessing=c("center","scale"), parameters=list(k=15), positive_class_cv="HP", seed=NULL, evaluation_method="Accuracy") {
	# Load the required libraries
	install_and_load_required_packages(c("caret", "pROC"))
	### PARALLEL BACKEND
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		install_and_load_required_packages("doMC")
		# Register the foreach backend
		registerDoMC(cores = cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		install_and_load_required_packages("doParallel")
		# Register the foreach backend
		cl <- makeCluster(cpu_thread_number, type='PSOCK')
		registerDoParallel(cl)
	}
	################# FEATURES
	if (autotuning == TRUE) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train and tune the model
		# Define the tune grid for the model
		knn_tune_grid <- expand.grid(k=tuning_parameters$k)#, l=tuning_parameters$l)
		# Training and tuning
		knn_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="knn", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=knn_tune_grid)
		# Extract the cross-validation information
		cv_knn_model <- confusionMatrix(knn_tuning)
		# Select automatically the best model after the tuning
		knn_model <- knn_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(k=knn_model$k)
		# Plots
		tuning_plot <- plot(knn_tuning)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train the model
		# Define the tune grid for the model
		knn_tune_grid <- expand.grid(k=parameters$k)
		# Training and tuning
		knn_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="knn", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=knn_tune_grid)
		# Extract the cross-validation information
		cv_knn_model <- confusionMatrix(knn_tuning)
		# Select automatically the best model after the tuning
		knn_model <- knn_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(k=parameters$k)
	}
	######. Close the cluster for parallel computation
	if (Sys.info()[1] == "Windows") {
		stopCluster(cl)
	}
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_knn <- predict(knn_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type="class")
		# Create the outcomes dataframe
		classification_results_knn <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_knn, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_knn <- confusionMatrix(data = predicted_classes_knn, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		knn_roc <- list()
		roc_curve <- roc(response=as.numeric(classification_results_knn$true), predictor=as.numeric(classification_results_knn$predicted))
		knn_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
		legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		knn_roc[[2]] <- recordPlot()
		###################### Pie chart classification
		correctly_classified <- 0
		misclassified <- 0
		for (i in 1:nrow(classification_results_knn)) {
			if (as.character(classification_results_knn$predicted[i]) == as.character(classification_results_knn$true[i])) {
				correctly_classified <- correctly_classified + 1
			} else {
				misclassified <- misclassified + 1
			}
		}
		classification_pie <- c(correctly_classified, misclassified)
		pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
		pie_chart_classification <- recordPlot()
		# Output the results
		return (list(model=knn_model, knn_features=knn_model$xNames, cross_validation=cv_knn_model, classification_results=classification_results_knn, parameters=parameters_output, tuning_plot=tuning_plot, performances=test_performances_knn, roc=knn_roc, pie_chart_classification=pie_chart_classification))
	} else {return (list(model=knn_model, knn_features=knn_model$xNames, parameters=parameters_output, tuning_plot=tuning_plot, cross_validation=cv_knn_model))}
}





################################################################################





###################################################### LDA TUNING AND VALIDATION
# This function operates the tuning of the Linear Discriminant Analysis (LDA), by testing all the parameters of the LDA (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
lda_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(dimen=seq(1,15, by=1)), k_fold_cv=10, repeats_cv=2, preprocessing=c("center","scale"), parameters=list(dimen=15), positive_class_cv="HP", seed=NULL, evaluation_method="Accuracy") {
	# Load the required libraries
	install_and_load_required_packages(c("caret", "pROC"))
	### PARALLEL BACKEND
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		install_and_load_required_packages("doMC")
		# Register the foreach backend
		registerDoMC(cores = cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		install_and_load_required_packages("doParallel")
		# Register the foreach backend
		cl <- makeCluster(cpu_thread_number, type='PSOCK')
		registerDoParallel(cl)
	}
	################# FEATURES
	if (autotuning == TRUE) {
		# Find the best tuning parameters for the LDA
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train and tune the model
		# Define the tune grid for the model
		lda_tune_grid <- expand.grid(dimen=tuning_parameters$dimen)
		# Training and tuning
		lda_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="lda2", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=lda_tune_grid)
		# Extract the cross-validation information
		cv_lda_model <- confusionMatrix(lda_tuning)
		# Select automatically the best model after the tuning
		lda_model <- lda_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(dimen=lda_model$tuneValue$dimen)
		# Plots
		tuning_plot <- plot(lda_tuning)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train the model
		# Define the tune grid for the model
		lda_tune_grid <- expand.grid(dimen=parameters$dimen)
		# Training and tuning
		lda_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="lda2", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=lda_tune_grid)
		# Extract the cross-validation information
		cv_lda_model <- confusionMatrix(lda_tuning)
		# Select automatically the best model after the tuning
		lda_model <- lda_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(dimen=parameters$dimen)
	}
	######. Close the cluster for parallel computation
	if (Sys.info()[1] == "Windows") {
		stopCluster(cl)
	}
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_lda <- predict(lda_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type="class")$class
		# Create the outcomes dataframe
		classification_results_lda <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_lda, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_lda <- confusionMatrix(data = predicted_classes_lda, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		lda_roc <- list()
		roc_curve <- roc(response=as.numeric(classification_results_lda$true), predictor=as.numeric(classification_results_lda$predicted))
		lda_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
		legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		lda_roc[[2]] <- recordPlot()
		###################### Pie chart classification
		correctly_classified <- 0
		misclassified <- 0
		for (i in 1:nrow(classification_results_lda)) {
			if (as.character(classification_results_lda$predicted[i]) == as.character(classification_results_lda$true[i])) {
				correctly_classified <- correctly_classified + 1
			} else {
				misclassified <- misclassified + 1
			}
		}
		classification_pie <- c(correctly_classified, misclassified)
		pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
		pie_chart_classification <- recordPlot()
		# Output the results
		return (list(model=lda_model, lda_features=lda_model$xNames, cross_validation=cv_lda_model, classification_results=classification_results_lda, parameters=parameters_output, tuning_plot=tuning_plot, performances=test_performances_lda, roc=lda_roc, pie_chart_classification=pie_chart_classification))
	} else {return (list(model=lda_model, lda_features=lda_model$xNames, parameters=parameters_output, tuning_plot=tuning_plot, cross_validation=cv_lda_model))}
}





################################################################################





###################################################### RF TUNING AND VALIDATION
# This function operates the tuning of the Random Forest (RF), by testing all the parameters of the RF (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
rf_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(mtry=seq(1,5, by=1)), k_fold_cv=10, repeats_cv=2, preprocessing=c("center","scale"), parameters=list(mtry=round(sqrt(ncol(peaklist_training)))), positive_class_cv="HP", seed=NULL, evaluation_method="Accuracy") {
	# Load the required libraries
	install_and_load_required_packages(c("caret", "pROC","randomForest"))
	### PARALLEL BACKEND
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		install_and_load_required_packages("doMC")
		# Register the foreach backend
		registerDoMC(cores = cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		install_and_load_required_packages("doParallel")
		# Register the foreach backend
		cl <- makeCluster(cpu_thread_number, type='PSOCK')
		registerDoParallel(cl)
	}
	################# FEATURES
	if (autotuning == TRUE) {
		# Find the best tuning parameters for the RF
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train and tune the model
		# Define the tune grid for the model
		rf_tune_grid <- expand.grid(mtry=tuning_parameters$mtry)
		# Training and tuning
		rf_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="rf", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=rf_tune_grid)
		# Extract the cross-validation information
		cv_rf_model <- confusionMatrix(rf_tuning)
		# Select automatically the best model after the tuning
		rf_model <- rf_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(mtry=rf_model$mtry)
		# Plots
		tuning_plot <- plot(rf_tuning)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train the model
		# Define the tune grid for the model
		rf_tune_grid <- expand.grid(mtry=parameters$mtry)
		# Training and tuning
		rf_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="rf", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=rf_tune_grid)
		# Extract the cross-validation information
		cv_rf_model <- confusionMatrix(rf_tuning)
		# Select automatically the best model after the tuning
		rf_model <- rf_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(mtry=parameters$mtry)
	}
	######. Close the cluster for parallel computation
	if (Sys.info()[1] == "Windows") {
		stopCluster(cl)
	}
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_rf <- predict(rf_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type="class")
		# Create the outcomes dataframe
		classification_results_rf <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_rf, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_rf <- confusionMatrix(data = predicted_classes_rf, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		rf_roc <- list()
		roc_curve <- roc(response=as.numeric(classification_results_rf$true), predictor=as.numeric(classification_results_rf$predicted))
		rf_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
		legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		rf_roc[[2]] <- recordPlot()
		###################### Pie chart classification
		correctly_classified <- 0
		misclassified <- 0
		for (i in 1:nrow(classification_results_rf)) {
			if (as.character(classification_results_rf$predicted[i]) == as.character(classification_results_rf$true[i])) {
				correctly_classified <- correctly_classified + 1
			} else {
				misclassified <- misclassified + 1
			}
		}
		classification_pie <- c(correctly_classified, misclassified)
		pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
		pie_chart_classification <- recordPlot()
		# Output the results
		return (list(model=rf_model, rf_features=rf_model$xNames, cross_validation=cv_rf_model, classification_results=classification_results_rf, parameters=parameters_output, tuning_plot=tuning_plot, performances=test_performances_rf, roc=rf_roc, pie_chart_classification=pie_chart_classification))
	} else {return (list(model=rf_model, rf_features=rf_model$xNames, parameters=parameters_output, tuning_plot=tuning_plot, cross_validation=cv_rf_model))}
}





################################################################################





###################################################### NN TUNING AND VALIDATION
# This function operates the tuning of the Neural Network (NN), by testing all the parameters of the NN (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
nn_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), autotuning=TRUE, tuning_parameters=list(size=seq(1,5, by=1),decay=seq(0,2,by=1)), k_fold_cv=10, repeats_cv=2, preprocessing=c("center","scale"), parameters=list(size=5, decay=0), positive_class_cv="HP", seed=NULL, evaluation_method="Accuracy") {
	# Load the required libraries
	install_and_load_required_packages(c("caret", "pROC","nnet"))
	### PARALLEL BACKEND
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		install_and_load_required_packages("doMC")
		# Register the foreach backend
		registerDoMC(cores = cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		install_and_load_required_packages("doParallel")
		# Register the foreach backend
		cl <- makeCluster(cpu_thread_number, type='PSOCK')
		registerDoParallel(cl)
	}
	################# FEATURES
	if (autotuning == TRUE) {
		# Find the best tuning parameters for the NN
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train and tune the model
		# Define the tune grid for the model
		nn_tune_grid <- expand.grid(size=tuning_parameters$size, decay=tuning_parameters$decay)
		# Training and tuning
		nn_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="nnet", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=nn_tune_grid)
		# Extract the cross-validation information
		cv_nn_model <- confusionMatrix(nn_tuning)
		# Select automatically the best model after the tuning
		nn_model <- nn_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(size=nn_model$tuneValue$size, decay=nn_model$decay)
		# Plots
		tuning_plot <- plot(nn_tuning)
	}
	#################### SVM with defined parameters
	if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
		# Find the best tuning parameters for the SVM
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		# Define the control function for training
		training_ctrl <- trainControl(method="repeatedcv", number=k_fold_cv, repeats=repeats_cv, classProbs=TRUE)#, seeds=seeds) #, summaryFunction=twoClassSummary)
		# Train the model
		# Define the tune grid for the model
		nn_tune_grid <- expand.grid(size=parameters$size, decay=parameters$decay)
		# Training and tuning
		nn_tuning <- train(x=peaklist_training [,!(names(peaklist_training) %in% non_features)], y=factor(peaklist_training$Class), method="nnet", preProcess = preprocessing, metric=evaluation_method, trControl=training_ctrl, tuneGrid=nn_tune_grid)
		# Extract the cross-validation information
		cv_nn_model <- confusionMatrix(nn_tuning)
		# Select automatically the best model after the tuning
		nn_model <- nn_tuning$finalModel
		# Parameters output
		parameters_output <- data.frame(size=parameters$size, decay=parameters$decay)
	}
	######. Close the cluster for parallel computation
	if (Sys.info()[1] == "Windows") {
		stopCluster(cl)
	}
	#################### EXTERNAL VALIDATION (If a dataset is provided)
	if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
		#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
		# Plant the seed only if a specified value is entered
		if (!is.null(seed)) {
			# Make the randomness reproducible
			set.seed(seed)
		}
		predicted_classes_nn <- predict(nn_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type="class")
		# Create the outcomes dataframe
		classification_results_nn <- data.frame(sample=peaklist_test$Sample, predicted=predicted_classes_nn, true=peaklist_test$Class)
		### Generate the confusion matrix to evaluate the performances
		test_performances_nn <- confusionMatrix(data = predicted_classes_nn, peaklist_test$Class, positive=positive_class_cv)
		#### ROC analysis
		nn_roc <- list()
		roc_curve <- roc(response=as.numeric(classification_results_nn$true), predictor=as.numeric(classification_results_nn$predicted))
		nn_roc[[1]] <- roc_curve$auc
		plot (roc_curve)
		roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
		legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
		nn_roc[[2]] <- recordPlot()
		###################### Pie chart classification
		correctly_classified <- 0
		misclassified <- 0
		for (i in 1:nrow(classification_results_nn)) {
			if (as.character(classification_results_nn$predicted[i]) == as.character(classification_results_nn$true[i])) {
				correctly_classified <- correctly_classified + 1
			} else {
				misclassified <- misclassified + 1
			}
		}
		classification_pie <- c(correctly_classified, misclassified)
		pie (x=classification_pie, labels=c("Correctly classified", "Misclassified"), col=c("green","blue"))
		pie_chart_classification <- recordPlot()
		# Output the results
		return (list(model=nn_model, nn_features=nn_model$xNames, cross_validation=cv_nn_model, classification_results=classification_results_nn, parameters=parameters_output, tuning_plot=tuning_plot, performances=test_performances_nn, roc=nn_roc, pie_chart_classification=pie_chart_classification))
	} else {return (list(model=nn_model, nn_features=nn_model$xNames, parameters=parameters_output, tuning_plot=tuning_plot, cross_validation=cv_nn_model))}
}





################################################################################






################################################# ENSEMBLE TUNING AND VALIDATION
# This function operates the tuning of the ensemble classifier, by relying upon other functions to train, tune and validate the individual classifiers.
# It returns the best models in terms of classification performances, along with their parameters and performances (cross-validation or external validation, according to if an external dataset is provided).
ensemble_tuning_and_validation <- function(peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), classifiers=c("svm","pls","bayes"), autotuning=TRUE,  classifier_parameters=list(svm=list(gamma=10^(-5:5), cost=10^(-5:5), epsilon=seq(1,2,by=1), degree=1:5, kernel="radial"), pls=data.frame(ncomp=1:5), bayes=NULL), k_fold_cv=10, repeats_cv=2, positive_class_cv="HP", seed=NULL, preprocessing=c("center","scale"), selection_criteria=data.frame(pls="Accuracy",bayes="Accuracy"), maximise_selection_criteria_values=data.frame(pls=TRUE,bayes=TRUE)) {
	if ("svm" %in% classifiers || "SVM" %in% classifiers) {
		svm_classifier <- svm_tuning_and_validation(peaklist_training=peaklist_training, peaklist_test=peaklist_test, non_features=non_features, autotuning=autotuning, tuning_parameters=classifier_parameters$svm, k_fold_cv=k_fold_cv, repeats_cv=repeats_cv, parameters=classifier_parameters$svm, positive_class_cv=positive_class_cv, seed=seed)
	}
	if ("pls" %in% classifiers || "PLS" %in% classifiers) {
		pls_classifier <- pls_tuning_and_validation(peaklist_training=peaklist_training, peaklist_test=peaklist_test, non_features=non_features, tuning_parameters=classifier_parameters$pls, k_fold_cv=k_fold_cv, repeats_cv=repeats_cv, positive_class_cv="HP", seed=seed, preprocessing=preprocessing, selection_criteria=as.character(selection_criteria$pls), maximise_selection_criteria_values=maximise_selection_criteria_values$pls)
	}
	if ("bayes" %in% classifiers || "Bayes" %in% classifiers) {
		bayes_classifier <- nbc_tuning_and_validation(peaklist_training=peaklist_training, peaklist_test=peaklist_test, non_features=non_features, tuning_parameters=classifier_parameters$bayes, k_fold_cv=k_fold_cv, repeats_cv=repeats_cv, positive_class_cv=positive_class_cv, seed=seed, preprocessing=preprocessing, selection_criteria=as.character(selection_criteria$bayes), maximise_selection_criteria_values=maximise_selection_criteria_values$bayes)
	}
	return(list(svm_classifier=svm_classifier, pls_classifier=pls_classifier, bayes_classifier=bayes_classifier))
}





################################################ ROUND NUMERIC FEATURES PEAKLIST
# This function rounds the numeric features to a certain decimal digit.
# It returns the same matrix with the rounded features, along with the original peaklist matrix.
round_features_peaklist <- function(peaklist, decimal_digits=3, non_features=c("Sample","Class","THY")) {
# Take the features
feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Round the features
for (i in 1:length(feature_list)) {
	feature_list[i] <- round(feature_list[i], digits=decimal_digits)
}
# Regather the features
features <- c(feature_list, non_features)
peaklist2 <- peaklist
names(peaklist2) <- features
#
return (list(original_peaklist=peaklist, peaklist_rounded=peaklist2))
}





################################################################################





#################### SPECTRAL TYPER SCORE ACCORDING TO THE HIERARCHICAL DISTANCE
# This function computes the Spectral Typer score by comparing the test spectra with the library spectra, determining the similarity (through the euclidean distance) and assigning a category according to the distance.
# Each sample gets compared with all the entries in the database, simultaneously.
spectral_typer_score_hierarchical_distance <- function (spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=NULL, peaks_filtering=TRUE, peaks_filtering_percentage_threshold=25, low_intensity_peaks_removal=FALSE, low_intensity_percentage_threshold=0.1, low_intensity_threshold_method="element-wise", tolerance_ppm=2000, spectra_path_output=TRUE, score_only=TRUE, spectra_format="brukerflex", normalise_distances=TRUE, normalisation_method="sum") {
# Load the required libraries
install_and_load_required_packages(c("MALDIquant", "stats", "ggplot2", "ggdendro"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Sample and Library size
if (isMassPeaksList(peaks_test)) {
    number_of_samples <- length(peaks_test)
} else if (isMassPeaks(peaks_test)) {
    number_of_samples <- 1
}
if (isMassPeaksList(peaks_database)) {
    database_size <- length(peaks_database)
} else if (isMassPeaks(peaks_database)) {
    database_size <- 1
}
####### Peak alignment
# Merge the peaklists and the spectra
peaks_all <- append(peaks_database, peaks_test)
spectra_all <- append(spectra_database, spectra_test)
# Align
peaks_all <- align_and_filter_peaks(peaks_all, tof_mode=tof_mode, peaks_filtering=peaks_filtering, frequency_threshold_percent=peaks_filtering_percentage_threshold, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=low_intensity_percentage_threshold, intensity_threshold_method=low_intensity_threshold_method)
# Restore the lists
peaks_database <- peaks_all [1:database_size]
peaks_test <- peaks_all [(database_size+1):length(peaks_all)]
#### Replace the sample name, both in the library and in the test set
peaks_test <- replace_sample_name(peaks_test, spectra_format=spectra_format)
peaks_database <- replace_class_name(peaks_database,  class_list=class_list_library, spectra_format=spectra_format)
####### Create the sample vector
sample_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	# If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
####### Create the library vector
database_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
	# If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
# Generate the path vector
spectra_path_vector <- character()
for (spectrum in spectra_test) {
	spectra_path_vector <- append(spectra_path_vector, spectrum@metaData$file[1])
}
# Generate the matrix (for hca)
peaklist_matrix <- intensityMatrix(peaks_all, spectra_all)
# Add additional info to the matrix
peaklist_matrix <- matrix_add_class_and_sample(peaklist_matrix, peaks=peaks_all, class_list=list(), spectra_format=spectra_format, sample_output=TRUE, class_output=FALSE)
rownames(peaklist_matrix) <- make.names(peaklist_matrix[,"Sample"], unique=TRUE)
# Compute the hca
distance_matrix <- dist(peaklist_matrix[,1:(ncol(peaklist_matrix)-1)], method="euclidean")
hierarchical_clustering <- hclust(distance_matrix)
#plot(hierarchical_clustering, main="Hierarchical clustering analysis - Spectral Typer, xlab="Samples", ylab="Tree height")
hca_dendrogram <- ggdendrogram(hierarchical_clustering, segments=TRUE, labels=TRUE, leaf_labels=TRUE, rotate=TRUE, theme_dendro=TRUE)#, main="Hierarchical clustering analysis - Spectral Typer", xlab="Samples", ylab="Tree height")
#hca_dendrogram <- recordPlot()
#
distance_matrix <- as.matrix(distance_matrix)
# The distance matrix displays the distance between the spectra
colnames(distance_matrix) <- peaklist_matrix[,"Sample"]
rownames(distance_matrix) <- peaklist_matrix[,"Sample"]
# Remove the first rows (the spectra from the database)
distance_matrix <- distance_matrix[(database_size+1):nrow(distance_matrix),]
# Keep only the first columns (the spectra from the database)
distance_matrix <- distance_matrix[,1:database_size]
### Normalise the euclidean distances
if (normalise_distances == TRUE) {
	# TIC (SUM)
	if (normalisation_method == "sum") {
		# Compute the sum of the rows
		row_sums <- apply(distance_matrix, MARGIN=1, FUN=sum)
		# Divide each element of the matrix by the sum of the row
		for (r in 1:nrow(distance_matrix)) {
			distance_matrix [r,] <- distance_matrix[r,] / row_sums[r]
		}
		# Multiply everything by 10, to have more readable results
		distance_matrix <- distance_matrix * 10
		# The classification is made by comparing the single sample spectrum with the spectrum of the database class (the distance is displayed in the distance matrix): the closer the better
		# Scroll the rows, assign the class based upon the distance, create the output matrix for results (create a function to apply to each matrix row)
		scoring_function <- function (x) {
			if (x < 1) {
		    	x <- paste("YES (", round(as.numeric(x),3), ")", sep="")
			} else if (x >= 1 && x < 1.2) {
				x <- paste("NI (", round(as.numeric(x),3), ")", sep="")
			} else if (x >= 1.2) {
		        x <- paste("NO (", round(as.numeric(x),3), ")", sep="")
		    }
		    return (x)
		}
		result_matrix <- apply(distance_matrix, MARGIN=c(1,2), FUN=function(x) scoring_function(x))
	}
	# MAX
	if (normalisation_method == "max") {
		# Divide each element of the matrix by the maximum of the row
		for (r in 1:nrow(distance_matrix)) {
			distance_matrix [r,] <- distance_matrix[r,] / max(distance_matrix[r,])
		}
		# Multiply everything by 100, to have more readable results (percentage of the max)
		distance_matrix <- distance_matrix * 100
		# The classification is made by comparing the single sample spectrum with the spectrum of the database class (the distance is displayed in the distance matrix): the closer the better
		# Scroll the rows, assign the class based upon the distance, create the output matrix for results (create a function to apply to each matrix row)
		scoring_function <- function (x) {
			if (x < 50) {
		    	x <- paste("YES (", round(as.numeric(x),3), ")", sep="")
			} else if (x >= 50 && x < 75) {
				x <- paste("NI (", round(as.numeric(x),3), ")", sep="")
			} else if (x >= 75) {
		        x <- paste("NO (", round(as.numeric(x),3), ")", sep="")
		    }
		    return (x)
		}
		result_matrix <- apply(distance_matrix, MARGIN=c(1,2), FUN=function(x) scoring_function(x))
	}
}
# Spectra path
if (spectra_path_output == TRUE) {
	result_matrix <- cbind(result_matrix, sample_vector)
}
return (list(result_matrix=result_matrix, hca_dendrogram=hca_dendrogram))
}





################################################################################





####################################### SPECTRAL TYPER SCORE: CORRELATION MATRIX
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity symmetry via the correlation matrix.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_correlation_matrix <- function(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library=NULL, peaks_filtering=TRUE, peaks_filtering_percentage_threshold=25, low_intensity_peaks_removal=FALSE, low_intensity_percentage_threshold=0.1, low_intensity_threshold_method="element-wise", tolerance_ppm=2000, intensity_correction_coefficient=1, spectra_format="brukerflex", spectra_path_output=TRUE, score_only=FALSE, multicore_processing=TRUE) {
install_and_load_required_packages(c("MALDIquant", "corrplot", "weights", "stats", "parallel"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
# Rename the trim function to avoid conflicts
trim_weights <- get(x="trim", pos="package:weights")
### Folder lists
#database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
#test_folder_list <- dir(filepath_test, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
# Sample and Library size
if (isMassPeaksList(peaks_test)) {
    number_of_samples <- length(peaks_test)
} else if (isMassPeaks(peaks_test)) {
    number_of_samples <- 1
}
if (isMassPeaksList(peaks_database)) {
    database_size <- length(peaks_database)
} else if (isMassPeaks(peaks_database)) {
    database_size <- 1
}
#### Replace the sample name, both in the library and in the test set
peaks_test <- replace_sample_name(peaks_test, spectra_format=spectra_format)
peaks_database <- replace_class_name(peaks_database, class_list=class_list_library, class_in_file_name=TRUE, spectra_format=spectra_format)
####### Create the sample vector
sample_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
####### Create the library vector
database_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
# Generate the path vector
spectra_path_vector <- character()
for (sp in 1:length(spectra_test)) {
	spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
}
############################################################ SCORE (FRI)
number_of_signals_database <- numeric(length=database_size)
for (d in 1:database_size) {
	number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
}
################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
# Create a list to be used for lapply
global_list <- list()
for (spl in 1:number_of_samples) {
	# Extract the peaklist and the spectrum
	peaks_database_temp <- peaks_database
	spectra_database_temp <- spectra_database
	peaks_sample <- peaks_test[[spl]]
	spectrum_sample <- spectra_test[[spl]]
	# Generate the entry of the global list
	global_list_entry <- append(peaks_sample, peaks_database_temp)
	global_list_entry <- append(global_list_entry, spectrum_sample)
	global_list_entry <- append(global_list_entry, spectra_database_temp)
	global_list[[spl]] <- global_list_entry
}
############################################## Define the function for parLapply
# x = each element of the global list
comparison_sample_db_subfunction_correlation <- function(x) {
	# Generate the matrix rows for the output
	matching_signals_matrix <- matrix(0, nrow=1, ncol=database_size)
	number_of_signals_database_matrix <- matrix(0, nrow=1, ncol=database_size)
	fit_matrix <- matrix(0, ncol=database_size, nrow=1)
	retrofit_matrix <- matrix(0, ncol=database_size, nrow=1)
	intensity_correlation_matrix <- matrix(0, ncol=database_size, nrow=1)
	pvalue_matrix <- matrix(0, ncol=database_size, nrow=1)
	slope_matrix <- matrix(0, ncol=database_size, nrow=1)
	###### Compare with all the elements in the library
	### For each entry in the library...
	for (db in 1:database_size) {
		# Extract the peaklist and the spectrum
		peaks_sample <- x[[1]]
		spectrum_sample <- x[[(database_size+2)]]
		peaks_database_temp_all <- x[2:(database_size+1)]
		spectra_database_temp_all <- x[(database_size+3):length(x)]
		peaks_database_temp <- peaks_database_temp_all[[db]]
		spectra_database_temp <- spectra_database_temp_all[[db]]
		####### Peak alignment
		# Merge the peaklists
		peaks_all <- append(peaks_database_temp, peaks_sample)
		spectra_all <- append(spectra_database_temp, spectrum_sample)
		# Align the peaks
		peaks_all <- align_and_filter_peaks(peaks_all, tof_mode=tof_mode, peaks_filtering=peaks_filtering, frequency_threshold_percent=peaks_filtering_percentage_threshold, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=low_intensity_percentage_threshold, intensity_threshold_method=low_intensity_threshold_method)
		# Restore the lists
		peaks_database_temp <- peaks_all[[1]]
		peaks_sample <- peaks_all[[2]]
		#################### Number of signals
		number_of_signals_samples <- length(peaks_sample@mass)
		number_of_signals_database <- length(peaks_database_temp@mass)
		number_of_signals_database_matrix[1,db] <- number_of_signals_database
		###### COUNTER 0 - MATCHING SIGNALS
		# Create a counter, symmetrical to the database Peaklist
		matching_signals_number <- 0
		matching_signals <- numeric()
		# For each peaklist in the Library
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_sample@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_database_temp@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
						matching_signals_number <- matching_signals_number + 1
						matching_signals <- append(matching_signals, peaks_sample@mass[i])
						break
					}
				}
			}
		} else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
				matching_signals_number <- 0
		}
		# Append this row to the global matrix
		matching_signals_matrix [1,db] <- matching_signals_number
		###### COUNTER 1 - FIT
		fit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			fit_sample <- matching_signals_number / length(peaks_sample@mass)
		}
		# Append this row to the global matrix
		fit_matrix[1,db] <- fit_sample
		###### COUNTER 2 - RETRO FIT
		retrofit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			retrofit_sample <- matching_signals_number / length(peaks_database_temp@mass)
		}
		# Append this row to the global matrix
		retrofit_matrix[1,db] <- retrofit_sample
		###### COUNTER 3
		# Symmetry -> comparison between intensities
		# Compute the correlation matrix with the library
		# Intensity matrix
		intensity_matrix_global <- intensityMatrix(peaks_all, spectra_all)
		# Keep only the matching signals
		#columns_to_keep <- as.character(matching_signals)
		#intensity_matrix_global <- intensity_matrix_global[,columns_to_keep]
		# Weighted correlation between samples (library + test samples) (samples must be as columns and features as test) - With weights
		if (intensity_correction_coefficient != 0 && intensity_correction_coefficient != 1) {
			# Compute the vector of weights
			weights_vector <- c(rep(1, length(database_vector)), rep(intensity_correction_coefficient, nrow(t(intensity_matrix_global))))
			correlation_sample <- wtd.cors(x=t(intensity_matrix_global), weight=weights_vector)
			intensity_correlation_sample <- as.matrix(correlation_matrix [(database_size+1):nrow(correlation_matrix), 1:database_size])
		} else if (intensity_correction_coefficient == 1) {
			t_intensity_matrix_global <- t(intensity_matrix_global)
			correlation_sample <- cor.test(t_intensity_matrix_global[,1], t_intensity_matrix_global[,2], method="pearson")
			intensity_correlation_sample <- correlation_sample$estimate
			# pvalue
			pvalue <- correlation_sample$p.value
			pvalue_replacement_function <- function(x, number_of_digits) {
				if (is.na(x)) {
					x <- "Not available"
				} else if (x < 0.00001) {
					x <- "< 0.00001"
				} else {
					x <- as.character(round(x, digits=number_of_digits))
				}
				return (x)
			}
			pvalue <- pvalue_replacement_function(pvalue, number_of_digits=6)
			# Append this row to the global matrix
			pvalue_matrix[1,db] <- pvalue
		} else if (intensity_correction_coefficient == 0) {
			intensity_correlation_sample <- 1
		}
		# Extract the absolute values and fix the NAs
		intensity_correlation_sample <- abs(intensity_correlation_sample)
		if (is.na(intensity_correlation_sample)) {
			intensity_correlation_sample <- 0
		}
		# Append this row to the global matrix
		intensity_correlation_matrix[1,db] <- intensity_correlation_sample
		###### COUNTER 4 - REGRESSION CURVE
		t_intensity_matrix_global <- t(intensity_matrix_global)
		t_intensity_matrix_database <- rbind(as.matrix(t_intensity_matrix_global[,1]))
		t_intensity_matrix_test <- rbind(as.matrix(t_intensity_matrix_global[,2]))
		linear_regression <- lm(t_intensity_matrix_database[,1] ~ t_intensity_matrix_test[,1])
		regression_slope <- linear_regression$coefficients[2]
		regression_intercept <- linear_regression$coefficients[1]
		slope_sample <- round(regression_slope, digits=3)
		# Append this row to the global matrix
		slope_matrix[1,db] <- slope_sample
	}
	# Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
	return(list(number_of_signals_samples=number_of_signals_samples, number_of_signals_database_matrix=number_of_signals_database_matrix, matching_signals_matrix=matching_signals_matrix, fit_matrix=fit_matrix, retrofit_matrix=retrofit_matrix, intensity_correlation_matrix=intensity_correlation_matrix, pvalue_matrix=pvalue_matrix, slope_matrix=slope_matrix))
}
if (multicore_processing == TRUE) {
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		output_list <- mclapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_correlation(global_list), mc.cores=cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		# Make the CPU cluster for parallelisation
		cls <- makeCluster(cpu_thread_number)
		# Make the cluster use the custom functions and the package functions along with their parameters
		clusterEvalQ(cls, {library(MALDIquant)})
		# Pass the variables to the cluster for running the function
		clusterExport(cl=cls, varlist=c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_correlation", "database_size", "tof_mode", "peaks_filtering", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_peaks_removal", "low_intensity_percentage_threshold", "low_intensity_threshold_method", "intensity_correction_coefficient", "correlation_pvalue"), envir=environment())
		output_list <- parLapply(cls, global_list, fun=function(global_list) comparison_sample_db_subfunction_correlation(global_list))
		stopCluster(cls)
	}
} else {
	output_list <- lapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_correlation(global_list))
}
############################ Merge the matrix pieces together
matching_signals_matrix_all <- NULL
number_of_signals_database_matrix_all <- NULL
fit_matrix_all <- NULL
retrofit_matrix_all <- NULL
intensity_correlation_matrix_all <- NULL
pvalue_matrix_all <- NULL
slope_matrix_all <- NULL
for (ns in 1:number_of_samples) {
	# Matching signals
	if (is.null(matching_signals_matrix_all)) {
		matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
	} else {
		matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
	}
	# Number of signals database
	if (is.null(number_of_signals_database_matrix_all)) {
		number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
	} else {
		number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
	}
	# Fit
	if (is.null(fit_matrix_all)) {
		fit_matrix_all <- output_list[[ns]]$fit_matrix
	} else {
		fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
	}
	# Retrofit
	if (is.null(retrofit_matrix_all)) {
		retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
	} else {
		retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
	}
	# Intensity correlation
	if (is.null(intensity_correlation_matrix_all)) {
		intensity_correlation_matrix_all <- output_list[[ns]]$intensity_correlation_matrix
	} else {
		intensity_correlation_matrix_all <- rbind(intensity_correlation_matrix_all, output_list[[ns]]$intensity_correlation_matrix)
	}
	# pvalue
	if (is.null(pvalue_matrix_all)) {
		pvalue_matrix_all <- output_list[[ns]]$pvalue_matrix
	} else {
		pvalue_matrix_all <- rbind(pvalue_matrix_all, output_list[[ns]]$pvalue_matrix)
	}
	# Slope
	if (is.null(slope_matrix_all)) {
		slope_matrix_all <- output_list[[ns]]$slope_matrix
	} else {
		slope_matrix_all <- rbind(slope_matrix_all, output_list[[ns]]$slope_matrix)
	}
}
######################################
# Fix the rownames and colnames
rownames(intensity_correlation_matrix_all) <- sample_vector
colnames(intensity_correlation_matrix_all) <- database_vector
################### Score calculation
if (intensity_correction_coefficient != 0) {
score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_correlation_matrix_all*1000)
} else {
	score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_correlation_matrix_all*100)
}
colnames(score) <- database_vector
rownames(score) <- sample_vector
#### Output the classification
output <- matrix ("", nrow=number_of_samples, ncol=database_size)
colnames(output) <- database_vector
rownames(output) <- sample_vector
if (spectra_path_output == TRUE) {
	output <- cbind(output, spectra_path_vector)
    colnames(output) <- c(database_vector, "Spectrum path")
}
if (score_only == TRUE) {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
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
} else {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
			if (score[r,w] >= 2) {
				output[r,w] <- paste("YES","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits=3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
			if (score[r,w] < 1.5) {
				output[r,w] <- paste("NO","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits=3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
			if (score[r,w] >= 1.5 && score[r,w] < 2) {
				output[r,w] <- paste("NI","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits=3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
		}
	}
}
return (output)
}





################################################################################





######################################### SPECTRAL TYPER SCORE: SIGNAL INTENSITY
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity comparison.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_signal_intensity <- function(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=NULL, comparison=c("intensity percentage", "standard deviation"), peaks_filtering=TRUE, peaks_filtering_percentage_threshold=25, low_intensity_peaks_removal=FALSE, low_intensity_percentage_threshold=0.1, low_intensity_threshold_method="element-wise", tolerance_ppm=2000, intensity_tolerance_percent_threshold=50, spectra_format="brukerflex", spectra_path_output=TRUE, score_only=TRUE, number_of_st_dev=1, multicore_processing=TRUE) {
# Load the required libraries
install_and_load_required_packages("MALDIquant","parallel")
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
### Folder lists
#database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
#test_folder_list <- dir(filepath_test, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
# Sample and Library size
if (isMassPeaksList(peaks_test)) {
    number_of_samples <- length(peaks_test)
} else if (isMassPeaks(peaks_test)) {
    number_of_samples <- 1
}
if (isMassPeaksList(peaks_database)) {
    database_size <- length(peaks_database)
} else if (isMassPeaks(peaks_database)) {
    database_size <- 1
}
#### Replace the sample name, both in the library and in the test set
peaks_test <- replace_sample_name(peaks_test, spectra_format=spectra_format)
peaks_database <- replace_class_name(peaks_database, class_list=class_list_library, class_in_file_name=TRUE, spectra_format=spectra_format)
####### Create the sample vector
sample_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
####### Create the library vector
database_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
# Generate the path vector
spectra_path_vector <- character()
for (sp in 1:length(spectra_test)) {
	spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
}
############################################################ SCORE (FRI)
number_of_signals_database <- numeric(length=database_size)
for (d in 1:database_size) {
	number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
}
################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
# Create a list to be used for lapply
global_list <- list()
for (spl in 1:number_of_samples) {
	# Extract the peaklist and the spectrum
	peaks_database_temp <- peaks_database
	spectra_database_temp <- spectra_database
	peaks_sample <- peaks_test[[spl]]
	spectrum_sample <- spectra_test[[spl]]
	# Generate the entry of the global list
	global_list_entry <- append(peaks_sample, peaks_database_temp)
	global_list_entry <- append(global_list_entry, spectrum_sample)
	global_list_entry <- append(global_list_entry, spectra_database_temp)
	global_list[[spl]] <- global_list_entry
}
############################################## Define the function for parLapply
# x = each element of the global list
comparison_sample_db_subfunction_intensity <- function(x) {
	# Generate the matrix rows for the output
	matching_signals_matrix <- matrix(0, nrow=1, ncol=database_size)
	number_of_signals_database_matrix <- matrix(0, nrow=1, ncol=database_size)
	fit_matrix <- matrix(0, ncol=database_size, nrow=1)
	retrofit_matrix <- matrix(0, ncol=database_size, nrow=1)
	intensity_matching_matrix <- matrix(0, ncol=database_size, nrow=1)
	###### Compare with all the elements in the library
	### For each entry in the database...
	for (db in 1:database_size) {
		# Extract the peaklist and the spectrum
		peaks_sample <- x[[1]]
		spectrum_sample <- x[[(database_size+2)]]
		peaks_database_temp_all <- x[2:(database_size+1)]
		spectra_database_temp_all <- x[(database_size+3):length(x)]
		peaks_database_temp <- peaks_database_temp_all[[db]]
		spectra_database_temp <- spectra_database_temp_all[[db]]
		####### Peak alignment
		# Merge the peaklists
		peaks_all <- append(peaks_database_temp, peaks_sample)
		spectra_all <- append(spectra_database_temp, spectrum_sample)
		# Align the peaks
		peaks_all <- align_and_filter_peaks(peaks_all, tof_mode=tof_mode, peaks_filtering=peaks_filtering, frequency_threshold_percent=peaks_filtering_percentage_threshold, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=low_intensity_percentage_threshold, intensity_threshold_method=low_intensity_threshold_method)
		# Restore the lists
		peaks_database_temp <- peaks_all[[1]]
		peaks_sample <- peaks_all[[2]]
		#################### Number of signals
		number_of_signals_samples <- length(peaks_sample@mass)
		number_of_signals_database <- length(peaks_database_temp@mass)
		number_of_signals_database_matrix[1,db] <- number_of_signals_database
		###### COUNTER 0 - MATCHING SIGNALS
		# Create a counter, symmetrical to the database Peaklist
		matching_signals_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_sample@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_database_temp@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
						matching_signals_sample <- matching_signals_sample + 1
						break
					}
				}
			}
		} else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
			matching_signals_sample <- 0
		}
		# Append this row to the global matrix
		matching_signals_matrix[1,db] <- matching_signals_sample
		###### COUNTER 1 - FIT
		fit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			fit_sample <- matching_signals_sample / length(peaks_sample@mass)
		}
		# Append this row to the global matrix
		fit_matrix[1,db] <- fit_sample
		###### COUNTER 2 - RETRO FIT
		retrofit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			retrofit_sample <- matching_signals_sample / length(peaks_database_temp@mass)
		}
		# Append this row to the global matrix
		retrofit_matrix[1,db] <- retrofit_sample
		###### COUNTER 3 (INTENSITY MATCHING)
		if ("intensity percentage" %in% comparison && !("standard deviation" %in% comparison)) {
			# Create a counter, symmetrical to the database Peaklist
			intensity_matching_sample <- 0
			if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
				# Scroll the peaks in the sample
				for (i in 1:length(peaks_sample@mass)) {
					# Scroll the peaklist in the library
					for (j in 1:length(peaks_database_temp@mass)) {
						# Count the number of closely matching signals with the reference in the sample peaklist
						if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
							# Evaluate the difference in intensity
							if ((abs(peaks_sample@intensity[i] - peaks_database_temp@intensity[j])*100/peaks_database_temp@intensity[j]) <= intensity_tolerance_percent_threshold) {
								intensity_matching_sample <- intensity_matching_sample + 1
							}
							break
						}
					}
				}
			} else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
				intensity_matching_sample <- 0
			}
			if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
				intensity_matching_sample <- intensity_matching_sample / matching_signals_sample
				if (is.na(intensity_matching_sample)) {
					intensity_matching_sample <- 0
				}
			}
			# Append this row to the global matrix
			intensity_matching_matrix[1,db] <- intensity_matching_sample
		} else if (!("intensity percentage" %in% comparison) && "standard deviation" %in% comparison) {
		    ########### To be implemented
			intensity_matching_sample <- NULL
		}
	}
	# Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
	return(list(number_of_signals_samples=number_of_signals_samples, number_of_signals_database_matrix=number_of_signals_database_matrix, matching_signals_matrix=matching_signals_matrix, fit_matrix=fit_matrix, retrofit_matrix=retrofit_matrix, intensity_matching_matrix=intensity_matching_matrix))
}
if (multicore_processing == TRUE) {
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		output_list <- mclapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_intensity(global_list), mc.cores=cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		# Make the CPU cluster for parallelisation
		cls <- makeCluster(cpu_thread_number)
		# Make the cluster use the custom functions and the package functions along with their parameters
		clusterEvalQ(cls, {library(MALDIquant)})
		clusterExport(cl=cls, varlist=c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_intensity", "database_size", "tof_mode", "peaks_filtering", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_peaks_removal", "low_intensity_percentage_threshold", "low_intensity_threshold_method", "intensity_correction_coefficient", "correlation_pvalue", "intensity_tolerance_percent_threshold", "comparison"), envir=environment())
		output_list <- parLapply(cls, global_list, fun=function(global_list) comparison_sample_db_subfunction_intensity(global_list))
		stopCluster(cls)
	}
} else {
	output_list <- lapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_intensity(global_list))
}
############################ Merge the matrix pieces together
matching_signals_matrix_all <- NULL
number_of_signals_database_matrix_all <- NULL
fit_matrix_all <- NULL
retrofit_matrix_all <- NULL
intensity_matching_matrix_all <- NULL
pvalue_matrix_all <- NULL
slope_matrix_all <- NULL
for (ns in 1:number_of_samples) {
	# Matching signals
	if (is.null(matching_signals_matrix_all)) {
		matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
	} else {
		matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
	}
	# Number of signals database
	if (is.null(number_of_signals_database_matrix_all)) {
		number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
	} else {
		number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
	}
	# Fit
	if (is.null(fit_matrix_all)) {
		fit_matrix_all <- output_list[[ns]]$fit_matrix
	} else {
		fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
	}
	# Retrofit
	if (is.null(retrofit_matrix_all)) {
		retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
	} else {
		retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
	}
	# Intensity correlation
	if (is.null(intensity_matching_matrix_all)) {
		intensity_matching_matrix_all <- output_list[[ns]]$intensity_matching_matrix
	} else {
		intensity_matching_matrix_all <- rbind(intensity_matching_matrix_all, output_list[[ns]]$intensity_matching_matrix)
	}
}
######################################
# Fix the rownames and colnames
colnames(intensity_matching_matrix_all) <- database_vector
rownames(intensity_matching_matrix_all) <- sample_vector
### Score calculation
score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_matching_matrix_all*1000)
colnames(score) <- database_vector
rownames(score) <- sample_vector
#### Output the classification
output <- matrix ("NO", nrow=number_of_samples, ncol=database_size)
colnames(output) <- database_vector
rownames(output) <- sample_vector
if (spectra_path_output == TRUE) {
	output <- cbind(output, spectra_path_vector)
    colnames(output) <- c(database_vector, "Spectrum path")
}
if (score_only == TRUE) {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
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
} else {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
			if (score[r,w]>=2) {
				output[r,w] <- paste("YES","(Score:", round(score[r,w], digits=3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits=3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
			if (score[r,w]<1.5) {
				output[r,w] <- paste("NO","(Score:", round(score[r,w], digits=3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits=3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
			if (score[r,w]>=1.5 && score[r,w]<2) {
				output[r,w] <- paste("NI","(Score:", round(score[r,w], digits=3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits=3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
			}
		}
	}
}
return(output)
}





######################################### SPECTRAL TYPER SCORE: SIMILARITY INDEX
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity symmetry via the similarity index computation.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_similarity_index <- function(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library=NULL, peaks_filtering=TRUE, peaks_filtering_percentage_threshold=25, low_intensity_peaks_removal=FALSE, low_intensity_percentage_threshold=0.1, low_intensity_threshold_method="element-wise", tolerance_ppm=2000, intensity_correction_coefficient=1, spectra_format="brukerflex", spectra_path_output=TRUE, score_only=FALSE, multicore_processing=TRUE) {
install_and_load_required_packages(c("MALDIquant", "stats", "parallel"))
# Rename the trim function
trim_spectra <- get(x="trim", pos="package:MALDIquant")
### Folder lists
#database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
#test_folder_list <- dir(filepath_test, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
# Sample and Library size
if (isMassPeaksList(peaks_test)) {
    number_of_samples <- length(peaks_test)
} else if (isMassPeaks(peaks_test)) {
    number_of_samples <- 1
}
if (isMassPeaksList(peaks_database)) {
    database_size <- length(peaks_database)
} else if (isMassPeaks(peaks_database)) {
    database_size <- 1
}
#### Replace the sample name, both in the library and in the test set
peaks_test <- replace_sample_name(peaks_test, spectra_format=spectra_format)
peaks_database <- replace_class_name(peaks_database, class_list=class_list_library, class_in_file_name=TRUE, spectra_format=spectra_format)
####### Create the sample vector
sample_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:number_of_samples) {
		sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
	}
}
####### Create the library vector
database_vector <- character()
if (spectra_format == "brukerflex" || spectra_format == "xmass") {
    # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
if (spectra_format == "imzml" | spectra_format == "imzML") {
	for (s in 1:database_size) {
		database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
	}
}
# Generate the path vector
spectra_path_vector <- character()
for (sp in 1:length(spectra_test)) {
	spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
}
############################################################ SCORE (FRI)
number_of_signals_database <- numeric(length=database_size)
for (d in 1:database_size) {
	number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
}
################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
# Create a list to be used for lapply
global_list <- list()
for (spl in 1:number_of_samples) {
	# Extract the peaklist and the spectrum
	peaks_database_temp <- peaks_database
	spectra_database_temp <- spectra_database
	peaks_sample <- peaks_test[[spl]]
	spectrum_sample <- spectra_test[[spl]]
	# Generate the entry of the global list
	global_list_entry <- append(peaks_sample, peaks_database_temp)
	global_list_entry <- append(global_list_entry, spectrum_sample)
	global_list_entry <- append(global_list_entry, spectra_database_temp)
	global_list[[spl]] <- global_list_entry
}
############################################## Define the function for parLapply
# x = each element of the global list
comparison_sample_db_subfunction_similarity_index <- function(x) {
	# Generate the matrix rows for the output
	matching_signals_matrix <- matrix(0, nrow=1, ncol=database_size)
	number_of_signals_database_matrix <- matrix(0, nrow=1, ncol=database_size)
	fit_matrix <- matrix(0, ncol=database_size, nrow=1)
	retrofit_matrix <- matrix(0, ncol=database_size, nrow=1)
	similarity_index_matrix <- matrix(0, ncol=database_size, nrow=1)
	###### Compare with all the elements in the library
	### For each entry in the library...
	for (db in 1:database_size) {
		# Extract the peaklist and the spectrum
		peaks_sample <- x[[1]]
		spectrum_sample <- x[[(database_size+2)]]
		peaks_database_temp_all <- x[2:(database_size+1)]
		spectra_database_temp_all <- x[(database_size+3):length(x)]
		peaks_database_temp <- peaks_database_temp_all[[db]]
		spectra_database_temp <- spectra_database_temp_all[[db]]
		####### Peak alignment
		# Merge the peaklists
		peaks_all <- append(peaks_database_temp, peaks_sample)
		spectra_all <- append(spectra_database_temp, spectrum_sample)
		# Align the peaks
		peaks_all <- align_and_filter_peaks(peaks_all, tof_mode=tof_mode, peaks_filtering=peaks_filtering, frequency_threshold_percent=peaks_filtering_percentage_threshold, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=low_intensity_percentage_threshold, intensity_threshold_method=low_intensity_threshold_method)
		# Restore the lists
		peaks_database_temp <- peaks_all[[1]]
		peaks_sample <- peaks_all[[2]]
		#################### Number of signals
		number_of_signals_samples <- length(peaks_sample@mass)
		number_of_signals_database <- length(peaks_database_temp@mass)
		number_of_signals_database_matrix[1,db] <- number_of_signals_database
		###### COUNTER 0 - MATCHING SIGNALS
		# Create a counter, symmetrical to the database Peaklist
		matching_signals_number <- 0
		matching_signals <- numeric()
		# For each peaklist in the Library
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			# Scroll the peaks in the sample
			for (i in 1:length(peaks_sample@mass)) {
				# Scroll the peaklist in the library
				for (j in 1:length(peaks_database_temp@mass)) {
					# Count the number of closely matching signals with the reference in the sample peaklist
					if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
						matching_signals_number <- matching_signals_number + 1
						matching_signals <- append(matching_signals, peaks_sample@mass[i])
						break
					}
				}
			}
		} else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
				matching_signals_number <- 0
		}
		# Append this row to the global matrix
		matching_signals_matrix [1,db] <- matching_signals_number
		###### COUNTER 1 - FIT
		fit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			fit_sample <- matching_signals_number / length(peaks_sample@mass)
		}
		# Append this row to the global matrix
		fit_matrix[1,db] <- fit_sample
		###### COUNTER 2 - RETRO FIT
		retrofit_sample <- 0
		if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
			retrofit_sample <- matching_signals_number / length(peaks_database_temp@mass)
		}
		# Append this row to the global matrix
		retrofit_matrix[1,db] <- retrofit_sample
		###### COUNTER 3
		# Symmetry -> comparison between intensities
		# Compute the similaroty index with the library
		similarity_index_matrix_global <- intensityMatrix(peaks_all, spectra_all)
		# Similarity index (E Id * Ix / sqrt(E Id^2 * E Ix^2)) = A / sqrt (B * E)
		A <- 0
		for (z in 1:ncol(similarity_index_matrix_global)) {
			A <- A + (similarity_index_matrix_global[1,z]*similarity_index_matrix_global[2,z])
		}
		B <- 0
		for (z in 1:ncol(similarity_index_matrix_global)) {
			B <- B + (similarity_index_matrix_global[1,z])^2
		}
		E <- 0
		for (z in 1:ncol(similarity_index_matrix_global)) {
			E <- E + (similarity_index_matrix_global[2,z])^2
		}
		similarity_index_sample <- A / sqrt (B * E)
		# Append this row to the global matrix
		similarity_index_matrix[1,db] <- similarity_index_sample
	}
	# Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
	return(list(number_of_signals_samples=number_of_signals_samples, number_of_signals_database_matrix=number_of_signals_database_matrix, matching_signals_matrix=matching_signals_matrix, fit_matrix=fit_matrix, retrofit_matrix=retrofit_matrix, similarity_index_matrix=similarity_index_matrix))
}
if (multicore_processing == TRUE) {
	# Detect the number of cores
	cpu_thread_number <- detectCores(logical=TRUE) - 1
	if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
		output_list <- mclapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_similarity_index(global_list), mc.cores=cpu_thread_number)
	} else if (Sys.info()[1] == "Windows") {
		# Make the CPU cluster for parallelisation
		cls <- makeCluster(cpu_thread_number)
		# Make the cluster use the custom functions and the package functions along with their parameters
		clusterEvalQ(cls, {library(MALDIquant)})
		clusterExport(cl=cls, varlist=c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_correlation", "database_size", "tof_mode", "peaks_filtering", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_peaks_removal", "low_intensity_percentage_threshold", "low_intensity_threshold_method"), envir=environment())
		output_list <- parLapply(cls, global_list, fun=function(global_list) comparison_sample_db_subfunction_similarity_index(global_list))
		stopCluster(cls)
	}
} else {
	output_list <- lapply(global_list, FUN=function(global_list) comparison_sample_db_subfunction_similarity_index(global_list))
}
############################ Merge the matrix pieces together
matching_signals_matrix_all <- NULL
number_of_signals_database_matrix_all <- NULL
fit_matrix_all <- NULL
retrofit_matrix_all <- NULL
similarity_index_matrix_all <- NULL
for (ns in 1:number_of_samples) {
	# Matching signals
	if (is.null(matching_signals_matrix_all)) {
		matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
	} else {
		matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
	}
	# Number of signals database
	if (is.null(number_of_signals_database_matrix_all)) {
		number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
	} else {
		number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
	}
	# Fit
	if (is.null(fit_matrix_all)) {
		fit_matrix_all <- output_list[[ns]]$fit_matrix
	} else {
		fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
	}
	# Retrofit
	if (is.null(retrofit_matrix_all)) {
		retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
	} else {
		retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
	}
	# Similarity index
	if (is.null(similarity_index_matrix_all)) {
		similarity_index_matrix_all <- output_list[[ns]]$similarity_index_matrix
	} else {
		similarity_index_matrix_all <- rbind(similarity_index_matrix_all, output_list[[ns]]$similarity_index_matrix)
	}
}
######################################
# Fix the rownames and colnames
rownames(similarity_index_matrix_all) <- sample_vector
colnames(similarity_index_matrix_all) <- database_vector
################### Score calculation
score <- log10(fit_matrix_all*retrofit_matrix_all*similarity_index_matrix_all*1000)
colnames(score) <- database_vector
rownames(score) <- sample_vector
#### Output the classification
output <- matrix ("", nrow=number_of_samples, ncol=database_size)
colnames(output) <- database_vector
rownames(output) <- sample_vector
if (spectra_path_output == TRUE) {
	output <- cbind(output, spectra_path_vector)
    colnames(output) <- c(database_vector, "Spectrum path")
}
if (score_only == TRUE) {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
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
} else {
	for (r in 1:number_of_samples) {
		for (w in 1:database_size) {
			if (score[r,w] >= 2) {
				output[r,w] <- paste("YES","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits=3), ")")
			}
			if (score[r,w] < 1.5) {
				output[r,w] <- paste("NO","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits=3), ")")
			}
			if (score[r,w] >= 1.5 && score[r,w] < 2) {
				output[r,w] <- paste("NI","(Score:", round(score[r,w], digits=3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, "=", round(fit_matrix_all[r,w], digits=3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], "=", round(retrofit_matrix_all[r,w], digits=3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits=3), ")")
			}
		}
	}
}
return (output)
}





################################################################################





#################################################### ADJACENCY MATRIX GENERATION
# The function generates an adjacency matrix from a peak list matrix, in order to generate a graph. The matrix is computed first by generating a correlation matrix and then replacing the correlation coefficients with 0 or 1 according to a threshold value.
generate_adjacency_matrix <- function(peaklist_matrix, correlation_method = "pearson", correlation_threshold=0.8, pvalue_threshold=0.05) {
##### Install the required packages
install_and_load_required_packages("stats")
##### Transpose the peaklist matrix to compute the correlation between observations and not features
peaklist_matrix_t <- t(peaklist_matrix)
##### Generate the correlation matrix
correlation_matrix <- cor(peaklist_matrix_t, method=correlation_method)
### If the p-value is not considered...
if (pvalue_threshold == 0 || is.null(pvalue_threshold)) {
	##### Generate the function to apply to the matrix (x = matrix entry)
	matrix_replacement_subfunction <- function(x, threshold) {
		if (abs(x) >= threshold) {
			x <- 1
		} else {
			x <- 0
		}
	}
	##### Generate the final adjacency matrix
	adjacency_matrix <- apply(correlation_matrix, MARGIN=c(1,2), FUN=function(x) matrix_replacement_subfunction(x, threshold=correlation_threshold))
} else {
	##### Install the required packages
	install_and_load_required_packages("psych")
	##### Generate the correlation matrix
	# Generate unique row and column names for the corr.test function
	rownames(peaklist_matrix_t) <- seq(1:nrow(peaklist_matrix_t))
	colnames(peaklist_matrix_t) <- seq(1:ncol(peaklist_matrix_t))
	correlation_matrix_pvalue <- corr.test(peaklist_matrix_t, adjust="none", ci=FALSE)$p
	##### Generate a global matrix (each entry displays "coefficient, pvalue")
	global_matrix <- matrix("", nrow=nrow(correlation_matrix_pvalue), ncol=ncol(correlation_matrix_pvalue))
	for (rw in 1:nrow(correlation_matrix_pvalue)) {
		for (cl in 1:ncol(correlation_matrix_pvalue)) {
			global_matrix[rw,cl] <- paste(correlation_matrix[rw,cl], correlation_matrix_pvalue[rw,cl], sep=",")
		}
	}
	##### Generate the function to apply to the matrix (x = matrix entry)
	matrix_replacement_subfunction2 <- function(x, coeff_threshold, p_threshold) {
		# Split the entry
		splitted_x <- as.numeric(unlist(strsplit(x, ",")))
		if (abs(splitted_x[1]) >= coeff_threshold && abs(splitted_x[2]) <= p_threshold) {
			x <- 1
		} else {
			x <- 0
		}
	}
	##### Generate the final adjacency matrix
	adjacency_matrix <- apply(global_matrix, MARGIN=c(1,2), FUN=function(x) matrix_replacement_subfunction2(x, coeff_threshold=correlation_threshold, p_threshold=pvalue_threshold))
}
##### Return the matrix
return(adjacency_matrix)
}





################################################################################





############################### GENERATE A CUSTOM INTENSITY MATRIX
# The function takes a list of spectra and a vector of custom features to be included in the generation of the final peaklist intensity matrix. The functions takes the spectra, preprocesses the spectra according to the specified parameters, performs the peak picking and outputs the intensity matrix only for the peaks specified as input (not all of those custom peaks if they are outside of the spectral mass range).
# If the range provided is too large, the function will return a NULL value, since some custom features cannot be found.
# This function is suited for aligning the spectral features (of an unknown dataset) with the model features.
generate_custom_intensity_matrix <- function(spectra, custom_feature_vector = NULL, tof_mode = "linear", spectra_preprocessing = TRUE, preprocessing_parameters = list(crop_spectra = TRUE, mass_range = c(800,3000), data_transformation = FALSE, transformation_algorithm = "sqrt", smoothing_algorithm = NULL, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 100, normalisation_algorithm = "TIC", normalisation_mass_range = NULL), peak_picking_algorithm = "SuperSmoother", peak_picking_SNR = 5, peaks_filtering = TRUE, frequency_threshold_percent = 10, low_intensity_peaks_removal = FALSE, intensity_threshold_percent = 1, intensity_threshold_method = "element-wise", process_in_packages_of = 0, multicore_processing = FALSE) {
	### Install the required packages
	install_and_load_required_packages("MALDIquant")
	# Rename the trim function
	trim_spectra <- get(x="trim", pos="package:MALDIquant")
	### Define the tolerance
	if (tof_mode == "linear" || tof_mode == "L") {
		tolerance_ppm <- 2000
	} else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
		tolerance_ppm <- 200
	}
	### Preprocessing
	if (spectra_preprocessing == TRUE) {
		spectra <- preprocess_spectra(spectra, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, process_in_packages_of = process_in_packages_of, align_spectra = FALSE, spectra_alignment_method = "cubic", multicore_processing = multicore_processing)
	}
	### Peak picking and alignment
	peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = peak_picking_SNR, multicore_processing = multicore_processing)
	peaks <- align_and_filter_peaks(peaks, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, peaks_filtering = peaks_filtering, frequency_threshold_percent = frequency_threshold_percent, low_intensity_peaks_removal = low_intensity_peaks_removal, intensity_threshold_percent = intensity_threshold_percent, intensity_threshold_method = intensity_threshold_method, reference_peaklist = NULL, spectra = spectra, alignment_iterations = 5, multicore_processing = multicore_processing)
	### Peaklist matrix
	# If there are more spectra...
	if (isMassSpectrumList(spectra) && isMassPeaksList(peaks)) {
		peaklist_matrix <- intensityMatrix(peaks, spectra)
	} else if (isMassSpectrum(spectra) && isMassPeaks(peaks)) {
		# If there is only one spectrum...
		peaklist_matrix <- as.matrix(rbind(peaks@intensity))
		colnames(peaklist_matrix) <- peaks@mass
	}
	### Run the alignment only if the vector of custom features is not null
	if (!is.null(custom_feature_vector)) {
		# Convert the custom feature vector in numeric
		custom_feature_vector <- as.numeric(custom_feature_vector)
		# Retrieve the peaks in the spectral dataset
		spectra_peaks <- as.numeric(colnames(peaklist_matrix))
		# Generate a temporary vector for the cistom peaks
		custom_feature_vector_final <- vector()
		# Check if every element of the custom vector is within the mass range (and create another vector with the compatible custom features)
		for (f in 1:length(custom_feature_vector)) {
			if (custom_feature_vector[f] >= spectra_peaks[1] && custom_feature_vector[f] <= spectra_peaks[length(spectra_peaks)]) {
				custom_feature_vector_final <- append(custom_feature_vector_final, custom_feature_vector[f])
			}
		}
		### Check the compatibility between the spectra and the provided mass list
		# Isolate the dataset features (sorted)
		spectral_dataset_features <- sort(as.numeric(colnames(peaklist_matrix)))
		# Sort the custom features
		custom_feature_vector_final <- sort(as.numeric(custom_feature_vector_final))
		# Define the compatibility
		if (spectral_dataset_features[1] <= custom_feature_vector_final[1] && spectral_dataset_features[length(spectral_dataset_features)] >= custom_feature_vector_final[length(custom_feature_vector_final)]) {
			feature_compatibility <- TRUE
		} else {
			feature_compatibility <- FALSE
		}
		### If there is feature compatibility...
		if (feature_compatibility == TRUE) {
			# Fix the original custom feature vector
			custom_feature_vector <- as.character(custom_feature_vector_final)
			### Determine the columns to keep and the column to add
			features_to_keep <- numeric()
			features_to_add <- numeric()
			adjusted_features_to_keep <- numeric()
			# For each feature in the custom feature vector
			for (csft in custom_feature_vector) {
				# Set the default presence of the signal in the sample to FALSE
				presence <- FALSE
				# Scroll the sample features
				for (ft in colnames(peaklist_matrix)) {
					# If there is a match
					if (abs((as.numeric(csft)-as.numeric(ft))*10^6/as.numeric(csft)) <= tolerance_ppm) {
						# Add it to the features to keep
						features_to_keep <- append(features_to_keep, ft)
						# Align the feature in the sample with the one in the custom feature vector
						ft <- csft
						# Add it to another list (it will be used to adjust the column names in the final sample peaklist)
						adjusted_features_to_keep <- append(adjusted_features_to_keep, ft)
						# Set the presence of the signal in the sample to TRUE
						presence <- TRUE
						# Avoid consecutive duplicates (once it is found there is no point in keep going)
						break
					}
				}
				# If after all the signal in the custom feature vector is not found in the sample
				if (presence == FALSE) {
					# Add this to the features to be added
					features_to_add <- append(features_to_add, csft)
				}
			}
			### Generate the final sample matrix (with the right column names)
			if (length(features_to_keep) > 0) {
				final_peaklist_matrix <- as.matrix(rbind(peaklist_matrix [,features_to_keep]))
				colnames(final_peaklist_matrix) <- adjusted_features_to_keep
			} else {
				final_peaklist_matrix <- NULL
			}
			### Add the missing features
			## Multiple spectra
			if (isMassSpectrumList(spectra)) {
				# If there are features to add...
				if (length(features_to_add > 0)) {
					# Generate a fake spectrum and a fake peaklist with the features to add
					fake_spectrum <- createMassSpectrum(mass=spectra[[1]]@mass, intensity=spectra[[1]]@intensity, metaData=list(name="Fake spectrum"))
					fake_peaks <- createMassPeaks(mass=as.numeric(features_to_add), intensity=rep(1, length(features_to_add)), snr=rep(3, length(features_to_add)), metaData=list(name="Fake peaklist"))
					# Detect the peaks in the spectra
					peaks <- detectPeaks(spectra, method="SuperSmoother", SNR=3)
					# Append the fake spectrum and the fake peaklist to he original lists
					spectra_all <- append(spectra, fake_spectrum)
					peaks_all <- append(peaks, fake_peaks)
					# Generate the intensity matrix
					intensity_matrix_all <- intensityMatrix(peaks_all, spectra_all)
					# Remove the last row (corresponding to the fake spectrum)
					intensity_matrix_all <- intensity_matrix_all[1:(nrow(intensity_matrix_all)-1),]
					# Keep only the columns that are corresponding to the desired features
					final_intensity_matrix <- intensity_matrix_all[,features_to_add]
					# If the final matrix does not exist yet and it is null, the final matrix becomes the feature column
					if (is.null(final_peaklist_matrix)) {
						final_peaklist_matrix <- final_peaklist_matrix
					} else {
						# If the final matrix exists, append the feature column to the matrix
						final_peaklist_matrix <- cbind(final_peaklist_matrix, final_intensity_matrix)
					}
				}
			} else if (isMassSpectrum(spectra)) {
				# If there are features to add...
				if (length(features_to_add > 0)) {
					# Scroll the features to add (in the model but not in the sample)
					for (f in 1:length(features_to_add)) {
						# Initialise the output
						x_intensity <- NA
						# Scroll the mass list of each spectrum
						for (m in 1:length(spectra@mass)) {
							# If there is a match
							if (abs(spectra@mass[m]-as.numeric(features_to_add[f]))*10^6/as.numeric(features_to_add[f]) <= tolerance_ppm) {
								# Store the corresponding intensity and generate the matrix column
								x_intensity <- as.matrix(cbind(spectra@intensity[m]))
								colnames(x_intensity) <- features_to_add[f]
								# If the final matrix does not exist yet and it is null, the final matrix becomes the feature column
								if (is.null(final_peaklist_matrix)) {
									final_peaklist_matrix <- cbind(x_intensity)
								} else {
									# If the final matrix exists, append the feature column to the matrix
									final_peaklist_matrix <- cbind(final_peaklist_matrix, x_intensity)
								}
								# Do not keep searching
								break
							}
						}
					}
				}
			}
			### Return the final matrix with the custom features
			return(final_peaklist_matrix)
		} else {
			### Return a NULL value if features are not compatible
			return(NULL)
		}
	} else {
		### Return the simple peaklist matrix if no custom vector is provided
		return(peaklist_matrix)
	}
}





################################################################################





##################################################### REARRANGE THE SPECTRAL DATASET
# The function takes a list of spectra and rearranges it in a certain way (defined by the user) on spatial, random or hierarchical basis. It returns a list of spectra, appropriately sorted.
# If space is selected, the spectral list is rearranged along the widest dimension: for example, if X is the widest dimension, the spectra are rearranged as X1,Y1 X2,Y1 X3,Y1 and so on... If the spectral list does not contain spatial coordinates (it is not an imaging dataset), nothing is performedin terms of space.
rearrange_spectral_dataset <- function(spectra, rearranging_method = c("space","random"), seed = NULL) {
	##### Install and load the required packages
	install_and_load_required_packages("MALDIquant")
	########## Do everything only if it is a list of spectra
	if (isMassSpectrumList(spectra)) {
		#################### SPACE
		if (rearranging_method == "space") {
			########## Do everything only if there are spatial coordinates
			if (!is.null(spectra[[1]]@metaData$imaging$pos)) {
				##### Generate the matrix of spectral coordinates
				spectral_coordinates <- matrix(0, ncol = 2, nrow = length(spectra))
				### Fill in the matrix
				for (s in 1:length(spectra)) {
					spectral_coordinates[s,] <- spectra[[s]]@metaData$imaging$pos
				}
				### Fix the matrix column and row names
				rownames(spectral_coordinates) <- seq(1:nrow(spectral_coordinates))
				colnames(spectral_coordinates) <- c("x", "y")
				### Convert it into a dataframe for sorting
				spectral_coordinates <- as.data.frame(spectral_coordinates)
				### Sort along the greater dimension
				## Define the rectangle containing the tissue section
				rectangle_coordinates <- matrix(0, nrow = 2, ncol = 2)
				rectangle_coordinates[1,1] <- min(spectral_coordinates[,1])
				rectangle_coordinates[2,1] <- max(spectral_coordinates[,1])
				rectangle_coordinates[1,2] <- min(spectral_coordinates[,2])
				rectangle_coordinates[2,2] <- max(spectral_coordinates[,2])
				rownames(rectangle_coordinates) <- c("min", "max")
				colnames(rectangle_coordinates) <- c("x", "y") 
				## Determine the WIDEST coordinate
				if (rectangle_coordinates[2,1] >= rectangle_coordinates[2,2]) {
					widest_coordinate <- "x"
				} else if (rectangle_coordinates[2,2] > rectangle_coordinates[2,1]) {
					widest_coordinate <- "y"
				}
				## X coordinate is the widest
				if (widest_coordinate == "x") {
					# Sort according the X
					spectral_coordinates_sorted <- spectral_coordinates[order(spectral_coordinates[,2], spectral_coordinates[,1], decreasing = FALSE), ]
					# Extract the row names, which are the ID numbers of the spectral list
					spectra_ID <- as.integer(rownames(spectral_coordinates_sorted))
					# Rearrange the spectral list
					spectra_riarranged <- spectra[spectra_ID]
				} else if (widest_coordinate == "y") {
					### Y coordinate is the widest
					# Sort according the Y
					spectral_coordinates_sorted <- spectral_coordinates[order(spectral_coordinates[,1], spectral_coordinates[,2], decreasing = FALSE), ]
					# Extract the row names, which are the ID numbers of the spectral list
					rearranged_spectra_IDs <- as.integer(rownames(spectral_coordinates_sorted))
					# Rearrange the spectral list
					spectra_riarranged <- spectra[rearranged_spectra_IDs]
				}
			} else {
				spectra_riarranged <- spectra
			}
		} else if (rearranging_method == "random") {
			#################### RANDOM
			##### Rearrange the spectra randomly in N parts
			##### Install and load the required packages
			install_and_load_required_packages("caret")
			# Set the seed (make randomness reproducible)
			if (!is.null(seed)) {
				set.seed(seed)
			}
			# Generate a random list of numbers for random rearrangement
			if (!is.null(seed)) {
				set.seed(seed)
			}
			rearranged_spectra_IDs <- sample(1:length(spectra), size = length(spectra))
			# Rearrange the spectral list
			spectra_riarranged <- spectra[rearranged_spectra_IDs]
		} else if (rearranging_method == "hca") {
			#################### HCA
			NULL
		}
		#################### Output
		return(spectra_riarranged)
	} else if (isMassSpectrum(spectra)) {
		########## Single spectrum
		return(spectra)
	}
}





################################################################################





##################################################### PARTITION THE SPECTRAL DATASET IN SUBSETS
# The function takes a list of spectra and partitions it in a certain number of subsets (defined by the user) on spatial, random or hierarchical basis. It returns a list of sub-lists of spectra.
partition_spectral_dataset <- function(spectra, partitioning_method = c("space","random", "hca"), number_of_partitions = 3, seed = NULL, tof_mode = "reflectron") {
	##### Install and load the required packages
	install_and_load_required_packages("MALDIquant")
	########## Do everything only if it is a list of spectra
	if (isMassSpectrumList(spectra)) {
		#################### SPACE
		if (partitioning_method == "space") {
			########## Do everything only if there are spatial coordinates
			if (!is.null(spectra[[1]]@metaData$imaging$pos)) {
				##### Generate the matrix of spectral coordinates
				spectral_coordinates <- matrix(0, ncol = 2, nrow = length(spectra))
				# Fill in the matrix
				for (s in 1:length(spectra)) {
					spectral_coordinates[s,] <- spectra[[s]]@metaData$imaging$pos
				}
				# Fix the matrix column and row names
				rownames(spectral_coordinates) <- seq(1:nrow(spectral_coordinates))
				colnames(spectral_coordinates) <- c("x", "y")
				# Define the rectangle containing the tissue section
				rectangle_coordinates <- matrix(0, nrow = 2, ncol = 2)
				rectangle_coordinates[1,1] <- min(spectral_coordinates[,1])
				rectangle_coordinates[2,1] <- max(spectral_coordinates[,1])
				rectangle_coordinates[1,2] <- min(spectral_coordinates[,2])
				rectangle_coordinates[2,2] <- max(spectral_coordinates[,2])
				rownames(rectangle_coordinates) <- c("min", "max")
				colnames(rectangle_coordinates) <- c("x", "y") 
				##### Split the dataset based upon the the WIDEST coordinate
				if (abs(rectangle_coordinates[2,1] - rectangle_coordinates[1,1]) >= abs(rectangle_coordinates[2,2] - rectangle_coordinates[1,2])) {
					widest_coordinate <- "x"
				} else {
					widest_coordinate <- "y"
				}
				### X coordinate is the widest
				if (widest_coordinate == "x") {
					# Define the intervals
					interval_width <- ceiling(abs(rectangle_coordinates[2,1] - rectangle_coordinates[1,1])/number_of_partitions)
					# Define the separating points
					interval_points <- seq.int(from = rectangle_coordinates[1,1], to = rectangle_coordinates[2,1], by = interval_width)
					# Add the maximum coordinate to the interval points (if not already present)
					if (interval_points[length(interval_points)] != rectangle_coordinates[2,1]) {
						interval_points <- c(interval_points, rectangle_coordinates[2,1])
					}
					## Generate the final list of partitioned spectra
					# Initialise the output
					spectra_partitioned <- list()
					# Populate the final list
					for (p in 1:(length(interval_points) - 1)) {
						# Generate the sublist of the spectra for that interval
						spectra_interval <- list()
						# Scroll the spectra...
						for (s in 1:length(spectra)) {
							# Populate the list...
							if (spectra[[s]]@metaData$imaging$pos[1] >= interval_points[p] && spectra[[s]]@metaData$imaging$pos[1] < interval_points[p + 1]) {
								spectra_interval <- append(spectra_interval, spectra[[s]])
							}
						}
						# Add the spectra with the max coordinate
						if (p == (length(interval_points) - 1)) {
							# Scroll the spectra...
							for (s in 1:length(spectra)) {
								# Populate the list...
								if (spectra[[s]]@metaData$imaging$pos[1] == interval_points[p + 1]) {
									spectra_interval <- append(spectra_interval, spectra[[s]])
								}
							}
						}
						spectra_partitioned[[p]] <- spectra_interval
					}
				} else if (widest_coordinate == "y") {
					### Y coordinate is the widest
					# Define the intervals
					interval_width <- ceiling((rectangle_coordinates[2,2] - rectangle_coordinates[1,2])/number_of_partitions)
					# Define the separating points
					interval_points <- seq.int(rectangle_coordinates[1,2], rectangle_coordinates[2,2], by = interval_width)
					# Add the maximum coordinate to the interval points (if not already present)
					if (interval_points[length(interval_points)] != rectangle_coordinates[2,2]) {
						interval_points <- c(interval_points, rectangle_coordinates[2,2])
					}
					## Generate the final list of partitioned spectra
					# Initialise the output
					spectra_partitioned <- list()
					# Populate the final list
					for (p in 1:(length(interval_points) - 1)) {
						# Generate the sublist of the spectra for that interval
						spectra_interval <- list()
						# Scroll the spectra...
						for (s in 1:length(spectra)) {
							# Populate the list...
							if (spectra[[s]]@metaData$imaging$pos[2] >= interval_points[p] && spectra[[s]]@metaData$imaging$pos[2] < interval_points[p + 1]) {
								spectra_interval <- append(spectra_interval, spectra[[s]])
							}
						}
						# Add the spectra with the max coordinate
						if (p == (length(interval_points) - 1)) {
							# Scroll the spectra...
							for (s in 1:length(spectra)) {
								# Populate the list...
								if (spectra[[s]]@metaData$imaging$pos[2] == interval_points[p + 1]) {
									spectra_interval <- append(spectra_interval, spectra[[s]])
								}
							}
						}
						spectra_partitioned[[p]] <- spectra_interval
					}
				}
			} else {
				spectra_partitioned <- spectra
			}
		} else if (partitioning_method == "random") {
			#################### RANDOM
			##### Split the spectra randomly in N parts
			##### Install and load the required packages
			install_and_load_required_packages("caret")
			# Set the seed (make randomness reproducible)
			if (!is.null(seed)) {
				set.seed(seed)
			}
			# Partition the spectral list
			spectra_partitioned_IDs <- createFolds(y = rep("spectra", length(spectra)), k = number_of_partitions, returnTrain = FALSE)
			# Generate the split sublist of spectra
			spectra_partitioned <- list()
			# For each element of the list containing the IDs of the fold...
			for (s in 1:length(spectra_partitioned_IDs)) {
				# Extract the corresponding spectra and put them in a final list
				spectra_partitioned[[s]] <- spectra[spectra_partitioned_IDs[[s]]]
			}
		} else if (partitioning_method == "hca") {
			#################### HCA
			### Detect and align peaks
			peaks <- peak_picking(spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), multicore_processing = FALSE)
			peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peaks_filtering = TRUE, frequency_threshold_percent = 10)
			### Generate the peaklist matrix
			peaklist <- intensityMatrix(peaks, spectra)
			### Compute the distance matrix
			distance_matrix <- dist(peaklist, method = "euclidean")
			# Generate the dendrogram
			hca <- hclust(distance_matrix)
			### Cut the tree to generate K number of sub-clusters
			hca_groups <- cutree(hca, k = number_of_partitions)
			### Initialise the output
			spectra_partitioned <- list()
			### For each subgroup to be isolated...
			for (p in 1:number_of_partitions) {
				# Index the spectra under in the selected subgroup of the HCA
				index <- which(hca_groups == p)
				spectra_hca <- spectra[index]
				# Store them in the final list
				spectra_partitioned[[p]] <- spectra_hca
			}
		}
		#################### Output
		return(spectra_partitioned)
	} else if (isMassSpectrum(spectra)) {
		########## Single spectrum
		return(spectra)
	}
}




































##########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

































################ SPECTRAL TYPER PROGRAM 2017.01.31

############## INSTALL AND LOAD THE REQUIRED PACKAGES
install_and_load_required_packages(c("tcltk", "XLConnect", "ggplot2"), repository="http://cran.mirror.garr.it/mirrors/CRAN/")


# Update packages and load the required packages
update.packages(repos="http://cran.mirror.garr.it/mirrors/CRAN/", ask=FALSE)









###################################### Initialise the variables (default values)
filepath_database <- NULL
filepath_test <- NULL
output_folder <- NULL
average_replicates_in_database <- FALSE
average_replicates_in_test <- FALSE
peaks_filtering <- FALSE
low_intensity_peaks_removal <- FALSE
intensity_threshold_method <- "element-wise"
tof_mode <- "linear"
spectra_format <- "brukerflex"
score_only <- FALSE
spectra_path_output <- TRUE
peak_picking_mode <- "all"
peak_picking_algorithm <- "SuperSmoother"
similarity_criteria <- "correlation"
signal_intensity_evaluation <- "intensity percentage"
file_type_export <- "xlsx"
spectra_database <- NULL
spectra_test <- NULL
peaks_database <- NULL
peaks_test <- NULL
multicore_processing <- TRUE
transform_data <- FALSE
transform_data_algorithm <- NULL
smoothing <- TRUE
smoothing_algorithm <- "SavitzkyGolay"
smoothing_strength <- "medium"
baseline_subtraction <- TRUE
baseline_subtraction_algorithm <- "SNIP"
baseline_subtraction_iterations <- 200
normalisation <- TRUE
normalisation_algorithm <- "TIC"
normalisation_mass_range <- NULL
preprocess_spectra_in_packages_of <- 200
mass_range <- c(3000,15000)
preprocessing_parameters_database <- list(crop_spectra=TRUE, mass_range=NULL, data_transformation=NULL, transformation_algorithm=NULL, smoothing_algorithm=NULL, smoothing_strength=NULL, baseline_subtraction_algorithm=NULL, baseline_subtraction_iterations=NULL, normalisation_algorithm=NULL, normalisation_mass_range=NULL)






################## Values of the variables (for displaying and dumping purposes)
tof_mode_value <- "linear"
average_replicates_in_database_value <- "NO"
average_replicates_in_test_value <- "NO"
peaks_filtering_value <- "NO"
low_intensity_peaks_removal_value <- "NO"
score_only_value <- "NO"
spectra_path_output_value <- "YES"
similarity_criteria_value <- "correlation"
peak_picking_algorithm_value <- "   SuperSmoother   "
peak_picking_mode_value <- "all"
intensity_threshold_method_value <- "element-wise"
spectra_format_value <- "Xmass"
multicore_processing_value <- "YES"
transform_data_value <- "    NO    "
smoothing_value <- "YES ( SavitzkyGolay , medium)"
baseline_subtraction_value <- "YES (SNIP , iterations: 200)"
normalisation_value <- "YES (TIC , mass range: )"
mass_range_value <- "3000 , 15000"
preprocess_spectra_in_packages_of_value <- "200"






##################################################### DEFINE WHAT THE BUTTONS DO

##### Preprocessing window
preprocessing_window_function <- function() {
	##### Functions
	# Transform the data
	transform_data_choice <- function() {
		# Catch the value from the menu
		transform_data <- select.list(c("YES","NO"), title="Choose")
		# Default
		if (transform_data == "YES") {
			transform_data <- TRUE
			# Ask for the algorithm
			transform_data_algorithm <- select.list(c("sqrt","log","log2","log10"), title="Choose")
			# Default
			if (transform_data_algorithm == "") {
				transform_data_algorithm <- "sqrt"
			}
		} else if (transform_data == "NO" || transform_data == "") {
			transform_data <- FALSE
		}
		# Set the value of the displaying label
		if (transform_data == TRUE) {
			transform_data_value <- paste("YES", "(", transform_data_algorithm, ")")
		} else {
			transform_data_value <- "    NO    "
		}
		transform_data_value_label <- tklabel(preproc_window, text=transform_data_value, font = label_font)
		tkgrid(transform_data_value_label, row=3, column=2)
		# Escape the function
		.GlobalEnv$transform_data <- transform_data
		.GlobalEnv$transform_data_algorithm <- transform_data_algorithm
		.GlobalEnv$transform_data_value <- transform_data_value
		# Add the parameters to the list of preprocessing parameters (to be dumped)
		.GlobalEnv$preprocessing_parameters_database$data_transformation <- transform_data
		.GlobalEnv$preprocessing_parameters_database$transformation_algorithm <- transform_data_algorithm
	}
	# Smoothing
	smoothing_choice <- function() {
		# Catch the value from the menu
		smoothing <- select.list(c("YES","NO"), title="Choose")
		# Default
		if (smoothing == "YES" || smoothing == "") {
			smoothing <- TRUE
			# Ask for the algorithm
			smoothing_algorithm <- select.list(c("SavitzkyGolay","MovingAverage"), title="Choose")
			# Default
			if (smoothing_algorithm == "") {
				smoothing_algorithm <- "SavitzkyGolay"
			}
			# Strength
			smoothing_strength <- select.list(c("medium","strong","stronger"), title="Choose")
			if (smoothing_strength == "") {
				smoothing_strength <- "medium"
			}
		} else if (smoothing == "NO") {
			smoothing <- FALSE
		}
		# Set the value of the displaying label
		if (smoothing == TRUE) {
			smoothing_value <- paste("YES", "(", smoothing_algorithm, "," , smoothing_strength, ")")
		} else {
			smoothing_value <- "    NO    "
		}
		smoothing_value_label <- tklabel(preproc_window, text=smoothing_value, font = label_font)
		tkgrid(smoothing_value_label, row=4, column=2)
		# Escape the function
		.GlobalEnv$smoothing <- smoothing
		.GlobalEnv$smoothing_strength <- smoothing_strength
		.GlobalEnv$smoothing_algorithm <- smoothing_algorithm
		.GlobalEnv$smoothing_value <- smoothing_value
		# Add the parameters to the list of preprocessing parameters (to be dumped)
		.GlobalEnv$preprocessing_parameters_database$smoothing_strength <- smoothing_strength
		.GlobalEnv$preprocessing_parameters_database$smoothing_algorithm <- smoothing_algorithm
	}
	# Baseline subtraction
	baseline_subtraction_choice <- function() {
		# Catch the value from the menu
		baseline_subtraction <- select.list(c("YES","NO"), title="Choose")
		# Default
		if (baseline_subtraction == "YES" || baseline_subtraction == "") {
			baseline_subtraction <- TRUE
			# Ask for the algorithm
			baseline_subtraction_algorithm <- select.list(c("SNIP","TopHat","ConvexHull","median"), title="Choose")
			# SNIP
			if (baseline_subtraction_algorithm == "SNIP") {
				baseline_subtraction_iterations <- tclvalue(baseline_subtraction_iterations2)
				baseline_subtraction_iterations_value <- as.character(baseline_subtraction_iterations)
				baseline_subtraction_iterations <- as.integer(baseline_subtraction_iterations)
			}
			# Default
			if (baseline_subtraction_algorithm == "") {
				baseline_subtraction_algorithm <- "SNIP"
				baseline_subtraction_iterations <- 200
			}
		} else if (baseline_subtraction == "NO") {
			baseline_subtraction <- FALSE
		}
		# Set the value of the displaying label
		if (baseline_subtraction == TRUE && baseline_subtraction_algorithm != "SNIP") {
			baseline_subtraction_value <- paste("YES", "(", baseline_subtraction_algorithm, ")")
		} else if (baseline_subtraction == TRUE && baseline_subtraction_algorithm == "SNIP") {
			baseline_subtraction_value <- paste("YES", "(", baseline_subtraction_algorithm, ", iterations:", baseline_subtraction_iterations, ")")
		} else {
			baseline_subtraction_value <- "    NO    "
		}
		baseline_subtraction_value_label <- tklabel(preproc_window, text=baseline_subtraction_value, font = label_font)
		tkgrid(baseline_subtraction_value_label, row=5, column=3)
		# Escape the function
		.GlobalEnv$baseline_subtraction <- baseline_subtraction
		.GlobalEnv$baseline_subtraction_iterations <- baseline_subtraction_iterations
		.GlobalEnv$baseline_subtraction_algorithm <- baseline_subtraction_algorithm
		.GlobalEnv$baseline_subtraction_value <- baseline_subtraction_value
		# Add the parameters to the list of preprocessing parameters (to be dumped)
		.GlobalEnv$preprocessing_parameters_database$baseline_subtraction_iterations <- baseline_subtraction_iterations
		.GlobalEnv$preprocessing_parameters_database$baseline_subtraction_algorithm <- baseline_subtraction_algorithm
	}
	# Normalisation
	normalisation_choice <- function() {
		# Catch the value from the menu
		normalisation <- select.list(c("YES","NO"), title="Choose")
		# Default
		if (normalisation == "YES" || normalisation == "") {
			normalisation <- TRUE
			# Ask for the algorithm
			normalisation_algorithm <- select.list(c("TIC","PQN","median"), title="Choose")
			# TIC
			if (normalisation_algorithm == "TIC") {
				normalisation_mass_range <- tclvalue(normalisation_mass_range2)
				normalisation_mass_range_value <- as.character(normalisation_mass_range)
				if (normalisation_mass_range != 0 && normalisation_mass_range != "") {
					normalisation_mass_range <- unlist(strsplit(normalisation_mass_range, ","))
					normalisation_mass_range <- as.numeric(normalisation_mass_range)
				} else if (normalisation_mass_range == 0 || normalisation_mass_range == "") {
					normalisation_mass_range <- NULL
				}
			}
			# Default
			if (normalisation_algorithm == "") {
				normalisation_algorithm <- "TIC"
			}
		} else if (normalisation == "NO") {
			normalisation <- FALSE
		}
		# Set the value of the displaying label
		if (normalisation == TRUE && normalisation_algorithm != "TIC") {
			normalisation_value <- paste("YES", "(", normalisation_algorithm, ")")
		} else if (normalisation == TRUE && normalisation_algorithm == "TIC") {
			normalisation_value <- paste("YES", "(", normalisation_algorithm, ", range:", normalisation_mass_range_value, ")")
		} else {
			normalisation_value <- "    NO    "
		}
		normalisation_value_label <- tklabel(preproc_window, text=normalisation_value, font = label_font)
		tkgrid(normalisation_value_label, row=6, column=3)
		# Escape the function
		.GlobalEnv$normalisation <- normalisation
		.GlobalEnv$normalisation_mass_range <- normalisation_mass_range
		.GlobalEnv$normalisation_algorithm <- normalisation_algorithm
		.GlobalEnv$normalisation_value <- normalisation_value
		# Add the parameters to the list of preprocessing parameters (to be dumped)
		.GlobalEnv$preprocessing_parameters_database$normalisation_mass_range <- normalisation_mass_range
		.GlobalEnv$preprocessing_parameters_database$normalisation_algorithm <- normalisation_algorithm
	}
	# TOF mode
	tof_mode_choice <- function() {
		# Catch the value from the menu
		tof_mode <- select.list(c("Linear","Reflector"), title="Choose")
		# Default
		if (tof_mode == "" || tof_mode == "Linear") {
			tof_mode <- "linear"
		}
		if (tof_mode == "Reflector") {
			tof_mode <- "reflector"
		}
		# Set the value of the displaying label
		tof_mode_value <- tof_mode
		if (tof_mode_value == "linear") {
			tof_mode_value <- "   linear   "
		}
		tof_mode_value_label <- tklabel(preproc_window, text=tof_mode_value, font = label_font)
		tkgrid(tof_mode_value_label, row=2, column=3)
		# Escape the function
		.GlobalEnv$tof_mode <- tof_mode
		.GlobalEnv$tof_mode_value <- tof_mode_value
	}
	# Commit preprocessing
	commit_preprocessing_function <- function() {
		# Get the values (they are filled with the default anyway)
		# Mass range
		mass_range <- tclvalue(mass_range2)
		mass_range_value <- as.character(mass_range)
		mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
		# Preprocessing
		preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of2)
		preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)
		preprocess_spectra_in_packages_of_value <- as.character(preprocess_spectra_in_packages_of)
		# Escape the function
		.GlobalEnv$mass_range <- mass_range
		.GlobalEnv$mass_range_value <- mass_range_value
		.GlobalEnv$preprocess_spectra_in_packages_of <- preprocess_spectra_in_packages_of
		.GlobalEnv$preprocess_spectra_in_packages_of_value <- preprocess_spectra_in_packages_of_value
		# Add the parameters to the list of preprocessing parameters (to be dumped)
		.GlobalEnv$preprocessing_parameters_database$mass_range <- mass_range
		# Destroy the window upon committing
		tkdestroy(preproc_window)
	}
	##### List of variables, whose values are taken from the entries in the GUI (create new variables for the sub window, that will replace the ones in the global environment, only if the default are changed)
	mass_range2 <- tclVar("")
	preprocess_spectra_in_packages_of2 <- tclVar("")
	baseline_subtraction_iterations2 <- tclVar("")
	normalisation_mass_range2 <- tclVar("")
	##### Window
	preproc_window <- tktoplevel()
	tktitle(preproc_window) <- "Spectra preprocessing parameters"
	# Mass range
	mass_range_label <- tklabel(preproc_window, text="Mass range", font = label_font)
	mass_range_entry <- tkentry(preproc_window, width=15, textvariable=mass_range2, font = entry_font)
	tkinsert(mass_range_entry, "end", as.character(paste(mass_range[1],",",mass_range[2])))
	# Preprocessing (in packages of)
	preprocess_spectra_in_packages_of_label <- tklabel(preproc_window, text="Preprocess spectra\nin packages of", font = label_font)
	preprocess_spectra_in_packages_of_entry <- tkentry(preproc_window, width=10, textvariable=preprocess_spectra_in_packages_of2, font = entry_font)
	tkinsert(preprocess_spectra_in_packages_of_entry, "end", as.character(preprocess_spectra_in_packages_of))
	# Tof mode
	tof_mode_label <- tklabel(preproc_window, text="Select the TOF mode", font = label_font)
	tof_mode_entry <- tkbutton(preproc_window, text="Choose the TOF mode", command=tof_mode_choice, font = button_font)
	# Transform the data
	transform_data_button <- tkbutton(preproc_window, text="Transform the data", command=transform_data_choice, font = button_font)
	# Smoothing
	smoothing_button <- tkbutton(preproc_window, text="Smoothing", command=smoothing_choice, font = button_font)
	# Baseline subtraction
	baseline_subtraction_button <- tkbutton(preproc_window, text="Baseline subtraction", command=baseline_subtraction_choice, font = button_font)
	baseline_subtraction_iterations_entry <- tkentry(preproc_window, width=15, textvariable=baseline_subtraction_iterations2, font = entry_font)
	tkinsert(baseline_subtraction_iterations_entry, "end", as.character(baseline_subtraction_iterations))
	# normalisation
	normalisation_button <- tkbutton(preproc_window, text="Normalisation", command=normalisation_choice, font = button_font)
	normalisation_mass_range_entry <- tkentry(preproc_window, width=15, textvariable=normalisation_mass_range2, font = entry_font)
	tkinsert(normalisation_mass_range_entry, "end", as.character(normalisation_mass_range))
	# Commit preprocessing
	commit_preprocessing_button <- tkbutton(preproc_window, text="Commit preprocessing", command=commit_preprocessing_function, font = button_font)
	##### Displaying labels
	tof_mode_value_label <- tklabel(preproc_window, text=tof_mode_value, font = label_font)
	transform_data_value_label <- tklabel(preproc_window, text=transform_data_value, font = label_font)
	smoothing_value_label <- tklabel(preproc_window, text=smoothing_value, font = label_font)
	baseline_subtraction_value_label <- tklabel(preproc_window, text=baseline_subtraction_value, font = label_font)
	normalisation_value_label <- tklabel(preproc_window, text=normalisation_value, font = label_font)
	#### Geometry manager
	tkgrid(mass_range_label, row=1, column=1)
	tkgrid(mass_range_entry, row=1, column=2)
	tkgrid(tof_mode_label, row=2, column=1)
	tkgrid(tof_mode_entry, row=2, column=2)
	tkgrid(tof_mode_value_label, row=2, column=3)
	tkgrid(transform_data_button, row=3, column=1)
	tkgrid(transform_data_value_label, row=3, column=2)
	tkgrid(smoothing_button, row=4, column=1)
	tkgrid(smoothing_value_label, row=4, column=2)
	tkgrid(baseline_subtraction_button, row=5, column=1)
	tkgrid(baseline_subtraction_iterations_entry, row=5, column=2)
	tkgrid(baseline_subtraction_value_label, row=5, column=3)
	tkgrid(normalisation_button, row=6, column=1)
	tkgrid(normalisation_mass_range_entry, row=6, column=2)
	tkgrid(normalisation_value_label, row=6, column=3)
	tkgrid(preprocess_spectra_in_packages_of_label, row=7, column=1)
	tkgrid(preprocess_spectra_in_packages_of_entry, row=7, column=2)
	tkgrid(commit_preprocessing_button, row=8, column=1)
}

##### File type (export)
file_type_export_choice <- function() {
	# Catch the value from the menu
	file_type_export <<- select.list(c("csv","xlsx","xls"), title="Choose")
	# Default
	if (file_type_export == "") {
		file_type_export <<- "xlsx"
	}
	# Escape the function
	#.GlobalEnv$file_type_export <- file_type_export
	# Set the value of the displaying label
	file_type_export_value_label <- tklabel(window, text=file_type_export, font = label_font)
	tkgrid(file_type_export_value_label, row=12, column=6)
}

##### Dump the peaklist of the database
database_dump_function <- function() {
	############### If there is a peaklist to be dumped
	if (!is.null(peaks_database) && !is.null(spectra_database)) {
		########## File name
		##### Catch the filename from the menu
		filename_peaklist <- tclvalue(file_name)
		filename_peaklist <- as.character(filename_peaklist)
		########## Dump the RData containing the list of the spectra and peaks in the database, along with the preprocessing parameters
		database_filename <- paste(filename_peaklist, " - Database.RData", sep = "")
		save(peaks_database, spectra_database, file = database_filename)
		##### Generate the output filename (based upon the filename)
		filename_peaklist <- paste(filename_peaklist, " - ", "Database peaklist", sep="")
		##### Add the extension if it is not present in the filename
		if (file_type_export == "csv") {
			if (length(grep(".csv", filename_peaklist, fixed=TRUE)) == 1) {
				filename_peaklist <- filename_peaklist
			}	else {filename_peaklist <- paste(filename_peaklist, ".csv", sep="")}
		} else if (file_type_export == "xlsx") {
			if (length(grep(".xlsx", filename_peaklist, fixed=TRUE)) == 1) {
				filename_peaklist <- filename_peaklist
			}	else {filename_peaklist <- paste(filename_peaklist, ".xlsx", sep="")}
		} else if (file_type_export == "xls") {
			if (length(grep(".xls", filename_peaklist, fixed=TRUE)) == 1) {
				filename_peaklist <- filename_peaklist
			}	else {filename_peaklist <- paste(filename_peaklist, ".xls", sep="")}
		}
		##### Database size
		if (isMassPeaksList(peaks_database)) {
			database_size <- length(peaks_database)
		} else if (isMassPeaks(peaks_database)) {
			database_size <- 1
		}
		##### Peak vector
		peak_vector <- character()
		for (i in 1:database_size) {
			peak_vector <- append(peak_vector, peaks_database[[i]]@metaData$file[1])
		}
		########## Database peaklist matrix
		##### Find out the longest peaklist to define the final matrix boundary
		highest_peak_number <- NULL
		for (p in 1:database_size) {
			if (is.null(highest_peak_number) || (length(peaks_database[[p]]@mass) > highest_peak_number)) {
				highest_peak_number <- length(peaks_database[[p]]@mass)
			}
		}
		##### Generate the final matrix
		peaklist_database_matrix <- NULL
		## Fill in the matrix (for each entry)
		for (j in 1:database_size) {
			# Two-row matrix for the database entry (mass, intensity)
			peaklist_database_matrix_entry <- matrix("", ncol=highest_peak_number, nrow=2)
			# Rownames
			rownames(peaklist_database_matrix_entry) <- c(paste(peak_vector[j], "m/z"), paste(peak_vector[j], "intensity"))
			# Mass
			peaklist_database_matrix_entry [1,(1:length(peaks_database[[j]]@mass))] <- peaks_database[[j]]@mass
			# Intensity
			peaklist_database_matrix_entry [2,(1:length(peaks_database[[j]]@intensity))] <- peaks_database[[j]]@intensity
			# Append it to the final matrix
			if (is.null(peaklist_database_matrix)) {
				peaklist_database_matrix <- peaklist_database_matrix_entry
			} else {
				peaklist_database_matrix <- rbind(peaklist_database_matrix, peaklist_database_matrix_entry)
			}
		}
		########## Dump the peaklist matrix
		if (!is.null(peaklist_database_matrix)) {
			if (file_type_export == "csv") {
				write.csv(peaklist_database_matrix, file = filename_peaklist, col.names = FALSE, row.names = TRUE)
			}
			if (file_type_export == "xls" || file_type_export == "xlsx") {
				# Convert it to a data frame
				peaklist_database_matrix <- as.data.frame(peaklist_database_matrix)
				# Generate unique row names
				unique_row_names <- make.names(rownames(peaklist_database_matrix), unique=TRUE)
				rownames(peaklist_database_matrix) <- unique_row_names
				# Export
				writeWorksheetToFile(file = filename_peaklist, data = peaklist_database_matrix, sheet = "Database peaklist", clearSheets = TRUE, header = FALSE, rownames = rownames(peaklist_database_matrix))
			}
			##### Message box
			tkmessageBox(title = "Peaklist dumped", message = "The database peaklist file has been dumped.", icon = "info")
		}
	} else {
		############### If there is no peaklist and spectra to be dumped
		##### Messagebox
		tkmessageBox(title = "Missing peaklist", message = "The spectra in the database or the database peaklist seem to be missing. Run the spectra import and the peak picking and try again.", icon = "warning")
	}
}

##### File name (export)
set_file_name <- function() {
	filename <<- tclvalue(file_name)
	# Add the extension if it is not present in the filename
	if (file_type_export == "csv") {
		if (length(grep(".csv", filename, fixed=TRUE)) == 1) {
			filename <<- filename
		}	else {filename <<- paste (filename, ".csv", sep="")}
	}
	if (file_type_export == "xlsx") {
		if (length(grep(".xlsx", filename, fixed=TRUE)) == 1) {
			filename <<- filename
		}	else {filename <<- paste (filename, ".xlsx", sep="")}
	}
	if (file_type_export == "xls") {
		if (length(grep(".xls", filename, fixed=TRUE)) == 1) {
			filename <<- filename
		}	else {filename <<- paste (filename, ".xls", sep="")}
	}
	# Set the value for displaying purposes
	filename_value <<- filename
	#### Exit the function and put the variable into the R workspace
	#.GlobalEnv$filename <- filename
	return (filename)
}

##### Library
select_database_function <- function() {
	########## Prompt if a folder has to be selected or a single file
	# Catch the value from the popping out menu
	spectra_input_type <- select.list(c("file","folder"), title="Folder or file?")
	if (spectra_input_type == "") {
		spectra_input_type <- "folder"
	}
	if (spectra_input_type == "folder") {
		filepath_database_select <- tkmessageBox(title = "Library", message = "Select the folder for the spectra for the database.\nThe database should be structured like this:\nDatabase folder/Database entry samples/Treatments/Spectra_replicates/Spectrum_coordinates/1/1SLin/Spectrum_data\n\nor\n\nDatabase folder/Classes - Entries/Sample imzML files", icon = "info")
		filepath_database <<- tclvalue(tkchooseDirectory())
		if (!nchar(filepath_database)) {
			tkmessageBox(message = "No folder selected")
		}	else {
			tkmessageBox(message = paste("The directory selected for the database is", filepath_database))
		}
	} else if (spectra_input_type == "file") {
		filepath_database_select <- tkmessageBox(title = "Library", message = "Select a previously dumped database RData file", icon = "info")
		filepath_database <<- tclvalue(tkgetOpenFile(filetypes="{{RData database files} {.RData}}"))
		if (!nchar(filepath_database)) {
			tkmessageBox(message = "No file selected")
		} else {
			tkmessageBox(message = paste("The spectra for the database will be read from:", filepath_database))
		}
	}
	# Set the value for displaying purposes
	filepath_database_value <<- filepath_database
	# Exit the function and put the variable into the R workspace
	#.GlobalEnv$filepath_database <- filepath_database
}

##### Samples
select_samples_function <-function() {
	filepath_test_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be tested.\nThe files should be organised like this:\nSample folder/Samples/Treatments/Spectra_replicates/Spectrum_coordinates/1/1SLin/Spectrum_data\n\nor\n\nSample folder/Sample imzML files", icon = "info")
	filepath_test <<- tclvalue(tkchooseDirectory())
	if (!nchar(filepath_test)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("The sample spectra will be read from:", filepath_test))
	}
	# Set the value for displaying purposes
	filepath_test_value <<- filepath_test
	# Exit the function and put the variable into the R workspace
	#.GlobalEnv$filepath_test <- filepath_test
}

##### Output
browse_output_function <- function() {
	output_folder <<- tclvalue(tkchooseDirectory())
	if (!nchar(output_folder)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("Every file will be saved in", output_folder))
	}
	setwd(output_folder)
	# Exit the function and put the variable into the R workspace
	#.GlobalEnv$output_folder <- output_folder
}

##### Close
quit_function <- function() {
	tkdestroy(window)
}

##### Exit
end_session_function <- function () {
	q(save="no")
}

##### Import the spectra
import_spectra_function <- function() {
	if (!is.null(filepath_database) && !is.null(filepath_test)) {
		# Load the required libraries
		install_and_load_required_packages(c("MALDIquantForeign", "MALDIquant"))
		# Generate the list of spectra (library and test)
		if (spectra_format == "brukerflex" || spectra_format == "xmass") {
			### Load the spectra
			if (!is.null(mass_range)) {
				if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
					## LOAD THE R WORKSPACE (FOR DATABASE)
					# Create a temporary environment
					temporary_environment <- new.env()
					# Load the workspace
					load(filepath_database, envir=temporary_environment)
					# Get the spectra for the database from the workspace
					spectra_database <- get("spectra_database", pos=temporary_environment)
				} else {
					spectra_database <- importBrukerFlex(filepath_database, massRange=mass_range)
					# Write the path inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$path <- filepath_database
					}
					# Read the folder list (database class list)
					database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
					# Store the folder list inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$databaseFolders <- database_folder_list
					}
				}
				spectra_test <- importBrukerFlex(filepath_test, massRange=mass_range)
				# Write the path inside the list
				for (x in 1:length(spectra_test)) {
					spectra_test[[x]]@metaData$path <- filepath_test
				}
			} else {
				if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
					## LOAD THE R WORKSPACE (FOR DATABASE)
					# Create a temporary environment
					temporary_environment <- new.env()
					# Load the workspace
					load(filepath_database, envir=temporary_environment)
					# Get the spectra for the database from the workspace
					spectra_database <- get("spectra_database", pos=temporary_environment)
				} else {
					spectra_database <- importBrukerFlex(filepath_database)
					# Write the path inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$path <- filepath_database
					}
					# Read the folder list (database class list)
					database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
					# Store the folder list inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$databaseFolders <- database_folder_list
					}
				}
				spectra_test <- importBrukerFlex(filepath_test)
				# Write the path inside the list
				for (x in 1:length(spectra_test)) {
					spectra_test[[x]]@metaData$path <- filepath_test
				}
			}
		} else if (spectra_format == "imzml" | spectra_format == "imzML") {
			### Load the spectra
			if (!is.null(mass_range)) {
				if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
					## LOAD THE R WORKSPACE (FOR DATABASE)
					# Create a temporary environment
					temporary_environment <- new.env()
					# Load the workspace
					load(filepath_database, envir=temporary_environment)
					# Get the spectra for the database from the workspace
					spectra_database <- get("spectra_database", pos=temporary_environment)
				} else {
					spectra_database <- importImzMl(filepath_database, massRange=mass_range)
					# Write the path inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$path <- filepath_database
					}
					# Read the folder list (database class list)
					database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
					# Store the folder list inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$databaseFolders <- database_folder_list
					}
				}
				spectra_test <- importImzMl(filepath_test, massRange=mass_range)
				# Write the path inside the list
				for (x in 1:length(spectra_test)) {
					spectra_test[[x]]@metaData$path <- filepath_test
				}
			} else {
				if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
					## LOAD THE R WORKSPACE (FOR DATABASE)
					# Create a temporary environment
					temporary_environment <- new.env()
					# Load the workspace
					load(filepath_database, envir=temporary_environment)
					# Get the spectra for the database from the workspace
					spectra_database <- get("spectra_database", pos=temporary_environment)
				} else {
					spectra_database <- importImzMl(filepath_database)
					# Write the path inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$path <- filepath_database
					}
					# Read the folder list (database class list)
					database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
					# Store the folder list inside the list
					for (x in 1:length(spectra_database)) {
						spectra_database[[x]]@metaData$databaseFolders <- database_folder_list
					}
				}
				spectra_test <- importImzMl(filepath_test)
				# Write the path inside the list
				for (x in 1:length(spectra_test)) {
					spectra_test[[x]]@metaData$path <- filepath_test
				}
			}
		}
		### Average the replicates
		if (average_replicates_in_database == TRUE) {
			### If the database is a RData file, the spectra_database are just there
			if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
				spectra_database <- spectra_database
			} else {
				spectra_database <- average_replicates_by_folder(spectra_database, filepath_database, spectra_format=spectra_format)
			}
		}
		if (average_replicates_in_test == TRUE) {
			spectra_test <- average_replicates_by_folder(spectra_test, filepath_test, spectra_format=spectra_format)
		}
		### Folder lists
		if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
			# Retrieve it from the variable in the spectral list
			database_folder_list <- spectra_database[[1]]@metaData$databaseFolders
		} else {
			database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
		}
		test_folder_list <- dir(filepath_test, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
		### Spectra grouping (class for database)
		if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
			spectra_database <- spectra_database
		} else {
			spectra_database <- group_spectra_class(spectra_database, class_list=database_folder_list, spectra_format=spectra_format, class_in_file_name=TRUE)
		}
		### Preprocessing
		if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
			spectra_database <- spectra_database
		} else {
			spectra_database <- preprocess_spectra(spectra_database, tof_mode=tof_mode, preprocessing_parameters=list(crop_spectra=TRUE, mass_range=NULL, data_transformation=transform_data, transformation_algorithm=transform_data_algorithm, smoothing_algorithm=smoothing_algorithm, smoothing_strength=smoothing_strength, baseline_subtraction_algorithm=baseline_subtraction_algorithm, baseline_subtraction_iterations=baseline_subtraction_iterations, normalisation_algorithm=normalisation_algorithm, normalisation_mass_range=normalisation_mass_range), process_in_packages_of=preprocess_spectra_in_packages_of, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=multicore_processing)
		}
		spectra_test <- preprocess_spectra(spectra_test, tof_mode=tof_mode, preprocessing_parameters=list(crop_spectra=TRUE, mass_range=NULL, data_transformation=transform_data, transformation_algorithm=transform_data_algorithm, smoothing_algorithm=smoothing_algorithm, smoothing_strength=smoothing_strength, baseline_subtraction_algorithm=baseline_subtraction_algorithm, baseline_subtraction_iterations=baseline_subtraction_iterations, normalisation_algorithm=normalisation_algorithm, normalisation_mass_range=normalisation_mass_range), process_in_packages_of=preprocess_spectra_in_packages_of, align_spectra=FALSE, spectra_alignment_method="cubic", multicore_processing=multicore_processing)
		# Number of samples
		if (isMassSpectrumList(spectra_test)) {
			number_of_samples <- length(spectra_test)
		} else if (isMassSpectrum(spectra_test)) {
			number_of_samples <- 1
		}
		# Exit the function and put the variable into the R workspace
		.GlobalEnv$spectra_database <- spectra_database
		.GlobalEnv$spectra_test <- spectra_test
		.GlobalEnv$database_folder_list <- database_folder_list
		.GlobalEnv$test_folder_list <- test_folder_list
		.GlobalEnv$number_of_samples <- number_of_samples
		### Messagebox
		tkmessageBox(title = "Import successful", message = "The spectra have been successfully imported and preprocessed", icon = "info")
	} else if (is.null(filepath_database) || is.null(filepath_test)) {
		### Messagebox
		tkmessageBox(title = "Folder not set", message = "The spectra folder has not been set. Set if before importing the spectra.", icon = "warning")
	}
}

##### Peak picking mode
peak_picking_mode_choice <- function() {
	# Catch the value from the menu
	peak_picking_mode <- select.list(c("all","most intense"), title="Choose")
	# Default
	if (peak_picking_mode == "") {
		peak_picking_mode <- "all"
	}
	# Set the value of the displaying label
	peak_picking_mode_value <- peak_picking_mode
	if (peak_picking_mode_value == "all") {
		peak_picking_mode_value <- "        all        "
	}
	peak_picking_mode_value_label <- tklabel(window, text=peak_picking_mode_value, font = label_font)
	tkgrid(peak_picking_mode_value_label, row=2, column=5)
	# Escape the function
	.GlobalEnv$peak_picking_mode <- peak_picking_mode
	.GlobalEnv$peak_picking_mode_value <- peak_picking_mode_value
}

##### Peak picking function
peak_picking_function <- function() {
	########## RDATA
	if (length(grep(".RData", filepath_database, fixed = TRUE)) > 0) {
		############ Do not run if the spectra have not been imported
		if (!is.null(spectra_test)) {
			###### Get the values
			## Signals to take in most intense peaks
			signals_to_take <- tclvalue(signals_to_take)
			signals_to_take <- as.integer(signals_to_take)
			signals_to_take_value <- as.character(signals_to_take)
			## SNR
			SNR <- tclvalue(SNR)
			SNR <- as.numeric(SNR)
			SNR_value <- as.character(SNR)
			## LOAD THE R WORKSPACE (FOR DATABASE)
			# Create a temporary environment
			temporary_environment <- new.env()
			# Load the workspace
			load(filepath_database, envir=temporary_environment)
			# Get the spectra for the database from the workspace
			peaks_database <- get("peaks_database", pos=temporary_environment)
			if (peak_picking_mode == "most intense") {
				# Peak picking on test
				peaks_test <- most_intense_signals(spectra_test, signals_to_take=signals_to_take)
			} else if (peak_picking_mode == "all") {
				# Peak picking on test
				peaks_test <- peak_picking(spectra_test, peak_picking_algorithm=peak_picking_algorithm, SNR=SNR, tof_mode=tof_mode)
			}
			# Exit the function and put the variable into the R workspace
			.GlobalEnv$peaks_database <- peaks_database
			.GlobalEnv$peaks_test <- peaks_test
			.GlobalEnv$signals_to_take_value <- signals_to_take_value
			.GlobalEnv$SNR_value <- SNR_value
			### Messagebox
			tkmessageBox(title = "Peak picking successful", message = "The peak picking process has been successfully performed", icon = "info")
		} else if (is.null(spectra_test)) {
			### Messagebox
			tkmessageBox(title = "Spectra not imported", message = "The spectra have not been imported yet.\nImport them before performing the peak picking", icon = "warning")
		}
	} else {
		########## FOLDER
		############ Do not run if the spectra have not been imported
		if (!is.null(spectra_database) && !is.null(spectra_test)) {
		###### Get the values
		## Signals to take in most intense peaks
		signals_to_take <- tclvalue(signals_to_take)
		signals_to_take <- as.integer(signals_to_take)
		signals_to_take_value <- as.character(signals_to_take)
		## SNR
		SNR <- tclvalue(SNR)
		SNR <- as.numeric(SNR)
		SNR_value <- as.character(SNR)
		if (peak_picking_mode == "most intense") {
			peaks_database <- most_intense_signals(spectra_database, signals_to_take=signals_to_take, tof_mode=tof_mode)
			peaks_test <- most_intense_signals(spectra_test, signals_to_take=signals_to_take)
		} else if (peak_picking_mode == "all") {
			peaks_database <- peak_picking(spectra_database, peak_picking_algorithm=peak_picking_algorithm, SNR=SNR, tof_mode=tof_mode)
			peaks_test <- peak_picking(spectra_test, peak_picking_algorithm=peak_picking_algorithm, SNR=SNR, tof_mode=tof_mode)
		}
		# Exit the function and put the variable into the R workspace
		.GlobalEnv$peaks_database <- peaks_database
		.GlobalEnv$peaks_test <- peaks_test
		.GlobalEnv$signals_to_take_value <- signals_to_take_value
		.GlobalEnv$SNR_value <- SNR_value
		### Messagebox
		tkmessageBox(title = "Peak picking successful", message = "The peak picking process has been successfully performed", icon = "info")
		} else if (is.null(spectra_database) || is.null(spectra_test)) {
			### Messagebox
			tkmessageBox(title = "Spectra not imported", message = "The spectra have not been imported yet.\nImport them before performing the peak picking", icon = "warning")
		}
	}
}

##### Run the Spectral Typer
run_spectral_typer_function <- function() {
	############ Do not run if the spectra have not been imported or the peaks have not been picked
	if (!is.null(spectra_database) && !is.null(spectra_test) && !is.null(peaks_database) && !is.null(peaks_test)) {
		#### Get the values
		## Intensity correction coefficient
		intensity_correction_coefficient <- tclvalue(intensity_correction_coefficient)
		intensity_correction_coefficient <- as.numeric(intensity_correction_coefficient)
		intensity_correction_coefficient_value <- as.character(intensity_correction_coefficient)
		## Intensity tolerance percent
		intensity_tolerance_percent <- tclvalue(intensity_tolerance_percent)
		intensity_tolerance_percent <- as.numeric(intensity_tolerance_percent)
		intensity_tolerance_percent_value <- as.character(intensity_tolerance_percent)
		## Peaks filtering threshold
		peaks_filtering_threshold_percent <- tclvalue(peaks_filtering_threshold_percent)
		peaks_filtering_threshold_percent <- as.numeric(peaks_filtering_threshold_percent)
		peaks_filtering_threshold_percent_value <- as.character(peaks_filtering_threshold_percent)
		## Low intensity threshold
		intensity_percentage_threshold <- tclvalue(intensity_percentage_threshold)
		intensity_percentage_threshold <- as.numeric(intensity_percentage_threshold)
		intensity_percentage_threshold_value <- as.character(intensity_percentage_threshold)
		### Tolerance (PPM)
		if (tof_mode == "linear") {
			tolerance_ppm <- 2000
		} else if (tof_mode == "reflectron" || tof_mode == "reflector") {
			tolerance_ppm <- 200
		}
		#### Run the function Spectral Typer score calculation
		score_correlation_matrix <- NULL
		score_hca <- NULL
		score_intensity_matrix <- NULL
		score_si_matrix <- NULL
		############### CORRELATION
		if (similarity_criteria == "correlation") {
			score_correlation_matrix <- spectral_typer_score_correlation_matrix(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library=database_folder_list, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, intensity_correction_coefficient=intensity_correction_coefficient, spectra_format=spectra_format, spectra_path_output=spectra_path_output, score_only=score_only, multicore_processing=multicore_processing)
		} else if (similarity_criteria == "hca") {
			############### HIERARCHICAL CLUSTERING ANALYSIS
			score_hca <- spectral_typer_score_hierarchical_distance(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=database_folder_list, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, spectra_path_output=spectra_path_output, score_only=score_only, spectra_format=spectra_format, normalise_distances=TRUE, normalisation_method="sum")
		} else if (similarity_criteria == "similarity index") {
			############### SIMILARITY INDEX
			score_si_matrix <- spectral_typer_score_similarity_index(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library=database_folder_list, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, spectra_format=spectra_format, spectra_path_output=spectra_path_output, score_only=score_only, multicore_processing=multicore_processing)
		} else if (similarity_criteria == "signal intensity") {
			############### INTENSITY
			score_intensity_matrix <- spectral_typer_score_signal_intensity(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=database_folder_list, comparison=signal_intensity_evaluation, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, intensity_tolerance_percent_threshold=intensity_tolerance_percent, spectra_format=spectra_format, spectra_path_output=spectra_path_output, score_only=score_only, number_of_st_dev=1, multicore_processing=multicore_processing)
		}
		### Parameters matrices
		# Parameters vector
		parameters_vector <- c(file_type_export, filepath_database, filepath_test, mass_range_value, tof_mode, spectra_format, preprocess_spectra_in_packages_of_value, peak_picking_mode, signals_to_take_value, intensity_tolerance_percent_value, similarity_criteria, intensity_correction_coefficient_value, SNR_value, peaks_filtering_value, peaks_filtering_threshold_percent_value, low_intensity_peaks_removal_value, intensity_percentage_threshold_value, average_replicates_in_database_value, average_replicates_in_test_value, score_only_value, spectra_path_output_value)
		names(parameters_vector) <- c("File type", "Database folder", "Samples folder", "Mass range", "TOF mode", "Spectra format", "Preprocess spectra in packages of", "Peak picking mode", "Most intense signals taken", "Intensity tolerance percent", "Similarity criteria", "Intensity correction coefficient", "Signal-to-noise ratio", "Peaks filtering", "Filtering threshold percentage", "Low intensity peaks removal", "Intensity threshold percent", "Average replicates in the database", "Average replicates in the samples", "Score only", "Spectra path in the output")
		# Fill in the matrices (the number of columns must be the same for rbind)
		if (!is.null(score_correlation_matrix)) {
			parameters_matrix_correlation <- matrix("", nrow=length(parameters_vector), ncol=ncol(score_correlation_matrix))
			parameters_matrix_correlation[,1] <- cbind(parameters_vector)
			rownames(parameters_matrix_correlation) <- names(parameters_vector)
		}
		if (!is.null(score_hca)) {
			score_hca_matrix <- score_hca$result_matrix
			parameters_matrix_hca <- matrix("", nrow=length(parameters_vector), ncol=ncol(score_hca_matrix))
			parameters_matrix_hca[,1] <- cbind(parameters_vector)
			rownames(parameters_matrix_hca) <- names(parameters_vector)
		}
		if (!is.null(score_si_matrix)) {
			parameters_matrix_si <- matrix("", nrow=length(parameters_vector), ncol=ncol(score_si_matrix))
			parameters_matrix_si[,1] <- cbind(parameters_vector)
			rownames(parameters_matrix_si) <- names(parameters_vector)
		}
		if (!is.null(score_intensity_matrix)) {
			parameters_matrix_intensity <- matrix("", nrow=length(parameters_vector), ncol=ncol(score_intensity_matrix))
			parameters_matrix_intensity[,1] <- cbind(parameters_vector)
			rownames(parameters_matrix_intensity) <- names(parameters_vector)
		}
		#### Exit the function and put the variable into the R workspace
		if (!is.null(score_hca)) {
			score_hca_matrix <- score_hca$result_matrix
			.GlobalEnv$score_hca_matrix_results <- rbind(score_hca_matrix, parameters_matrix_hca)
		}
		if (!is.null(score_si_matrix)) {
			.GlobalEnv$score_si_matrix_results <- rbind(score_si_matrix, parameters_matrix_si)
		}
		if (!is.null(score_intensity_matrix)) {
			.GlobalEnv$score_intensity_matrix_results <- rbind(score_intensity_matrix, parameters_matrix_intensity)
		}
		if (!is.null(score_correlation_matrix)) {
			.GlobalEnv$score_correlation_matrix_results <- rbind(score_correlation_matrix, parameters_matrix_correlation)
		}
		# Save the files (CSV)
		if (file_type_export == "csv") {
			filename <- set_file_name()
			if (!is.null(score_intensity_matrix)) {
			    write.csv(score_intensity_matrix_results, file=paste("int_", filename, sep=""))
			}
			if (!is.null(score_hca)) {
				# Dump the hca plot
				#png(filename="hca.png", width=1900, height=1280)
				#score$score_hca$plots
				scaling_factor <- number_of_samples/25
				if (scaling_factor > 3) {
					scaling_factor <- 3
				}
				ggsave(plot=score_hca$hca_dendrogram, device="png", filename="hca.png", width=12.8, height=7.2, units="in", dpi=300, scale=scaling_factor)
				#savePlot(filename="hca.png", type="png")
				#dev.print(X11, file="hca.png", width=1900, height=1280)
				#dev.off()
			    write.csv(score_hca_matrix_results, file=paste("hca_", filename, sep=""))
			}
			if (!is.null(score_correlation_matrix)) {
			    write.csv(score_correlation_matrix_results, file=paste("corr_", filename, sep=""))
			}
			if (!is.null(score_si_matrix)) {
			    write.csv(score_si_matrix_results, file=paste("si_", filename, sep=""))
			}
		}
		# Save the files (Excel)
		if (file_type_export == "xlsx" || file_type_export == "xls") {
			filename <- set_file_name()
			# Check if PERL is installed by using the built in function, for WriteXLS package
			#perl_installed <- testPerl(verbose=FALSE)
			#if (perl_installed == FALSE) {
				#tkmessageBox(title = "Perl", message = "The software PERL has not been found on the computer. Please install the Perl software from\nhttps://www.perl.org/get.html\nand restart the R session", icon = "warning")
			#}
			if (!is.null(score_intensity_matrix)) {
				# Convert it to a data frame
				score_intensity_matrix_results <- as.data.frame(score_intensity_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_intensity_matrix_results), unique=TRUE)
				rownames(score_intensity_matrix_results) <- unique_row_names
				# Export
				writeWorksheetToFile(file = paste("int_", filename, sep=""), data = score_intensity_matrix_results, sheet = "Scores - Intensity", clearSheets = TRUE, rownames = rownames(score_intensity_matrix_results))

			}
			if (!is.null(score_hca)) {
				# Dump the hca plot
				#png(filename="hca.png", width=1900, height=1280)
				#score$score_hca$plots
				scaling_factor <- number_of_samples/25
				if (scaling_factor > 3) {
					scaling_factor <- 3
				}
				ggsave(plot=score_hca$hca_dendrogram, device="png", filename="hca.png", width=12.8, height=7.2, units="in", dpi=300, scale=scaling_factor)
				#savePlot(filename="hca.png", type="png")
				#dev.print(X11, file="hca.png", width=1900, height=1280)
				#dev.off()
				# Convert it to a data frame
				score_hca_matrix_results <- as.data.frame(score_hca_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_hca_matrix_results), unique=TRUE)
				rownames(score_hca_matrix_results) <- unique_row_names
				# Export
				writeWorksheetToFile(file = paste("hca_", filename, sep=""), data = score_hca_matrix_results, sheet = "Scores - HCA", clearSheets = TRUE, rownames = rownames(score_hca_matrix_results))
			}
			if (!is.null(score_correlation_matrix)) {
				# Convert it to a data frame
				score_correlation_matrix_results <- as.data.frame(score_correlation_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_correlation_matrix_results), unique=TRUE)
				rownames(score_correlation_matrix_results) <- unique_row_names
				# Export
				writeWorksheetToFile(file = paste("corr_", filename, sep=""), data = score_correlation_matrix_results, sheet = "Scores - Correlation", clearSheets = TRUE, rownames = rownames(score_correlation_matrix_results))
			}
			if (!is.null(score_si_matrix)) {
				# Convert it to a data frame
				score_si_matrix_results <- as.data.frame(score_si_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_si_matrix_results), unique=TRUE)
				rownames(score_si_matrix_results) <- unique_row_names
				# Export
				writeWorksheetToFile(file = paste("si_", filename, sep=""), data = score_si_matrix_results, sheet = "Scores - Similarity Index", clearSheets = TRUE, rownames = rownames(score_si_matrix_results))
			}
		}
		### Messagebox
		tkmessageBox(title = "Done!", message = "The file(s) have been dumped\n\nLegend:\nF: Fit\nRF: Retrofit\nCorr: intensity Pearson's correlation coefficient\nIntMtch: signal intensity matching\nsl: slope of the regression curve\nns: number of signals\nSI: Similarity Index\n\n\nFit = number of sample-database matching signals / number of signals in the sample\nRetrofit = number of database-sample matching signals / number of signals in the database entry", icon = "info")
	} else if (is.null(spectra_database) || is.null(spectra_test) || is.null(peaks_database) || is.null(peaks_test)) {
		### Messagebox
		tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the spectra have been imported and the peak picking process has been performed", icon = "warning")
	}
}

##### Similarity criteria
similarity_criteria_choice <- function() {
	# Catch the value from the menu
	similarity_criteria <- select.list(c("correlation","hca","signal intensity","similarity index"), title="Choose")
	# Default
	if (similarity_criteria == "") {
		similarity_criteria <- "correlation"
	}
	# Set the value of the displaying label
	similarity_criteria_value <- similarity_criteria
	if (similarity_criteria_value == "hca") {
		similarity_criteria_value <- "         hca         "
	} else if (similarity_criteria_value == "correlation") {
		similarity_criteria_value <- "   correlation   "
	}
	similarity_criteria_value_label <- tklabel(window, text=similarity_criteria_value, font = label_font)
	tkgrid(similarity_criteria_value_label, row=3, column=3)
	# Escape the function
	.GlobalEnv$similarity_criteria <- similarity_criteria
	.GlobalEnv$similarity_criteria_value <- similarity_criteria_value
}

##### Signal intensity evaluation
signal_intensity_evaluation_choice <- function() {
	# Catch the value from the menu
	signal_intensity_evaluation <- select.list(c("intensity percentage"), title="Choose")
	# Default
	if (signal_intensity_evaluation == "") {
		signal_intensity_evaluation <- "intensity percentage"
	}
	# Escape the function
	.GlobalEnv$signal_intensity_evaluation <- signal_intensity_evaluation
	# Set the value of the displaying label
	signal_intensity_evaluation_value_label <- tklabel(window, text=signal_intensity_evaluation, font = label_font)
	tkgrid(signal_intensity_evaluation_value_label, row=4, column=6)
}

##### Peak picking algorithm
peak_picking_algorithm_choice <- function() {
	# Catch the value from the menu
	peak_picking_algorithm <- select.list(c("MAD","SuperSmoother"), title="Choose")
	# Default
	if (peak_picking_algorithm == "") {
		peak_picking_algorithm <- "MAD"
	}
	# Set the value of the displaying label
	peak_picking_algorithm_value <- peak_picking_algorithm
	if (peak_picking_algorithm_value == "MAD") {
		peak_picking_algorithm_value <- "          MAD          "
	} else if (peak_picking_algorithm_value == "SuperSmoother") {
		peak_picking_algorithm_value <- "Super Smoother"
	}
	peak_picking_algorithm_value_label <- tklabel(window, text=peak_picking_algorithm_value, font = label_font)
	tkgrid(peak_picking_algorithm_value_label, row=3, column=5)
	# Escape the function
	.GlobalEnv$peak_picking_algorithm <- peak_picking_algorithm
	.GlobalEnv$peak_picking_algorithm_value <- peak_picking_algorithm_value
}

##### Peaks filtering
peaks_filtering_choice <- function() {
	# Catch the value from the menu
	peaks_filtering <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (peaks_filtering == "YES") {
		peaks_filtering <- TRUE
	}
	if (peaks_filtering == "NO" || peaks_filtering == "") {
		peaks_filtering <- FALSE
	}
	# Set the value of the displaying label
	if (peaks_filtering == TRUE) {
		peaks_filtering_value <- "YES"
	} else {
		peaks_filtering_value <- "NO"
	}
	peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value, font = label_font)
	tkgrid(peaks_filtering_value_label, row=7, column=3)
	# Escape the function
	.GlobalEnv$peaks_filtering <- peaks_filtering
	.GlobalEnv$peaks_filtering_value <- peaks_filtering_value
}

##### Low intensity peaks removal
low_intensity_peaks_removal_choice <- function() {
	# Catch the value from the menu
	low_intensity_peaks_removal <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (low_intensity_peaks_removal == "" || low_intensity_peaks_removal == "NO") {
		low_intensity_peaks_removal <- FALSE
	}
	if (low_intensity_peaks_removal == "YES") {
		low_intensity_peaks_removal <- TRUE
	}
	# Set the value of the displaying label
	if (low_intensity_peaks_removal == TRUE) {
		low_intensity_peaks_removal_value <- "YES"
	} else {
		low_intensity_peaks_removal_value <- "NO"
	}
	low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value, font = label_font)
	tkgrid(low_intensity_peaks_removal_value_label, row=8, column=3)
	# Escape the function
	.GlobalEnv$low_intensity_peaks_removal <- low_intensity_peaks_removal
	.GlobalEnv$low_intensity_peaks_removal_value <- low_intensity_peaks_removal_value
}

##### Low intensity peaks removal Method
intensity_threshold_method_choice <- function() {
	# Catch the value from the menu
	intensity_threshold_method <- select.list(c("whole","element-wise"), title="Choose")
	# Default
	if (intensity_threshold_method == "") {
		intensity_threshold_method <- "element-wise"
	}
	# Set the value of the displaying label
	intensity_threshold_method_value <- intensity_threshold_method
	if (intensity_threshold_method_value == "whole") {
		intensity_threshold_method_value <- "     whole     "
	}
	intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value, font = label_font)
	tkgrid(intensity_threshold_method_value_label, row=9, column=3)
	# Escape the function
	.GlobalEnv$intensity_threshold_method <- intensity_threshold_method
	.GlobalEnv$intensity_threshold_method_value <- intensity_threshold_method_value
}

##### Multicore processing
multicore_processing_choice <- function() {
	# Catch the value from the menu
	multicore_processing <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (multicore_processing == "YES" || multicore_processing == "") {
		multicore_processing <- TRUE
	}
	if (multicore_processing == "NO") {
		multicore_processing <- FALSE
	}
	# Set the value of the displaying label
	if (multicore_processing == TRUE) {
		multicore_processing_value <- "YES"
	} else {
		multicore_processing_value <- "NO"
	}
	multicore_processing_value_label <- tklabel(window, text=multicore_processing_value, font = label_font)
	tkgrid(multicore_processing_value_label, row=9, column=5)
	# Escape the function
	.GlobalEnv$multicore_processing <- multicore_processing
	.GlobalEnv$multicore_processing_value <- multicore_processing_value
}


##### Average replicates in database
average_replicates_in_database_choice <- function() {
	# Catch the value from the menu
	average_replicates_in_database <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (average_replicates_in_database == "" || average_replicates_in_database == "NO") {
		average_replicates_in_database <- FALSE
	}
	if (average_replicates_in_database == "YES") {
		average_replicates_in_database <- TRUE
	}
	# Set the value of the displaying label
	if (average_replicates_in_database == TRUE) {
		average_replicates_in_database_value <- "YES"
	} else {
		average_replicates_in_database_value <- "NO"
	}
	average_replicates_in_database_value_label <- tklabel(window, text=average_replicates_in_database_value, font = label_font)
	tkgrid(average_replicates_in_database_value_label, row=10, column=3)
	# Escape the function
	.GlobalEnv$average_replicates_in_database <- average_replicates_in_database
	.GlobalEnv$average_replicates_in_database_value <- average_replicates_in_database_value
}

##### Average replicates in test
average_replicates_in_test_choice <- function() {
	# Catch the value from the menu
	average_replicates_in_test <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (average_replicates_in_test == "" || average_replicates_in_test == "NO") {
		average_replicates_in_test <- FALSE
	}
	if (average_replicates_in_test == "YES") {
		average_replicates_in_test <- TRUE
	}
	# Set the value of the displaying label
	if (average_replicates_in_test == TRUE) {
		average_replicates_in_test_value <- "YES"
	} else {
		average_replicates_in_test_value <- "NO"
	}
	average_replicates_in_test_value_label <- tklabel(window, text=average_replicates_in_test_value, font = label_font)
	tkgrid(average_replicates_in_test_value_label, row=10, column=6)
	# Escape the function
	.GlobalEnv$average_replicates_in_test <- average_replicates_in_test
	.GlobalEnv$average_replicates_in_test_value <- average_replicates_in_test_value
}

##### Score only (or all the score components)
score_only_choice <- function() {
	# Catch the value from the menu
	score_only <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (score_only == "" || score_only == "NO") {
		score_only <- FALSE
	}
	if (score_only == "YES") {
		score_only <- TRUE
	}
	# Set the value of the displaying label
	if (score_only == TRUE) {
		score_only_value <- "YES"
	} else {
		score_only_value <- "NO"
	}
	score_only_value_label <- tklabel(window, text=score_only_value, font = label_font)
	tkgrid(score_only_value_label, row=11, column=3)
	# Escape the function
	.GlobalEnv$score_only <- score_only
	.GlobalEnv$score_only_value <- score_only_value
}

##### Spectra path in the output
spectra_path_output_choice <- function() {
	# Catch the value from the menu
	spectra_path_output <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (spectra_path_output == "" || spectra_path_output == "YES") {
		spectra_path_output <- TRUE
	}
	if (spectra_path_output == "NO") {
		spectra_path_output <- FALSE
	}
	# Set the value of the displaying label
	if (spectra_path_output == TRUE) {
		spectra_path_output_value <- "YES"
	} else {
		spectra_path_output_value <- "NO"
	}
	spectra_path_output_value_label <- tklabel(window, text=spectra_path_output_value, font = label_font)
	tkgrid(spectra_path_output_value_label, row=11, column=6)
	# Escape the function
	.GlobalEnv$spectra_path_output <- spectra_path_output
	.GlobalEnv$spectra_path_output_value <- spectra_path_output_value
}

##### File format
spectra_format_choice <- function() {
	# Catch the value from the menu
	spectra_format <- select.list(c("imzML","Xmass"), title="Choose")
	# Default
	if (spectra_format == "" || spectra_format == "Xmass") {
		spectra_format <- "xmass"
		spectra_format_value <- "Xmass"
	}
	if (spectra_format == "imzML") {
		spectra_format <- "imzml"
		spectra_format_value <- "imzML"
	}
	# Escape the function
	.GlobalEnv$spectra_format <- spectra_format
	.GlobalEnv$spectra_format_value <- spectra_format_value
	# Set the value of the displaying label
	spectra_format_value_label <- tklabel(window, text=spectra_format_value, font = label_font)
	tkgrid(spectra_format_value_label, row=12, column=3)
}










##################################################################### WINDOW GUI

########## List of variables, whose values are taken from the entries in the GUI
SNR <- tclVar("")
intensity_correction_coefficient <- tclVar("")
peaks_filtering_threshold_percent <- tclVar("")
intensity_percentage_threshold <- tclVar("")
signals_to_take <- tclVar("")
file_name <- tclVar("")
intensity_tolerance_percent <- tclVar("")



######################## GUI

### FONTS
# Get system info (Platform - Release - Version (- Linux Distro))
system_os = Sys.info()[1]
os_release = Sys.info()[2]
os_version = Sys.info()[3]
# Windows
if (system_os == "Windows") {
	# Define the fonts
	garamond_24_bold = tkfont.create(family = "Garamond", size = 24, weight = "bold")
	garamond_12_normal = tkfont.create(family = "Garamond", size = 12, weight = "normal")
	arial_24_bold = tkfont.create(family = "Arial", size = 24, weight = "bold")
	arial_12_normal = tkfont.create(family = "Arial", size = 12, weight = "normal")
	trebuchet_24_bold = tkfont.create(family = "Trebuchet MS", size = 24, weight = "bold")
	trebuchet_11_normal = tkfont.create(family = "Trebuchet MS", size = 11, weight = "normal")
	trebuchet_11_bold = tkfont.create(family = "Trebuchet MS", size = 11, weight = "bold")
	# Use them in the GUI
	title_font = trebuchet_24_bold
	label_font = trebuchet_11_normal
	entry_font = trebuchet_11_normal
	button_font = trebuchet_11_bold
} else if (system_os == "Linux") {
	# Linux
	# Ubuntu
	if (os_release == "Ubuntu") {
		# Define the fonts
		ubuntu_24_bold = tkfont.create(family = "Ubuntu", size = 24, weight = "bold")
		ubuntu_12_normal = tkfont.create(family = "Ubuntu", size = 12, weight = "normal")
		ubuntu_12_bold = tkfont.create(family = "Ubuntu", size = 12, weight = "bold")
		# Use them in the GUI
		title_font = ubuntu_24_bold
		label_font = ubuntu_12_normal
		entry_font = ubuntu_12_normal
		button_font = ubuntu_12_bold
	} else if (os_release == "Fedora") {
		# Fedora
		# Define the fonts
		cantarell_24_bold = tkfont.create(family = "Cantarell", size = 24, weight = "bold")
		cantarell_12_normal = tkfont.create(family = "Cantarell", size = 12, weight = "normal")
		cantarell_12_bold = tkfont.create(family = "Cantarell", size = 12, weight = "bold")
		# Use them in the GUI
		title_font = cantarell_24_bold
		label_font = cantarell_12_normal
		entry_font = cantarell_12_normal
		button_font = cantarell_12_bold
	} else {
		# Other linux distros
		# Define the fonts
		liberation_24_bold = tkfont.create(family = "Liberation Sans", size = 24, weight = "bold")
		liberation_12_normal = tkfont.create(family = "Liberation Sans", size = 12, weight = "normal")
		liberation_12_bold = tkfont.create(family = "Liberation Sans", size = 12, weight = "bold")
		# Use them in the GUI
		title_font = liberation_24_bold
		label_font = liberation_12_normal
		entry_font = liberation_12_normal
		button_font = liberation_12_bold
	}
} else if (system_os == "Darwin") {
	# macOS
	# Define the fonts
	helvetica_24_bold = tkfont.create(family = "Helvetica", size = 24, weight = "bold")
	helvetica_16_normal = tkfont.create(family = "Helvetica", size = 16, weight = "normal") 
	helvetica_16_bold = tkfont.create(family = "Helvetica", size = 16, weight = "bold")
	# Use them in the GUI
	title_font = helvetica_24_bold
	label_font = helvetica_16_normal
	entry_font = helvetica_16_normal
	button_font = helvetica_16_bold
}



### Initial Messagebox
#tkmessageBox(title = "Before starting", message = "The library should be structured like this: Main folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/Spectrum_data", icon = "info")

# The "area" where we will put our input lines
window <- tktoplevel()
tktitle(window) <- "Spectral Typer"
#### Browse
# Library
select_database_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data", font = label_font)
select_database_button <- tkbutton(window, text="Browse database folder", command=select_database_function, font = button_font)
# Samples
select_samples_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data", font = label_font)
select_samples_button <- tkbutton(window, text="Browse samples folder", command=select_samples_function, font = button_font)
# Output
select_output_label <- tklabel(window, text="Select the folder where to save all the outputs", font = label_font)
browse_output_button <- tkbutton(window, text="Browse output folder", command=browse_output_function, font = button_font)
#### Entries
# Similarity criteria
similarity_criteria_label <- tklabel(window, text="Similarity criteria", font = label_font)
similarity_criteria_entry <- tkbutton(window, text="Choose similarity\ncriteria", command=similarity_criteria_choice, font = button_font)
# Intensity correction coefficient
intensity_correction_coefficient_label <- tklabel(window, text="Intensity correction coefficient\n(0: discard the intensities,\n1: unweighted correlation)", font = label_font)
intensity_correction_coefficient_entry <- tkentry(window, width=10, textvariable=intensity_correction_coefficient, font = entry_font)
tkinsert(intensity_correction_coefficient_entry, "end", "1")
# Intensty tolerance percent
intensity_tolerance_percent_label <- tklabel(window, text="Intensity tolerance percent\n(if 'signal intensity' is selected)", font = label_font)
intensity_tolerance_percent_entry <- tkentry(window, width=10, textvariable=intensity_tolerance_percent, font = entry_font)
tkinsert(intensity_tolerance_percent_entry, "end", "80")
# Signal intensity evaluation
signal_intensity_evaluation_label <- tklabel(window, text="Signal intensity evaluation", font = label_font)
signal_intensity_evaluation_entry <- tkbutton(window, text="Choose signal intensity\nevaluation method", command=signal_intensity_evaluation_choice, font = button_font)
# Peak picking mode
peak_picking_mode_label <- tklabel(window, text="Peak picking mode", font = label_font)
peak_picking_mode_entry <- tkbutton(window, text="Choose peak picking\nmode", command=peak_picking_mode_choice, font = button_font)
# Peak picking mode
peak_picking_algorithm_button <- tkbutton(window, text="Choose peak picking\nalgorithm", command=peak_picking_algorithm_choice, font = button_font)
peak_picking_algorithm_value_label <- tklabel(window, text=peak_picking_algorithm_value, font = label_font)
# Signals to take
signals_to_take_label <- tklabel(window, text="Most intense signals to take\n(if 'most intense' is selected)", font = label_font)
signals_to_take_entry <- tkentry(window, width=10, textvariable=signals_to_take, font = entry_font)
tkinsert(signals_to_take_entry, "end", "25")
# SNR
SNR_label <- tklabel(window, text="Signal-to-noise ratio", font = label_font)
SNR_entry <- tkentry(window, width=10, textvariable=SNR, font = entry_font)
tkinsert(SNR_entry, "end", "5")
# Peaks filtering
peaks_filtering_label <- tklabel(window, text="Peaks filtering", font = label_font)
peaks_filtering_entry <- tkbutton(window, text="Choose peak filtering", command=peaks_filtering_choice, font = button_font)
# Peaks filtering threshold
peaks_filtering_threshold_percent_label <- tklabel(window, text="Peaks filtering threshold frequency percentage", font = label_font)
peaks_filtering_threshold_percent_entry <- tkentry(window, width=10, textvariable=peaks_filtering_threshold_percent, font = entry_font)
tkinsert(peaks_filtering_threshold_percent_entry, "end", "25")
# Low intensity peaks removal
low_intensity_peaks_removal_label <- tklabel(window, text="Low intensity peaks removal", font = label_font)
low_intensity_peaks_removal_entry <- tkbutton(window, text="Choose low intensity\npeaks removal", command=low_intensity_peaks_removal_choice, font = button_font)
# Intensity percentage threshold
intensity_percentage_threshold_label <- tklabel(window, text="Intensity percentage threshold", font = label_font)
intensity_percentage_threshold_entry <- tkentry(window, width=10, textvariable=intensity_percentage_threshold, font = entry_font)
tkinsert(intensity_percentage_threshold_entry, "end", "0.1")
# Intensiry percentage theshold method
intensity_threshold_method_label <- tklabel(window, text="Intensity threshold method", font = label_font)
intensity_threshold_method_entry <- tkbutton(window, text="Choose the method for\nthe intensity threshold", command=intensity_threshold_method_choice, font = button_font)
# Average replicates in database
average_replicates_in_database_label <- tklabel(window, text="Average replicates in the database", font = label_font)
average_replicates_in_database_entry <- tkbutton(window, text="Choose average replicates\nin the database", command=average_replicates_in_database_choice, font = button_font)
# Average replicates in samples
average_replicates_in_test_label <- tklabel(window, text="Average replicates in the samples", font = label_font)
average_replicates_in_test_entry <- tkbutton(window, text="Choose average replicates\nin the samples", command=average_replicates_in_test_choice, font = button_font)
# Score only
score_only_label <- tklabel(window, text="Score only\n('NO', all the score\ncomponents are displayed)", font = label_font)
score_only_entry <- tkbutton(window, text="Choose", command=score_only_choice, font = button_font)
# Spectra path output
spectra_path_output_label <- tklabel(window, text="Spectra path in the output", font = label_font)
spectra_path_output_entry <- tkbutton(window, text="Choose to display\nthe spectra path", command=spectra_path_output_choice, font = button_font)
# File format
spectra_format_label <- tklabel(window, text="Select the spectra format", font = label_font)
spectra_format_entry <- tkbutton(window, text="Choose the spectra format", command=spectra_format_choice, font = button_font)
# File type export
file_type_export_label <- tklabel(window, text="Select the format\nof the exported file", font = label_font)
file_type_export_entry <- tkbutton(window, text="Choose the file type", command=file_type_export_choice, font = button_font)
#### Close
exit_label <- tklabel(window, text="Exit", font = label_font)
quit_button <- tkbutton(window, text="Quit", command=quit_function, font = button_font)
# End session
#end_session_label <- tklabel(window, text="Quit")
end_session_button <- tkbutton(window, text="QUIT", command=end_session_function, font = button_font)
# Multicore
multicore_processing_button <- tkbutton(window, text="ALLOW PARALLEL\nPROCESSING", command=multicore_processing_choice, font = button_font)
# Import the spectra
import_spectra_button <- tkbutton(window, text="SPECTRA IMPORT AND\nPREPROCESSING", command=import_spectra_function, font = button_font)
# Peak picking
peak_picking_button <- tkbutton(window, text="PEAK PICKING", command=peak_picking_function, font = button_font)
# Run the Spectral typer!
run_spectral_typer_button <- tkbutton(window, text="RUN THE SPECTRAL TYPER", command=run_spectral_typer_function, font = button_font)
# Spectra preprocessing button
spectra_preprocessing_button <- tkbutton(window, text="SPECTRA PREPROCESSING\nPARAMETERS", command=preprocessing_window_function, font = button_font)
# Set the file name
set_file_name_label <- tklabel(window, text="<-- Set the file name", font = label_font)
set_file_name_entry <- tkentry(window, width=30, textvariable=file_name, font = entry_font)
tkinsert(set_file_name_entry, "end", "Score")
# Dump the database peaklist
database_peaklist_dump_button <- tkbutton(window, text="Dump the database", command=database_dump_function, font = button_font)




#### Displaying labels
file_type_export_value_label <- tklabel(window, text=file_type_export, font = label_font)
similarity_criteria_value_label <- tklabel(window, text=similarity_criteria_value, font = label_font)
signal_intensity_evaluation_value_label <- tklabel(window, text=signal_intensity_evaluation, font = label_font)
peak_picking_mode_value_label <- tklabel(window, text=peak_picking_mode_value, font = label_font)
peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value, font = label_font)
low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value, font = label_font)
intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value, font = label_font)
average_replicates_in_database_value_label <- tklabel(window, text=average_replicates_in_database_value, font = label_font)
average_replicates_in_test_value_label <- tklabel(window, text=average_replicates_in_test_value, font = label_font)
score_only_value_label <- tklabel(window, text=score_only_value, font = label_font)
spectra_path_output_value_label <- tklabel(window, text=spectra_path_output_value, font = label_font)
spectra_format_value_label <- tklabel(window, text=spectra_format_value, font = label_font)
multicore_processing_value_label <- tklabel(window, text=multicore_processing_value, font = label_font)




#### Geometry manager
# Scrollbar
#window_scrollbar <- tkscrollbar(window, command=function(...)tkyview(window,...))
# tkgrid
tkgrid(set_file_name_label, row=1, column=5)
tkgrid(set_file_name_entry, row=1, column=4)
#tkgrid(select_database_label, row=1, column=1)
tkgrid(select_database_button, row=1, column=1)
#tkgrid(select_samples_label, row=1, column=3)
tkgrid(select_samples_button, row=1, column=2)
#tkgrid(select_output_label, row=2, column=3)
tkgrid(browse_output_button, row=1, column=3)
tkgrid(similarity_criteria_label, row=3, column=1)
tkgrid(similarity_criteria_entry, row=3, column=2)
tkgrid(similarity_criteria_value_label, row=3, column=3)
tkgrid(peak_picking_algorithm_button, row=3, column=4)
tkgrid(peak_picking_algorithm_value_label, row=3, column=5)
tkgrid(intensity_correction_coefficient_label, row=4, column=1)
tkgrid(intensity_correction_coefficient_entry, row=4, column=2)
tkgrid(signal_intensity_evaluation_label, row=4, column=4)
tkgrid(signal_intensity_evaluation_entry, row=4, column=5)
tkgrid(signal_intensity_evaluation_value_label, row=4, column=6)
tkgrid(peak_picking_mode_label, row=5, column=1)
tkgrid(peak_picking_mode_entry, row=5, column=2)
tkgrid(peak_picking_mode_value_label, row=5, column=3)
tkgrid(signals_to_take_label, row=5, column=4)
tkgrid(signals_to_take_entry, row=5, column=5)
tkgrid(SNR_label, row=6, column=2)
tkgrid(SNR_entry, row=6, column=3)
tkgrid(peaks_filtering_label, row=7, column=1)
tkgrid(peaks_filtering_entry, row=7, column=2)
tkgrid(peaks_filtering_value_label, row=7, column=3)
tkgrid(peaks_filtering_threshold_percent_label, row=7, column=4)
tkgrid(peaks_filtering_threshold_percent_entry, row=7, column=5)
tkgrid(low_intensity_peaks_removal_label, row=8, column=1)
tkgrid(low_intensity_peaks_removal_entry, row=8, column=2)
tkgrid(low_intensity_peaks_removal_value_label, row=8, column=3)
tkgrid(intensity_threshold_method_label, row=9, column=1)
tkgrid(intensity_threshold_method_entry, row=9, column=2)
tkgrid(intensity_threshold_method_value_label, row=9, column=3)
tkgrid(intensity_percentage_threshold_label, row=8, column=4)
tkgrid(intensity_percentage_threshold_entry, row=8, column=5)
tkgrid(average_replicates_in_database_label, row=10, column=1)
tkgrid(average_replicates_in_database_entry, row=10, column=2)
tkgrid(average_replicates_in_database_value_label, row=10, column=3)
tkgrid(average_replicates_in_test_label, row=10, column=4)
tkgrid(average_replicates_in_test_entry, row=10, column=5)
tkgrid(average_replicates_in_test_value_label, row=10, column=6)
tkgrid(multicore_processing_button, row=9, column=4)
tkgrid(multicore_processing_value_label, row=9, column=5)
tkgrid(score_only_label, row=11, column=1)
tkgrid(score_only_entry, row=11, column=2)
tkgrid(score_only_value_label, row=11, column=3)
tkgrid(spectra_path_output_label, row=11, column=4)
tkgrid(spectra_path_output_entry, row=11, column=5)
tkgrid(spectra_path_output_value_label, row=11, column=6)
tkgrid(spectra_format_label, row=12, column=1)
tkgrid(spectra_format_entry, row=12, column=2)
tkgrid(spectra_format_value_label, row=12, column=3)
tkgrid(intensity_tolerance_percent_label, row=6, column=4)
tkgrid(intensity_tolerance_percent_entry, row=6, column=5)
tkgrid(file_type_export_label, row=12, column=4)
tkgrid(file_type_export_entry, row=12, column=5)
tkgrid(file_type_export_value_label, row=12, column=6)
tkgrid(import_spectra_button, row=14, column=2)
tkgrid(peak_picking_button, row=14, column=3)
tkgrid(run_spectral_typer_button, row=14, column=4)
tkgrid(spectra_preprocessing_button, row=14, column=5)
tkgrid(database_peaklist_dump_button, row=15, column=3)
#tkgrid(exit_label, row=15, column=1)
#tkgrid(quit_button, row=15, column=2)
tkgrid(end_session_button, row=15, column=4)
#window_scrollbar
