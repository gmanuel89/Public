###INSTALL THE REQUIRED PACKAGES
##### MEMORY LIMIT
# Set the memory (RAM) limit to 100GB (100000MB)
memory.limit (100000)

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
	library (required_packages[i], character.only=TRUE)
}
}

install_and_load_required_packages("beepr")




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






















######################################### AVERAGE SIGNAL NUMBER
# Average number of signals with that S/N in the peaklist set
average_signal_number_with_selected_SNR <- function (spectra, SNR=5) {
# Load the required libraries
install_and_load_required_packages ("MALDIquant")
#
peaks <- detectPeaks(spectra, method="MAD", SNR=SNR)
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
# Number of spectra with more than a defined number of signals with that SNR
number_of_representative_spectra <- function (spectra, SNR=15, signal_number_threshold=20) {
# Load the required libraries
install_and_load_required_packages ("MALDIquant")
#
counter <- 0
## Peak picking
peaks <- detectPeaks(spectra, method="MAD", SNR=SNR)
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
# Average number of spectra with more than a defined number of signals with that SNR
average_number_of_representative_spectra <- function (filepath, SNR=15, signal_number_threshold=20, tof_mode="linear", file_format="imzml") {
# Load the required libraries
install_and_load_required_packages (c("MALDIquant", "MALDIquantForeign"))
#
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
	peaks <- detectPeaks(spectra, method="MAD", SNR=SNR)
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
# Remove all the bad spectra, keep only the best ones (Number of signals with a certain SNR)
filter_spectra_signal_number <- function (spectra, SNR_filter=15, signal_number_threshold=15) {
# Load the required libraries
install_and_load_required_packages ("parallel")
#
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
##################################################### FILTERING FUNCTION
spectra_filtering_function <- function (spectra) {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR_filter)
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
		skyline_spectrum  <- generate_skyline_spectrum(spectra_list)
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
		skyline_spectrum  <- generate_skyline_spectrum(spectra_list)
		# Add the average spectrum to another final list of average spectra
		spectra_grouped <- append(spectra_grouped, skyline_spectrum )
	}
}
return (spectra_grouped)
}





########################################################################










########################################################################





########################################## SPECTRA GROUPING (CLASSES) - SKYLINE
# Obtain one spectrum per class (average)
group_spectra_class_skyline <- function (spectra, class_list, file_format="imzml") {
class_list <- sort(class_list)
if (file_format == "imzml" | file_format == "imzML") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name(spectra, class_list=class_list)
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
		skyline_spectrum  <- generate_skyline_spectrum(spectra_list)
		# Add the average spectrum to another final list of average spectra
		class_spectra_grouped <- append(class_spectra_grouped, skyline_spectrum )
	}
}
if (file_format == "brukerflex") {
	## REPLACE THE filepath PARAMETER FOR EACH SPECTRUM WITH THE CLASS
	spectra <- replace_class_name(spectra, class_list=class_list)
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
		skyline_spectrum  <- generate_skyline_spectrum(spectra_list)
		# Add the average spectrum to another final list of average spectra
		class_spectra_grouped <- append(class_spectra_grouped, skyline_spectrum )
	}
}
return (class_spectra_grouped)
}





################################################################################











########################################################################





###################################### NO HEMOGLOBIN
no_hemoglobin <- function (spectra, hemoglobin_mass=15400, tolerance_ppm=2000, average_spectrum_evaluation=FALSE) {
# Load the required libraries
install_and_load_required_packages (c("parallel", "MALDIquant"))
#
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
#################################################### FILTERING FUNCTION
hemoglobin_removal_function <- function (spectra, hemoglobin_mass, tolerance_ppm) {
## Peak picking
peaks <- detectPeaks(spectra, method="MAD", SNR=3)
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
	average_spectrum <- averageMassSpectra(spectra_no_hemoglobin, method="mean")
	## Peak picking on the AVG
	peaks_average <- detectPeaks(average_spectrum, method="MAD", SNR=3)
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










########################################################################










########################################################################











#################### SIGNAL COUNTER IN THE AVERAGE SPECTRUM
##### Count signals in a mass range
countSignalsRangeAverageSpectrum <- function (spectra, SNR=5, lower=4000, interval_width=1000, max=20000) {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
average_spectrum <- averageMassSpectra (spectra, method="mean")
average_spectrum <- removeBaseline (average_spectrum, method="TopHat")
# Peak picking on the average_spectrum
peaks_average <- detectPeaks (average_spectrum, method="MAD", SNR=SNR)
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
count_signals_range_spectra <- function (spectra, SNR=5, lower=4000, interval_width=1000, max=20000) {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
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
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
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





########################################################################










########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N MOST INTENSE
most_intense_signals <- function (spectra, signals_to_take=20) {
# Load the required libraries
install_and_load_required_packages(c("parallel", "MALDIquant"))
#
####################################################### PICKING FUNCTION
picking_function <- function (peaks, signals_to_take) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, SNR = peaks@SNR)
	# Sort the dataframe according to the SNR
	peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$SNR),]
	# Select only the first most intense signals
	selected_signals <- peaks_data_frame [1:signals_to_take,]
	# Sort the dataframe back according to mass
	selected_signals <- selected_signals [order(selected_signals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selected_signals$mass
	peaks@intensity <- selected_signals$intensity
	peaks@SNR <- selected_signals$SNR
	return (peaks)
}
########################################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=3)
most_intense_peaks <- mclapply(peaks, FUN= function (peaks) picking_function (peaks, signals_to_take=signals_to_take))
return (most_intense_peaks)
}






########################################################################
###### SORT THE SIGNALS AND TAKE ONLY THE N LESS INTENSE
less_intense_signals <- function (spectra, base_SNR=2, signals_to_take=20) {
# Load the required libraries
install_and_load_required_packages(c("parallel", "MALDIquant"))
#
############################################################
picking_function <- function (peaks, signals_to_take) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, SNR = peaks@SNR)
	# Sort the dataframe according to the SNR
	peaks_data_frame <- peaks_data_frame [order(peaks_data_frame$SNR),]
	# Select only the first most intense signals
	selected_signals <- peaks_data_frame [1:signals_to_take,]
	# Sort the dataframe back according to mass
	selected_signals <- selected_signals [order(selected_signals$mass),]
	# Put these signals back into the peaklist
	peaks@mass <- selected_signals$mass
	peaks@intensity <- selected_signals$intensity
	peaks@SNR <- selected_signals$SNR
	return (peaks)
}
###########################################################
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=base_SNR)
less_intense_peaks <- mclapply(peaks, FUN= function (peaks) picking_function (peaks, signals_to_take=signals_to_take))
return (less_intense_peaks)
}





########################################################################
####### LESS AND MOST INTENSE SIGNALS
highest_and_lowest_intensity_signals <- function (spectra, base_SNR=3, signals_to_take) {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=base_SNR)
# For each peaklist...
for (p in 1:length(peaks)) {
	# Create a dataframe with mass and intensity
	peaks_data_frame <- data.frame (mass = peaks[[p]]@mass, intensity = peaks[[p]]@intensity, SNR = peaks[[p]]@SNR)
	# Sort the dataframe according to the SNR
	most_intense_peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$SNR),]
	less_intense_peaks_data_frame <- peaks_data_frame [order(peaks_data_frame$SNR),]
	# Select only the first most intense signals
	selected_signals_most_intense <- most_intense_peaks_data_frame [1:(signals_to_take/2),]
	selected_signals_less_intense <- less_intense_peaks_data_frame [1:(signals_to_take/2),]
	# Sort the dataframe back according to mass
	selected_signals_most_intense <- selected_signals_most_intense [order(selected_signals_most_intense$mass),]
	selected_signals_less_intense <- selected_signals_less_intense [order(selected_signals_less_intense$mass),]
	# Put these signals back into the peaklist
	peaks[[p]]@mass <- selected_signals_less_intense$mass
	peaks[[p]]@intensity <- selected_signals_less_intense$intensity
	peaks[[p]]@SNR <- selected_signals_less_intense$SNR
	peaks[[p]]@mass <- append(peaks[[p]]@mass, selected_signals_most_intense$mass)
	peaks[[p]]@intensity <- append(peaks[[p]]@intensity, selected_signals_most_intense$intensity)
	peaks[[p]]@SNR <- append(peaks[[p]]@SNR, selected_signals_most_intense$SNR)
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
# Load the required libraries
install_and_load_required_packages("parallel")
#
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
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
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
############################ peaks_library SHULD HAVE SNR REPLACED BY SD!
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
						if ((peaks_test[[s]]@intensity[i] <= (peaks_library[[l]]@intensity[j] + (number_of_st_dev * peaks_library[[l]]@SNR[j]))) && (peaks_test[[s]]@intensity[i] >= (peaks_library[[l]]@intensity[j] - (number_of_st_dev * peaks_library[[l]]@SNR[j])))) {
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










########################################################################









########################################################################





################### AVERAGE THE REPLICATES (BY NAME)
average_replicates_by_name <- function (spectra, pattern_for_replicates="name_of_the_file", number_of_characters_for_similarity=10, file_format="brukerflex") {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
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
filter_peaklist <- function (peaks, SNR_filter=15, signal_number_threshold=15) {
# Load the required libraries
install_and_load_required_packages("parallel")
#
if (SNR_filter != 0 && signal_number_threshold != 0) {
	# Create the clusters
	cl <- makeCluster(cpu_core_number)
	############################################# FILTERING FUNCTION
	peak_filtering_function <- function (peaks, SNR_filter, signal_number_threshold) {
		### Create a dataframe with the SNR
		peaks_data_frame <- data.frame (SNR = peaks@SNR)
		### elect only the SNR more than the threshold
		peaks_data_frame <- peaks_data_frame [peaks_data_frame$SNR >= SNR_filter,]
		### Select only the desired peaklists, based upon the threshold characteristics
			if (length(peaks_data_frame) < signal_number_threshold) {
				# Remove the bad spectrum from the list
				peaks <- NULL
			}
		return (peaks)
	}
	################################################################
	peaks_purified <- parLapply(cl, peaks, fun=function(peaks) peak_filtering_function(peaks, SNR_filter=SNR_filter, signal_number_threshold=signal_number_threshold))
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
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
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
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
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
	peaks <- replace_class_name (peaks, class_list=class_list)
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
	peaks <- replace_class_name (peaks, class_list=class_list)
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















########################################################################










########################################################################





############################ PEAK FILTERING
discard_peaks <- function (peaks, list_of_peaks=list(), mass_type="monoisotopic", cluster_size=4, tof_mode="reflectron") {
# Load the required libraries
install_and_load_required_packages("parallel")
#
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
	peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, SNR = peaks@SNR)
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
	peaks@SNR <- peaks_data_frame$SNR
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
filter_spectra_noise <- function (spectra, SNR_filter=20, percentage_above_threshold=10) {
# Load the required libraries
install_and_load_required_packages(c("parallel", "MALDIquant"))
#
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
################################################## NOISY SPECTRA REMOVAL
noisy_spectra_removal_function <- function (spectra, SNR_filter, percentage_above_threshold) {
	# Peak picking
	peaks <- detectPeaks (spectra, method="MAD", SNR=3)
	# Keep only the spectra with a certain percentage of signals above the threshold
	peaks_new <- detectPeaks (spectra, method="MAD", SNR=SNR_filter)
	if ((length(peaks_new@mass) / length(peaks@mass))*100 < percentage_above_threshold) {
		spectra <- NULL
	}
	return (spectra)
}
########################################################################
spectra <- parLapply(cl, spectra, fun=function (spectra) noisy_spectra_removal_function (spectra, SNR_filter, percentage_above_threshold))
stopCluster(cl)
### Keep only the elements that are different from NULL
spectra <- spectra [!sapply(spectra, is.null)]
return (spectra)
}





########################################################################









################################################################################










################################################################################





######## RICHEST PEAKLIST AS REFERENCE FOR THE SPECTRA ALIGNMENT/WARPING
richest_peaklist_for_alignment <- function (spectra, SNR=5, tolerance_ppm=2000) {
# Load the required libraries
install_and_load_required_packages("MALDIquant")
#
# Peak picking
peaks <- detectPeaks (spectra, method="MAD", SNR=SNR)
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
spectra <- alignSpectra(spectra, reference=reference_peaklist, tolerance=(tolerance_ppm/10^6), noiseMethod="MAD", SNR=SNR, warpingMethod="quadratic")
return (spectra)
}





################################################################################





############################# MATRIX ZERO FILLING
# Replace the NA values with a zero
matrix_zero_filling <- function (input_matrix) {
# Load the required libraries
install_and_load_required_packages("parallel")
#
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










################################################################################





############################## PEAK STATISTICS (on processed Spectra)
peak_statistics <- function (spectra, SNR=3, class_list=list(), tof_mode="linear", file_format="imzml", tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25, remove_outliers=TRUE) {
# Load the required libraries
install_and_load_required_packages(c("MALDIquant", "stats"))
#
# Determine the number of classes
if (length(class_list) == 0 | length(class_list) == 1) {
number_of_classes <- 1
}
if (length(class_list) > 1) {
	number_of_classes <- length(class_list)
}
# Detect and Align Peaks
if (tof_mode == "linear") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR, halfWindowSize=20)
}
if (tof_mode == "reflector" | tof_mode == "reflectron") {
	peaks <- detectPeaks (spectra, method="MAD", SNR=SNR, halfWindowSize=5)
}
peaks <- align_and_filter_peaks(peaks, tolerance_ppm=tolerance_ppm, peaks_filtering=peaks_filtering, frequency_threshold=frequency_threshold)
# Generate the matrix (and convert it into a data frame)
signal_matrix <- intensityMatrix(peaks, spectra)
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
		st_dev_intensity <- sd(intensity_vector)
		summary_intensity_vector <- summary (intensity_vector)
		mean_intensity <- summary_intensity_vector [4]
		coeff_variation <- (st_dev_intensity / mean_intensity) *100
		median_intensityensity <- summary_intensity_vector [3]
		first_quartile <- summary_intensity_vector [2]
		third_quartile <- summary_intensity_vector [5]
		inter_quartile_range <- third_quartile - first_quartile
		spectra_counter <- length(intensity_vector)
		# Fill the matrix with the values
		peak_stat_matrix [p,1] <- distribution_type
		peak_stat_matrix [p,2] <- mean_intensity
		peak_stat_matrix [p,3] <- st_dev_intensity
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
			summary_intensity_vector [[l]] <- summary (intensity_vector[[l]])
			mean_intensity[[l]] <- summary_intensity_vector[[l]] [4]
			coeff_variation[[l]] <- (st_dev_intensity[[l]] / mean_intensity[[l]]) *100
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
			st_dev_name <- paste(st_dev_name, " ", st_dev_intensity[[l]], " - ", class_list[l], sep="")
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

















################################# REMOVE UNDIGESTED SPECTRA
remove_undigested_spectra <- function (spectra, mass_range_of_interest=c(0,2500), percentage_of_spectrum_in_the_range=75) {
# Load the required libraries
install_and_load_required_packages("parallel")
#
# Make the cluster (one for each core/thread)
cl <- makeCluster(cpu_core_number)
clusterEvalQ (cl, {library(MALDIquant)})
################################################# UNDIGESTED FILTERING FUNCTION
undigested_filtering_function <- function (spectra, mass_range_of_interest, percentage_of_spectrum_in_the_range) {
	## Peak picking
	peaks <- detectPeaks (spectra, method="MAD", SNR=3, halfWindowSize=20)
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







################################################################################





############################################## DEISOTOPING
de_isotope_peaks <- function (peaks) {
# Load the required libraries
install_and_load_required_packages("parallel")
#
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




###################### BIOTYPER SCORE ACCORDING TO THE HIERARCHICAL DISTANCE
biotyper_score_hierarchical_distance <- function (library_list, test_list) {
# Merge the library and test peaklist together (for alignment)
global_peaklist <- c(library_list$peaks, test_list$peaks)
# Align peaks
global_peaklist <- align_and_filter_peaks(global_peaklist, tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25)
# Merge spectra
global_spectralist <- c(library_list$spectra, test_list$spectra)
# Generate the matrix (for hca)
peaklist_matrix <- intensityMatrix(global_peaklist, global_spectralist)
# Add additional info to the matrix
peaklist_matrix <- matrix_add_class_and_sample(peaklist_matrix, peaks=global_peaklist, class_list=list(), file_format="brukerflex", sample_output=TRUE, class_output=FALSE)
rownames(peaklist_matrix) <- make.names(peaklist_matrix[,"Sample"], unique=TRUE)
# Compute the hca
distance_matrix <- dist(peaklist_matrix[,1:(ncol(peaklist_matrix)-1)], method="euclidean")
hierarchical_clustering <- hclust(distance_matrix)
plot(hierarchical_clustering, main="Hierarchical clustering analysis - Biotyper-like", xlab="Samples", ylab="Tree height")
hca_dendrogram <- recordPlot()
#
distance_matrix <- as.matrix(distance_matrix)
# The distance matrix displays the distance between the spectra
colnames(distance_matrix) <- peaklist_matrix[,"Sample"]
rownames(distance_matrix) <- peaklist_matrix[,"Sample"]
# Remove the first rows (the spectra from the database)
distance_matrix <- distance_matrix[(length(library_list$spectra)+1):nrow(distance_matrix),]
# Keep only the first columns (the spectra from the database)
distance_matrix <- distance_matrix[,1:length(library_list$spectra)]
# Normalise the euclidean distances, by multiplying the value by 100
distance_matrix <- distance_matrix * 100
# The classification is made by comparing the single sample spectrum with the spectrum of the database class (the distance is displayed in the distance matrix): the closer the better
result_matrix <- distance_matrix
# Scroll the rows, assign the class based upon the distance (the lowest is the class), create the output matrix for results (create a function to apply to each matrix row) (also check the absolute distance, not only the relative!! A sample might be far from all the elements in the database!)
scoring_function <- function (matrix_row) {
    minimum_value_position <- which(matrix_row == min(matrix_row))
    other_positions <- which(matrix_row != min(matrix_row))
    matrix_row [minimum_value_position] <- paste("YES (", round(as.numeric(matrix_row [minimum_value_position]),3), ")", sep="")
    for (i in other_positions) {
        matrix_row[i] <- paste("NO (", round(as.numeric(matrix_row[i]),3), ")", sep="")
    }
    return (matrix_row)
}
result_matrix <- apply(result_matrix, MARGIN=1, FUN=function(matrix_row) scoring_function(matrix_row))
result_matrix <- t(result_matrix)
return (list(result_matrix=result_matrix, plots=hca_dendrogram))
}





###############################################################################










################################################################################







































































################## FEATURE SELECTION FUNCTIONS
################################################ INSTALL REQUIRED PACKAGES
install_and_load_required_packages <- function (required_packages) {
installed_packages <- installed.packages()[,1]
missing_packages <- character()
for (p in 1:length(required_packages)) {
	if ((required_packages[p] %in% installed_packages) == FALSE) {
		missing_packages <- append(missing_packages, required_packages[p])
	}
}
if (length(missing_packages) > 0) {
	install.packages(missing_packages)
}
for (p in required_packages) {
	library(p, character.only=TRUE)
}
}





####################################### READ (SOURCE) SCRIPT FROM URL
source_github <- function(url, ...) {
  # Load required package (install it if not present)
  if ("RCurl" %in% installed.packages()[,1] == FALSE) {
	  install.packages("RCurl")
  }
  library(RCurl)
  # Parse and evaluate each R script (in the global environement)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}





###############################################################################










################################################################################





################################################## NO HEMOGLOBIN IN THE PEAKLIST
no_hemoglobin_peaklist <- function (peaklist, hemoglobin_mass=15400, tolerance_ppm=2000, non_features=c("Sample", "Class")) {
# Identify the monocharged and the bicharged ion
hemoglobin <- c(hemoglobin_mass, hemoglobin_mass/2)
# Extract the feature names (as numeric)
feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Create the index for the hemoglobin mass(es)
index <- integer()
# Scroll the feature list
for (f in 1:length(feature_list)) {
	for (h in 1:length(hemoglobin)) {
		# Identify the hemoglobin
		if (abs(hemoglobin[h] - feature_list[f])*10^6/hemoglobin_mass <= tolerance_ppm) {
			index <- append(index, f)
		}
	}
}
# Keep only the part of the dataframe without the hemoglobin column(s)
if (length(index) > 0) {
	peaklist_no_hemoglobin <- peaklist [,-index]
}	else {peaklist_no_hemoglobin <- peaklist}
return (peaklist_no_hemoglobin)
}





################################################################################





























################################################## DEISOTOPING IN PEAKLIST
deisotope_peaklist <- function (peaklist, non_features=c("Sample", "Class")) {
ions <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Deisotoped Ions
deisotoped_ions <- numeric()
# Take each ion and evaluate the consecutive ions
isotope_bin <- numeric ()
# Scroll the ion list
for (i in 1:(length(ions)-1)) {
	# Create the isotope bin, if the consecutive peaks belong to a isotope cluster
	if (abs(ions[i+1] - ions[i]) < 1.5) {
		isotope_bin <- append(isotope_bin, ions[i])
	}
	# When a major distance is found... We are out of the cluster
	if (abs(ions[i+1] - ions[i]) >= 1.5) {
		# Add this mass to the isotope bin (it will be the last belonging to a this isotope cluster
		isotope_bin <- append(isotope_bin, ions[i])
		# Store the first isotope of the isotope Bin created so far in a new vector
		deisotoped_ions <- append(deisotoped_ions, isotope_bin[1])
		# Empty the isotope Bin in order for it to become available for the future cycle
		isotope_bin <- numeric()

	}
}
# Keep only the part of the matrix with the selected ions
deisotoped_peaklist <- peaklist [,as.character(deisotoped_ions)]
# Add the non features back
column_names <- names(deisotoped_peaklist)
for (i in 1:length(non_features)) {
	deisotoped_peaklist <- cbind(deisotoped_peaklist, peaklist[,non_features[i]])
}
names(deisotoped_peaklist) <- c(column_names, non_features)
return (deisotoped_peaklist)
}





################################################################################





################################################## NO TRYPSIN IN THE PEAKLIST
no_trypsin_peaklist <- function (peaklist, trypsin_mass=c(841.50, 905.50, 1005.48, 1044.56, 1468.72, 1735.84, 1767.79, 2157.02, 2210.10, 2282.17, 3012.32, 4488.11), tolerance_ppm=200, non_features=c("Sample", "Class")) {
# Extract the feature names (as numeric)
feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Create the index for the hemoglobin mass(es)
index <- integer()
# Scroll the feature list
for (f in 1:length(feature_list)) {
	for (t in 1:length(trypsin_mass)) {
		# Identify the hemoglobin
		if (abs(trypsin_mass[t] - feature_list[f])*10^6/trypsin_mass[t] <= tolerance_ppm) {
			index <- append(index, f)
		}
	}
}
# Keep only the part of the dataframe without the hemoglobin column(s)
if (length(index) > 0) {
	peaklist_no_trypsin <- peaklist [,-index]
}	else {peaklist_no_trypsin <- peaklist}
return (peaklist_no_trypsin)
}





################################################################################

















############################################ OUTLIERS REMOVAL
outliers_removal <- function (vector, replace_with="") {
summary_vector <- summary(vector)
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
		# Replace the outliers with the vector mean (without outliers)
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





################################################################################





############################## PEAK STATISTICS (on the final peaklist)
peak_statistics_peaklist <- function (peaklist, non_features=c("Sample","Class","THY"), remove_outliers=TRUE, replace_outliers_with="") {
class_list <- levels(peaklist$Class)
# Determine the number of classes
if (length(class_list) == 0 | length(class_list) == 1) {
number_of_classes <- 1
}
if (length(class_list) > 1) {
	number_of_classes <- length(class_list)
}
# Peak vector
#peakVector <- as.numeric(names(peaklist))
############################################################## ONE CLASS
if (number_of_classes == 1) {
	# Output matrix
	peak_statistics_matrix <- matrix (0, nrow=(ncol(peaklist)-length(non_features)), ncol=8)
	rownames (peak_statistics_matrix) <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
	colnames (peak_statistics_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter")
	# For each peak
	for (p in 1:(ncol(peaklist)-length(non_features))) {
		intensity_vector <- as.numeric(peaklist[,p])
		if (remove_outliers == TRUE) {
			intensity_vector <- outliers_removal (intensity_vector, replace_with=replace_outliers_with)$vector
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
		#spectraNames <- levels(peaklist$Sample)
		st_dev_intensity <- sd (intensity_vector)
		summary_intensity_vector <- summary (intensity_vector)
		mean_intensity <- summary_intensity_vector [4]
		coefficient_of_variation <- (st_dev_intensity / mean_intensity) *100
		median_intensity <- summary_intensity_vector [3]
		first_quartile <- summary_intensity_vector [2]
		third_quartile <- summary_intensity_vector [5]
		inter_quartile_range <- third_quartile - first_quartile
		spectra_counter <- length (intensity_vector)
		# Fill the matrix with the values
		peak_statistics_matrix [p,1] <- distribution_type
		peak_statistics_matrix [p,2] <- mean_intensity
		peak_statistics_matrix [p,3] <- st_dev_intensity
		peak_statistics_matrix [p,4] <- coefficient_of_variation
		peak_statistics_matrix [p,5] <- median_intensity
		peak_statistics_matrix [p,6] <- inter_quartile_range
		peak_statistics_matrix [p,7] <- paste ("1st quartile", first_quartile, "; 3rd quartile", third_quartile)
		peak_statistics_matrix [p,8] <- spectra_counter
		#peak_statistics_matrix [p,9] <- spectraNames
	}
}
############################################################ TWO OR MORE CLASSES
# Every variable now is a list, each element of which corresponds to a certain value from a class
# So every variable is a list with the same length of the class list (each element of the list
# is referred to a class
if (number_of_classes > 1) {
	# Output matrix
	peak_statistics_matrix <- matrix (0, nrow=(ncol(peaklist)-length(non_features)), ncol=14)
	rownames (peak_statistics_matrix) <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
	colnames (peak_statistics_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %",
		"Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)",
		"t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
	# For each peak
	for (p in 1:(ncol(peaklist)-length(non_features))) {
		# Put the intensity of that peak into one vector per class (in a global list)
		intensity_vector <- list()
		# Scroll the peaklists and Add the peak intensity to a vector (one for each class)
		for (l in 1:length(class_list)) {
			# Allocate in the intensity vector the rows for that peak belonging to the certain class
			intensity_vector [[l]] <- peaklist [peaklist$Class == class_list[l],p]
		}
		if (remove_outliers == TRUE) {
			for (i in 1:length(intensity_vector)) {
				intensity_vector[[l]] <- outliers_removal (intensity_vector[[l]], replace_with=replace_outliers_with)$vector
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
			variance_test_parametric  <- var.test(intensity_vector[[1]], intensity_vector[[2]])
		}
		if (length(class_list) >= 2) {
			variance_test_non_parametric  <- bartlett.test(as.numeric(peaklist[,p]), g=peaklist$Class)
		}
		########################################### Other parameters (per class)
		st_dev_intensity <- list()
		summary_intensity_vector <- list()
		mean_intensity <- list()
		coefficient_of_variation <- list()
		median_intensity <- list()
		first_quartile <- list()
		third_quartile <- list()
		inter_quartile_range <- list()
		spectra_counter <- list()
		variance <- list()
		for (l in 1:length(class_list)) {
			st_dev_intensity[[l]] <- sd (intensity_vector[[l]])
			summary_intensity_vector [[l]] <- summary (intensity_vector[[l]])
			mean_intensity[[l]] <- summary_intensity_vector[[l]] [4]
			coefficient_of_variation[[l]] <- (st_dev_intensity[[l]] / mean_intensity[[l]]) *100
			median_intensity[[l]] <- summary_intensity_vector[[l]] [3]
			first_quartile[[l]] <- summary_intensity_vector[[l]] [2]
			third_quartile[[l]] <- summary_intensity_vector[[l]] [5]
			inter_quartile_range[[l]] <- third_quartile[[l]] - first_quartile[[l]]
			spectra_counter[[l]] <- length (intensity_vector[[l]])
			variance[[l]] <- var (intensity_vector[[l]])
		}
		############################################# Parameters between classes
		# T-test (Two classes, parametric)
		if (length(class_list) == 2) {
			t_test <- t.test (intensity_vector[[1]], intensity_vector[[2]])
		}
		# ANOVA TEST (More than two classes, parametric)
		if (length(class_list) >= 2) {
		anova_test <- aov (peaklist[,p] ~ peaklist$Class)
		}
		# WILCOXON - MANN-WHITNEY TEST (Two classes, non parametric)
		if (length(class_list) == 2) {
			wilcoxon_test <- wilcox.test (intensity_vector[[1]], intensity_vector[[2]])
		}
		# KRUSKAL-WALLIS TEST (more than two classes, non parametric)
		if (length(class_list) >= 2) {
			kruskal_wallis_test <- kruskal.test (peaklist[,p], g=peaklist$Class)
		}
		######################################## Fill the matrix with the values
		# Distribution Type
		distribution_type_name <- character()
		for (l in length(class_list)) {
			distribution_type_name <- paste (distribution_type_name, " ", distribution_type[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,1] <- paste (distribution_type_name)
		# Mean
		mean_intensity_name <- character()
		for (l in length(class_list)) {
			mean_intensity_name <- paste (mean_intensity_name, " ", mean_intensity[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,2] <- mean_intensity_name
		# Standard Deviation
		st_dev_name <- character()
		for (l in length(class_list)) {
			st_dev_name <- paste (st_dev_name, " ", st_dev_intensity[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,3] <- st_dev_name
		# Coefficient of Variation
		coefficient_of_variation_name <- character()
		for (l in length(class_list)) {
			coefficient_of_variation_name <- paste (coefficient_of_variation_name, " ", coefficient_of_variation[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,4] <- coefficient_of_variation_name
		# Median
		median_intensity_name <- character()
		for (l in length(class_list)) {
			median_intensity_name <- paste (median_intensity_name, " ", median_intensity[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,5] <- median_intensity_name
		# Interquartile Range (IQR)
		inter_quartile_range_name <- character()
		for (l in length(class_list)) {
			inter_quartile_range_name <- paste (inter_quartile_range_name, " ", inter_quartile_range[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,6] <- inter_quartile_range_name
		# Spectra counter
		spectra_counter_name <- character()
		for (l in length(class_list)) {
			spectra_counter_name <- paste (inter_quartile_range_name, " ", spectra_counter[[l]], " - ", class_list[l], sep="")
		}
		peak_statistics_matrix [p,7] <- spectra_counter_name
		# Class
		#className <- class_list
		#peak_statistics_matrix [p,8] <- className
		# Homoscedasticity (Parametric)
		if (variance_test_parametric $p.value < 0.05) {
		homoscedasticity_parametric <- paste ("Non homoscedastic data", "(p-value:", variance_test_parametric $p.value, ")")
		}
		if (variance_test_parametric $p.value >= 0.05) {
		homoscedasticity_parametric <- paste ("Homoscedastic data", "(p-value:", variance_test_parametric $p.value, ")")
		}
		if (variance_test_non_parametric $p.value < 0.05) {
		homoscedasticity_non_parametric <- paste ("Non homoscedastic data", "(p-value:", variance_test_non_parametric $p.value, ")")
		}
		if (variance_test_non_parametric $p.value >= 0.05) {
		homoscedasticity_non_parametric <- paste ("Homoscedastic data", "(p-value:", variance_test_non_parametric $p.value, ")")
		}
		peak_statistics_matrix [p,9] <- homoscedasticity_parametric
		peak_statistics_matrix [p,10] <- homoscedasticity_non_parametric
		# t-Test
		peak_statistics_matrix [p,11] <- t_test$p.value
		# ANOVA
		peak_statistics_matrix [p,12] <- summary(anova_test)[[1]]$"Pr(>F)"[1]
		# Wilcoxon / Mann-Whitney test
		peak_statistics_matrix [p,13] <- wilcoxon_test$p.value
		# Kruskal-Wallis test
		peak_statistics_matrix [p,14] <- kruskal_wallis_test$p.value
	}
}
return (peak_statistics_matrix)
}



################################################################################


























beep(5)
