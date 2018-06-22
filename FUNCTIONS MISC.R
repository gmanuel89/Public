#################### FUNCTIONS - MISC AND DEPRECATED 2017.03.22 ####################




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
deisotope_peaklist <- function (peaklist, non_features=c("Sample", "Class"), tolerance_da=1.5) {
# Extract the ion masses (numbers)
ions <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Extract the peaklist matrix with only the features
peaklist_features <- peaklist [,!(names(peaklist) %in% non_features)]
# Deisotoped Ions
deisotoped_ions <- numeric()
# Scroll the ion list (For each ion...)
for (i in 1:length(ions)) {
	isotope_bin <- numeric ()
	# Create the isotope bin, if the consecutive peaks belong to a isotope cluster
	# Scroll the peaks
	for (p in 1:length(ions)) {
		# Evaluate the difference between the masses
		if (abs(ions[i] - ions[p]) <= tolerance_da) {
			# Add the isotope cluster masses to the bin
			isotope_bin <- append(isotope_bin, ions[p])
		}
	}
	# Sort the bin elements
	isotope_bin <- sort(isotope_bin)
	# Store the first isotope of the isotope Bin created so far in a new vector
	deisotoped_ions <- append(deisotoped_ions, isotope_bin[1])
}
# Discard the duplicates
deisotoped_ions <- unique(deisotoped_ions)
# Convert the list of ions into characters
deisotoped_ions <- as.character(deisotoped_ions)
# Keep only the part of the matrix with the selected ions
deisotoped_peaklist <- peaklist [,deisotoped_ions]
# Add the non features back
column_names <- names(deisotoped_peaklist)
for (n in 1:length(non_features)) {
	deisotoped_peaklist <- cbind(deisotoped_peaklist, peaklist[,non_features[n]])
}
names(deisotoped_peaklist) <- as.character(c(column_names, non_features))
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





############################### ADD CUSTOM FEATURES TO THE PEAKLIST INTENSITY MATRIX
# The function takes a list of spectra and a vector of custom features to be included in the generation of the final peaklist intensity matrix. If a matrix is specified as input, the columns corresponding to the features to be searched for are appended to the matrix itself. The input matrix must have the number of spectra as the number of rows.
custom_peaklist_intensity_matrix <- function (spectra, features_to_add = numeric(), final_sample_matrix = NULL, allow_parallelization = FALSE, tolerance_ppm = 2000) {
    ## Install the required packages
    install_and_load_required_packages("MALDIquant")
    #################### Multiple spectra
    if (isMassSpectrumList(spectra)) {
        ### If there are features to add...
        if (length(features_to_add > 0)) {
            ## Generate a fake spectrum and a fake peaklist with the features to add
            fake_spectrum <- createMassSpectrum(mass = spectra[[1]]@mass, intensity = spectra[[1]]@intensity, metaData = list(name = "Fake spectrum"))
            fake_peaks <- createMassPeaks(mass = as.numeric(features_to_add), intensity = rep(1, length(features_to_add)), snr = rep(3, length(features_to_add)), metaData = list(name = "Fake peaklist"))
            ## Detect the peaks in the spectra
            peaks <- detectPeaks(spectra, method = "SuperSmoother", SNR = 3)
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
                # Initialize the output
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





################################################################################





################################################################ SIGNAL FOLLOWER
# This function takes a folder containing imzML files in it and a list of signals of interest (along with their possible names). It imports the spectra and returns a matrix listing the signal of interests along with their statistics.
signal_follower_statistics <- function (filepath, signal_list, mass_labels = list(), SNR = 5, spectra_format = "imzml", tof_mode = "linear", smoothing_strength_preprocessing = "medium", process_in_packages_of = 0, tolerance_ppm = 2000, peak_picking_algorithm = "SuperSmoother") {
    # Load the required libraries
    install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Define the column and row headers (if mass labels are provided or not)
    if (length(mass_labels) != 0) {
        column_names_st_dev <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_st_dev[i] <- paste(signal_list[i], "-", mass_labels[i], "- StDev")
        }
        column_names_coeff_var <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_coeff_var[i] <- paste(signal_list[i], "-", mass_labels[i], "- CV")
        }
        column_names_mean_abs_dev <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_mean_abs_dev[i] <- paste(signal_list[i], "-", mass_labels[i], "- MeanAbsDev")
        }
    } else {
        column_names_st_dev <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_st_dev[i] <- paste(signal_list[i], "- StDev")
        }
        column_names_coeff_var <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_coeff_var[i] <- paste(signal_list[i], "- CV")
        }
        column_names_mean_abs_dev <- vector(length = length(signal_list))
        for (i in 1:length(signal_list)) {
            column_names_mean_abs_dev[i] <- paste(signal_list[i], "- MeanAbsDev")
        }
    }
    sample_names <- filepath
    ### Create the result matrices
    st_dev_result_matrix <- matrix (ncol = length(signal_list), nrow = length(filepath))
    colnames(st_dev_result_matrix) <- column_names_st_dev
    rownames(st_dev_result_matrix) <- sample_names
    coeff_var_result_matrix <- matrix (ncol = length(signal_list), nrow = length(filepath))
    colnames(coeff_var_result_matrix) <- column_names_coeff_var
    rownames(coeff_var_result_matrix) <- sample_names
    mean_abs_dev_result_matrix <- matrix (ncol = length(signal_list), nrow = length(filepath))
    colnames(mean_abs_dev_result_matrix) <- column_names_mean_abs_dev
    rownames(mean_abs_dev_result_matrix) <- sample_names
    # The script is run for every imzML file
    for (lib in 1:length(filepath)) {
        # Spectra import and preprocessing
        if (!is.null(mass_range)) {
            spectra <- importImzMl(filepath[lib], massRange = mass_range)
        } else {
            spectra <- importImzMl(filepath[lib])
        }
        spectra <- preprocess_spectra(spectra, tof_mode = tof_mode, smoothing_strength = smoothing_strength_preprocessing, process_in_packages_of = process_in_packages_of)
        ### Peak Picking and alignment
        peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
        peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
        # Generate the intensity matrix
        intensity_matrix <- intensityMatrix(peaks, spectra)
        #### Create the partial result matrices
        st_dev_matrix <- matrix (ncol = length(signal_list), nrow = 1)
        colnames(st_dev_matrix) <- column_names_st_dev
        rownames(st_dev_matrix) <- spectra[[1]]@metaData$file[1]
        coeff_var_matrix <- matrix (0, ncol = length(signal_list), nrow = 1)
        colnames(coeff_var_matrix) <- column_names_coeff_var
        rownames(coeff_var_matrix) <- spectra[[1]]@metaData$file[1]
        mean_abs_dev_matrix <- matrix (0, ncol = length(signal_list), nrow = 1)
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
                st_dev_intensity <- sd(intensity_vector, na.rm = TRUE)
                mean_abs_dev_int <- mad(intensity_vector, na.rm = TRUE)
                mean_intensity <- mean(intensity_vector, na.rm = TRUE)
                median_intensity <- median(intensity_vector, na.rm = FALSE)
                coeff_var_intensity <- (st_dev_intensity / mean_intensity)*100
                # Fill in the partial result matrices with the values
                st_dev_matrix [1,s] <- paste(mean_intensity, "+/-", st_dev_intensity)
                coeff_var_matrix [1,s] <- paste(coeff_var_intensity, "%", "( mean  = ", mean_intensity, ")")
                mean_abs_dev_matrix [1,s] <- paste(median_intensity, "+/-", mean_abs_dev_int)
            }
        }
        # Put the partial result matrices together in one matrix
        st_dev_result_matrix [lib,] <- st_dev_matrix
        coeff_var_result_matrix [lib,] <- coeff_var_matrix
        mean_abs_dev_result_matrix [lib,] <- mean_abs_dev_matrix
    }
    result_matrix <- cbind(st_dev_result_matrix, coeff_var_result_matrix, mean_abs_dev_result_matrix)
    return(result_matrix)
}





################################################################################





###################### SPECTRA PURIFICATION BASED ON THE TIC (Before processing)
# This function discards spectra with a total ion current (TIC) below the selected threshold.
# A specific threshold can be specified; if just the TOF MS mode is specified, a default TIC threshold value is used.
spectra_tic_purification <- function (spectra, tof_mode = "reflectron", absolute_tic_threshold = 0) {
    # If not specified, set it based on the TOF mode
    if (absolute_tic_threshold == 0) {
        if (tof_mode == "linear") {
            absolute_tic_threshold <- 2000
        }
        if (tof_mode == "reflectron" | tof_mode == "reflector") {
            absolute_tic_threshold <- 200
        }
    }    else {absolute_tic_threshold <- absolute_tic_threshold}
    # Before preprocessing (and thus normalization), evaluate the TIC
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
    }    else {spectra <- spectra}
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
    }    else {spectra <- spectra}
    #
    return (spectra)
}





################################################################################





################## PLOT THE SIGNALS OF INTEREST WITH THE SD BARS ON THE AVERAGE
average_spectrum_bars_signals_of_interest <- function (spectra, SNR = 5, signals_of_interest = peaks_average@mass, tolerance_ppm = 2000, tof_mode = "linear", half_window_plot = 1000, graph_title = "Spectrum", average_spectrum_color = "black", peak_points = TRUE, points_color = "red", bar_width = 40, bar_color = "blue", peak_picking_algorithm = "SuperSmoother") {
    # Load the required libraries
    install_and_load_required_packages("MALDIquant")
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Outputs
    spectrum_images <- list()
    # Generate the average spectrum
    average_spectrum <- averageMassSpectra(spectra, method = "mean")
    average_spectrum <- removeBaseline(average_spectrum, method = "TopHat")
    # Peak picking on the average spectrum (for plotting)
    peaks_average <- peak_picking(average_spectrum, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
    # Peak picking on the dataset
    peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
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
        intensity_vector <- vector(length = 0)
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
        plot(average_spectrum, main = graph_title, col.main = average_spectrum_color, xlab = "m/z", ylab = "Intensity (a.i.)", xlim = c(low_mass_plot,up_mass_plot), col = average_spectrum_color)
        # Peaks
        if (peak_points == TRUE) {
            points(peaks_average[peak_index], pch = 4, col = points_color)
        }
        # Bars
        # epsilon: lengthof the horizontal segment
        epsilon = bar_width
        # Define the upper and the lower limit of the vertical segment
        up <- peaks_average@intensity[peak_index] + st_dev_intensity
        low <- peaks_average@intensity[peak_index] - st_dev_intensity
        # Vertical bar (x,y x,y)
        segments(signals_of_interest[a], low, signals_of_interest[a], up, col = bar_color)
        # Horizontal segments(x,y , x,y)
        segments(signals_of_interest[a]-epsilon, low, signals_of_interest[a]+epsilon, low, col = bar_color)
        segments(signals_of_interest[a]-epsilon, up, signals_of_interest[a]+epsilon, up, col = bar_color)
        # Store the zoomed part of the spectrum with the signal
        avg_spectrum_zoom <- recordPlot()
        # Add this to a final list of images
        spectrum_images <- append(spectrum_images, avg_spectrum_zoom)
    }
    return (spectrum_images)
}





################################################################################





################################################### SPECTRA GROUPING (PATIENTS)
# The functions takes a list of already preprocessed spectra (MALDIquant) and generates a list of representative spectra, averaging spectra randomly or according to a cluster algorithm.
# It is advisable to use a list of spectra coming from one patient (one imzML)
# Spectra_per_patient = 1 returns the average spectrum of the spectra dataset.
# The function returns a list of elements: representative spectra, MS images of pixels under the same nodes colored the same way, list of spectra under the discarded nodes (with their average spectrum), list of spectra grouped according to the node they belong to.
group_spectra <- function(spectra, spectra_per_patient = 1, spectra_format = "imzml", tof_mode = "linear", seed = NULL, algorithm = "random", clustering_method = "agglomerative", discarded_nodes = 1, balanced = TRUE, method = "mean", peak_picking_algorithm = "SuperSmoother") {
    # Load the required libraries
    install_and_load_required_packages (c("MALDIquant", "caret", "stats"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
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
            patient_spectra_final <- averageMassSpectra(spectra, labels = file_vector, method = "mean")
        }
        if (method == "skyline") {
            patient_spectra_final <- group_spectra_skyline(spectra, spectra_format = spectra_format)
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
                    index <- createFolds(file_vector[file_vector == patient_vector[p]], k = spectra_per_patient)
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
                            temporary_patient_average <- averageMassSpectra (patient_spectra_temp, method = "mean")
                        }
                        if (method == "skyline") {
                            temporary_patient_average <- generate_skyline_spectrum (patient_spectra_temp)
                        }
                        # Store them in the final spectra list
                        patient_spectra_final <- append(patient_spectra_final, temporary_patient_average)
                    }
                    # Store the plot into the list
                    plotMsiSlice(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[1], tolerance = 1, legend = FALSE)
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
                    peaks <- peak_picking(patient_spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
                    peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
                    # Generate the peaklist matrix
                    peaklist <- intensityMatrix(peaks, patient_spectra)
                    # Compute the distance matrix
                    distance_matrix <- dist(peaklist, method = "euclidean")
                    hca <- hclust(distance_matrix)
                    # Store the plot
                    plot(hca, xlab = "MSI elements", main = "Hierarchical clustering analysis", sub = "")
                    legend_text <- paste("Cluster method:", hca$method, "\nDistance method: ", hca$dist.method, "\n")
                    legend("topright", title = "Hierarchical clustering", legend = legend_text, xjust = 0.5, border = "black")
                    plots[[p]] <- recordPlot()
                    # Associate to each row/spectrum the subgroup of the tree to which it belongs
                    if (spectra_per_patient > 56) {
                        spectra_per_patient <- 56
                    }
                    hca_groups <- cutree(hca, k = spectra_per_patient)
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
                            }
                            if (method == "skyline") {
                                spectra_average <- generate_skyline_spectrum(spectra_hca)
                            }
                            # Add the average to the final list
                            patient_spectra_final <- append(patient_spectra_final, spectra_average)
                        }
                    }
                    # Store the plot into the list
                    plotMsiSlice(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[1], tolerance = 1, legend = FALSE)
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
                    peaks <- peak_picking(patient_spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
                    peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
                    # Generate the peaklist matrix
                    peaklist <- intensityMatrix(peaks, patient_spectra)
                    # Compute the k-Means clustering
                    hca <- kmeans(peaklist, centers = spectra_per_patient)
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
                            }
                            if (method == "skyline") {
                                spectra_average <- generate_skyline_spectrum(spectra_hca, method = "mean")
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
                                spectra_average <- averageMassSpectra(spectra_hca, method = "mean")
                            }
                            if (method == "skyline") {
                                spectra_average <- generate_skyline_spectrum(spectra_hca, method = "mean")
                            }
                            # Add the average to the final list
                            patient_spectra_final <- append(patient_spectra_final, spectra_average)
                        }
                    }
                    # Store the plot into the list
                    plotMsiSlice(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[1], tolerance = 1, legend = FALSE)
                    msi_plots [[p]] <- recordPlot()
                    # Clustered spectra
                    spectra_hca_grouped <- spectra_hca_grouped_patient
                }
                #####################
            }
        }
    }
    # Processing before returning
    patient_spectra_final <- removeBaseline(patient_spectra_final, method = "TopHat")
    patient_spectra_final <- calibrateIntensity(patient_spectra_final, method = "TIC")
    #
    return (list(spectra = patient_spectra_final, ms_images = msi_plots, discarded_spectra = discarded_spectra, discarded_average_spectra = discarded_spectra_average, plots = plots, spectra_hca_grouped = spectra_hca_grouped))
}





################################################################################











###############################################################################





################################################ LIBRARY CREATION
# This function reads the files contained in a provided folder (no memory efficient importing), it can average them according to the class they belong to or according to the folder they are into (average replicate/patient). It also computes the peak picking and can replace the SNR field in the peaklist with the standard deviation of the peaks or the coefficient of variation.
# It returns a list containing: the spectra and the peaks (NOT aligned).
# It is used to create the library/database.
library_creation <- function (filepath_database, peak_picking_algorithm = "SuperSmoother", class_grouping = TRUE, mass_range = c(3000,15000), spectra_preprocessing = list(smoothing_strength = "medium", preprocess_in_packages_of = length(spectra)), average_replicates = FALSE, average_patients = FALSE, SNR = 5, most_intense_peaks = FALSE, signals_to_take = 20, reference_peaklist_for_alignment = NULL, tof_mode = "linear", spectra_format = "brukerflex") {
    # Load the required libraries
    install_and_load_required_packages(c("MALDIquantForeign", "MALDIquant"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ### Define the classes (the classes are the folders in the library directory)
    class_list <- dir(filepath_database, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
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
            spectra <- importBrukerFlex(filepath_database, massRange = mass_range)
        } else {
            spectra <- importBrukerFlex(filepath_database)
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        ### Load the spectra
        if (!is.null(mass_range)) {
            spectra <- importImzMl(filepath_database, massRange = mass_range)
        } else {
            spectra <- importImzMl(filepath_database)
        }
    }
    ### Average the replicates
    if (average_replicates == TRUE) {
        spectra <- average_replicates_by_folder(spectra, filepath_database, spectra_format = spectra_format)
    }
    # Average the patients
    if (average_patients == TRUE) {
        spectra <- group_spectra(spectra, spectra_per_patient = 1, spectra_format = spectra_format, tof_mode = tof_mode)
    }
    ### Class grouping
    if (class_grouping == TRUE) {
        ### Spectra grouping (class)
        spectra <- group_spectra_class(spectra, class_list = class_list, spectra_format = spectra_format, class_in_file_name = TRUE)
    } else {spectra <- spectra}
    ### Preprocessing
    spectra <- preprocess_spectra(spectra, tof_mode = tof_mode, smoothing_strength = spectra_preprocessing$smoothing_strength, process_in_packages_of = spectra_preprocessing$preprocess_in_packages_of, sampling_method = "sequential")
    ##########################
    ### Peak picking on the individual spectra
    if (most_intense_peaks == TRUE) {
        peaks_database <- most_intense_signals(spectra, signals_to_take = signals_to_take)
    } else {
        peaks_database <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
        #peaks_database <- align_and_filter_peaks(peaks_database, tof_mode = tof_mode, low_intensity_peaks_removal = FALSE, reference_peaklist = reference_peaklist_for_alignment)
    }
    ####
    library_list <- list(spectra = spectra, peaks = peaks_database)
    return(library_list)
}





################################################################################





#################################### CLASSIFICATION VIA HIERARCHICAL CLUSTERING
# The input is a S4 list of spectra (MALDIquant), both for the classification and for the database, the class names (if in the spectrum file name, the filenames will be replaced with the class name and averaged to create representative spectra for the classes; if not in the files, one representative spectrum per class should be provided in the same order as the class list), the TOF-MS mode, the algorithm for the clustering and the nodes to discard.
# The script returns a list containing: the hierarchical analysis dendrograms, a matrix with the classification (of pixels/spectra and samples) and the MS image with pixel or area classification, along with the classification of the patients.
hierarchical_clustering_classification <- function(spectra_to_be_classified, spectra_database, peak_picking_algorithm = "SuperSmoother", class_list = c("HP","PTC"), class_in_file_name = TRUE, tof_mode = "linear", clustering_method = "agglomerative", nodes = 4, discarded_nodes = 0, seed = NULL, classification_of = c("average","subareas","pixels"), true_class_in_file_name = FALSE, true_class_list = NULL) {
    # Do everything if the user wants some outputs
    if (!is.null(classification_of) && length(classification_of) > 0) {
        ################################################## Load the required libraries
        install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign", "stats", "caret", "pROC"))
        # Rename the trim function
        trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
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
        spectra_database <- group_spectra_class(spectra_database, class_list = class_list, spectra_format = "imzML")
        ########################## PATIENT SPECTRA
        # Fix the filenames
        spectra_to_be_classified <- replace_sample_name(spectra_to_be_classified, spectra_format = "imzML")
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
                patient_spectra_grouped <- group_spectra(patient_spectra, spectra_per_patient = nodes, spectra_format = "imzml", tof_mode = tof_mode, seed = seed, algorithm = "hca", clustering_method = clustering_method, discarded_nodes = discarded_nodes, balanced = TRUE, method = "mean")
                patient_spectra_clustered_average <- patient_spectra_grouped$spectra
                patient_spectra_clustered <- patient_spectra_grouped$spectra_hca_grouped[[1]]
                # Put the database and the spectra together
                global_spectra <- append(spectra_database, patient_spectra_clustered_average)
                global_peaks <- peak_picking(global_spectra, SNR = SNR, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode)
                global_peaks <- align_and_filter_peaks(global_peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
                # Compute the intensity matrix
                intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
                intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks = global_peaks, class_list = class_list, spectra_format = "imzml", sample_output = TRUE, class_output = FALSE)
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
                    classification_hca  <- cbind(apply(distance_matrix, MARGIN = 1, FUN = function(x) classification_hca_function(x)))
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
                    plotMsiSlice(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[1], tolerance = 1, legend = FALSE)
                    legend(x = "bottomright", legend = class_list, fill = c("green","red"), xjust = 0.5, yjust = 0.5)
                    legend(x = "topright", legend = spectra_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                    msi_plots [[p]] <- recordPlot()
                }
            }
            ######################### PIXEL-BY-PIXEL CLASSIFICATION
            if ("pixels" %in% classification_of) {
                # Put the database and the spectra together
                global_spectra <- append(spectra_database, patient_spectra)
                global_peaks <- peak_picking(global_spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
                global_peaks <- align_and_filter_peaks(global_peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
                # Compute the intensity matrix
                intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
                intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks = global_peaks, class_list = class_list, spectra_format = "imzml", sample_output = TRUE, class_output = FALSE)
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
                    classification_hca  <- cbind(apply(distance_matrix, MARGIN = 1, FUN = function(x) classification_hca_function(x)))
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
                    plotMsiSlice(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[1], tolerance = 1, legend = FALSE)
                    legend(x = "bottomright", legend = class_list, fill = c("green","red"), xjust = 0.5, yjust = 0.5)
                    legend(x = "topright", legend = spectra_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                    msi_plots_pixel_by_pixel [[p]] <- recordPlot()
                }
            }
            ############################# AVERAGE SPECTRUM
            if ("average" %in% classification_of) {
                # If the spectra per patient are more, average them for patient classification...
                if (isMassSpectrumList(patient_spectra)) {
                    patient_spectra <- averageMassSpectra(patient_spectra, method = "mean")
                }
                # Put the database and the spectra together
                global_spectra <- append(spectra_database, patient_spectra)
                global_peaks <- peak_picking(global_spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
                global_peaks <- align_and_filter_peaks(global_peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
                # Compute the intensity matrix
                intensity_matrix <- intensityMatrix(global_peaks, global_spectra)
                intensity_matrix <- matrix_add_class_and_sample(intensity_matrix, peaks = global_peaks, class_list = class_list, spectra_format = "imzml", sample_output = TRUE, class_output = FALSE)
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
            spectra_to_be_classified <- replace_class_name(spectra_to_be_classified, class_list = class_list)
            for (sp in spectra_to_be_classified) {
                true_class_list <- append(true_class_list, sp@metaData$file[1])
            }
            classification_hca_results_avg <- cbind(classification_hca_results_avg, true_class_list)
            hca_classification_performances <- confusionMatrix(classification_hca_results_avg[,1], classification_hca_results_avg[,2], positive = "HP")
            #### ROC analysis
            #hca_roc <- roc(response = as.numeric(classification_hca_results_avg[,2]), predictor = as.numeric(classification_hca_results_avg[,1]))
        }
        # Output returning
        return (list(classification_hca_results_avg = classification_hca_results_avg, classification_hca_results = classification_hca_results, classification_hca_results_pixel_by_pixel = classification_hca_results_pixel_by_pixel, hca_graph = hca_graph, ms_images = msi_plots, ms_images_pixel_by_pixel = msi_plots_pixel_by_pixel, hca_classification_performances = hca_classification_performances, hca_roc = hca_roc))
    }
}





################################################################################





#################################################### PEARSON CORRELATION p-value
# This function returns a p-value for the significance of a correlation coefficient.
correlation_pvalue <- function(correlation_coefficient, number_of_samples, correlation_type = "pearson", tails = 2) {
    degrees_of_freedom <- number_of_samples - 2
    if (correlation_type == "pearson") {
        t_value = correlation_coefficient * sqrt((number_of_samples-2)/(1-(correlation_coefficient)^2))
        p_value <- pt(t_value, df = degrees_of_freedom, lower.tail = FALSE, log.p = FALSE)
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
group_peaklist <- function (peaklist, rows_per_patient = 100, seed = NULL, balanced = TRUE, discard_poor_samples = FALSE, discard_if_lower_than = 100, non_features = c("Sample","Class"), algorithm = "random", clustering_method = "agglomerative", grouping_variable = "Sample") {
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
                    distance_matrix <- dist(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], method = "euclidean")
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
                    hca_groups <- cutree(hca, k = rows_per_patient)
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
                    k_means <- kmeans(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], centers = rows_per_patient)
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
    return(result_data_frame)
}





################################################################################





####################################################### TRUNCATE PEAKLIST MATRIX
# This function operates the truncation at the peaklist level, by removing peaks (columns) out of a certain range.
truncate_peaklist <- function (peaklist, range = c(4000, 15000), non_features = c("Sample", "Class")){
    features_to_keep <- numeric()
    # Do not do anything if the values are set to zero
    if (range[1] == 0 & range[2] == 0) {
        # Keep all the features
        features_to_keep <- as.numeric(names(peaklist[,!(names(peaklist) %in% non_features)]))
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
cross_validation_svm <- function (training_dataset, seed = NULL, k_fold = 10, repeats = 5, svm_model, non_features = c("Sample","Class","THY"), positive_class = levels(training_dataset$Class)[1]) {
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
    result_matrix <- matrix (NA, ncol = 14, nrow = 0)
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
            model <- svm (x = training_predictors, y = training_outcomes, scale = TRUE, kernel = svm_kernel, degree = svm_degree, gamma = svm_gamma, cost = svm_cost)
            # Use the model to predict the testing subset
            # Plant the seed only if a specified value is entered
            if (!is.null(seed)) {
                # Make the randomness reproducible
                set.seed(seed)
            }
            predicted_classes <- predict(model, newdata = test_predictors)
            # Compute the performances
            performances <- confusionMatrix(data = predicted_classes, test_outcomes, positive = positive_class)
            confusion_matrix_list [[k*i]] <- performances
            #### ROC analysis
            svm_roc <- roc(response = test_outcomes, predictor = as.numeric(predicted_classes))
            # Create a matrix row with the partial results
            result_matrix_row <- matrix (0, ncol = ncol(result_matrix), nrow = 1)
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
    result_matrix_average <- apply(result_matrix, 2, mean, na.rm = TRUE)
    #return (list(confusion_matrix = confusion_matrix_list, result_matrix = result_matrix, features = names(training_dataset)[!(names(training_dataset) %in% non_features)]))
    return (list(result_matrix_average = result_matrix_average))
}





################################################################################





###################################################### SVM TUNING AND VALIDATION
# This function operates the tuning of the Support Vector Machine (SVM), by testing all the parameters of the SVM (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
svm_tuning_and_validation2 <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(gamma = 10^(-5:5), cost = 10^(-5:5), epsilon = seq(1,2,by = 1), degree = 1:5, kernel = "radial"), k_fold_cv = 10, repeats_cv = 2, parameters = list(gamma = 0.1, cost = 10, epsilon = 0.1, degree = 3, kernel = "radial"), positive_class_cv = "HP", seed = NULL, pca = FALSE, number_of_components = 3) {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "kernlab", "e1071", "pROC"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
            svm_tuning <- tune.svm(as.data.frame(pca_training$x[,1:number_of_components]), factor(peaklist_training$Class), gamma = tuning_parameters$gamma, cost = tuning_parameters$cost, kernel = tuning_parameters$kernel, epsilon = tuning_parameters$epsilon)
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
            parameters_output <- data.frame (kernel = svm_kernel, cost = svm_tuning$best.model$cost, degree = svm_tuning$best.model$degree, epsilon = svm_tuning$best.model$epsilon, gamma = svm_tuning$best.model$gamma)
        }
        #################### SVM with defined parameters
        if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
            svm_model <- svm(pca_training, factor(peaklist_training$Class), kernel = parameters$kernel, cost = parameters$cost, epsilon = parameters$epsilon, gamma = parameters$gamma)
            parameters_output <- list (kernel = parameters$kernel, cost = parameters$cost, degree = parameters$degree, epsilon = parameters$epsilon, gamma = parameters$gamma)
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
            classification_results_svm <- data.frame (sample = peaklist_test$Sample, predicted = predicted_classes_svm, true = peaklist_test$Class)
            ### Generate the confusion matrix to evaluate the performances
            test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive = positive_class_cv)
            #### ROC analysis
            svm_roc <- list()
            roc_curve <- roc(response = classification_results_svm$true, predictor = as.numeric(classification_results_svm$predicted))
            svm_roc[[1]] <- roc_curve$auc
            plot (roc_curve)
            roc_legend <- paste ("ROC area under the curve:", roc_curve$auc)
            legend ("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
            svm_roc[[2]] <- recordPlot()
            # Output the results
            return (list (model = svm_model, svm_features = colnames(svm_model$SV), svm_parameters = parameters_output, classification_results = classification_results_svm, performances = test_performances_svm, roc = svm_roc))
        }    else {return (list (model = svm_model, svm_features = colnames(svm_model$SV), svm_parameters = parameters_output, cross_validation = cv_svm_model))}
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
            svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma = tuning_parameters$gamma, cost = tuning_parameters$cost, kernel = tuning_parameters$kernel, epsilon = tuning_parameters$epsilon, degree = tuning_parameters$degree)
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
            parameters_output <- data.frame (kernel = svm_kernel, cost = svm_tuning$best.model$cost, degree = svm_tuning$best.model$degree, epsilon = svm_tuning$best.model$epsilon, gamma = svm_tuning$best.model$gamma)
        }
        #################### SVM with defined parameters
        if (autotuning == FALSE || is.null(tuning_parameters) || length(tuning_parameters) == 0) {
            svm_model <- svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), kernel = parameters$kernel, cost = parameters$cost, epsilon = parameters$epsilon, gamma = parameters$gamma)
            parameters_output <- data.frame (kernel = parameters$kernel, cost = parameters$cost, degree = parameters$degree, epsilon = parameters$epsilon, gamma = parameters$gamma)
        }
        #################### CROSS-VALIDATION
        cv_svm_model <- cross_validation_svm(peaklist_training, k_fold = k_fold_cv, repeats = repeats_cv, svm_model = svm_model, non_features = non_features, positive_class = positive_class_cv, seed = seed)
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
            classification_results_svm <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_svm, true = peaklist_test$Class)
            ### Generate the confusion matrix to evaluate the performances
            test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive = positive_class_cv)
            #### ROC analysis
            svm_roc <- list()
            roc_curve <- roc(response = as.numeric(classification_results_svm$true), predictor = as.numeric(classification_results_svm$predicted))
            svm_roc[[1]] <- roc_curve$auc
            plot (roc_curve)
            roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
            legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
            pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
            pie_chart_classification <- recordPlot()
            # Output the results
            return (list(model = svm_model, svm_features = colnames(svm_model$SV), cross_validation = cv_svm_model, classification_results = classification_results_svm, parameters = parameters_output, performances = test_performances_svm, roc = svm_roc, pie_chart_classification = pie_chart_classification))
        } else {return (list(model = svm_model, svm_features = colnames(svm_model$SV), parameters = parameters_output, cross_validation = cv_svm_model))}
    }
}





################################################################################





###################################################### SVM TUNING AND VALIDATION
# This function operates the tuning of the Support Vector Machine (SVM), by testing all the parameters of the SVM (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
svm_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(sigma = 10^(-5:5), cost = 10^(-5:5), epsilon = seq(1,2,by = 1), degree = 1:5, scale = 1, kernel = "radial"), k_fold_cv = 10, repeats_cv = 2, preprocessing = c("center","scale"), parameters = list(sigma = 0.001, scale = 1, gamma = 0.1, cost = 10, epsilon = 0.1, degree = 3, kernel = "radial"), positive_class_cv = "HP", seed = NULL, evaluation_method = "Accuracy") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "pROC", "kernlab"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train and tune the model
        if (tuning_parameters$kernel == "radial") {
            if (!is.null(tuning_parameters$sigma)) {
                # Define the tune grid for the model
                svm_tune_grid <- expand.grid(sigma = tuning_parameters$sigma, C = tuning_parameters$cost)
                # Training and tuning
                svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmRadial", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
            } else {
                # Define the tune grid for the model
                svm_tune_grid <- expand.grid(C = tuning_parameters$cost)
                # Training and tuning
                svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmRadialCost", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
            }
        } else if (tuning_parameters$kernel == "polynomial") {
            # Define the tune grid for the model
            svm_tune_grid <- expand.grid(scale = tuning_parameters$scale, C = tuning_parameters$cost, degree = tuning_parameters$degree)
            # Training and tuning
            svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmPoly", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
        } else if (tuning_parameters$kernel == "linear") {
            # Define the tune grid for the model
            svm_tune_grid <- expand.grid(cost = tuning_parameters$cost, gamma = tuning_parameters$gamma)
            # Training and tuning
            svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmLinear2", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
        }
        # Extract the cross-validation information
        cv_svm_model <- confusionMatrix(svm_tuning)
        # Select automatically the best model after the tuning
        svm_model <- svm_tuning$finalModel
        # Parameters output
        if (tuning_parameters$kernel == "radial") {
            parameters_output <- data.frame(kernel = tuning_parameters$kernel, cost = svm_model@param$C, sigma = svm_model@kernelf@kpar[[1]])
        } else if (tuning_parameters$kernel == "polynomial") {
            parameters_output <- data.frame(kernel = tuning_parameters$kernel, cost = svm_model@param$C, degree = svm_model@kernelf@kpar$degree, scale = svm_model@kernelf@kpar$scale)
        } else if (tuning_parameters$kernel == "linear") {
            parameters_output <- data.frame(kernel = tuning_parameters$kernel, cost = svm_model$cost, gamma = svm_model$gamma)
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train the model
        if (tuning_parameters$kernel == "radial") {
            if (!is.null(tuning_parameters$sigma)) {
                # Define the tune grid for the model
                svm_tune_grid <- expand.grid(sigma = parameters$sigma, C = parameters$cost)
                # Training and tuning
                svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmRadial", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
            } else {
                # Define the tune grid for the model
                svm_tune_grid <- expand.grid(C = parameters$cost)
                # Training and tuning
                svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmRadialCost", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
            }
        } else if (tuning_parameters$kernel == "polynomial") {
            # Define the tune grid for the model
            svm_tune_grid <- expand.grid(scale = parameters$scale, C = parameters$cost, degree = parameters$degree)
            # Training and tuning
            svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmPoly", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
        } else if (tuning_parameters$kernel == "linear") {
            # Define the tune grid for the model
            svm_tune_grid <- expand.grid(cost = parameters$cost, gamma = parameters$gamma)
            # Training and tuning
            svm_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "svmLinear2", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = svm_tune_grid)
        }
        # Extract the cross-validation information
        cv_svm_model <- confusionMatrix(svm_tuning)
        # Select automatically the best model after the tuning
        svm_model <- svm_tuning$finalModel
        # Parameters output
        if (tuning_parameters$kernel == "radial") {
            parameters_output <- data.frame(kernel = parameters$kernel, cost = svm_model@param$C, sigma = svm_model@kernelf@kpar[[1]])
        } else if (tuning_parameters$kernel == "polynomial") {
            parameters_output <- data.frame(kernel = parameters$kernel, cost = svm_model@param$C, degree = svm_model@kernelf@kpar$degree, scale = svm_model@kernelf@kpar$scale)
        } else if (tuning_parameters$kernel == "linear") {
            parameters_output <- data.frame(kernel = parameters$kernel, cost = svm_model$cost, gamma = svm_model$gamma)
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
        classification_results_svm <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_svm, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        svm_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_svm$true), predictor = as.numeric(classification_results_svm$predicted))
        svm_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        # Output the results
        return (list(model = svm_model, svm_features = colnames(svm_model@xmatrix[[1]]), cross_validation = cv_svm_model, classification_results = classification_results_svm, parameters = parameters_output, tuning_plot = tuning_plot, performances = test_performances_svm, roc = svm_roc, pie_chart_classification = pie_chart_classification))
    } else {return (list(model = svm_model, svm_features = colnames(svm_model@xmatrix[[1]]), parameters = parameters_output, tuning_plot = tuning_plot, cross_validation = cv_svm_model))}
}





################################################################################





###################################################### PLS TUNING AND VALIDATION
# This function operates the training of a Partial Least Squares (PLS) model, by testing all the parameters of the PLS and selecting the best (in terms of accuracy). The tuning is performed via cross-validation onto the training dataset.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
pls_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), tuning_parameters = data.frame(ncomp = 1:5), k_fold_cv = 10, repeats_cv = 2, positive_class_cv = "HP", seed = NULL, preprocessing = c("center","scale"), selection_criteria = "Accuracy", maximise_selection_criteria_values = TRUE) {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "e1071"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
    train_control_pls <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv)
    if (!is.null(seed)) {
        # Make the randomness reproducible
        set.seed(seed)
    }
    pls_model <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "pls", trControl = train_control_pls, preProcess = preprocessing, metric = selection_criteria, maximize = maximise_selection_criteria_values, tuneGrid = tune_grid)
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
        classification_results_pls <- data.frame (sample = peaklist_test$Sample, predicted = predicted_classes_pls, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_pls <- confusionMatrix(data = predicted_classes_pls, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        pls_roc <- list()
        roc_curve <- roc(response = classification_results_pls$true, predictor = as.numeric(classification_results_pls$predicted))
        pls_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie(x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        plots <- append(plots, pie_chart_classification)
    }
    # Output the results
    return (list(model = pls_model, classification_results = classification_results_pls, performances_pls_cv = pls_performances, pls_performances = test_performances_pls, roc = pls_roc, plots = plots, pls_features = pls_model$finalModel$xNames))
}





################################################################################





################################### NAIVE BAYES CLASSIFIER TUNING AND VALIDATION
# This function operates the training of a Naive Bayes Classifier (NBC) model, by testing all the parameters of the NBC and selecting the best (in terms of accuracy). The tuning is performed via cross-validation onto the training dataset.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
nbc_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), tuning_parameters = data.frame(fL = NULL,usekernel = NULL), k_fold_cv = 10, repeats_cv = 2, positive_class_cv = "HP", seed = NULL, preprocessing = c("center","scale"), selection_criteria = "Accuracy", maximise_selection_criteria_values = TRUE) {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "e1071", "klaR", "MASS"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
    train_control_nbc <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv)
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
    nbc_model <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "nb", trControl = train_control_nbc, preProcess = preprocessing, metric = selection_criteria, maximize = maximise_selection_criteria_values, tuneGrid = tune_grid)
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
        classification_results_nbc <- data.frame (sample = peaklist_test$Sample, predicted = predicted_classes_nbc, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_nbc <- confusionMatrix(data = predicted_classes_nbc, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        nbc_roc <- list()
        roc_curve <- roc(response = classification_results_nbc$true, predictor = as.numeric(classification_results_nbc$predicted))
        nbc_roc[[1]] <- roc_curve$auc
        plot(roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie(x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        plots <- append(plots, pie_chart_classification)
    }
    # Output the results
    return (list(model = nbc_model, classification_results = classification_results_nbc, performances_nbc_cv = nbc_performances, nbc_performances = test_performances_nbc, roc = nbc_roc, plots = plots, nbc_features = nbc_model$finalModel$xNames))
}





################################################################################





###################################################### KNN TUNING AND VALIDATION
# This function operates the tuning of the k-Nearest Neighbour (KNN), by testing all the parameters of the KNN (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
knn_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(k = seq(1,15, by = 1)), k_fold_cv = 10, repeats_cv = 2, preprocessing = c("center","scale"), parameters = list(k = 15), positive_class_cv = "HP", seed = NULL, evaluation_method = "Accuracy") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "pROC"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train and tune the model
        # Define the tune grid for the model
        knn_tune_grid <- expand.grid(k = tuning_parameters$k)#, l = tuning_parameters$l)
        # Training and tuning
        knn_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "knn", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = knn_tune_grid)
        # Extract the cross-validation information
        cv_knn_model <- confusionMatrix(knn_tuning)
        # Select automatically the best model after the tuning
        knn_model <- knn_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(k = knn_model$k)
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train the model
        # Define the tune grid for the model
        knn_tune_grid <- expand.grid(k = parameters$k)
        # Training and tuning
        knn_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "knn", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = knn_tune_grid)
        # Extract the cross-validation information
        cv_knn_model <- confusionMatrix(knn_tuning)
        # Select automatically the best model after the tuning
        knn_model <- knn_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(k = parameters$k)
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
        predicted_classes_knn <- predict(knn_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type = "class")
        # Create the outcomes dataframe
        classification_results_knn <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_knn, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_knn <- confusionMatrix(data = predicted_classes_knn, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        knn_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_knn$true), predictor = as.numeric(classification_results_knn$predicted))
        knn_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        # Output the results
        return (list(model = knn_model, knn_features = knn_model$xNames, cross_validation = cv_knn_model, classification_results = classification_results_knn, parameters = parameters_output, tuning_plot = tuning_plot, performances = test_performances_knn, roc = knn_roc, pie_chart_classification = pie_chart_classification))
    } else {return (list(model = knn_model, knn_features = knn_model$xNames, parameters = parameters_output, tuning_plot = tuning_plot, cross_validation = cv_knn_model))}
}





################################################################################





###################################################### LDA TUNING AND VALIDATION
# This function operates the tuning of the Linear Discriminant Analysis (LDA), by testing all the parameters of the LDA (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
lda_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(dimen = seq(1,15, by = 1)), k_fold_cv = 10, repeats_cv = 2, preprocessing = c("center","scale"), parameters = list(dimen = 15), positive_class_cv = "HP", seed = NULL, evaluation_method = "Accuracy") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "pROC"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train and tune the model
        # Define the tune grid for the model
        lda_tune_grid <- expand.grid(dimen = tuning_parameters$dimen)
        # Training and tuning
        lda_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "lda2", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = lda_tune_grid)
        # Extract the cross-validation information
        cv_lda_model <- confusionMatrix(lda_tuning)
        # Select automatically the best model after the tuning
        lda_model <- lda_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(dimen = lda_model$tuneValue$dimen)
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train the model
        # Define the tune grid for the model
        lda_tune_grid <- expand.grid(dimen = parameters$dimen)
        # Training and tuning
        lda_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "lda2", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = lda_tune_grid)
        # Extract the cross-validation information
        cv_lda_model <- confusionMatrix(lda_tuning)
        # Select automatically the best model after the tuning
        lda_model <- lda_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(dimen = parameters$dimen)
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
        predicted_classes_lda <- predict(lda_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type = "class")$class
        # Create the outcomes dataframe
        classification_results_lda <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_lda, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_lda <- confusionMatrix(data = predicted_classes_lda, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        lda_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_lda$true), predictor = as.numeric(classification_results_lda$predicted))
        lda_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        # Output the results
        return (list(model = lda_model, lda_features = lda_model$xNames, cross_validation = cv_lda_model, classification_results = classification_results_lda, parameters = parameters_output, tuning_plot = tuning_plot, performances = test_performances_lda, roc = lda_roc, pie_chart_classification = pie_chart_classification))
    } else {return (list(model = lda_model, lda_features = lda_model$xNames, parameters = parameters_output, tuning_plot = tuning_plot, cross_validation = cv_lda_model))}
}





################################################################################





###################################################### RF TUNING AND VALIDATION
# This function operates the tuning of the Random Forest (RF), by testing all the parameters of the RF (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
rf_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(mtry = seq(1,5, by = 1)), k_fold_cv = 10, repeats_cv = 2, preprocessing = c("center","scale"), parameters = list(mtry = round(sqrt(ncol(peaklist_training)))), positive_class_cv = "HP", seed = NULL, evaluation_method = "Accuracy") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "pROC","randomForest"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train and tune the model
        # Define the tune grid for the model
        rf_tune_grid <- expand.grid(mtry = tuning_parameters$mtry)
        # Training and tuning
        rf_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "rf", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = rf_tune_grid)
        # Extract the cross-validation information
        cv_rf_model <- confusionMatrix(rf_tuning)
        # Select automatically the best model after the tuning
        rf_model <- rf_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(mtry = rf_model$mtry)
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train the model
        # Define the tune grid for the model
        rf_tune_grid <- expand.grid(mtry = parameters$mtry)
        # Training and tuning
        rf_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "rf", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = rf_tune_grid)
        # Extract the cross-validation information
        cv_rf_model <- confusionMatrix(rf_tuning)
        # Select automatically the best model after the tuning
        rf_model <- rf_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(mtry = parameters$mtry)
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
        predicted_classes_rf <- predict(rf_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type = "class")
        # Create the outcomes dataframe
        classification_results_rf <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_rf, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_rf <- confusionMatrix(data = predicted_classes_rf, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        rf_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_rf$true), predictor = as.numeric(classification_results_rf$predicted))
        rf_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        # Output the results
        return (list(model = rf_model, rf_features = rf_model$xNames, cross_validation = cv_rf_model, classification_results = classification_results_rf, parameters = parameters_output, tuning_plot = tuning_plot, performances = test_performances_rf, roc = rf_roc, pie_chart_classification = pie_chart_classification))
    } else {return (list(model = rf_model, rf_features = rf_model$xNames, parameters = parameters_output, tuning_plot = tuning_plot, cross_validation = cv_rf_model))}
}





################################################################################





###################################################### NN TUNING AND VALIDATION
# This function operates the tuning of the Neural Network (NN), by testing all the parameters of the NN (choosing them from a list provided by the user) and selecting the best.
# It returns the best model in terms of classification performances, along with its parameters and its performances (cross-validation or external validation, according to if an external dataset is provided).
#This function uses the functions from CARET only
nn_tuning_and_validation <- function (peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), autotuning = TRUE, tuning_parameters = list(size = seq(1,5, by = 1),decay = seq(0,2,by = 1)), k_fold_cv = 10, repeats_cv = 2, preprocessing = c("center","scale"), parameters = list(size = 5, decay = 0), positive_class_cv = "HP", seed = NULL, evaluation_method = "Accuracy") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "pROC","nnet"))
    ### PARALLEL BACKEND
    # Detect the number of cores
    cpu_thread_number <- detectCores(logical = TRUE)
    cpu_thread_number <- cpu_thread_number / 2
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train and tune the model
        # Define the tune grid for the model
        nn_tune_grid <- expand.grid(size = tuning_parameters$size, decay = tuning_parameters$decay)
        # Training and tuning
        nn_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "nnet", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = nn_tune_grid)
        # Extract the cross-validation information
        cv_nn_model <- confusionMatrix(nn_tuning)
        # Select automatically the best model after the tuning
        nn_model <- nn_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(size = nn_model$tuneValue$size, decay = nn_model$decay)
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
        training_ctrl <- trainControl(method = "repeatedcv", number = k_fold_cv, repeats = repeats_cv, classProbs = TRUE)#, seeds = seeds) #, summaryFunction = twoClassSummary)
        # Train the model
        # Define the tune grid for the model
        nn_tune_grid <- expand.grid(size = parameters$size, decay = parameters$decay)
        # Training and tuning
        nn_tuning <- train(x = peaklist_training [,!(names(peaklist_training) %in% non_features)], y = factor(peaklist_training$Class), method = "nnet", preProcess = preprocessing, metric = evaluation_method, trControl = training_ctrl, tuneGrid = nn_tune_grid)
        # Extract the cross-validation information
        cv_nn_model <- confusionMatrix(nn_tuning)
        # Select automatically the best model after the tuning
        nn_model <- nn_tuning$finalModel
        # Parameters output
        parameters_output <- data.frame(size = parameters$size, decay = parameters$decay)
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
        predicted_classes_nn <- predict(nn_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)], type = "class")
        # Create the outcomes dataframe
        classification_results_nn <- data.frame(sample = peaklist_test$Sample, predicted = predicted_classes_nn, true = peaklist_test$Class)
        ### Generate the confusion matrix to evaluate the performances
        test_performances_nn <- confusionMatrix(data = predicted_classes_nn, peaklist_test$Class, positive = positive_class_cv)
        #### ROC analysis
        nn_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_nn$true), predictor = as.numeric(classification_results_nn$predicted))
        nn_roc[[1]] <- roc_curve$auc
        plot (roc_curve)
        roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
        legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
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
        pie (x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
        pie_chart_classification <- recordPlot()
        # Output the results
        return (list(model = nn_model, nn_features = nn_model$xNames, cross_validation = cv_nn_model, classification_results = classification_results_nn, parameters = parameters_output, tuning_plot = tuning_plot, performances = test_performances_nn, roc = nn_roc, pie_chart_classification = pie_chart_classification))
    } else {return (list(model = nn_model, nn_features = nn_model$xNames, parameters = parameters_output, tuning_plot = tuning_plot, cross_validation = cv_nn_model))}
}





################################################################################






################################################# ENSEMBLE TUNING AND VALIDATION
# This function operates the tuning of the ensemble classifier, by relying upon other functions to train, tune and validate the individual classifiers.
# It returns the best models in terms of classification performances, along with their parameters and performances (cross-validation or external validation, according to if an external dataset is provided).
ensemble_tuning_and_validation <- function(peaklist_training, peaklist_test = NULL, non_features = c("Sample","Class","THY"), classifiers = c("svm","pls","bayes"), autotuning = TRUE,  classifier_parameters = list(svm = list(gamma = 10^(-5:5), cost = 10^(-5:5), epsilon = seq(1,2,by = 1), degree = 1:5, kernel = "radial"), pls = data.frame(ncomp = 1:5), bayes = NULL), k_fold_cv = 10, repeats_cv = 2, positive_class_cv = "HP", seed = NULL, preprocessing = c("center","scale"), selection_criteria = data.frame(pls = "Accuracy",bayes = "Accuracy"), maximise_selection_criteria_values = data.frame(pls = TRUE,bayes = TRUE)) {
    if ("svm" %in% classifiers || "SVM" %in% classifiers) {
        svm_classifier <- svm_tuning_and_validation(peaklist_training = peaklist_training, peaklist_test = peaklist_test, non_features = non_features, autotuning = autotuning, tuning_parameters = classifier_parameters$svm, k_fold_cv = k_fold_cv, repeats_cv = repeats_cv, parameters = classifier_parameters$svm, positive_class_cv = positive_class_cv, seed = seed)
    }
    if ("pls" %in% classifiers || "PLS" %in% classifiers) {
        pls_classifier <- pls_tuning_and_validation(peaklist_training = peaklist_training, peaklist_test = peaklist_test, non_features = non_features, tuning_parameters = classifier_parameters$pls, k_fold_cv = k_fold_cv, repeats_cv = repeats_cv, positive_class_cv = "HP", seed = seed, preprocessing = preprocessing, selection_criteria = as.character(selection_criteria$pls), maximise_selection_criteria_values = maximise_selection_criteria_values$pls)
    }
    if ("bayes" %in% classifiers || "Bayes" %in% classifiers) {
        bayes_classifier <- nbc_tuning_and_validation(peaklist_training = peaklist_training, peaklist_test = peaklist_test, non_features = non_features, tuning_parameters = classifier_parameters$bayes, k_fold_cv = k_fold_cv, repeats_cv = repeats_cv, positive_class_cv = positive_class_cv, seed = seed, preprocessing = preprocessing, selection_criteria = as.character(selection_criteria$bayes), maximise_selection_criteria_values = maximise_selection_criteria_values$bayes)
    }
    return(list(svm_classifier = svm_classifier, pls_classifier = pls_classifier, bayes_classifier = bayes_classifier))
}





################################################################################





################################################ ROUND NUMERIC FEATURES PEAKLIST
# This function rounds the numeric features to a certain decimal digit.
# It returns the same matrix with the rounded features, along with the original peaklist matrix.
round_features_peaklist <- function(peaklist, decimal_digits = 3, non_features = c("Sample","Class","THY")) {
    # Take the features
    feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
    # Round the features
    for (i in 1:length(feature_list)) {
        feature_list[i] <- round(feature_list[i], digits = decimal_digits)
    }
    # Regather the features
    features <- c(feature_list, non_features)
    peaklist2 <- peaklist
    names(peaklist2) <- features
    #
    return (list(original_peaklist = peaklist, peaklist_rounded = peaklist2))
}





################################################################################





