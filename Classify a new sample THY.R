###INSTALL THE REQUIRED PACKAGES
required_packages <- c ("parallel", "MALDIquant", "MALDIquantForeign", "caret", "tcltk", "beepr", "kernlab", "e1071", "doMC")

install_and_load_required_packages (required_packages)



tkmessageBox (title = "Patient Classification", message = "Classify unknown patients", icon = "warning")

# Path where the imzML file is located (GUI)
filepath_test_select <- tkmessageBox (title = "Path", message = "Select the folder containing the imzML files to be imported", icon = "info")
filepath_test <- tclvalue(tkchooseDirectory())
if (!nchar(filepath_test)) {
    tkmessageBox(message = "No folder selected!")
}	else {
    tkmessageBox(message = paste("The directory selected is", filepath_test))
}

# Path where the R workspace file is located (GUI) (the one with the model)
tkmessageBox (title = "R workspace", message = "Select the R workspace with the model to be imported", icon = "info")
filepath_R <- tclvalue(tkgetOpenFile())
if (!nchar(filepath_R)) {
	tkmessageBox(message = "No R workspace file was selected!")
}	else {
	tkmessageBox(message = paste("The R workspace file selected is", filepath_R))
}
## Define the path where to save output files
filepath_output_select <- tkmessageBox (title = "Path", message = "Select the folder where to save all the outputs", icon = "info")
filepath_output <- tclvalue(tkchooseDirectory())
if (!nchar(filepath_output)) {
    tkmessageBox(message = "No folder selected!")
}	else {
    tkmessageBox(message = paste("All the output files will be saved in", filepath_output))
}



### LOAD THE R WORKSPACE
# Create a temporary environment
temporary_environment <- new.env()
load (filepath_R, envir=temporary_environment)
svm_model <- get ("SVMModel", pos=temporary_environment)


##### MULTICORE
# Detect the number of cores
cpu_thread_number <- detectCores (logical=TRUE)
cpu_core_number <- cpu_thread_number/2
# Use only the 3/4 of this power
#CPUcoreNumber <- floor (CPUcoreNumber *3/4)
# Register the foreach backend
registerDoMC (cores = cpu_core_number)








#### CLASSIFICATION FUNCTION
classify_patients <- function (spectra_folder, svm_model, smoothing_strength_preprocessing="medium", preprocess_spectra_in_packages_of=length(sample_spectra), mass_range=c(4000,15000)) {
# Class list
class_list <- levels(factor(svm_model$levels))
## List the imzML files in the selected folder
filepath_test_imzml <- read_spectra_files (spectra_folder, file_format="imzml", full_path=TRUE)
########################################## Outputs
#classification_output <- list()
final_result_matrix <- NULL
final_result_matrix_avg <- NULL
classification_msi <- list()
average_spectra_with_bars <- list()
############################ The sample peaklist must have the exact same features that are used for the model, so we need to align the sample features to the one of the model, discard the features that are in the sample but not in the model and add the features that are in the model but not in the sample.
##################################### Isolate the peaks used to create the model
features_model <- colnames(svm_model$SV)
# Remove the X
for (f in 1:length(features_model)) {
	name_splitted <- unlist(strsplit(features_model[f],""))
	feature_def <- name_splitted [2]
	for (i in 3:length(name_splitted)) {
		feature_def <- paste(feature_def, name_splitted[i], sep="")
	}
	features_model[f] <- feature_def
}
# Process the sample to be classified
###### For each imzML file...
for (p in 1:length(filepath_test_imzml)) {
################ SPECTRA IMPORT AND PROCESSING
# Import the spectra (one imzML at a time)
sample_spectra <- importImzMl(filepath_test_imzml[p])
# Trim the spectra
if (mass_range[1] == 0 && mass_range[2] == 0) {
	sample_spectra <- trim(sample_spectra)
}
if (mass_range[1] == 0 && mass_range[2] != 0) {
	sample_spectra <- trim(sample_spectra, range=mass_range)
}
if (mass_range[1] != 0 && mass_range[2] == 0) {
	mass_range[2] <- Inf
	sample_spectra <- trim(sample_spectra, range=mass_range)
}
if (mass_range[1] != 0 && mass_range[2] != 0) {
	sample_spectra <- trim(sample_spectra, range=mass_range)
}
################################################# PIXEL BY PIXEL CLASSIFICATION
sample_spectra <- preprocess_spectra(sample_spectra, tof_mode="linear", smoothing_strength=smoothing_strength_preprocessing, process_in_packages_of=preprocess_spectra_in_packages_of)
# Peak picking and alignment
sample_peaks <- detectPeaks(sample_spectra, method="MAD", SNR=5)
sample_peaks <- align_and_filter_peaks(sample_peaks, tolerance_ppm=2000, peaks_filtering=TRUE, frequency_threshold=0.25, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1)
# Generate the intensity matrix for classification
sample_matrix <- intensityMatrix(sample_peaks, sample_spectra)
########################### Determine the columns to keep and the column to add
features_to_keep <- numeric()
features_to_add <- numeric()
adjusted_features_to_keep <- numeric()
# For each feature in the model
for (feature_model in features_model) {
	# Set the default presence of the signal in the sample to FALSE
	presence <- FALSE
	# Scroll the sample features
	for (feature_sample in colnames(sample_matrix)) {
		# If there is a match
		if (abs((as.numeric(feature_model)-as.numeric(feature_sample))*10^6/as.numeric(feature_model)) <= 2000) {
			# Add it to the features to keep
			features_to_keep <- append(features_to_keep, feature_sample)
			# Align the feature in the sample with the one in the model
			feature_sample <- feature_model
			# Add it to another list (it will be used to adjust the column names in the final sample peaklist)
			adjusted_features_to_keep <- append(adjusted_features_to_keep, feature_sample)
			# Set the presence of the signal in the sample to TRUE
			presence <- TRUE
			# Avoid consecutive duplicates (once it is found there is no point in keep going)
			break
		}
	}
	# If after all the signal in the model is not found in the sample
	if (presence == FALSE) {
		# Add this to the features to be added
		features_to_add <- append(features_to_add, feature_model)
	}
}
# Generate the final sample matrix (with the right column names)
final_sample_matrix <- sample_matrix [,features_to_keep]
colnames(final_sample_matrix) <- adjusted_features_to_keep
# Scroll the features to add (in the model but not in the sample)
for (f in features_to_add) {
	# Create the empty vector in which the intensity of that peak in the sample spectra will be stored
	feature_intensity_vector <- numeric()
	# Scroll the spectra
	for (s in 1:length(sample_spectra)) {
		# Scroll the mass list of each spectrum
		for (m in 1:length(sample_spectra[[s]]@mass)) {
			# If there is a match
			if (abs(sample_spectra[[s]]@mass[m]-as.numeric(f))*10^6/as.numeric(f) <= 200) {
				# Add the corresponding intensity to the vector
				feature_intensity_vector <- append(feature_intensity_vector, sample_spectra[[s]]@intensity[m])
				# Break the for cycle to avoid dupicates and avoid going further once a match has been found
				break
			}
		}
	}
	# Generate the matrix column
	feature_intensity_column <- cbind(feature_intensity_vector)
	colnames(feature_intensity_column) <- f
	# Attach this to the final sample matrix
	final_sample_matrix <- cbind(final_sample_matrix, feature_intensity_column)
}
# Put the X at the beginning of the peak names
for (n in 1:length(colnames(final_sample_matrix))) {
	name <- paste("X", colnames(final_sample_matrix)[n], sep="")
	colnames(final_sample_matrix)[n] <- name
}
######################################################## Classify pixel by pixel
# Predictions (spectra by spectra)
predicted_classes <- predict(svm_model, newdata = final_sample_matrix)
# Generate a matrix with the results
result_matrix <- matrix (nrow=length(predicted_classes), ncol=2)
sample_name <- unlist(strsplit(filepath_test_imzml[p],"/"))
sample_name <- sample_name[length(sample_name)]
sample_name <- unlist(strsplit(sample_name, ".imzML"))
sample_name <- sample_name[1]
result_matrix [,1] <- cbind(rep(sample_name, length(sample_spectra)))
result_matrix [,2] <- cbind(predicted_classes)
colnames(result_matrix) <- c("Sample", "Predicted Class")
# Turn the numbers in the class column to the corresponding class name
for (l in 1:length(class_list)) {
	for (z in 1:length(result_matrix[,2])) {
		if (result_matrix[,2][z] == l) {
			result_matrix[,2][z] <- class_list[l]
		}
	}
}
#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
if (is.null(final_result_matrix)) {
	final_result_matrix <- result_matrix
} else {
	final_result_matrix <- rbind(final_result_matrix, result_matrix)
}
# Add the list of output for this patient to the final list of outputs
#classification_output <- append(classification_output, result_matrix)
############################### Generate a molecular image of the classification
# Replace the spectra intensities with the class number for plotting purposes
class_as_number <- as.numeric(predicted_classes)
spectra_for_plotting <- sample_spectra
for (s in 1:length(spectra_for_plotting)) {
	spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s], length(spectra_for_plotting[[s]]@intensity))
}
plotMsiSlice (spectra_for_plotting, center=spectra_for_plotting[[1]]@mass[1], tolerance=4, legend=FALSE)
legend(x="bottomright", legend=class_list, fill=c("green","red"))
legend(x="topright", legend=sample_name)
# Store the plot into the list of images
classification_msi[[p]] <- recordPlot()
#################################################### CLASSIFICATION OF AVERAGES
sample_spectra_avg <- averageMassSpectra(sample_spectra, method="mean")
sample_spectra_avg <- smoothIntensity(sample_spectra_avg, method="SavitzkyGolay")
sample_spectra_avg <- removeBaseline(sample_spectra_avg, method="TopHat")
sample_spectra_avg <- calibrateIntensity(sample_spectra_avg, method="TIC")
# Peak picking
sample_peaks_avg <- detectPeaks(sample_spectra_avg, method="MAD", SNR=5)
# Generate the intensity matrix for classification
sample_matrix_avg <- matrix (nrow=1, ncol=length(sample_peaks_avg@mass))
colnames(sample_matrix_avg) <- sample_peaks_avg@mass
sample_matrix_avg[1,] <- sample_peaks_avg@intensity
########################### Determine the columns to keep and the column to add
features_to_keep_avg <- numeric()
features_to_add_avg <- numeric()
adjusted_features_to_keep_avg <- numeric()
# For each feature in the model
for (feature_model in features_model) {
	# Set the default presence of the signal in the sample to FALSE
	presence <- FALSE
	# Scroll the sample features
	for (feature_sample_avg in colnames(sample_matrix_avg)) {
		# If there is a match
		if (abs((as.numeric(feature_model)-as.numeric(feature_sample_avg))*10^6/as.numeric(feature_model)) <= 2000) {
			# Add it to the features to keep
			features_to_keep_avg <- append(features_to_keep_avg, feature_sample_avg)
			# Align the feature in the sample with the one in the model
			feature_sample_avg <- feature_model
			# Add it to another list (it will be used to adjust the column names in the final sample peaklist)
			adjusted_features_to_keep_avg <- append(adjusted_features_to_keep_avg, feature_sample_avg)
			# Set the presence of the signal in the sample to TRUE
			presence <- TRUE
			# Avoid consecutive duplicates (once it is found there is no point in keep going)
			break
		}
	}
	# If after all the signal in the model is not found in the sample
	if (presence == FALSE) {
		# Add this to the features to be added
		features_to_add_avg <- append(features_to_add_avg, feature_model)
	}
}
# Generate the final sample matrix (with the right column names)
final_sample_matrix_avg <- as.matrix(rbind(sample_matrix_avg [,features_to_keep_avg]))
colnames(final_sample_matrix_avg) <- adjusted_features_to_keep_avg
# Scroll the features to add (in the model but not in the sample)
for (f in features_to_add_avg) {
	# Create the empty vector in which the intensity of that peak in the sample spectra will be stored
	feature_intensity_vector_avg <- numeric()
	# Scroll the mass list of the mean spectrum
	for (m in 1:length(sample_spectra_avg@mass)) {
		# If there is a match
		if (abs(sample_spectra_avg@mass[m]-as.numeric(f))*10^6/as.numeric(f) <= 2000) {
			# Add the corresponding intensity to the vector
			feature_intensity_vector_avg <- append(feature_intensity_vector_avg, sample_spectra_avg@intensity[m])
			# Break the for cycle to avoid dupicates and avoid going further once a match has been found
			break
		}
	}
	# Generate the matrix column
	feature_intensity_column_avg <- as.matrix(cbind(feature_intensity_vector_avg))
	colnames(feature_intensity_column_avg) <- f
	# Attach this to the final sample matrix
	final_sample_matrix_avg <- cbind(final_sample_matrix_avg, feature_intensity_column_avg)
}
# Put the X at the beginning of the peak names
for (n in 1:length(colnames(final_sample_matrix_avg))) {
	name <- paste("X", colnames(final_sample_matrix_avg)[n], sep="")
	colnames(final_sample_matrix_avg)[n] <- name
}
# Predictions (spectra by spectra)
predicted_classes_avg <- predict(svm_model, newdata = final_sample_matrix_avg)
# Generate a matrix with the results
result_matrix_avg <- matrix (nrow=1, ncol=2)
result_matrix_avg [,1] <- sample_name
result_matrix_avg [,2] <- predicted_classes_avg
colnames(result_matrix_avg) <- c("Sample", "Predicted Class")
# Turn the numbers in the class column to the corresponding class name
for (l in 1:length(class_list)) {
	for (z in 1:length(result_matrix_avg[,2])) {
		if (result_matrix_avg[,2][z] == l) {
			result_matrix_avg[,2][z] <- class_list[l]
		}
	}
}
#### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
if (is.null(final_result_matrix_avg)) {
	final_result_matrix_avg <- result_matrix_avg
} else {
	final_result_matrix_avg <- rbind(final_result_matrix_avg, result_matrix_avg)
}
################# Average spectrum with bars onto the signals used by the model
# Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
# Detect peaks in the avg (SNR=1)
sample_peaks_avg_for_bars <- detectPeaks(sample_spectra_avg, method="MAD", SNR=3)
# Determine the coordinates of the bars
coordinates_of_bars <- list(x=numeric(), y=numeric())
# Check if the features used for the model are in the spectrum
for (f in features_model) {
	presence_in_the_avg <- FALSE
	# Scroll the peaks
	for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
		# If there is a match...
		if ((abs(sample_peaks_avg_for_bars@mass[z]-as.numeric(f))*10^6/as.numeric(f)) <= 2000) {
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
		for (j in 1:length(sample_spectra_avg@mass)) {
			# If there is a match...
			if ((abs(sample_spectra_avg@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= 2000) {
				# Add the intensity of this peak to the y coordinates of the bars
				coordinates_of_bars$x = append(coordinates_of_bars$x,sample_spectra_avg@mass[j])
				coordinates_of_bars$y = append(coordinates_of_bars$y,sample_spectra_avg@intensity[j])
				# Break the for cycle to avoid duplicates and to continue
				break
			}
		}
	}
}
plot(sample_spectra_avg, xlab="m/z", ylab="Intensity (a.i.)")
legend(x="topright", legend=sample_name)
# Draw the bars
for (s in 1:length(coordinates_of_bars$x)) {
	# Vertical bars (x,y x,y)
	segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col="red", lwd=2)
	# Horizontal segments(x,y , x,y)
	segments(coordinates_of_bars$x[s]-20, 0, coordinates_of_bars$x[s]+20, 0, col="red", lwd=2)
	segments(coordinates_of_bars$x[s]-20, coordinates_of_bars$y[s], coordinates_of_bars$x[s]+20, coordinates_of_bars$y[s], col="red", lwd=2)
}
average_spectra_with_bars[[p]] <- recordPlot()
#
}
return (list(pixel_by_pixel_classification=final_result_matrix, patient_classification_matrix=final_result_matrix_avg, pixel_by_pixel_classification_images=classification_msi, average_spectra_with_bars=average_spectra_with_bars))
}













classification <- classify_patients (filepath_test, svm_model, smoothing_strength_preprocessing="medium", preprocess_spectra_in_packages_of=100, mass_range=c(4000,15000))





############ EXPORT THE CSV
#write.csv (final_result_matrix, file=filename, row.names=FALSE)

beep (3)
