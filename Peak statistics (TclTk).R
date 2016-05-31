############################ PEAK STATISTICS 2016.05.31

############## INSTALL AND LOAD THE REQUIRED PACKAGES

install_and_load_required_packages(c("tcltk", "xlsx", "parallel"), repository="http://cran.mirror.garr.it/mirrors/CRAN/")



# Update packages and load the required packages
update.packages(repos="http://cran.mirror.garr.it/mirrors/CRAN/", ask=FALSE)







###################################### Initialize the variables (default values)
filepath_import <- NULL
tof_mode <- "linear"
filepath_export <- NULL
spectra_format <- "imzml"
peaks_filtering <- TRUE
low_intensity_peaks_removal <- FALSE
intensity_threshold_method <- "element-wise"
peak_picking_mode <- "all"
peak_picking_algorithm <- "SuperSmoother"
file_type_export <- "xlsx"
spectra <- NULL
peaks <- NULL
remove_outliers <- FALSE
exclude_spectra_without_peak <- FALSE
multicore_processing <- TRUE
transform_data <- FALSE
transform_data_algorithm <- NULL




################## Values of the variables (for displaying and dumping purposes)
tof_mode_value <- "linear"
filepath_import_value <- NULL
filepath_export_value <- NULL
spectra_format_value <- "imzML"
peaks_filtering_value <- "YES"
low_intensity_peaks_removal_value <- "NO"
peak_picking_mode_value <- "all"
peak_picking_algorithm_value <- "Super Smoother"
intensity_threshold_method_value <- "element-wise"
spectra_format_value <- "imzML"
remove_outliers_value <- "NO"
exclude_spectra_without_peak_value <- "NO"
multicore_processing_value <- "YES"
transform_data_value <- "    NO    "






##################################################### DEFINE WHAT THE BUTTONS DO

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
	file_type_export_value_label <- tklabel(window, text=file_type_export)
	tkgrid(file_type_export_value_label, row=10, column=4)
}

##### File name (export)
set_file_name <- function() {
	filename <- tclvalue(file_name)
	# Add the extension if it is not present in the filename
	if (file_type_export == "csv") {
		if (length(grep(".csv", filename, fixed=TRUE)) == 1) {
			filename <- filename
		}	else {filename <- paste (filename, ".csv", sep="")}
	}
	if (file_type_export == "xlsx") {
		if (length(grep(".xlsx", filename, fixed=TRUE)) == 1) {
			filename <- filename
		}	else {filename <- paste (filename, ".xlsx", sep="")}
	}
	if (file_type_export == "xls") {
		if (length(grep(".xls", filename, fixed=TRUE)) == 1) {
			filename <- filename
		}	else {filename <- paste (filename, ".xls", sep="")}
	}
	# Set the value for displaying purposes
	filename_value <- filename
	#### Exit the function and put the variable into the R workspace
	.GlobalEnv$filename <- filename
}

##### Samples
select_samples_function <-function() {
	########## Prompt if a folder has to be selected or a single file
	# Catch the value from the popping out menu
	spectra_input_type <- select.list(c("file","folder"), title="Folder or file?")
	if (spectra_input_type == "") {
		spectra_input_type <- "file"
	}
	if (spectra_input_type == "folder") {
		filepath_import_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be imported.\nIf there are imzML files in the folder, each one of them should contain spectra from the same class.\nIf there are subfolders, each one of them should contain spectra from the same class as imzML files.", icon = "info")
		filepath_import <- tclvalue(tkchooseDirectory())
		if (!nchar(filepath_import)) {
		    tkmessageBox(message = "No folder selected")
		}	else {
		    tkmessageBox(message = paste("The sample spectra will be read from:", filepath_import))
		}
	} else if (spectra_input_type == "file") {
		filepath_import_select <- tkmessageBox(title = "Samples", message = "Select the file for the spectra to be imported.", icon = "info")
		filepath_import <- tclvalue(tkgetOpenFile())
		if (!nchar(filepath_import)) {
		    tkmessageBox(message = "No file selected")
		}	else {
		    tkmessageBox(message = paste("The sample spectra will be read from:", filepath_import))
		}
	}
	# Set the value for displaying purposes
	filepath_import_value <- filepath_import
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$filepath_import <- filepath_import
}

##### Output
browse_output_function <- function() {
	output_folder <- tclvalue(tkchooseDirectory())
	if (!nchar(output_folder)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("Every file will be saved in", output_folder))
	}
	setwd(output_folder)
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$output_folder <- output_folder
}

##### Exit
end_session_function <- function () {
	q(save="no")
}

##### Import the spectra
import_spectra_function <- function() {
	# Load the required libraries
	install_and_load_required_packages(c("MALDIquantForeign", "MALDIquant"))
	# Rename the trim function
	trim_spectra <- get(x="trim", pos="package:MALDIquant")
	###### Get the values
	## Mass range
	mass_range <- tclvalue(mass_range)
	mass_range_value <- as.character(mass_range)
	mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
	# Preprocessing
	preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of)
	preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)
	preprocess_spectra_in_packages_of_value <- as.character(preprocess_spectra_in_packages_of)
	# Generate the list of spectra (library and test)
	if (spectra_format == "brukerflex" || spectra_format == "xmass") {
		### Load the spectra
		spectra <- importBrukerFlex(filepath_import)
	}
	if (spectra_format == "imzml" | spectra_format == "imzML") {
		### Fix the name, if the ibd file was selected instead of the imzML
		if (length(grep(".ibd",filepath_import)) != 0) {
			# Split the filepath
			filepath_splitted <- unlist(strsplit(filepath_import, ".ibd"))
			# Rebuild the imzML filename
			filepath_import <- paste(filepath_splitted, ".imzML", sep="")
		}
		### Load the spectra
		spectra <- importImzMl(filepath_import)
	}
	### Preprocessing
	spectra <- preprocess_spectra(spectra, tof_mode=tof_mode, smoothing_strength="medium", process_in_packages_of=preprocess_spectra_in_packages_of, multicore_processing=multicore_processing, align_spectra=TRUE, data_transformation=transform_data, transformation_algorithm=transform_data_algorithm, crop_spectra=TRUE, mass_range=mass_range)
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$spectra <- spectra
	.GlobalEnv$mass_range_value <- mass_range_value
	.GlobalEnv$preprocess_spectra_in_packages_of_value <- preprocess_spectra_in_packages_of_value
	### Messagebox
	tkmessageBox(title = "Import successful", message = "The spectra have been successfully imported and preprocessed", icon = "info")
}

##### Peak picking function
peak_picking_function <- function() {
	############ Do not run if the spectra have not been imported
	if (!is.null(spectra)) {
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
			peaks <- most_intense_signals(spectra, signals_to_take=signals_to_take, tof_mode=tof_mode, multicore_processing=multicore_processing)
		} else if (peak_picking_mode == "all"){
			peaks <- peak_picking(spectra, peak_picking_algorithm=peak_picking_algorithm, SNR=SNR, tof_mode=tof_mode, multicore_processing=multicore_processing)
		}
		# Exit the function and put the variable into the R workspace
		.GlobalEnv$peaks <- peaks
		.GlobalEnv$signals_to_take_value <- signals_to_take_value
		.GlobalEnv$SNR_value <- SNR_value
		### Messagebox
		tkmessageBox(title = "Peak picking successful", message = "The peak picking process has been successfully performed", icon = "info")
	} else if (is.null(spectra)) {
		### Messagebox
		tkmessageBox(title = "Spectra not imported", message = "The spectra have not been imported yet.\nImport them before performing the peak picking", icon = "warning")
	}
}

##### Output the average number of signals with the SD
signals_avg_and_sd_function <- function() {
	############ Do not run if the peaks have not been picked
	if (!is.null(peaks)) {
		# Generate the vector recording the number of signals
		number_of_signals_vector <- numeric()
		for (p in 1:length(peaks)) {
			number_of_signals_vector <- append(number_of_signals_vector, length(peaks[[p]]@mass))
		}
		# Compute the mean and the standard deviation
		mean_signal_number <- mean(number_of_signals_vector, na.rm=TRUE)
		sd_signal_number <- sd(number_of_signals_vector, na.rm=TRUE)
		cv_signal_number <- (sd_signal_number / mean_signal_number) *100
		### Messagebox
		# Message
		message_avg_sd <- paste("The mean number of signals in the spectral dataset is:", mean_signal_number, ",\n\nthe standard deviation is:", sd_signal_number, ",\n\nthe coefficient of variation is:", cv_signal_number, "%")
		tkmessageBox(title = "Mean and SD of the number of signals", message = message_avg_sd, icon = "info")
	} else if (is.null(peaks)) {
		### Messagebox
		tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the peak picking process has been performed", icon = "warning")
	}
}

##### Run the Peak Statistics
run_peak_statistics_function <- function() {
	############ Do not run if the spectra have not been imported or the peaks have not been picked
	if (!is.null(spectra) && !is.null(peaks)) {
		### List the directories in the filepath_import folder
		folder_list <- list.dirs(filepath_import, full.names=FALSE, recursive=FALSE)
		# If there are only imzML files <- one for class
		if (length(folder_list) == 0) {
			# Each imzML file is a class
			class_list <- read_spectra_files(filepath_import, spectra_format=spectra_format)
		} else if ((length(folder_list) == 1 && folder_list != "") || (length(folder_list) >= 1)) {
			class_list <- folder_list
		}
		## Peaks filtering threshold
		peaks_filtering_threshold_percent <- tclvalue(peaks_filtering_threshold_percent)
		peaks_filtering_threshold_percent <- as.numeric(peaks_filtering_threshold_percent)
		peaks_filtering_threshold_percent_value <- as.character(peaks_filtering_threshold_percent)
		## Low intensity threshold
		intensity_percentage_threshold <- tclvalue(intensity_percentage_threshold)
		intensity_percentage_threshold <- as.numeric(intensity_percentage_threshold)
		intensity_percentage_threshold_value <- as.character(intensity_percentage_threshold)
		### Run the peak statistics function
		peak_statistics_results <- peak_statistics(spectra, peaks, class_list=class_list, class_in_file_name=TRUE, tof_mode=tof_mode, spectra_format=spectra_format, peaks_filtering=peaks_filtering, peak_picking_algorithm=peak_picking_algorithm, frequency_threshold_percent=peaks_filtering_threshold_percent, remove_outliers=remove_outliers, low_intensity_peaks_removal=low_intensity_peaks_removal, intensity_threshold_percent=intensity_percentage_threshold, intensity_threshold_method=intensity_threshold_method, alignment_iterations=5)
		# Save the files (CSV)
		if (file_type_export == "csv") {
			filename <- set_file_name()
			    write.csv(peak_statistics_results, file=filename)
		} else if (file_type_export == "xlsx" || file_type_export == "xls") {
		# Save the files (Excel)
			filename <- set_file_name()
			peak_statistics_results <- as.data.frame(peak_statistics_results)
			# Generate unique row names
			unique_row_names <- make.names(rownames(peak_statistics_results), unique=TRUE)
			rownames(peak_statistics_results) <- unique_row_names
			# Export
			write.xlsx(x=peak_statistics_results, file=filename, sheetName="Peak statistics", row.names=TRUE)
		}
		### Messagebox
		tkmessageBox(title = "Done!", message = "The peak statistics file has been dumped!", icon = "info")
	} else if (is.null(spectra) || is.null(peaks)) {
		### Messagebox
		tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the spectra have been imported and the peak picking process has been performed", icon = "warning")
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
	peak_picking_mode_value_label <- tklabel(window, text=peak_picking_mode_value)
	tkgrid(peak_picking_mode_value_label, row=2, column=5)
	# Escape the function
	.GlobalEnv$peak_picking_mode <- peak_picking_mode
	.GlobalEnv$peak_picking_mode_value <- peak_picking_mode_value
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
	peak_picking_algorithm_value_label <- tklabel(window, text=peak_picking_algorithm_value)
	tkgrid(peak_picking_algorithm_value_label, row=6, column=6)
	# Escape the function
	.GlobalEnv$peak_picking_algorithm <- peak_picking_algorithm
	.GlobalEnv$peak_picking_algorithm_value <- peak_picking_algorithm_value
}

##### Peaks filtering
peaks_filtering_choice <- function() {
	# Catch the value from the menu
	peaks_filtering <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (peaks_filtering == "YES" || peaks_filtering == "") {
		peaks_filtering <- TRUE
	}
	if (peaks_filtering == "NO") {
		peaks_filtering <- FALSE
	}
	# Set the value of the displaying label
	if (peaks_filtering == TRUE) {
		peaks_filtering_value <- "YES"
	} else {
		peaks_filtering_value <- "NO"
	}
	peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value)
	tkgrid(peaks_filtering_value_label, row=4, column=3)
	# Escape the function
	.GlobalEnv$peaks_filtering <- peaks_filtering
	.GlobalEnv$peaks_filtering_value <- peaks_filtering_value
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
	multicore_processing_value_label <- tklabel(window, text=multicore_processing_value)
	tkgrid(multicore_processing_value_label, row=12, column=3)
	# Escape the function
	.GlobalEnv$multicore_processing <- multicore_processing
	.GlobalEnv$multicore_processing_value <- multicore_processing_value
}

##### Transform the data
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
	transform_data_value_label <- tklabel(window, text=transform_data_value)
	tkgrid(transform_data_value_label, row=12, column=5)
	# Escape the function
	.GlobalEnv$transform_data <- transform_data
	.GlobalEnv$transform_data_value <- transform_data_value
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
	low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value)
	tkgrid(low_intensity_peaks_removal_value_label, row=5, column=3)
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
	intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value)
	tkgrid(intensity_threshold_method_value_label, row=6, column=3)
	# Escape the function
	.GlobalEnv$intensity_threshold_method <- intensity_threshold_method
	.GlobalEnv$intensity_threshold_method_value <- intensity_threshold_method_value
}

##### Remove outliers from the peak intensity evaluation
remove_outliers_choice <- function() {
	# Catch the value from the menu
	remove_outliers <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (remove_outliers == "" || remove_outliers == "NO") {
		remove_outliers <- FALSE
	}
	if (remove_outliers == "YES") {
		remove_outliers <- TRUE
	}
	# Set the value of the displaying label
	if (remove_outliers == TRUE) {
		remove_outliers_value <- "YES"
	} else {
		remove_outliers_value <- "NO"
	}
	remove_outliers_value_label <- tklabel(window, text=remove_outliers_value)
	tkgrid(remove_outliers_value_label, row=7, column=3)
	# Escape the function
	.GlobalEnv$remove_outliers <- remove_outliers
	.GlobalEnv$remove_outliers_value <- remove_outliers_value
}

##### Exclude spectra without the peak
exclude_spectra_without_peak_function <- function() {
	# Catch the value from the menu
	exclude_spectra_without_peak <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (exclude_spectra_without_peak == "" || exclude_spectra_without_peak == "NO") {
		exclude_spectra_without_peak <- FALSE
	} else if (exclude_spectra_without_peak == "YES") {
		exclude_spectra_without_peak <- TRUE
	}
	# Set the value of the displaying label
	if (exclude_spectra_without_peak == TRUE) {
		exclude_spectra_without_peak_value <- "YES"
	} else {
		exclude_spectra_without_peak_value <- "NO"
	}
	exclude_spectra_without_peak_value_label <- tklabel(window, text=exclude_spectra_without_peak_value)
	tkgrid(exclude_spectra_without_peak_value_label, row=7, column=5)
	# Escape the function
	.GlobalEnv$exclude_spectra_without_peak <- exclude_spectra_without_peak
	.GlobalEnv$exclude_spectra_without_peak_value <- exclude_spectra_without_peak_value
}

##### TOF mode
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
	tof_mode_value_label <- tklabel(window, text=tof_mode_value)
	tkgrid(tof_mode_value_label, row=8, column=3)
	# Escape the function
	.GlobalEnv$tof_mode <- tof_mode
	.GlobalEnv$tof_mode_value <- tof_mode_value
}

##### File format
spectra_format_choice <- function() {
	# Catch the value from the menu
	spectra_format <- select.list(c("imzML","Xmass"), title="Choose")
	# Default
	if (spectra_format == "Xmass") {
		spectra_format <- "xmass"
		spectra_format_value <- "Xmass"
	}
	if (spectra_format == "" || spectra_format == "imzML") {
		spectra_format <- "imzml"
		spectra_format_value <- "imzML"
	}
	# Escape the function
	.GlobalEnv$spectra_format <- spectra_format
	.GlobalEnv$spectra_format_value <- spectra_format_value
	# Set the value of the displaying label
	spectra_format_value_label <- tklabel(window, text=spectra_format_value)
	tkgrid(spectra_format_value_label, row=9, column=3)
}










##################################################################### WINDOW GUI

########## List of variables, whose values are taken from the entries in the GUI
mass_range <- tclVar("")
SNR <- tclVar("")
peaks_filtering_threshold_percent <- tclVar("")
intensity_percentage_threshold <- tclVar("")
signals_to_take <- tclVar("")
file_name <- tclVar("")
preprocess_spectra_in_packages_of <- tclVar("")



######################## GUI

# The "area" where we will put our input lines
window <- tktoplevel()
tktitle(window) <- "Peak statistics"
#### Browse
# Library
select_samples_label <- tklabel(window, text="Select the file/folder containing the spectra")
select_samples_button <- tkbutton(window, text="Browse spectra...", command=select_samples_function)
# Output
select_output_label <- tklabel(window, text="Select the folder where to save all the outputs")
browse_output_button <- tkbutton(window, text="Browse output folder", command=browse_output_function)
#### Entries
# Mass range
mass_range_label <- tklabel(window, text="Mass range")
mass_range_entry <- tkentry(window, width=15, textvariable=mass_range)
tkinsert(mass_range_entry, "end", "3000, 15000")
# Preprocessing (in packages of)
preprocess_spectra_in_packages_of_label <- tklabel(window, text="Preprocess spectra\nin packages of")
preprocess_spectra_in_packages_of_entry <- tkentry(window, width=10, textvariable=preprocess_spectra_in_packages_of)
tkinsert(preprocess_spectra_in_packages_of_entry, "end", "200")
# Peak picking mode
peak_picking_mode_label <- tklabel(window, text="Peak picking mode")
peak_picking_mode_entry <- tkbutton(window, text="Choose peak picking\nmode", command=peak_picking_mode_choice)
# Peak picking method
peak_picking_algorithm_label <- tklabel(window, text="Peak picking method")
peak_picking_algorithm_entry <- tkbutton(window, text="Choose peak picking\nmethod", command=peak_picking_algorithm_choice)
# Signals to take
signals_to_take_label <- tklabel(window, text="Most intense signals to take\n(if 'most intense' is selected)")
signals_to_take_entry <- tkentry(window, width=10, textvariable=signals_to_take)
tkinsert(signals_to_take_entry, "end", "25")
# SNR
SNR_label <- tklabel(window, text="Signal-to-noise ratio")
SNR_entry <- tkentry(window, width=10, textvariable=SNR)
tkinsert(SNR_entry, "end", "5")
# Peaks filtering
peaks_filtering_label <- tklabel(window, text="Peaks filtering")
peaks_filtering_entry <- tkbutton(window, text="Choose peak filtering", command=peaks_filtering_choice)
# Peaks filtering threshold
peaks_filtering_threshold_percent_label <- tklabel(window, text="Peaks filtering threshold frequency percentage")
peaks_filtering_threshold_percent_entry <- tkentry(window, width=10, textvariable=peaks_filtering_threshold_percent)
tkinsert(peaks_filtering_threshold_percent_entry, "end", "25")
# Low intensity peaks removal
low_intensity_peaks_removal_label <- tklabel(window, text="Low intensity peaks removal")
low_intensity_peaks_removal_entry <- tkbutton(window, text="Choose low intensity\npeaks removal", command=low_intensity_peaks_removal_choice)
# Intensity percentage threshold
intensity_percentage_threshold_label <- tklabel(window, text="Intensity percentage threshold")
intensity_percentage_threshold_entry <- tkentry(window, width=10, textvariable=intensity_percentage_threshold)
tkinsert(intensity_percentage_threshold_entry, "end", "0.1")
# Intensiry percentage theshold method
intensity_threshold_method_label <- tklabel(window, text="Intensity threshold method")
intensity_threshold_method_entry <- tkbutton(window, text="Choose the method for\nthe intensity threshold", command=intensity_threshold_method_choice)
# Average replicates in database
remove_outliers_label <- tklabel(window, text="Discard outliers in intensity evaluation")
remove_outliers_entry <- tkbutton(window, text="Choose if outliers are to be discarded", command=remove_outliers_choice)
# Tof mode
tof_mode_label <- tklabel(window, text="Select the TOF mode")
tof_mode_entry <- tkbutton(window, text="Choose the TOF mode", command=tof_mode_choice)
# Exclude spectra without the peak
exclude_spectra_without_peak_button <- tkbutton(window, text="Exclude spectra without\nthe peak", command=exclude_spectra_without_peak_function)
# File format
spectra_format_label <- tklabel(window, text="Select the spectra format")
spectra_format_entry <- tkbutton(window, text="Choose the spectra format", command=spectra_format_choice)
# File type export
file_type_export_label <- tklabel(window, text="Select the format\nof the exported file")
file_type_export_entry <- tkbutton(window, text="Choose the file type", command=file_type_export_choice)
# End session
#end_session_label <- tklabel(window, text="Quit")
end_session_button <- tkbutton(window, text="QUIT", command=end_session_function)
# Import the spectra
import_spectra_button <- tkbutton(window, text="SPECTRA IMPORT AND\nPREPROCESSING", command=import_spectra_function)
# Peak picking
peak_picking_button <- tkbutton(window, text="PEAK PICKING", command=peak_picking_function)
# Run the Peak Statistics!!
run_peak_statistics_button <- tkbutton(window, text="COMPUTE THE PEAK STATISTICS", command=run_peak_statistics_function)
# Average number of signals and standard deviation
signals_avg_and_sd_button <- tkbutton(window, text="MEAN +/- SD of number of signals", command=signals_avg_and_sd_function)
# Multicore
multicore_processing_button <- tkbutton(window, text="ALLOW PARALLEL\nPROCESSING", command=multicore_processing_choice)
# Transform the data
transform_data_button <- tkbutton(window, text="TRANSFORM THE DATA", command=transform_data_choice)
# Set the file name
set_file_name_label <- tklabel(window, text="<-- Set the file name")
set_file_name_entry <- tkentry(window, width=30, textvariable=file_name)
tkinsert(set_file_name_entry, "end", "Peak statistics")

#### Displaying labels
file_type_export_value_label <- tklabel(window, text=file_type_export)
peak_picking_mode_value_label <- tklabel(window, text=peak_picking_mode_value)
peak_picking_algorithm_value_label <- tklabel(window, text=peak_picking_algorithm_value)
peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value)
low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value)
intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value)
remove_outliers_value_label <- tklabel(window, text=remove_outliers_value)
tof_mode_value_label <- tklabel(window, text=tof_mode_value)
spectra_format_value_label <- tklabel(window, text=spectra_format_value)
exclude_spectra_without_peak_value_label <- tklabel(window, text=exclude_spectra_without_peak_value)
multicore_processing_value_label <- tklabel(window, text=multicore_processing_value)
transform_data_value_label <- tklabel(window, text=transform_data_value)

#### Geometry manager
# Scrollbar
#window_scrollbar <- tkscrollbar(window, command=function(...)tkyview(window,...))
# tkgrid
tkgrid(select_samples_button, row=1, column=1)
tkgrid(browse_output_button, row=1, column=2)
tkgrid(set_file_name_entry, row=1, column=3)
tkgrid(set_file_name_label, row=1, column=4)
tkgrid(mass_range_label, row=2, column=1)
tkgrid(mass_range_entry, row=2, column=2)
tkgrid(peak_picking_mode_label, row=2, column=3)
tkgrid(peak_picking_mode_entry, row=2, column=4)
tkgrid(peak_picking_mode_value_label, row=2, column=5)
tkgrid(signals_to_take_label, row=3, column=1)
tkgrid(signals_to_take_entry, row=3, column=2)
tkgrid(SNR_label, row=3, column=3)
tkgrid(SNR_entry, row=3, column=4)
tkgrid(peaks_filtering_label, row=4, column=1)
tkgrid(peaks_filtering_entry, row=4, column=2)
tkgrid(peaks_filtering_value_label, row=4, column=3)
tkgrid(peaks_filtering_threshold_percent_label, row=4, column=4)
tkgrid(peaks_filtering_threshold_percent_entry, row=4, column=5)
tkgrid(low_intensity_peaks_removal_label, row=5, column=1)
tkgrid(low_intensity_peaks_removal_entry, row=5, column=2)
tkgrid(low_intensity_peaks_removal_value_label, row=5, column=3)
tkgrid(intensity_percentage_threshold_label, row=5, column=4)
tkgrid(intensity_percentage_threshold_entry, row=5, column=5)
tkgrid(intensity_threshold_method_label, row=6, column=1)
tkgrid(intensity_threshold_method_entry, row=6, column=2)
tkgrid(intensity_threshold_method_value_label, row=6, column=3)
tkgrid(peak_picking_algorithm_label, row=6, column=4)
tkgrid(peak_picking_algorithm_entry, row=6, column=5)
tkgrid(peak_picking_algorithm_value_label, row=6, column=6)
tkgrid(remove_outliers_label, row=7, column=1)
tkgrid(remove_outliers_entry, row=7, column=2)
tkgrid(remove_outliers_value_label, row=7, column=3)
tkgrid(exclude_spectra_without_peak_button, row=7, column=4)
tkgrid(exclude_spectra_without_peak_value_label, row=7, column=5)
tkgrid(tof_mode_label, row=8, column=1)
tkgrid(tof_mode_entry, row=8, column=2)
tkgrid(tof_mode_value_label, row=8, column=3)
tkgrid(spectra_format_label, row=9, column=1)
tkgrid(spectra_format_entry, row=9, column=2)
tkgrid(spectra_format_value_label, row=9, column=3)
tkgrid(preprocess_spectra_in_packages_of_label, row=9, column=4)
tkgrid(preprocess_spectra_in_packages_of_entry, row=9, column=5)
tkgrid(file_type_export_label, row=10, column=2)
tkgrid(file_type_export_entry, row=10, column=3)
tkgrid(file_type_export_value_label, row=10, column=4)
tkgrid(multicore_processing_button, row=12, column=2)
tkgrid(multicore_processing_value_label, row=12, column=3)
tkgrid(transform_data_button, row=12, column=4)
tkgrid(transform_data_value_label, row=12, column=5)
tkgrid(import_spectra_button, row=13, column=2)
tkgrid(peak_picking_button, row=13, column=3)
tkgrid(run_peak_statistics_button, row=13, column=4)
tkgrid(signals_avg_and_sd_button, row=13, column=5)
tkgrid(end_session_button, row=14, column=3)
