################ BIOTYPER-LIKE PROGRAM 2015.12.14

install_and_load_required_packages("tcltk")







############################# DEFINE WHAT THE BUTTONS DO

# Library
select_library_function <- function() {
	#filepath_library_select <- tkmessageBox(title = "Library", message = "Select the folder for the spectra for the library.\nThe library should be structured like this: Database folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/1/1SLin/Spectrum_data", icon = "info")
	filepath_library <- tclvalue(tkchooseDirectory())
	if (!nchar(filepath_library)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("The directory selected is", filepath_library))
	}
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$filepath_library <- filepath_library
}
# Samples
select_samples_function <-function() {
	#filepath_test_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be tested.\nThe library should be structured like this: Sample folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/1/1SLin/Spectrum_data", icon = "info")
	filepath_test <- tclvalue(tkchooseDirectory())
	if (!nchar(filepath_test)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("The directory selected is", filepath_test))
	}
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$filepath_test <- filepath_test
}

# Output
browse_output_function <- function() {
	output_folder <- tclvalue(tkchooseDirectory())
	if (!nchar(output_folder)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("Every file will be saved in", output_folder))
	}
	setwd(output_folder)
}

# Close
quit_function <-function() {
	tkdestroy(window)
}

# Exit
end_session_function <- function () {
	q(save="no")
}

# Import the spectra
import_spectra_function <- function() {
	###### Get the values
	## Peak picking mode
	peak_picking_mode <- tclvalue(peak_picking_mode)
	## Signals to take in most intense peaks
	signals_to_take <- tclvalue(signals_to_take)
	signals_to_take <- as.integer(signals_to_take)
	## Mass range
	mass_range <- tclvalue(mass_range)
	mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
	# Preprocessing
	preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of)
	if (preprocess_spectra_in_packages_of == "Preprocess spectra in packages of (default 200)") {
		preprocess_spectra_in_packages_of <- 200
	} else {preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)}
	## Average replicates database
	average_replicates_in_database <- tclvalue(average_replicates_in_database)
	if (average_replicates_in_database == "Y" || average_replicates_in_database == "y" || average_replicates_in_database == "YES" || average_replicates_in_database == "yes") {
		average_replicates_in_database <- TRUE
	}
	if (average_replicates_in_database == "N" || average_replicates_in_database == "n" || average_replicates_in_database == "NO" || average_replicates_in_database == "no") {
		average_replicates_in_database <- FALSE
	}
	#tkconfigure(average_replicates_in_database_yes, variable=average_replicates_in_database, value=TRUE)
	#tkconfigure(average_replicates_in_database_no, variable=average_replicates_in_database, value=FALSE)
	## Average replicates samples
	average_replicates_in_test <- tclvalue(average_replicates_in_test)
	if (average_replicates_in_test == "Y" || average_replicates_in_test == "y" || average_replicates_in_test == "YES" || average_replicates_in_test == "yes") {
		average_replicates_in_test <- TRUE
	}
	if (average_replicates_in_test == "N" || average_replicates_in_test == "n" || average_replicates_in_test == "NO" || average_replicates_in_test == "no") {
		average_replicates_in_test <- FALSE
	}
	#tkconfigure(average_replicates_in_test_yes, variable=average_replicates_in_test, value=TRUE)
	#tkconfigure(average_replicates_in_test_no, variable=average_replicates_in_test, value=FALSE)
	## SNR
	SNR <- tclvalue(SNR)
	SNR <- as.numeric(SNR)
	## TOF mode
	tof_mode <- tclvalue(tof_mode)
	if (tof_mode == "L" || tof_mode == "l" || tof_mode == "linear" || tof_mode == "Linear") {
		tof_mode <- "linear"
	}
	if (tof_mode == "R" || tof_mode == "r" || tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "Reflectron" || tof_mode == "Reflector") {
		tof_mode <- "reflectron"
	}
	#tkconfigure(linear_mode, variable=file_format, value="linear")
	#tkconfigure(reflectron_mode, variable=file_format, value="reflectron")
	## File format
	file_format <- tclvalue(file_format)
	if (file_format == "imzml" || file_format == "imzML") {
		file_format <- "imzml"
	}
	if (file_format == "brukerflex" || file_format == "xmass" || file_format == "Xmass" || file_format == "XMass") {
		file_format <- "brukerflex"
	}
	#tkconfigure(imzml_format, variable=file_format, value="imzml")
	#tkconfigure(brukerflex_format, variable=file_format, value="brukerflex")
	if ("most intense" %in% peak_picking_mode) {
		########## DATABASE
		# Generate the library
		library_list <- library_creation(filepath_library, class_grouping=TRUE, spectra_preprocessing=list(smoothing_strength="medium", preprocess_in_packages_of=preprocess_spectra_in_packages_of), mass_range=mass_range, average_replicates=average_replicates_in_database, SNR=SNR, most_intense_peaks=TRUE, signals_to_take=signals_to_take, tof_mode=tof_mode, file_format=file_format, average_patients=FALSE)
		# Isolate the peaks and the spectra
		peaks_library <- library_list$peaks
		spectra_library <- library_list$spectra
		########## SAMPLES
		# Generate the library
		test_list <- library_creation(filepath_test, class_grouping=FALSE, spectra_preprocessing=list(smoothing_strength="medium", preprocess_in_packages_of=200), mass_range=mass_range, average_replicates=average_replicates_in_test, SNR=SNR, most_intense_peaks=TRUE, signals_to_take=signals_to_take, tof_mode=tof_mode, file_format=file_format, average_patients=FALSE)
		# Isolate the peaks and the spectra
		peaks_test <- test_list$peaks
		spectra_test <- test_list$spectra
	} else {
		########## DATABASE
		# Generate the library
		library_list <- library_creation(filepath_library, class_grouping=TRUE, spectra_preprocessing=list(smoothing_strength="medium", preprocess_in_packages_of=200), mass_range=mass_range, average_replicates=average_replicates_in_database, SNR=SNR, most_intense_peaks=FALSE, signals_to_take=25, tof_mode=tof_mode, file_format=file_format, average_patients=FALSE)
		# Isolate the peaks and the spectra
		peaks_library <- library_list$peaks
		spectra_library <- library_list$spectra
		########## SAMPLES
		# Generate the library
		test_list <- library_creation(filepath_test, class_grouping=FALSE, spectra_preprocessing=list(smoothing_strength="medium", preprocess_in_packages_of=200), mass_range=mass_range, average_replicates=average_replicates_in_test, SNR=SNR, most_intense_peaks=FALSE, signals_to_take=25, tof_mode=tof_mode, file_format=file_format, average_patients=FALSE)
		# Isolate the peaks and the spectra
		peaks_test <- test_list$peaks
		spectra_test <- test_list$spectra
	}
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$library_list <- library_list
	.GlobalEnv$peaks_library <- peaks_library
	.GlobalEnv$spectra_library <- spectra_library
	.GlobalEnv$test_list <- test_list
	.GlobalEnv$peaks_test <- peaks_test
	.GlobalEnv$spectra_test <- spectra_test
	### Messagebox
	tkmessageBox(title = "Import successful", message = "The spectra have been successfully imported", icon = "info")
}

# Run the biotyper_like
run_biotyper_like_function <- function() {
	#### Get the values
	## Similarity criteria
	similarity_criteria <- tclvalue(similarity_criteria)
	## Intensity correction coefficient
	intensity_correction_coefficient <- tclvalue(intensity_correction_coefficient)
	intensity_correction_coefficient <- as.numeric(intensity_correction_coefficient)
	## Signal intensity evaluation
	signal_intensity_evaluation <- tclvalue(signal_intensity_evaluation)
	## Peak picking mode
	peak_picking_mode <- tclvalue(peak_picking_mode)
	## Peaks filtering
	peaks_filtering <- tclvalue(peaks_filtering)
	if (peaks_filtering == "Y" || peaks_filtering == "y" || peaks_filtering == "YES" || peaks_filtering == "yes") {
		peaks_filtering <- TRUE
	}
	if (peaks_filtering == "N" || peaks_filtering == "n" || peaks_filtering == "NO" || peaks_filtering == "no") {
		peaks_filtering == "N" <- FALSE
	}
	## Peaks filtering threshold
	peaks_filtering_threshold_percent <- tclvalue(peaks_filtering_threshold_percent)
	peaks_filtering_threshold_percent <- as.numeric(peaks_filtering_threshold_percent)
	## Low intensity peaks removal
	low_intensity_peaks_removal <- tclvalue(low_intensity_peaks_removal)
	if (low_intensity_peaks_removal == "Y" || low_intensity_peaks_removal == "y" || low_intensity_peaks_removal == "YES" || low_intensity_peaks_removal == "yes") {
		low_intensity_peaks_removal <- TRUE
	}
	if (low_intensity_peaks_removal == "N" || low_intensity_peaks_removal == "n" || low_intensity_peaks_removal == "NO" || low_intensity_peaks_removal == "no") {
		low_intensity_peaks_removal <- FALSE
	}
	#tkconfigure(low_intensity_peaks_removal_yes, variable=low_intensity_peaks_removal, value=TRUE)
	#tkconfigure(low_intensity_peaks_removal_no, variable=low_intensity_peaks_removal, value=FALSE)
	## Low intensity threshold
	intensity_percentage_threshold <- tclvalue(intensity_percentage_threshold)
	intensity_percentage_threshold <- as.numeric(intensity_percentage_threshold)
	## TOF mode
	tof_mode <- tclvalue(tof_mode)
	if (tof_mode == "L" || tof_mode == "l" || tof_mode == "linear" || tof_mode == "Linear") {
		tof_mode <- "linear"
	}
	if (tof_mode == "R" || tof_mode == "r" || tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "Reflectron" || tof_mode == "Reflector") {
		tof_mode <- "refletron"
	}
	#tkconfigure(linear_mode, variable=file_format, value="linear")
	#tkconfigure(reflectron_mode, variable=file_format, value="reflectron")
	## File format
	file_format <- tclvalue(file_format)
	if (file_format == "imzml" || file_format == "imzML") {
		file_format <- "imzml"
	}
	if (file_format == "brukerflex" || file_format == "xmass" || file_format == "Xmass" || file_format == "XMass") {
		file_format <- "brukerflex"
	}
	#tkconfigure(imzml_format, variable=file_format, value="imzml")
	#tkconfigure(brukerflex_format, variable=file_format, value="brukerflex")
	## Score only
	score_only <- tclvalue(score_only)
	if (score_only == "Y" || score_only == "y" || score_only == "YES" || score_only == "yes") {
		score_only <- TRUE
	}
	if (score_only == "N" || score_only == "n" || score_only == "NO" || score_only == "no") {
		score_only <- FALSE
	}
	#tkconfigure(score_only_yes, variable=score_only, value=TRUE)
	#tkconfigure(score_only_no, variable=score_only, value=FALSE)
	## Spectra path output
	spectra_path_output <- tclvalue(spectra_path_output)
	if (spectra_path_output == "Y" || spectra_path_output == "y" || spectra_path_output == "YES" || spectra_path_output == "yes") {
		spectra_path_output <- TRUE
	}
	if (spectra_path_output == "N" || spectra_path_output == "n" || spectra_path_output == "NO" || spectra_path_output == "no") {
		spectra_path_output <- FALSE
	}
	#tkconfigure(spectra_path_output_yes, variable=spectra_path_output, value=TRUE)
	#tkconfigure(spectra_path_output_no, variable=spectra_path_output, value=FALSE)
	#### Run the function
	score <- biotyper_like(library_list, test_list, similarity_criteria=similarity_criteria, intensity_correction_coefficient=intensity_correction_coefficient, signal_intensity_evaluation=signal_intensity_evaluation, intensity_tolerance_percent=70, peaks_filtering=peaks_filtering, peaks_filtering_threshold_percent=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, tof_mode=tof_mode, file_format=file_format, score_only=score_only, spectra_path_output=spectra_path_output)
	# Matrices
	score_hca_matrix <- score$score_hca$result_matrix$output
	score_intensity_matrix <- score$score_intensity
	score_correlation_matrix <- score$score_correlation_matrix$output
	#### Exit the function and put the variable into the R workspace
	.GlobalEnv$score <- score
	# Matrices
	.GlobalEnv$score_hca_matrix <- score$score_hca$result_matrix$output
	.GlobalEnv$score_intensity_matrix <- score$score_intensity
	.GlobalEnv$score_correlation_matrix <- score$score_correlation_matrix$output
	# Save the files
	filename <- set_file_name()
	if (!is.null(score_intensity_matrix)) {
	    write.csv(score_intensity_matrix, file=paste("1_", filename, sep=""))
	}
	if (!is.null(score_hca_matrix)) {
	    write.csv(score_hca_matrix, file=paste("2_", filename, sep=""))
	}
	if (!is.null(score_correlation_matrix)) {
	    write.csv(score_correlation_matrix, file=paste("3_", filename, sep=""))
	}
	### Messagebox
	tkmessageBox(title = "Done!", message = "The file(s) have been dumped", icon = "info")
}

set_file_name <- function() {
	filename <- tclvalue(file_name)
	# Add the extension if it is not present in the filename
	if (length(grep(".csv", filename, fixed=TRUE)) == 1) {
		filename <- filename
	}	else {filename <- paste (filename, ".csv", sep="")}
	#### Exit the function and put the variable into the R workspace
	.GlobalEnv$filename <- filename
	return (filename)
}





######################################### WINDOW GUI
# List of variables, whose values are taken from the entries in the GUI
mass_range <- tclVar("")
SNR <- tclVar("")
average_replicates_in_database <- tclVar("")
average_replicates_in_test <- tclVar("")
tof_mode <- tclVar("")
file_format <- tclVar("")
similarity_criteria <- tclVar("")
intensity_correction_coefficient <- tclVar("")
signal_intensity_evaluation <- tclVar("")
peak_picking_mode <- tclVar("")
peaks_filtering <- tclVar("")
peaks_filtering_threshold_percent <- tclVar("")
low_intensity_peaks_removal <- tclVar("")
intensity_percentage_threshold <- tclVar("")
score_only <- tclVar("")
spectra_path_output <- tclVar("")
signals_to_take <- tclVar("")
file_name <- tclVar("")
preprocess_spectra_in_packages_of <- tclVar("")



####### GUI
# The "area" where we will put our input lines
window <- tktoplevel()
tktitle(window) <- "Biotyper-like"
#### Browse
# Library
select_library_label <- tklabel(window, text="Select the folder for the spectra for the library.\nThe library should be structured like this:\nDatabase folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/1/1SLin/Spectrum_data")
select_library_button <- tkbutton(window, text="Browse database", command=select_library_function)
# Samples
select_samples_label <- tklabel(window, text="Select the folder for the spectra for the library.\nThe library should be structured like this:\nDatabase folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/1/1SLin/Spectrum_data")
select_samples_button <- tkbutton(window, text="Browse samples", command=select_samples_function)
# Output
select_output_label <- tklabel(window, text="Select the folder where to save all the outputs")
browse_output_button <- tkbutton(window, text="Browse", command=browse_output_function)
#### Entries
# Mass range
mass_range_label <- tklabel(window, text="Mass range")
mass_range_entry <- tkentry(window, width=30, textvariable=mass_range)
tkinsert(mass_range_entry, "end", "3000, 15000")
# Preprocessing (in packages of)
preprocess_spectra_in_packages_of_entry <- tkentry(window, width=30, textvariable=preprocess_spectra_in_packages_of)
tkinsert(preprocess_spectra_in_packages_of_entry, "end", "Preprocess spectra in packages of (default 200)")
# Similarity criteria
similarity_criteria_label <- tklabel(window, text="Similarity criteria")
similarity_criteria_entry <- tkentry(window, width=30, textvariable=similarity_criteria)
tkinsert(similarity_criteria_entry, "end", "correlation")
# Intensity correction coefficient
intensity_correction_coefficient_label <- tklabel(window, text="Intensity correction coefficient")
intensity_correction_coefficient_entry <- tkentry(window, width=30, textvariable=intensity_correction_coefficient)
tkinsert(intensity_correction_coefficient_entry, "end", "1")
# Signal intensity evaluation
signal_intensity_evaluation_label <- tklabel(window, text="Signal intensity evaluation")
signal_intensity_evaluation_entry <- tkentry(window, width=30, textvariable=signal_intensity_evaluation)
tkinsert(signal_intensity_evaluation_entry, "end", "intensity percentage")
# Peak picking mode
peak_picking_mode_label <- tklabel(window, text="Peak picking mode")
peak_picking_mode_entry <- tkentry(window, width=30, textvariable=peak_picking_mode)
tkinsert(peak_picking_mode_entry, "end", "all")
# Signals to take
signals_to_take_label <- tklabel(window, text="Most intense signals to take")
signals_to_take_entry <- tkentry(window, width=30, textvariable=signals_to_take)
tkinsert(signals_to_take_entry, "end", "25")
# SNR
SNR_label <- tklabel(window, text="Signal-to-noise ratio")
SNR_entry <- tkentry(window, width=30, textvariable=SNR)
tkinsert(SNR_entry, "end", "5")
# Peaks filtering
peaks_filtering_label <- tklabel(window, text="Peaks filtering")
peaks_filtering_entry <- tkentry(window, width=30, textvariable=peaks_filtering)
tkinsert(peaks_filtering_entry, "end", "yes")
#peaks_filtering_yes <- tkradiobutton(window)
#peaks_filtering_no <- tkradiobutton(window)
#tkconfigure(peaks_filtering_yes, variable=peaks_filtering, value=TRUE)
#tkconfigure(peaks_filtering_no, variable=peaks_filtering, value=FALSE)
# Peaks filtering threshold
peaks_filtering_threshold_percent_label <- tklabel(window, text="Peaks filtering threshold frequency")
peaks_filtering_threshold_percent_entry <- tkentry(window, width=30, textvariable=peaks_filtering_threshold_percent)
tkinsert(peaks_filtering_threshold_percent_entry, "end", "25")
# Low intensity peaks removal
low_intensity_peaks_removal_label <- tklabel(window, text="Low intensity peaks removal")
low_intensity_peaks_removal_entry <- tkentry(window, width=30, textvariable=low_intensity_peaks_removal)
tkinsert(low_intensity_peaks_removal_entry, "end", "no")
#low_intensity_peaks_removal_yes <- tkradiobutton(window)
#low_intensity_peaks_removal_no <- tkradiobutton(window)
# Intensity percentage threshold
intensity_percentage_threshold_label <- tklabel(window, text="Intensity percentage threshold")
intensity_percentage_threshold_entry <- tkentry(window, width=30, textvariable=intensity_percentage_threshold)
tkinsert(intensity_percentage_threshold_entry, "end", "0.1")
# Average replicates in database
average_replicates_in_database_label <- tklabel(window, text="Average replicates in the database")
average_replicates_in_database_entry <- tkentry(window, width=30, textvariable=average_replicates_in_database)
tkinsert(average_replicates_in_database_entry, "end", "no")
#average_replicates_in_database_yes <- tkradiobutton(window)
#average_replicates_in_database_no <- tkradiobutton(window)
# Average replicates in samples
average_replicates_in_test_label <- tklabel(window, text="Average replicates in the samples")
average_replicates_in_test_entry <- tkentry(window, width=30, textvariable=average_replicates_in_test)
tkinsert(average_replicates_in_test_entry, "end", "no")
#average_replicates_in_test_yes <- tkradiobutton(window)
#average_replicates_in_test_no <- tkradiobutton(window)
# Score only
score_only_label <- tklabel(window, text="Score only")
score_only_entry <- tkentry(window, width=30, textvariable=score_only)
tkinsert(score_only_entry, "end", "no")
#score_only_yes <- tkradiobutton(window)
#score_only_no <- tkradiobutton(window)
# Spectra path output
spectra_path_output_label <- tklabel(window, text="Spectra path output")
spectra_path_output_entry <- tkentry(window, width=30, textvariable=spectra_path_output)
tkinsert(spectra_path_output_entry, "end", "yes")
#spectra_path_output_yes <- tkradiobutton(window)
#spectra_path_output_no <- tkradiobutton(window)
# Tof mode
tof_mode_label <- tklabel(window, text="Select the TOF mode")
tof_mode_entry <- tkentry(window, width=30, textvariable=tof_mode)
tkinsert(tof_mode_entry, "end", "linear")
#linear_mode <- tkradiobutton(window)
#reflectron_mode <- tkradiobutton(window)
# File format
file_format_label <- tklabel(window, text="Select the file format")
file_format_entry <- tkentry(window, width=30, textvariable=file_format)
tkinsert(file_format_entry, "end", "xmass")
#imzml_format <- tkradiobutton(window)
#brukerflex_format <- tkradiobutton(window)
#### Close
exit_label <- tklabel(window, text="Exit")
quit_button <- tkbutton(window, text="Quit", command=quit_function)
# End session
#end_session_label <- tklabel(window, text="Quit")
end_session_button <- tkbutton(window, text="End R session", command=end_session_function)
# Import the spectra
import_spectra_button <- tkbutton(window, text="Import spectra, preprocessing and peak picking", command=import_spectra_function)
# Run the Biotyper-like!
run_biotyper_like_button <- tkbutton(window, text="Run the Biotyper-like!", command=run_biotyper_like_function)
# Set the file name
set_file_name_label <- tklabel(window, text="Set the file name")
set_file_name_entry <- tkentry(window, width=30, textvariable=file_name)
tkinsert(set_file_name_entry, "end", "Biotyper-like score output")

#### Geometry manager
# Scrollbar
#window_scrollbar <- tkscrollbar(window, command=function(...)tkyview(window,...))
# tkgrid
tkgrid(set_file_name_label, row=1, column=4)
tkgrid(set_file_name_entry, row=2, column=4)
tkgrid(select_library_label, row=1, column=2)
tkgrid(select_library_button, row=1, column=3)
tkgrid(select_samples_label, row=2, column=2)
tkgrid(select_samples_button, row=2, column=3)
tkgrid(select_output_label, row=3, column=2)
tkgrid(browse_output_button, row=3, column=3)
tkgrid(mass_range_label, row=4, column=1)
tkgrid(mass_range_entry, row=4, column=2)
tkgrid(similarity_criteria_label, row=4, column=3)
tkgrid(similarity_criteria_entry, row=4, column=4)
tkgrid(intensity_correction_coefficient_label, row=5, column=1)
tkgrid(intensity_correction_coefficient_entry, row=5, column=2)
tkgrid(signal_intensity_evaluation_label, row=5, column=3)
tkgrid(signal_intensity_evaluation_entry, row=5, column=4)
tkgrid(peak_picking_mode_label, row=6, column=1)
tkgrid(peak_picking_mode_entry, row=6, column=2)
tkgrid(signals_to_take_label, row=6, column=3)
tkgrid(signals_to_take_entry, row=6, column=4)
tkgrid(SNR_label, row=7, column=2)
tkgrid(SNR_entry, row=7, column=3)
tkgrid(peaks_filtering_label, row=8, column=1)
tkgrid(peaks_filtering_entry, row=8, column=2)
tkgrid(peaks_filtering_threshold_percent_label, row=8, column=3)
tkgrid(peaks_filtering_threshold_percent_entry, row=8, column=4)
tkgrid(low_intensity_peaks_removal_label, row=9, column=1)
tkgrid(low_intensity_peaks_removal_entry, row=9, column=2)
tkgrid(intensity_percentage_threshold_label, row=9, column=3)
tkgrid(intensity_percentage_threshold_entry, row=9, column=4)
tkgrid(average_replicates_in_database_label, row=10, column=1)
tkgrid(average_replicates_in_database_entry, row=10, column=2)
tkgrid(average_replicates_in_test_label, row=10, column=3)
tkgrid(average_replicates_in_test_entry, row=10, column=4)
tkgrid(score_only_label, row=11, column=1)
tkgrid(score_only_entry, row=11, column=2)
tkgrid(spectra_path_output_label, row=11, column=3)
tkgrid(spectra_path_output_entry, row=11, column=4)
tkgrid(tof_mode_label, row=12, column=1)
tkgrid(tof_mode_entry, row=12, column=2)
tkgrid(file_format_label, row=12, column=3)
tkgrid(file_format_entry, row=12, column=4)
tkgrid(import_spectra_button, row=13, column=1)
tkgrid(preprocess_spectra_in_packages_of_entry, row=13, column=2)
tkgrid(run_biotyper_like_button, row=13, column=4)
tkgrid(exit_label, row=14, column=1)
tkgrid(quit_button, row=14, column=2)
tkgrid(end_session_button, row=14, column=4)
#window_scrollbar
