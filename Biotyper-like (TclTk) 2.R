################ BIOTYPER-LIKE PROGRAM 2016.01.22

# Update packages and load the required packages
update.packages(repos="http://cran.mirror.garr.it/mirrors/CRAN/")
install_and_load_required_packages(c("tcltk", "WriteXLS"), repository="http://cran.mirror.garr.it/mirrors/CRAN/")





###################################### Initialise the variables (default values)
average_replicates_in_database <- FALSE
average_replicates_in_test <- FALSE
peaks_filtering <- TRUE
low_intensity_peaks_removal <- FALSE
tof_mode <- "linear"
file_format <- "brukerflex"
score_only <- FALSE
spectra_path_output <- TRUE
peak_picking_mode <- "all"
similarity_criteria <- "correlation"
signal_intensity_evaluation <- "intensity percentage"
file_type_export <- "xlsx"




##################################################### DEFINE WHAT THE BUTTONS DO

##### File type (export)
file_type_export_choice <- function() {
	# Catch the value from the menu
	file_type_export <- select.list(c("csv","xlsx","xls"), title="Choose")
	# Default
	if (file_type_export == "") {
		file_type_export <- "xlsx"
	}
	# Escape the function
	.GlobalEnv$file_type_export <- file_type_export
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
	#### Exit the function and put the variable into the R workspace
	.GlobalEnv$filename <- filename
	return (filename)
}

##### Library
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

##### Samples
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

##### Output
browse_output_function <- function() {
	output_folder <- tclvalue(tkchooseDirectory())
	if (!nchar(output_folder)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("Every file will be saved in", output_folder))
	}
	setwd(output_folder)
}

##### Close
quit_function <-function() {
	tkdestroy(window)
}

##### Exit
end_session_function <- function () {
	q(save="no")
}

##### Import the spectra
import_spectra_function <- function() {
	###### Get the values
	## Signals to take in most intense peaks
	signals_to_take <- tclvalue(signals_to_take)
	signals_to_take <- as.integer(signals_to_take)
	## Mass range
	mass_range <- tclvalue(mass_range)
	mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
	# Preprocessing
	preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of)
	preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)
	## SNR
	SNR <- tclvalue(SNR)
	SNR <- as.numeric(SNR)
	# Peak picking mode
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

##### Run the biotyper_like
run_biotyper_like_function <- function() {
	#### Get the values
	## Intensity correction coefficient
	intensity_correction_coefficient <- tclvalue(intensity_correction_coefficient)
	intensity_correction_coefficient <- as.numeric(intensity_correction_coefficient)
	## Intensity tolerance percent
	intensity_tolerance_percent <- tclvalue(intensity_tolerance_percent)
	intensity_tolerance_percent <- as.numeric(intensity_tolerance_percent)
	## Peaks filtering threshold
	peaks_filtering_threshold_percent <- tclvalue(peaks_filtering_threshold_percent)
	peaks_filtering_threshold_percent <- as.numeric(peaks_filtering_threshold_percent)
	## Low intensity threshold
	intensity_percentage_threshold <- tclvalue(intensity_percentage_threshold)
	intensity_percentage_threshold <- as.numeric(intensity_percentage_threshold)
	#### Run the function
	score <- biotyper_like(library_list, test_list, similarity_criteria=similarity_criteria, intensity_correction_coefficient=intensity_correction_coefficient, signal_intensity_evaluation=signal_intensity_evaluation, intensity_tolerance_percent=intensity_tolerance_percent, peaks_filtering=peaks_filtering, peaks_filtering_threshold_percent=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, tof_mode=tof_mode, file_format=file_format, score_only=score_only, spectra_path_output=spectra_path_output)
	# Matrices
	score_hca_matrix <- score$score_hca$result_matrix
	score_intensity_matrix <- score$score_intensity$output
	score_correlation_matrix <- score$score_correlation_matrix$output
	#### Exit the function and put the variable into the R workspace
	.GlobalEnv$score <- score
	# Matrices
	.GlobalEnv$score_hca_matrix <- score$score_hca$result_matrix
	.GlobalEnv$score_intensity_matrix <- score$score_intensity$output
	.GlobalEnv$score_correlation_matrix <- score$score_correlation_matrix$output
	# Save the files (CSV)
	if (file_type_export == "csv") {
		filename <- set_file_name()
		if (!is.null(score_intensity_matrix)) {
		    write.csv(score_intensity_matrix, file=paste("int_", filename, sep=""))
		}
		if (!is.null(score_hca_matrix)) {
		    write.csv(score_hca_matrix, file=paste("hca_", filename, sep=""))
		}
		if (!is.null(score_correlation_matrix)) {
		    write.csv(score_correlation_matrix, file=paste("corr_", filename, sep=""))
		}
	}
	if (file_type_export == "xlsx" || file_type_export == "xls") {
		filename <- set_file_name()
		# Check if PERL is installed by using the built in function, for WriteXLS package
		perl_installed <- testPerl(verbose=FALSE)
		if (perl_installed == FALSE) {
			tkmessageBox(title = "Perl", message = "The software PERL has not been found on the computer. Please install the Perl software from\nhttps://www.perl.org/get.html\nand restart the R session", icon = "warning")
		}
		if (!is.null(score_intensity_matrix) && perl_installed == TRUE) {
			# Convert it to a data frame
			score_intensity_matrix <- as.data.frame(score_intensity_matrix)
			# Generate unique row names
			unique_row_names <- make.names(rownames(score_intensity_matrix), unique=TRUE)
			rownames(score_intensity_matrix) <- unique_row_names
			# Export
		    WriteXLS(x=score_intensity_matrix, ExcelFileName=paste("int_", filename, sep=""), SheetNames="Scores - Intensity", row.names=TRUE, AdjWidth=TRUE)
		}
		if (!is.null(score_hca_matrix) && perl_installed == TRUE) {
			# Convert it to a data frame
			score_hca_matrix <- as.data.frame(score_hca_matrix)
			# Generate unique row names
			unique_row_names <- make.names(rownames(score_hca_matrix), unique=TRUE)
			rownames(score_hca_matrix) <- unique_row_names
			# Export
		    WriteXLS(x=score_hca_matrix, ExcelFileName=paste("hca_", filename, sep=""), SheetNames="Scores - HCA", row.names=TRUE, AdjWidth=TRUE)
		}
		if (!is.null(score_correlation_matrix) && perl_installed == TRUE) {
			# Convert it to a data frame
			score_correlation_matrix <- as.data.frame(score_correlation_matrix)
			# Generate unique row names
			unique_row_names <- make.names(rownames(score_correlation_matrix), unique=TRUE)
			rownames(score_correlation_matrix) <- unique_row_names
			# Export
		    WriteXLS(x=score_correlation_matrix, ExcelFileName=paste("corr_", filename, sep=""), SheetNames="Scores - Correlation", row.names=TRUE, AdjWidth=TRUE)
		}
	}
	### Messagebox
	tkmessageBox(title = "Done!", message = "The file(s) have been dumped", icon = "info")
}

##### Similarity criteria
similarity_criteria_choice <- function() {
	# Catch the value from the menu
	similarity_criteria <- select.list(c("correlation","hca","signal intensity"), title="Choose")
	# Default
	if (similarity_criteria == "") {
		similarity_criteria <- "correlation"
	}
	# Escape the function
	.GlobalEnv$similarity_criteria <- similarity_criteria
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
}

##### Peak picking mode
peak_picking_mode_choice <- function() {
	# Catch the value from the menu
	peak_picking_mode <- select.list(c("all","most intense"), title="Choose")
	# Default
	if (peak_picking_mode == "") {
		peak_picking_mode <- "all"
	}
	# Escape the function
	.GlobalEnv$peak_picking_mode <- peak_picking_mode
}

##### Peaks filtering
peaks_filtering_choice <- function() {
	# Catch the value from the menu
	peaks_filtering <- select.list(c("YES","NO"), title="Choose")
	# Default
	if (peaks_filtering == "" || peaks_filtering == "YES") {
		peaks_filtering <- TRUE
	}
	if (peaks_filtering == "NO") {
		peaks_filtering <- FALSE
	}
	# Escape the function
	.GlobalEnv$peaks_filtering <- peaks_filtering
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
	# Escape the function
	.GlobalEnv$low_intensity_peaks_removal <- low_intensity_peaks_removal
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
	# Escape the function
	.GlobalEnv$average_replicates_in_database <- average_replicates_in_database
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
	# Escape the function
	.GlobalEnv$average_replicates_in_test <- average_replicates_in_test
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
	# Escape the function
	.GlobalEnv$score_only <- score_only
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
	# Escape the function
	.GlobalEnv$spectra_path_output <- spectra_path_output
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
	# Escape the function
	.GlobalEnv$tof_mode <- tof_mode
}

##### File format
file_format_choice <- function() {
	# Catch the value from the menu
	file_format <- select.list(c("imzML","Xmass"), title="Choose")
	# Default
	if (file_format == "" || file_format == "Xmass") {
		file_format <- "brukerflex"
	}
	if (file_format == "imzML") {
		file_format <- "imzml"
	}
	# Escape the function
	.GlobalEnv$file_format <- file_format
}










##################################################################### WINDOW GUI

########## List of variables, whose values are taken from the entries in the GUI
mass_range <- tclVar("")
SNR <- tclVar("")
intensity_correction_coefficient <- tclVar("")
peaks_filtering_threshold_percent <- tclVar("")
intensity_percentage_threshold <- tclVar("")
signals_to_take <- tclVar("")
file_name <- tclVar("")
preprocess_spectra_in_packages_of <- tclVar("")
intensity_tolerance_percent <- tclVar("")



######################## GUI
# The "area" where we will put our input lines
window <- tktoplevel()
tktitle(window) <- "Biotyper-like"
#### Browse
# Library
select_library_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data")
select_library_button <- tkbutton(window, text="Browse database folder", command=select_library_function)
# Samples
select_samples_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data")
select_samples_button <- tkbutton(window, text="Browse samples folder", command=select_samples_function)
# Output
select_output_label <- tklabel(window, text="Select the folder where to save all the outputs")
browse_output_button <- tkbutton(window, text="Browse", command=browse_output_function)
#### Entries
# Mass range
mass_range_label <- tklabel(window, text="Mass range")
mass_range_entry <- tkentry(window, width=30, textvariable=mass_range)
tkinsert(mass_range_entry, "end", "3000, 15000")
# Preprocessing (in packages of)
preprocess_spectra_in_packages_of_label <- tklabel(window, text="Preprocess spectra in packages of")
preprocess_spectra_in_packages_of_entry <- tkentry(window, width=30, textvariable=preprocess_spectra_in_packages_of)
tkinsert(preprocess_spectra_in_packages_of_entry, "end", "200")
# Similarity criteria
similarity_criteria_label <- tklabel(window, text="Similarity criteria\n(default: 'correlation')")
similarity_criteria_entry <- tkbutton(window, text="Choose similarity\ncriteria", command=similarity_criteria_choice)
# Intensity correction coefficient
intensity_correction_coefficient_label <- tklabel(window, text="Intensity correction coefficient")
intensity_correction_coefficient_entry <- tkentry(window, width=30, textvariable=intensity_correction_coefficient)
tkinsert(intensity_correction_coefficient_entry, "end", "1")
# Intensty tolerance percent
intensity_tolerance_percent_label <- tklabel(window, text="Intensity tolerance percent\n(if 'signal intensity' is selected)")
intensity_tolerance_percent_entry <- tkentry(window, width=30, textvariable=intensity_tolerance_percent)
tkinsert(intensity_tolerance_percent_entry, "end", "80")
# Signal intensity evaluation
signal_intensity_evaluation_label <- tklabel(window, text="Signal intensity evaluation\n(default: 'intensity percentage')")
signal_intensity_evaluation_entry <- tkbutton(window, text="Choose signal intensity\nevaluation method", command=signal_intensity_evaluation_choice)
# Peak picking mode
peak_picking_mode_label <- tklabel(window, text="Peak picking mode\n(default: 'all')")
peak_picking_mode_entry <- tkbutton(window, text="Choose pick picking\nmode", command=peak_picking_mode_choice)
# Signals to take
signals_to_take_label <- tklabel(window, text="Most intense signals to take\n(if 'most intense' is selected)")
signals_to_take_entry <- tkentry(window, width=30, textvariable=signals_to_take)
tkinsert(signals_to_take_entry, "end", "25")
# SNR
SNR_label <- tklabel(window, text="Signal-to-noise ratio")
SNR_entry <- tkentry(window, width=30, textvariable=SNR)
tkinsert(SNR_entry, "end", "5")
# Peaks filtering
peaks_filtering_label <- tklabel(window, text="Peaks filtering\n(default: 'YES')")
peaks_filtering_entry <- tkbutton(window, text="Choose peak filtering", command=peaks_filtering_choice)
# Peaks filtering threshold
peaks_filtering_threshold_percent_label <- tklabel(window, text="Peaks filtering threshold frequency percentage")
peaks_filtering_threshold_percent_entry <- tkentry(window, width=30, textvariable=peaks_filtering_threshold_percent)
tkinsert(peaks_filtering_threshold_percent_entry, "end", "25")
# Low intensity peaks removal
low_intensity_peaks_removal_label <- tklabel(window, text="Low intensity peaks removal\n(default: 'NO')")
low_intensity_peaks_removal_entry <- tkbutton(window, text="Choose low intensity\npeaks removal", command=low_intensity_peaks_removal_choice)
# Intensity percentage threshold
intensity_percentage_threshold_label <- tklabel(window, text="Intensity percentage threshold")
intensity_percentage_threshold_entry <- tkentry(window, width=30, textvariable=intensity_percentage_threshold)
tkinsert(intensity_percentage_threshold_entry, "end", "0.1")
# Average replicates in database
average_replicates_in_database_label <- tklabel(window, text="Average replicates in the database\n(default: 'NO')")
average_replicates_in_database_entry <- tkbutton(window, text="Choose average replicates\nin the database", command=average_replicates_in_database_choice)
# Average replicates in samples
average_replicates_in_test_label <- tklabel(window, text="Average replicates in the samples\n(default: 'NO')")
average_replicates_in_test_entry <- tkbutton(window, text="Choose average replicates\nin the samples", command=average_replicates_in_test_choice)
# Score only
score_only_label <- tklabel(window, text="Score only\n(default: 'NO', all the score components are displayed)")
score_only_entry <- tkbutton(window, text="Choose", command=score_only_choice)
# Spectra path output
spectra_path_output_label <- tklabel(window, text="Spectra path in the output\n(default: 'YES')")
spectra_path_output_entry <- tkbutton(window, text="Choose to display\nthe spectra path", command=spectra_path_output_choice)
# Tof mode
tof_mode_label <- tklabel(window, text="Select the TOF mode\n(default: 'Linear')")
tof_mode_entry <- tkbutton(window, text="Choose the TOF mode", command=tof_mode_choice)
# File format
file_format_label <- tklabel(window, text="Select the file format\n(default: 'Xmass')")
file_format_entry <- tkbutton(window, text="Choose the file format", command=file_format_choice)
# File type export
file_type_export_label <- tklabel(window, text="Select the format of the exported file\n(default: 'xlsx')")
file_type_export_entry <- tkbutton(window, text="Choose the file type", command=file_type_export_choice)
#### Close
exit_label <- tklabel(window, text="Exit")
quit_button <- tkbutton(window, text="Quit", command=quit_function)
# End session
#end_session_label <- tklabel(window, text="Quit")
end_session_button <- tkbutton(window, text="End R session", command=end_session_function)
# Import the spectra
import_spectra_button <- tkbutton(window, text="Import spectra, preprocessing\nand peak picking", command=import_spectra_function)
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
tkgrid(intensity_tolerance_percent_label, row=13, column=1)
tkgrid(intensity_tolerance_percent_entry, row=13, column=2)
tkgrid(preprocess_spectra_in_packages_of_label, row=13, column=3)
tkgrid(preprocess_spectra_in_packages_of_entry, row=13, column=4)
tkgrid(file_type_export_label, row=14, column=2)
tkgrid(file_type_export_entry, row=14, column=3)
tkgrid(import_spectra_button, row=14, column=1)
tkgrid(run_biotyper_like_button, row=14, column=4)
tkgrid(exit_label, row=15, column=1)
tkgrid(quit_button, row=15, column=2)
tkgrid(end_session_button, row=15, column=4)
#window_scrollbar
