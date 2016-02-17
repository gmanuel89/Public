################ SPECTRAL TYPER PROGRAM 2016.02.17

# Update packages and load the required packages
update.packages(repos="http://cran.mirror.garr.it/mirrors/CRAN/", ask=FALSE)

install_and_load_required_packages(c("tcltk", "xlsx", "ggplot2"), repository="http://cran.mirror.garr.it/mirrors/CRAN/")





###################################### Initialise the variables (default values)
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
similarity_criteria <- "correlation"
signal_intensity_evaluation <- "intensity percentage"
file_type_export <- "xlsx"
spectra_database <- NULL
spectra_test <- NULL
peaks_database <- NULL
peaks_test <- NULL





################## Values of the variables (for displaying and dumping purposes)
average_replicates_in_database_value <- "NO"
average_replicates_in_test_value <- "NO"
peaks_filtering_value <- "NO"
low_intensity_peaks_removal_value <- "NO"
score_only_value <- "NO"
spectra_path_output_value <- "YES"
similarity_criteria_value <- "correlation"
peak_picking_mode_value <- "all"
tof_mode_value <- "linear"
intensity_threshold_method_value <- "element-wise"
spectra_format_value <- "Xmass"




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
	# Set the value of the displaying label
	file_type_export_value_label <- tklabel(window, text=file_type_export)
	tkgrid(file_type_export_value_label, row=13, column=4)
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
	return (filename)
}

##### Library
select_database_function <- function() {
	filepath_database_select <- tkmessageBox(title = "Library", message = "Select the folder for the spectra for the database.\nThe database should be structured like this:\nDatabase folder/Database entry samples/Treatments/Spectra_replicates/Spectrum_coordinates/1/1SLin/Spectrum_data\n\nor\n\nDatabase folder/Classes - Entries/Sample imzML files", icon = "info")
	filepath_database <- tclvalue(tkchooseDirectory())
	if (!nchar(filepath_database)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("The directory selected for the database is", filepath_database))
	}
	# Set the value for displaying purposes
	filepath_database_value <- filepath_database
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$filepath_database <- filepath_database
}

##### Samples
select_samples_function <-function() {
	filepath_test_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be tested.\nThe files should be organised like this:\nSample folder/Samples/Treatments/Spectra_replicates/Spectrum_coordinates/1/1SLin/Spectrum_data\n\nor\n\nSample folder/Sample imzML files", icon = "info")
	filepath_test <- tclvalue(tkchooseDirectory())
	if (!nchar(filepath_test)) {
	    tkmessageBox(message = "No folder selected")
	}	else {
	    tkmessageBox(message = paste("The sample spectra will be read from:", filepath_test))
	}
	# Set the value for displaying purposes
	filepath_test_value <- filepath_test
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
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$output_folder <- output_folder
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
		spectra_database <- importBrukerFlex(filepath_database)
		spectra_test <- importBrukerFlex(filepath_test)
	}
	if (spectra_format == "imzml" | spectra_format == "imzML") {
		### Load the spectra
		spectra_database <- importImzMl(filepath_database)
		spectra_test <- importImzMl(filepath_test)
	}
	### Truncation
	spectra_database <- trim_spectra(spectra_database, range = mass_range)
	spectra_test <- trim_spectra(spectra_test, range = mass_range)
	### Average the replicates
	if (average_replicates_in_database == TRUE) {
		spectra_database <- average_replicates_by_folder(spectra_database, filepath_database, spectra_format=spectra_format)
	}
	if (average_replicates_in_test == TRUE) {
		spectra_test <- average_replicates_by_folder(spectra_test, filepath_test, spectra_format=spectra_format)
	}
	### Folder lists
	database_folder_list <- dir(filepath_database, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
	test_folder_list <- dir(filepath_test, ignore.case=TRUE, full.names=FALSE, recursive=FALSE, include.dirs=TRUE)
	### Spectra grouping (class for database)
	spectra_database <- group_spectra_class(spectra_database, class_list=database_folder_list, spectra_format=spectra_format, class_in_file_name=TRUE)
	### Preprocessing
	spectra_database <- preprocess_spectra(spectra_database, tof_mode=tof_mode, smoothing_strength="medium", process_in_packages_of=preprocess_spectra_in_packages_of)
	spectra_test <- preprocess_spectra(spectra_test, tof_mode=tof_mode, smoothing_strength="medium", process_in_packages_of=preprocess_spectra_in_packages_of)
	# Exit the function and put the variable into the R workspace
	.GlobalEnv$spectra_database <- spectra_database
	.GlobalEnv$spectra_test <- spectra_test
	.GlobalEnv$mass_range_value <- mass_range_value
	.GlobalEnv$preprocess_spectra_in_packages_of_value <- preprocess_spectra_in_packages_of_value
	.GlobalEnv$database_folder_list <- database_folder_list
	.GlobalEnv$test_folder_list <- test_folder_list
	### Messagebox
	tkmessageBox(title = "Import successful", message = "The spectra have been successfully imported and preprocessed", icon = "info")
}

##### Peak picking function
peak_picking_function <- function() {
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
		# Define the halfWindowSize for peak picking
		if (tof_mode == "linear") {
			halfWindowSize <- 20
		} else if (tof_mode == "reflector") {
			halfWindowSize <- 5
		}
		if (peak_picking_mode == "most intense") {
			peaks_database <- most_intense_signals(spectra_database, signals_to_take=signals_to_take, tof_mode=tof_mode)
			peaks_test <- most_intense_signals(spectra_test, signals_to_take=signals_to_take)
		} else if (peak_picking_mode == "all"){
			peaks_database <- detectPeaks(spectra_database, method="MAD", SNR=SNR, halfWindowSize=halfWindowSize)
			peaks_test <- detectPeaks(spectra_test, method="MAD", SNR=SNR, halfWindowSize=halfWindowSize)
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
		############### CORRELATION
		if (similarity_criteria == "correlation") {
			score_correlation_matrix <- spectral_typer_score_correlation_matrix(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library=database_folder_list, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, intensity_correction_coefficient=intensity_correction_coefficient, spectra_format=spectra_format, spectra_path_output=spectra_path_output, score_only=score_only)
		} else if (similarity_criteria == "hca") {
			score_hca <- spectral_typer_score_hierarchical_distance(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=database_folder_list, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, spectra_path_output=spectra_path_output, score_only=score_only, spectra_format=spectra_format, normalise_distances=TRUE, normalisation_method="sum")
		} else if (similarity_criteria == "signal intensity") {
			score_intensity_matrix <- spectral_typer_score_signal_intensity(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library=database_folder_list, comparison=signal_intensity_evaluation, peaks_filtering=peaks_filtering, peaks_filtering_percentage_threshold=peaks_filtering_threshold_percent, low_intensity_peaks_removal=low_intensity_peaks_removal, low_intensity_percentage_threshold=intensity_percentage_threshold, low_intensity_threshold_method=intensity_threshold_method, tolerance_ppm=tolerance_ppm, intensity_tolerance_percent_threshold=intensity_tolerance_percent, spectra_format=spectra_format, spectra_path_output=spectra_path_output, score_only=score_only, number_of_st_dev=1)
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
				ggsave(plot=score_hca$hca_dendrogram, device="png", filename="hca.png", width=6.33, height=4, dpi=300)
				#savePlot(filename="hca.png", type="png")
				#dev.print(X11, file="hca.png", width=1900, height=1280)
				#dev.off()
			    write.csv(score_hca_matrix_results, file=paste("hca_", filename, sep=""))
			}
			if (!is.null(score_correlation_matrix)) {
			    write.csv(score_correlation_matrix_results, file=paste("corr_", filename, sep=""))
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
			    #WriteXLS(x=score_intensity_matrix_results, ExcelFileName=paste("int_", filename, sep=""), SheetNames="Scores - Intensity", row.names=TRUE, AdjWidth=TRUE, verbose=TRUE)
				write.xlsx(x=score_intensity_matrix_results, file=paste("int_", filename, sep=""), sheetName="Scores - Intensity", row.names=TRUE)

			}
			if (!is.null(score_hca)) {
				# Dump the hca plot
				#png(filename="hca.png", width=1900, height=1280)
				#score$score_hca$plots
				ggsave(plot=score_hca$hca_dendrogram, device="png", filename="hca.png", width=6.33, height=4, dpi=300)
				#savePlot(filename="hca.png", type="png")
				#dev.print(X11, file="hca.png", width=1900, height=1280)
				#dev.off()
				# Convert it to a data frame
				score_hca_matrix_results <- as.data.frame(score_hca_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_hca_matrix_results), unique=TRUE)
				rownames(score_hca_matrix_results) <- unique_row_names
				# Export
			    #WriteXLS(x=score_hca_matrix_results, ExcelFileName=paste("hca_", filename, sep=""), SheetNames="Scores - HCA", row.names=TRUE, AdjWidth=TRUE, verbose=TRUE)
				write.xlsx(x=score_hca_matrix_results, file=paste("hca_", filename, sep=""), sheetName="Scores - HCA", row.names=TRUE)
			}
			if (!is.null(score_correlation_matrix)) {
				# Convert it to a data frame
				score_correlation_matrix_results <- as.data.frame(score_correlation_matrix_results)
				# Generate unique row names
				unique_row_names <- make.names(rownames(score_correlation_matrix_results), unique=TRUE)
				rownames(score_correlation_matrix_results) <- unique_row_names
				# Export
			    #WriteXLS(x=score_correlation_matrix_results, ExcelFileName=paste("corr_", filename, sep=""), SheetNames="Scores - Correlation", row.names=TRUE, AdjWidth=TRUE, verbose=TRUE)
				write.xlsx(x=score_correlation_matrix_results, file=paste("corr_", filename, sep=""), sheetName="Scores - Correlation", row.names=TRUE)
			}
		}
		### Messagebox
		tkmessageBox(title = "Done!", message = "The file(s) have been dumped\n\nLegend:\nF: Fit\nRF: Retrofit\nCorr: intensity Pearson's correlation coefficient\nIntMtch: signal intensity matching\nsl: slope of the regression curve\nns: number of signals\n\n\nFit = number of sample-database matching signals / number of signals in the sample\nRetrofit = number of database-sample matching signals / number of signals in the database entry", icon = "info")
	} else if (is.null(spectra_database) || is.null(spectra_test) || is.null(peaks_database) || is.null(peaks_test)) {
		### Messagebox
		tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the spectra have been imported and the peak picking process has been performed", icon = "warning")
	}
}

##### Similarity criteria
similarity_criteria_choice <- function() {
	# Catch the value from the menu
	similarity_criteria <- select.list(c("correlation","hca","signal intensity"), title="Choose")
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
	similarity_criteria_value_label <- tklabel(window, text=similarity_criteria_value)
	tkgrid(similarity_criteria_value_label, row=3, column=6)
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
	signal_intensity_evaluation_value_label <- tklabel(window, text=signal_intensity_evaluation)
	tkgrid(signal_intensity_evaluation_value_label, row=4, column=6)
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
	tkgrid(peak_picking_mode_value_label, row=5, column=3)
	# Escape the function
	.GlobalEnv$peak_picking_mode <- peak_picking_mode
	.GlobalEnv$peak_picking_mode_value <- peak_picking_mode_value
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
	peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value)
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
	low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value)
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
	intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value)
	tkgrid(intensity_threshold_method_value_label, row=9, column=3)
	# Escape the function
	.GlobalEnv$intensity_threshold_method <- intensity_threshold_method
	.GlobalEnv$intensity_threshold_method_value <- intensity_threshold_method_value
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
	average_replicates_in_database_value_label <- tklabel(window, text=average_replicates_in_database_value)
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
	average_replicates_in_test_value_label <- tklabel(window, text=average_replicates_in_test_value)
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
	score_only_value_label <- tklabel(window, text=score_only_value)
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
	spectra_path_output_value_label <- tklabel(window, text=spectra_path_output_value)
	tkgrid(spectra_path_output_value_label, row=11, column=6)
	# Escape the function
	.GlobalEnv$spectra_path_output <- spectra_path_output
	.GlobalEnv$spectra_path_output_value <- spectra_path_output_value
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
	tkgrid(tof_mode_value_label, row=12, column=3)
	# Escape the function
	.GlobalEnv$tof_mode <- tof_mode
	.GlobalEnv$tof_mode_value <- tof_mode_value
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
	spectra_format_value_label <- tklabel(window, text=spectra_format_value)
	tkgrid(spectra_format_value_label, row=12, column=6)
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

### Initial Messagebox
#tkmessageBox(title = "Before starting", message = "The library should be structured like this: Main folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/Spectrum_data", icon = "info")

# The "area" where we will put our input lines
window <- tktoplevel()
tktitle(window) <- "Spectral Typer"
#### Browse
# Library
select_database_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data")
select_database_button <- tkbutton(window, text="Browse database folder", command=select_database_function)
# Samples
select_samples_label <- tklabel(window, text="The library should be structured like this:\nMain folder/Classes/Samples/Replicates/Spectra/\nSpectrum_coordinates/Spectrum_data")
select_samples_button <- tkbutton(window, text="Browse samples folder", command=select_samples_function)
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
# Similarity criteria
similarity_criteria_label <- tklabel(window, text="Similarity criteria")
similarity_criteria_entry <- tkbutton(window, text="Choose similarity\ncriteria", command=similarity_criteria_choice)
# Intensity correction coefficient
intensity_correction_coefficient_label <- tklabel(window, text="Intensity correction coefficient\n(0: discard the intensities,\n1: unweighted correlation)")
intensity_correction_coefficient_entry <- tkentry(window, width=10, textvariable=intensity_correction_coefficient)
tkinsert(intensity_correction_coefficient_entry, "end", "1")
# Intensty tolerance percent
intensity_tolerance_percent_label <- tklabel(window, text="Intensity tolerance percent\n(if 'signal intensity' is selected)")
intensity_tolerance_percent_entry <- tkentry(window, width=10, textvariable=intensity_tolerance_percent)
tkinsert(intensity_tolerance_percent_entry, "end", "80")
# Signal intensity evaluation
signal_intensity_evaluation_label <- tklabel(window, text="Signal intensity evaluation")
signal_intensity_evaluation_entry <- tkbutton(window, text="Choose signal intensity\nevaluation method", command=signal_intensity_evaluation_choice)
# Peak picking mode
peak_picking_mode_label <- tklabel(window, text="Peak picking mode")
peak_picking_mode_entry <- tkbutton(window, text="Choose pick picking\nmode", command=peak_picking_mode_choice)
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
average_replicates_in_database_label <- tklabel(window, text="Average replicates in the database")
average_replicates_in_database_entry <- tkbutton(window, text="Choose average replicates\nin the database", command=average_replicates_in_database_choice)
# Average replicates in samples
average_replicates_in_test_label <- tklabel(window, text="Average replicates in the samples")
average_replicates_in_test_entry <- tkbutton(window, text="Choose average replicates\nin the samples", command=average_replicates_in_test_choice)
# Score only
score_only_label <- tklabel(window, text="Score only\n('NO', all the score\ncomponents are displayed)")
score_only_entry <- tkbutton(window, text="Choose", command=score_only_choice)
# Spectra path output
spectra_path_output_label <- tklabel(window, text="Spectra path in the output")
spectra_path_output_entry <- tkbutton(window, text="Choose to display\nthe spectra path", command=spectra_path_output_choice)
# Tof mode
tof_mode_label <- tklabel(window, text="Select the TOF mode")
tof_mode_entry <- tkbutton(window, text="Choose the TOF mode", command=tof_mode_choice)
# File format
spectra_format_label <- tklabel(window, text="Select the spectra format")
spectra_format_entry <- tkbutton(window, text="Choose the spectra format", command=spectra_format_choice)
# File type export
file_type_export_label <- tklabel(window, text="Select the format\nof the exported file")
file_type_export_entry <- tkbutton(window, text="Choose the file type", command=file_type_export_choice)
#### Close
exit_label <- tklabel(window, text="Exit")
quit_button <- tkbutton(window, text="Quit", command=quit_function)
# End session
#end_session_label <- tklabel(window, text="Quit")
end_session_button <- tkbutton(window, text="QUIT", command=end_session_function)
# Import the spectra
import_spectra_button <- tkbutton(window, text="SPECTRA IMPORT AND\nPREPROCESSING", command=import_spectra_function)
# Peak picking
peak_picking_button <- tkbutton(window, text="PEAK PICKING", command=peak_picking_function)
# Run the Spectral typer!
run_spectral_typer_button <- tkbutton(window, text="RUN THE SPECTRAL TYPER", command=run_spectral_typer_function)
# Set the file name
set_file_name_label <- tklabel(window, text="<-- Set the file name")
set_file_name_entry <- tkentry(window, width=30, textvariable=file_name)
tkinsert(set_file_name_entry, "end", "Score")

#### Displaying labels
file_type_export_value_label <- tklabel(window, text=file_type_export)
similarity_criteria_value_label <- tklabel(window, text=similarity_criteria_value)
signal_intensity_evaluation_value_label <- tklabel(window, text=signal_intensity_evaluation)
peak_picking_mode_value_label <- tklabel(window, text=peak_picking_mode_value)
peaks_filtering_value_label <- tklabel(window, text=peaks_filtering_value)
low_intensity_peaks_removal_value_label <- tklabel(window, text=low_intensity_peaks_removal_value)
intensity_threshold_method_value_label <- tklabel(window, text=intensity_threshold_method_value)
average_replicates_in_database_value_label <- tklabel(window, text=average_replicates_in_database_value)
average_replicates_in_test_value_label <- tklabel(window, text=average_replicates_in_test_value)
score_only_value_label <- tklabel(window, text=score_only_value)
spectra_path_output_value_label <- tklabel(window, text=spectra_path_output_value)
tof_mode_value_label <- tklabel(window, text=tof_mode_value)
spectra_format_value_label <- tklabel(window, text=spectra_format_value)

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
tkgrid(mass_range_label, row=3, column=1)
tkgrid(mass_range_entry, row=3, column=2)
tkgrid(similarity_criteria_label, row=3, column=4)
tkgrid(similarity_criteria_entry, row=3, column=5)
tkgrid(similarity_criteria_value_label, row=3, column=6)
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
tkgrid(score_only_label, row=11, column=1)
tkgrid(score_only_entry, row=11, column=2)
tkgrid(score_only_value_label, row=11, column=3)
tkgrid(spectra_path_output_label, row=11, column=4)
tkgrid(spectra_path_output_entry, row=11, column=5)
tkgrid(spectra_path_output_value_label, row=11, column=6)
tkgrid(tof_mode_label, row=12, column=1)
tkgrid(tof_mode_entry, row=12, column=2)
tkgrid(tof_mode_value_label, row=12, column=3)
tkgrid(spectra_format_label, row=12, column=4)
tkgrid(spectra_format_entry, row=12, column=5)
tkgrid(spectra_format_value_label, row=12, column=6)
tkgrid(intensity_tolerance_percent_label, row=6, column=4)
tkgrid(intensity_tolerance_percent_entry, row=6, column=5)
tkgrid(preprocess_spectra_in_packages_of_label, row=9, column=4)
tkgrid(preprocess_spectra_in_packages_of_entry, row=9, column=5)
tkgrid(file_type_export_label, row=13, column=2)
tkgrid(file_type_export_entry, row=13, column=3)
tkgrid(file_type_export_value_label, row=13, column=4)
tkgrid(import_spectra_button, row=14, column=2)
tkgrid(peak_picking_button, row=14, column=3)
tkgrid(run_spectral_typer_button, row=14, column=4)
#tkgrid(exit_label, row=15, column=1)
#tkgrid(quit_button, row=15, column=2)
tkgrid(end_session_button, row=14, column=5)
#window_scrollbar
