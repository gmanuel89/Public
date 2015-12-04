###INSTALL THE REQUIRED PACKAGES
install_and_load_required_packages(c("parallel", "MALDIquant", "MALDIquantForeign", "tcltk", "doMC"))


##### MULTICORE
# Detect the number of cores
cpu_thread_number <- detectCores(logical=TRUE)
cpu_core_number <- cpu_thread_number/2
# Use only the 3/4 of this power
#CPUcoreNumber <- floor (CPUcoreNumber *3/4)
# Register the foreach backend
registerDoMC(cores = cpu_core_number)





############################################ Define the parameters (GUI)
## Define the path from where to load the spectra
filepath_library_select <- tkmessageBox(title = "Library", message = "Select the folder for the spectra for the library.\nThe library should be structured like this: Database folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/1/1SLin/Spectrum_data", icon = "info")
filepath_library <- tclvalue(tkchooseDirectory())
if (!nchar(filepath_library)) {
    tkmessageBox(message = "No folder selected")
}	else {
    tkmessageBox(message = paste("The directory selected is", filepath_library))
}
##
filepath_test_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be tested.\nThe library should be structured like this: Database folder/Classes/Samples/Replicates/Spectra/Spectrum_coordinates/1/1SLin/Spectrum_data", icon = "info")
filepath_test <- tclvalue(tkchooseDirectory())
if (!nchar(filepath_test)) {
    tkmessageBox(message = "No folder selected")
}	else {
    tkmessageBox(message = paste("The directory selected is", filepath_test))
}
##
## Define the path where to save output files
output_folder_select <- tkmessageBox(title = "Output directory", message = "Select the folder where to save all the outputs", icon = "info")
output_folder <- tclvalue(tkchooseDirectory())
if (!nchar(output_folder)) {
    tkmessageBox(message = "No folder selected")
}	else {
    tkmessageBox(message = paste("Every file will be saved in", output_folder))
}
setwd(output_folder)






################## SCORE
score <- biotyper_like(filepath_library, filepath_test, mass_range=c(3000,15000), similarity_criteria=c("hca", "signal intensity"), signal_itensity_evaluation="intensity percentage", intensity_tolerance_percent=70, peak_picking_mode=c("all", "most intense"), SNR=5, number_of_peaks=20, low_intensity_peaks_removal=FALSE, intensity_threshold_percent=0.1, tof_mode="linear", file_format="brukerflex", average_replicates_in_database=FALSE, average_replicates_in_test=FALSE, score_only=TRUE, spectra_path_output=TRUE)



# Matrices
score_hca_matrix <- score$score_hca$result_matrix
score_intensity_matrix <- score$score_intensity

# Plots
score_hca_plots <- score$score_hca$plots







#testPerformance <- biotyperPerformanceWide (score)








######################## OUTPUTS
# Score
write.csv(score_intensity_matrix, file="Classification_results_Biotyper-like (intensity).csv")
write.csv(score_hca_matrix, file="Classification_results_Biotyper-like (hierarchical clustering).csv")

# HCA
png(file="Hierarchical clustering", width=1900, height=1200)
score_hca_plots
dev.off()
