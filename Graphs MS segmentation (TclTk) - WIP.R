#################### GRAPHS AND CLASSIFICATION (MODULAR) 2017.08.22 ####################

install_and_load_required_packages(c("parallel", "tcltk", "XML", "MALDIquantForeign", "MALDIquant", "parallel", "caret", "pls", "igraph", "GA", "psych", "doParallel"), repository = "http://cran.mirror.garr.it/mirrors/CRAN/", update_packages = TRUE)


# Path where the R workspace file is located (GUI) (the one with the models) (FILE RData)
#tkmessageBox(title = "Input R workspace", message = "Select the R workspace file containing all the models", icon = "info")
#filepath_R <- tclvalue(tkgetOpenFile(filetypes = "{{R workspace files} {.RData}}"))
#if (!nchar(filepath_R)) {
#    tkmessageBox(message = "No R workspace file was selected!")
#}    else {
#    tkmessageBox(message = paste("The R workspace file selected is", filepath_R))
#}



####################################### INPUT DATA
# Path where the imzML file is located (GUI)
tkmessageBox(title = "Input imzML", message = "Select the MSI imzML dataset", icon = "info")
filepath_imzml <- tclvalue(tkgetOpenFile(filetypes = "{{imzML files} {.imzML .imzml}}"))
if (!nchar(filepath_imzml)) {
    tkmessageBox(message = "No imzML file was selected!")
}    else {
    tkmessageBox(message = paste("The imzML file selected is", filepath_imzml))
}

## Define the path where to save output files
filepath_output_select <- tkmessageBox(title = "Path", message = "Select the folder where to save all the outputs", icon = "info")
filepath_output <- tclvalue(tkchooseDirectory())
if (!nchar(filepath_output)) {
    tkmessageBox(message = "No folder selected!")
}    else {
    tkmessageBox(message = paste("All the output files will be saved in", filepath_output))
}
setwd(filepath_output)















##### LOAD THE R WORKSPACE WITH THE MODEL LIST
# Create a temporary environment
#temporary_environment <- new.env()
# Load the workspace
#load(filepath_R, envir = temporary_environment)
# Get the models (R objects) from the workspace
#model_list <- get("model_list", pos = temporary_environment)
# Get the list of models
#list_of_models <- names(model_list)
# Get the list of features for all the models
#feature_list <- get("feature_list", pos = temporary_environment)






#custom_feature_vector <- model_list[[1]]$features_model
#custom_feature_vector <- feature_list$model_features
#custom_feature_vector <- NULL









########## RUN THE GRAPH SEGMENTATION FUNCTION
# Reflectron
graph_segm <- graph_MSI_segmentation(filepath_imzml = filepath_imzml, custom_feature_vector = NULL, spectra_format = "imzml", preprocessing_parameters = list(mass_range = c(800, 3000), transformation_algorithm = NULL, smoothing_algorithm = NULL, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 200, normalization_algorithm = "TIC", normalization_mass_range = NULL, spectral_alignment_algorithm = "cubic", spectral_alignment_reference = "auto", preprocess_spectra_in_packages_of = 0), tolerance_ppm = 200, allow_parallelization = T, peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = T, envelope_peaklist = F, SNR = 3, tof_mode = "reflectron", peak_filtering_frequency_threshold_percent = 0, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", correlation_method_for_adjacency_matrix = "pearson", correlation_threshold_for_adjacency_matrix = 0.90, pvalue_threshold_for_adjacency_matrix = 0.01, number_of_high_degree_vertices_for_subgraph = 0, vertices_not_in_induced_subgraph = "independent", max_GA_generations = 5, iterations_with_no_change = 5, plot_figures = TRUE, plot_graphs = TRUE, plot_legends = c("legend"), seed = 12345)

# Linear
graph_segm <- graph_MSI_segmentation(filepath_imzml = filepath_imzml, custom_feature_vector = NULL, spectra_format = "imzml", preprocessing_parameters = list(mass_range = c(4000, 15000), transformation_algorithm = NULL, smoothing_algorithm = "SavitzkyGolay", smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 200, normalization_algorithm = "TIC", normalization_mass_range = NULL, spectral_alignment_algorithm = "cubic", spectral_alignment_reference = "auto", preprocess_spectra_in_packages_of = 0), tolerance_ppm = NULL, allow_parallelization = FALSE, peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = FALSE, envelope_peaklist = FALSE, SNR = 3, tof_mode = "linear", peak_filtering_frequency_threshold_percent = 0, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", correlation_method_for_adjacency_matrix = "spearman", correlation_threshold_for_adjacency_matrix = 0.90, pvalue_threshold_for_adjacency_matrix = 0.05, number_of_high_degree_vertices_for_subgraph = 0, vertices_not_in_induced_subgraph = "independent", max_GA_generations = 15, iterations_with_no_change = 5, plot_figures = TRUE, plot_graphs = TRUE, plot_legends = c("legend"), seed = 12345)


#print(graph_segm$peaklist_matrix)
print(colnames(graph_segm$peaklist_matrix))


### Export the graph image PNG
png(filename = "Final_graph.png", width = 1920, height = 1080, res = 150, pointsize = 12)
graph_segm$final_graph_plot
dev.off()
### Export the graph image EPS
postscript(file = "Final_graph.eps", width = 1920, height = 1080, title = "Graph segmentation MSI", pointsize = 12, fonts = c("serif", "Palatino"))
graph_segm$final_graph_plot
dev.off()



### Export the graph image PNG
png(filename = "Graph_segmentation_MSI.png", width = 1920, height = 1080, res = 300, pointsize = 5)
graph_segm$msi_segmentation
dev.off()
### Export the graph image EPS
postscript(file = "Graph_segmentation_MSI.eps", width = 1920, height = 1080, title = "Graph segmentation MSI", pointsize = 12)
graph_segm$msi_segmentation
dev.off()


### Export the graph image PNG
png(filename = "Initial_graph.png", width = 1920, height = 1080, res = 150, pointsize = 12)
graph_segm$graph_spectra_plot
dev.off()
### Export the graph image EPS
postscript(file = "Initial_graph.eps", width = 1920, height = 1080, title = "Initial graph", pointsize = 12, fonts=c("serif", "Palatino"))
graph_segm$graph_spectra_plot
dev.off()



### Save the workspace
#save.image(file = "Graph R workspace.RData")
