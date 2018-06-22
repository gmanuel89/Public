# MS PEAKLIST EXPORT

## import_spectra

This function imports the spectra into R from files

#### Inputs

* filepath --> absolute path (folder) where the spectra are located

* spectra_format --> the format in which the spectra are stored. It can be one among "imzML", "fid", "txt", "csv" or "msd"

* mass_range --> the mass range to which the spectra are cut to (the first number is the lower boundary, while the second number is the upper boundary); if NULL, the spectra are brought to the same scale by using highest overlapping mass range

* allow_parallelization --> if the parallel computation should be exploited

* spectral_names --> how to set the list names: "name" uses the spectrum name (sample name) as the name for the items in the list, while "number" uses sequential integer numbers as list names

* replace_sample_name_field --> if the file path in the spectrum should be replaced by the sample name

* remove_empty_spectra --> if the empty spectra should be automatically discarded in the import phase

#### Returns

* a list of 'MALDIquant' spectra, in which each item corresponds to a spectrum; it returns NULL if no spectra were found in the specified filepath.

#### Calls

* functions from the 'MALDIquant' and 'MALDIquantForeign' packages

## preprocess_spectra

This function runs the preprocessing steps onto the spectra

#### Inputs

* spectra --> a list of 'MALDIquant' spectra

* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"

* preprocessing_parameters --> a named list of preprocessing parameters to be used:
  * mass_range --> the mass range to which the spectra are cut to (the first number is the lower boundary, while the second number is the upper boundary); if  NULL, the spectra are brought to the same scale by using highest overlapping mass range
  * transformation_algorithm --> the algorithm to be used for data transformation; it can be "sqrt", "log", "log2", "log10", "sin", "cos", or "exp"; if NULL, no  smoothing is performed
  * smoothing_algorithm --> the algorithm to be used for smoothing; it can be "SavitzkyGolay" or "MovingAverage"; if NULL, no smoothing is performed
  * smoothing_strength --> the width of the smoothing filter to use; it can be "medium", "strong" or "stronger"
  * baseline_subtraction_algorithm --> the algorithm to be used for baseline subtraction; it can be "SNIP", "TopHat", "ConvexHull", or "median"; if NULL, no  baseline subtraction is performed
  * baseline_subtraction_algorithm_parameter --> the parameter to be used for defining the filter width for baseline subtraction
  * normalization_algorithm --> the algorithm to be used for normalization; it can be "TIC", "RMS", "PQN", or "median"; if NULL, no normalization is performed
  * normalization_mass_range --> the mass range to use for normalization: the first number is the lower boundary, while the second number is the upper boundary;  if NULL, the entire spectrum is used as the window for normalization
  * preprocess_spectra_in_packages_of --> the number of spectra to be processed simultaneously, in packets, if NULL or less or equal than 0, all the spectra are  preprocessed at once
  * spectral_alignment_algorithm --> the algorithm to be used for peak-based spectral alignment; it can be "cubic", "quadratic" or "linear", if NULL, no spectral  alignment is performed
  * spectral_alignment_reference --> the peak list to be used as reference for peak-based spectral alignment; it can be "average spectrum", "skyline spectrum" or "auto": the first two options allow to use the peaks of the average spectrum or the skyline spectrum as the reference for alignment, while the third option automatically estimates the peaks to be used for alignment as the most frequent peaks in the spectral list

* allow_parallelization --> if the parallel computation should be exploited

* tolerance_ppm --> the tolerance (in ppm) for the peak-based spectral alignment; if NULL, the tolerance is set to 1000 ppm for "linear" TOF mode and to 100 ppm for "reflectron" TOF mode

#### Returns

* a list of 'MALDIquant' spectra, in which each item corresponds to a preprocessed spectrum

#### Calls

* functions from the 'MALDIquant' package

## align_spectra

This function runs the spectral alignment

#### Inputs

* spectra --> a list of 'MALDIquant' spectra

* spectral_alignment_algorithm --> the algorithm to be used for peak-based spectral alignment; it can be "cubic", "quadratic" or "linear", if NULL, no spectral alignment is performed

* spectral_alignment_reference --> the peak list to be used as reference for peak-based spectral alignment; it can be "average spectrum", "skyline spectrum" or "auto": the first two options allow to use the peaks of the average spectrum or the skyline spectrum as the reference for alignment, while the third option automatically estimates the peaks to be used for alignment as the most frequent peaks in the spectral list

* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"

* deisotope_peaklist --> if the peak deisotoping should be performed, by preserving only the monoisotopic peak out of the isotope cluster

* envelope_peaklist --> if the peak enveloping should be performed, by preserving only the most intense peak out of the isotope cluster

* tolerance_ppm --> the tolerance (in ppm) for the peak-based spectral alignment; if NULL, the tolerance is set to 1000 ppm for "linear" TOF mode and to 100 ppm for "reflectron" TOF mode

#### Returns

* a list of aligned 'MALDIquant' spectra

#### Calls

* functions from the 'MALDIquant' package

### peak_picking

This function performs the peak picking

#### Inputs

* spectra --> a list of aligned 'MALDIquant' spectra

* peak_picking_algorithm --> the algorithm to be used for peak picking; it can be either "SuperSmoother" or "MAD"

* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"

* SNR --> the signal-to-noise ratio to be used for peak picking

* allow_parallelization --> if the parallel computation should be exploited

* deisotope_peaklist --> if the peak deisotoping should be performed, by preserving only the monoisotopic peak out of the isotope cluster

* envelope_peaklist --> if the peak enveloping should be performed, by preserving only the most intense peak out of the isotope cluster

* signals_to_take --> number of maximum number of peaks to retain; if less or equal than 0, all the peaks are retained

#### Returns

* a list of 'MALDIquant' peaks, in which each item corresponds to the peaks of a spectrum

#### Calls

## align_and_filter_peaks

This function performs the peak alignment and filtering

#### Inputs

* peaks --> a list of 'MALDIquant' peaks

* peak_picking_algorithm --> the algorithm to be used for peak picking; it can be either "SuperSmoother" or "MAD"

* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"

* peak_filtering_frequency_threshold_percent --> the percentage of spectra in which a peak has to be present for being preserved; peaks under the percentage frequency are discarded

* class_vector_for_peak_filtering --> a vector of the same length of the peaks list, to be used for class-wise peak frequency estimation; if NULL, the peak frequency is estimated onto the entire peaks list

* low_intensity_peak_removal_threshold_percent --> the percentage of the most intense peak intensity to be used as threshold for peak filtering: peaks with an intensity below the percentage of intensity of the most intense peak are discarded

* low_intensity_peak_removal_threshold_method --> how to determine the most intense peak; it can be "element-wise", the most intense peak is identified for each spectrum, or "whole dataset", the most intense peak corresponds to the most intense peak of the entire spectral dataset

* reference_peaklist --> a numeric vector of numbers to be considered as absolute reference values to which the peaks should be aligned to; or a character, indicating a peaklist to be used as reference for alignment (it can be "average spectrum" or "skyline spectrum"); if NULL, no absolute reference values are taken, and peaks are simply aligned to each other

* spectra --> a list of 'MALDIquant' spectra

* alignment_iterations --> number of peak alignment iterations to be run

* allow_parallelization --> if the parallel computation should be exploited

* tolerance_ppm --> the tolerance (in ppm) for the peak-based spectral alignment; if NULL, the tolerance is set to 1000 ppm for "linear" TOF mode and to 100 ppm for "reflectron" TOF mode

#### Returns

* a list of aligned and filtered 'MALDIquant' peaks

#### Calls

* functions from the 'MALDIquant' package

## intensityMatrix

* function from the 'MALDIquant' package

# ENSEMBLE MS TUNER

## model_ensemble_embedded_fs

This function performs feature selection, model training and tuning onto the training dataset, using different methods

#### Inputs

* training_set --> a matrix or dataframe (yielded by the MS PEAKLIST EXPORT functions), in which row represent spectra (single spectra per patient or average spectra per patient) and columns represent (aligned) peaks in the spectral dataset

* features_to_select --> the (maximum) number of features to preserve

* common_features_to_select --> the number of common features among the different feature selection methods to preserve; if 0, all the features (common and uncommon) are preserved, up to a maximum number defined by 'features_to_select'

* model_tuning --> if and when to perform the tuning of the model hyperparameters; it can be "after" the feature selection phase, or "embedded" within the feature selection process

* selection_metric --> the metric for selecting the best discriminatory features; it can be "Accuracy" or "Kappa"

* discriminant_attribute --> the name of the column of the training set matrix that contains the outcome variable

* non_features --> a vector of character indicating the columns to be discarded from the feature selection step

* seed --> the seed for defining randomness

* automatically_select_features --> if to allow the algorithm to automatically determine the number of features to preserve, up to a maximum number defined by 'features_to_select'

* generate_plots --> if to generate plots during feature selection

* cv_repeats_control --> the number of cross-validation iterations during the training phase

* k_fold_cv_control --> the k of the k-fold cross-validation

* preprocessing --> the preprocessing to be applied to the feature values (i.e., "center" or "scale")

* allow_parallelization --> if the parallel computation should be exploited

* feature_reranking --> whether the predictor rankings are recomputed on the model on the reduced feature set

* try_combination_of_parameters --> if to allow to iteratively try all the combination of preprocessing and feature reranking parameters

* outcome_list --> a character vector, corresponding to the levels of the outcome variable, indicating if the corresponding outcome is bening or malignant (for coloring the pixels in red and green in the next phase)

* progress_bar --> if a progress bar should be displayed ("tcltk" or "text" style)

* test_set --> a matrix or dataframe corresponding to the test set, for the validation phase

#### Returns

* a list containing:
  * model_list --> a lit containing, for each classification model:
    * the model itself (object)
    * the features preserved by the model
    * the class list
    * an identification name
    * the model performances
    * the cross-validation and/or external validation confusion matrix
    * the model hyperparameters
  * feature_list --> the list of features preserved for each feature selection method
  * common_features_list --> list of the common features selected among the different feature selection methods
  * training_set_common_features --> the same training set but with only the common features selected among the different feature selection methods
  * common_features_matrix --> a matrix with a column listing the common features selected among the different feature selection methods
  * model_performance_matrix --> a matrix displaying the performances for each classification model after the feature selection step
  * model_performance_parameter_list --> a named list containing, for each classification model, the complete list of classification parameters, for each outcome class, e.g. sensitivity, specificity, etc...

#### Calls

* automated_embedded_rfe (for each model)

* extract_feature_list_from_model_list

* extract_common_features_from_model_list

## automated_embedded_rfe

This function performs feature selection, model training and tuning onto the training dataset, using different methods

#### Inputs

* training_set --> a matrix or dataframe (yielded by the MS PEAKLIST EXPORT functions), in which row represent spectra (single spectra per patient or average spectra per patient) and columns represent (aligned) peaks in the spectral dataset

* features_to_select --> the (maximum) number of features to preserve

* selection_method --> the selection method (i.e. classifier) to be used during the feature selection phase

* model_tuning --> if and when to perform the tuning of the model hyperparameters; it can be "after" the feature selection phase, or "embedded" within the feature selection process

* model_tune_grid --> the tune grid of the model

* selection_metric --> the metric for selecting the best discriminatory features; it can be "Accuracy" or "Kappa"

* cv_repeats_control --> the number of cross-validation iterations during the training phase

* k_fold_cv_control --> the k of the k-fold cross-validation

* discriminant_attribute --> the name of the column of the training set matrix that contains the outcome variable

* non_features --> a vector of character indicating the columns to be discarded from the feature selection step

* seed --> the seed for defining randomness

* automatically_select_features --> if to allow the algorithm to automatically determine the number of features to preserve, up to a maximum number defined by 'features_to_select'

* generate_plots --> if to generate plots during feature selection

* preprocessing --> the preprocessing to be applied to the feature values (i.e., "center" or "scale")

* allow_parallelization --> if the parallel computation should be exploited

* feature_reranking --> whether the predictor rankings are recomputed on the model on the reduced feature set

* test_set --> a matrix or dataframe corresponding to the test set, for the validation phase

* positive_class_cv --> the class that should be considered as the positive class for performance evaluation; if NULL, the first of the classes is taken

* try_combination_of_parameters --> if to allow to iteratively try all the combination of preprocessing and feature reranking parameters

#### Returns

* a list containing:
  * training_set_feature_selection --> the same training set but with only the features selected among the best feature selection method
  * predictors_feature_selection --> the list of predictors preserved
  * feature_selection_graphics --> plots obtained during the feature selection phase
  * model_tuning_graphics --> plots obtained during the tuning phase
  * feature_weights --> the weight for each feature
  * variable_importance --> the importance level of each variable in the feature selection phase
  * fs_model_performance --> the model performances
  * feature_selection_model --> the model itself (object)
  * class_list --> the levels of the 'discriminant_attribute' variable
  * pie_chart_classification --> a pie chart displaying the rate of hits and misclassification
  * model_roc --> the ROC curve of the model
  * external_validation_confusion_matrix --> the confusion matrix computed during the external validation phase
  * external_validation_confusion_matrix_df --> the confusion matrix computed during the external validation phase (as dataframe)
  * cross_validation_confusion_matrix --> the confusion matrix computed during the cross-validation phase
  * cross_validation_confusion_matrix_df --> the confusion matrix computed during the cross-validation phase (as dataframe)
  * model_cv_performance_parameter_list --> a named list containing, for each classification model, the complete list of classification parameters  (cross-validation), for each outcome class, e.g. sensitivity, specificity, etc...
  * model_external_performance_parameter_list --> a named list containing, for each classification model, the complete list of classification parameters (external validation), for each outcome class, e.g. sensitivity, specificity, etc...

#### Calls

* embedded_rfe

## embedded_rfe

This function performs feature selection, model training and tuning onto the training dataset, using a single methods

#### Inputs

* training_set --> a matrix or dataframe (yielded by the MS PEAKLIST EXPORT functions), in which row represent spectra (single spectra per patient or average spectra per patient) and columns represent (aligned) peaks in the spectral dataset

* features_to_select --> the (maximum) number of features to preserve

* selection_method --> the selection method (i.e. classifier) to be used during the feature selection phase

* model_tuning --> if and when to perform the tuning of the model hyperparameters; it can be "after" the feature selection phase, or "embedded" within the feature selection process

* model_tune_grid --> the tune grid of the model

* selection_metric --> the metric for selecting the best discriminatory features; it can be "Accuracy" or "Kappa"

* cv_repeats_control --> the number of cross-validation iterations during the training phase

* k_fold_cv_control --> the k of the k-fold cross-validation

* discriminant_attribute --> the name of the column of the training set matrix that contains the outcome variable

* non_features --> a vector of character indicating the columns to be discarded from the feature selection step

* seed --> the seed for defining randomness

* automatically_select_features --> if to allow the algorithm to automatically determine the number of features to preserve, up to a maximum number defined by 'features_to_select'

* generate_plots --> if to generate plots during feature selection

* preprocessing --> the preprocessing to be applied to the feature values (i.e., "center" or "scale")

* allow_parallelization --> if the parallel computation should be exploited

* feature_reranking --> whether the predictor rankings are recomputed on the model on the reduced feature set

* test_set --> a matrix or dataframe corresponding to the test set, for the validation phase

* positive_class_cv --> the class that should be considered as the positive class for performance evaluation; if NULL, the first of the classes is taken

#### Returns

* a list containing:
  * training_set_feature_selection --> the same training set but with only the features selected among the feature selection method
  * predictors_feature_selection --> the list of predictors preserved
  * feature_selection_graphics --> plots obtained during the feature selection phase
  * model_tuning_graphics --> plots obtained during the tuning phase
  * feature_weights --> the weight for each feature
  * variable_importance --> the importance level of each variable in the feature selection phase
  * fs_model_performance --> the model performances
  * feature_selection_model --> the model itself (object)
  * class_list --> the levels of the 'discriminant_attribute' variable
  * pie_chart_classification --> a pie chart displaying the rate of hits and misclassification
  * model_roc --> the ROC curve of the model
  * external_validation_confusion_matrix --> the confusion matrix computed during the external validation phase
  * external_validation_confusion_matrix_df --> the confusion matrix computed during the external validation phase (as dataframe)
  * cross_validation_confusion_matrix --> the confusion matrix computed during the cross-validation phase
  * cross_validation_confusion_matrix_df --> the confusion matrix computed during the cross-validation phase (as dataframe)
  * model_cv_performance_parameter_list --> a named list containing, for each classification model, the complete list of classification parameters  (cross-validation), for each outcome class, e.g. sensitivity, specificity, etc...
  * model_external_performance_parameter_list --> a named list containing, for each classification model, the complete list of classification parameters  (external validation), for each outcome class, e.g. sensitivity, specificity, etc...

#### Calls

* functions from the 'caret' package

# MS PIXEL TYPER

## spectral_classification

This function performs the ensemble classification of spectra

#### Inputs

* spectra_path --> the path (folder) in which the spectra are stored
* filepath_R --> the path of the .RData file, containing the model list from the 'ENSEMBLE MS TUNER' program
* model_list --> the 'model_list' extracted from the results of the 'ENSEMBLE MS TUNER' software
* model_performance_parameter_list -->
* classification_mode = c("pixel", "profile")
* peak_picking_algorithm --> the algorithm to be used for peak picking; it can be either "SuperSmoother" or "MAD"
* deisotope_peaklist --> if the peak deisotoping should be performed, by preserving only the monoisotopic peak out of the isotope cluster
* preprocessing_parameters --> a named list of preprocessing parameters to be used:
  * mass_range --> the mass range to which the spectra are cut to (the first number is the lower boundary, while the second number is the upper boundary); if  NULL, the spectra are brought to the same scale by using highest overlapping mass range
  * transformation_algorithm --> the algorithm to be used for data transformation; it can be "sqrt", "log", "log2", "log10", "sin", "cos", or "exp"; if NULL, no  smoothing is performed
  * smoothing_algorithm --> the algorithm to be used for smoothing; it can be "SavitzkyGolay" or "MovingAverage"; if NULL, no smoothing is performed
  * smoothing_strength --> the width of the smoothing filter to use; it can be "medium", "strong" or "stronger"
  * baseline_subtraction_algorithm --> the algorithm to be used for baseline subtraction; it can be "SNIP", "TopHat", "ConvexHull", or "median"; if NULL, no  baseline subtraction is performed
  * baseline_subtraction_algorithm_parameter --> the parameter to be used for defining the filter width for baseline subtraction
  * normalization_algorithm --> the algorithm to be used for normalization; it can be "TIC", "RMS", "PQN", or "median"; if NULL, no normalization is performed
  * normalization_mass_range --> the mass range to use for normalization: the first number is the lower boundary, while the second number is the upper boundary;  if NULL, the entire spectrum is used as the window for normalization
  * preprocess_spectra_in_packages_of --> the number of spectra to be processed simultaneously, in packets, if NULL or less or equal than 0, all the spectra are  preprocessed at once
  * spectral_alignment_algorithm --> the algorithm to be used for peak-based spectral alignment; it can be "cubic", "quadratic" or "linear", if NULL, no spectral  alignment is performed
  * spectral_alignment_reference --> the peak list to be used as reference for peak-based spectral alignment; it can be "average spectrum", "skyline spectrum" or "auto": the first two options allow to use the peaks of the average spectrum or the skyline spectrum as the reference for alignment, while the third option automatically estimates the peaks to be used for alignment as the most frequent peaks in the spectral list
* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"
* allow_parallelization --> if the parallel computation should be exploited
* decision_method_ensemble --> how the predictions of the individual classifiers should be combined: "unweighted majority" simply takes the mode among all votes, while "bayesian probabilities" uses the Bayesian framework to estimate the model's classification reliability and weigh the individual votes accordingly
* pixel_grouping --> how the pixels (i.e. spectra) of a patient's spectral dataset should be grouped: "single" takes individual spectra and predictions are performed spectrum-by-spectrum; "moving window average" slides a window, takes the average of spectra in the window, makes predictions onto the average spectrum and colors the pixel at the center of the window according to the prediction; "hca" runs a hierarchical clustering analysis to generate clusters of spectra, takes the average spectra from the clusters, makes predictions onto the average spectrum and colors the pixels belonging to the cluster according to the predictions
* moving_window_size --> the size (in pixels) of the sliding window ('moving window average')
* number_of_hca_nodes --> number of clusters to generate by the hierarchical clustering analysis ('hca')
* plot_figures --> if to generate the segmentation images
* plot_legends --> if to plot the figure legends and which legend ("sample name", color "legend", "plot name")
* progress_bar --> if a progress bar should be displayed ("tcltk" or "text" style)
* tolerance_ppm --> the tolerance (in ppm) for the peak-based spectral alignment; if NULL, the tolerance is set to 1000 ppm for "linear" TOF mode and to 100 ppm for "reflectron" TOF mode
* alpha --> in the weighing Bayesian system, if a probability of 0 is found, it is replaced by the alpha value, to avoid the generation of zeros that impairs the weighing system

#### Returns

* a list containing:
  * final_result_matrix_msi_list --> a named list containing, for each patient, a matrix with the predictions for all the models of the ensemble classifier for each spectrum (corresponding to a pixel) of a patient's dataset
  * final_result_matrix_profile_list --> a named list containing, for each patient, a matrix with the predictions for all the models of the ensemble classifier for each spectrum (corresponding to a pixel) of a patient's dataset
  * classification_ms_images_list --> a named list containing, for each patient, a red-green segmentation image, for each classification model, colored according to the predictions made onto the spectra of the patient's dataset
  * classification_ensemble_matrix_msi_all --> a matrix containing the ensemble predictions onto the spectra of the patient's dataset
  * classification_ensemble_matrix_profile_all --> a matrix containing the ensemble predictions onto the average spectra of the patients
  * classification_ensemble_ms_image_list --> a named list containing, for each patient, the red-green segmentation image derived from the predictions of the ensemble classifier onto the entire patient's spectral dataset
  * average_spectrum_with_bars_profile_list --> a named list of images of the average spectra of the patients, with overlay bars indicating the peaks chosen by each classifier along with the standard deviation of the peak intensities across the patient's dataset

#### Calls

* import_spectra

* preprocess_spectra

* single_model_classification_of_spectra

* ensemble_vote_classification

* outcome_and_class_to_MS

## single_model_classification_of_spectra

This function performs the classification of spectra using a single method

#### Inputs

* spectra --> a list of 'MALDIquant' spectra

* model_x

* model_name --> the name of the model

* preprocessing_parameters --> a named list of preprocessing parameters to be used:
  * mass_range --> the mass range to which the spectra are cut to (the first number is the lower boundary, while the second number is the upper boundary); if  NULL, the spectra are brought to the same scale by using highest overlapping mass range
  * transformation_algorithm --> the algorithm to be used for data transformation; it can be "sqrt", "log", "log2", "log10", "sin", "cos", or "exp"; if NULL, no  smoothing is performed
  * smoothing_algorithm --> the algorithm to be used for smoothing; it can be "SavitzkyGolay" or "MovingAverage"; if NULL, no smoothing is performed
  * smoothing_strength --> the width of the smoothing filter to use; it can be "medium", "strong" or "stronger"
  * baseline_subtraction_algorithm --> the algorithm to be used for baseline subtraction; it can be "SNIP", "TopHat", "ConvexHull", or "median"; if NULL, no  baseline subtraction is performed
  * baseline_subtraction_algorithm_parameter --> the parameter to be used for defining the filter width for baseline subtraction
  * normalization_algorithm --> the algorithm to be used for normalization; it can be "TIC", "RMS", "PQN", or "median"; if NULL, no normalization is performed
  * normalization_mass_range --> the mass range to use for normalization: the first number is the lower boundary, while the second number is the upper boundary;  if NULL, the entire spectrum is used as the window for normalization
  * preprocess_spectra_in_packages_of --> the number of spectra to be processed simultaneously, in packets, if NULL or less or equal than 0, all the spectra are  preprocessed at once
  * spectral_alignment_algorithm --> the algorithm to be used for peak-based spectral alignment; it can be "cubic", "quadratic" or "linear", if NULL, no spectral  alignment is performed
  * spectral_alignment_reference --> the peak list to be used as reference for peak-based spectral alignment; it can be "average spectrum", "skyline spectrum" or "auto": the first two options allow to use the peaks of the average spectrum or the skyline spectrum as the reference for alignment, while the third option automatically estimates the peaks to be used for alignment as the most frequent peaks in the spectral list

* peak_picking_algorithm --> the algorithm to be used for peak picking; it can be either "SuperSmoother" or "MAD"

* deisotope_peaklist --> if the peak deisotoping should be performed, by preserving only the monoisotopic peak out of the isotope cluster

* envelope_peaklist --> if the peak enveloping should be performed, by preserving only the most intense peak out of the isotope cluster

* peak_picking_SNR --> the signal-to-noise ratio to be used for peak picking

* peak_filtering_frequency_threshold_percent --> the percentage of spectra in which a peak has to be present for being preserved; peaks under the percentage frequency are discarded

* low_intensity_peak_removal_threshold_percent --> the percentage of the most intense peak intensity to be used as threshold for peak filtering: peaks with an intensity below the percentage of intensity of the most intense peak are discarded

* low_intensity_peak_removal_threshold_method --> how to determine the most intense peak; it can be "element-wise", the most intense peak is identified for each spectrum, or "whole dataset", the most intense peak corresponds to the most intense peak of the entire spectral dataset

* tof_mode --> the mode in which the Time of Flight (TOF) is used; it can be either "linear" or "reflectron"

* allow_parallelization --> if the parallel computation should be exploited

* pixel_grouping --> how the pixels (i.e. spectra) of a patient's spectral dataset should be grouped: "single" takes individual spectra and predictions are performed spectrum-by-spectrum; "moving window average" slides a window, takes the average of spectra in the window, makes predictions onto the average spectrum and colors the pixel at the center of the window according to the prediction; "hca" runs a hierarchical clustering analysis to generate clusters of spectra, takes the average spectra from the clusters, makes predictions onto the average spectrum and colors the pixels belonging to the cluster according to the predictions

* number_of_hca_nodes --> number of clusters to generate by the hierarchical clustering analysis ('hca')

* moving_window_size --> the size (in pixels) of the sliding window ('moving window average')

* plot_figures --> if to generate the segmentation images

* plot_legends --> if to plot the figure legends and which legend ("sample name", color "legend", "plot name")

* tolerance_ppm --> the tolerance (in ppm) for the peak-based spectral alignment; if NULL, the tolerance is set to 1000 ppm for "linear" TOF mode and to 100 ppm for "reflectron" TOF mode

#### Returns

* a list containing:
  * result_matrix_model
  * classification_msi_model
  * final_result_matrix
  * average_spectrum_with_bars

#### Calls

* preprocess_spectra

* outcome_and_class_to_MS

* functions from the 'caret' package

## ensemble_vote_classification

This function performs the combination of the votes coming from the individual classifiers

#### Inputs

* classification_matrix --> the result matrix of an ensemble classification, in which each row is an observation/spectrum (patient or pixel) and each column is the predicted class of that observation by one model

* class_list --> a vector of classes; if null, it is estimated from the input matrix

* weighted_decision_method --> how the predictions of the individual classifiers should be combined: "unweighted majority" simply takes the mode among all votes, while "bayesian probabilities" uses the Bayesian framework to estimate the model's classification reliability and weigh the individual votes accordingly

* performance_parameter_list --> a named list with the performance parameters for each classification model, exported from the 'ENSEMBLE MS TUNER' software

* type_of_validation_for_performance_estimation --> the parameters from which validation (cross-validation "cv" or external validation "ev") to be used

* alpha --> in the weighing Bayesian system, if a probability of 0 is found, it is replaced by the alpha value, to avoid the generation of zeros that impairs the weighing system

#### Returns

* classification_ensemble_matrix --> a single column matrix with the ensemble classification results

## outcome_and_class_to_MS

This function performs the association of the levels of the response variable with the outcome ('benign' and 'malignant') to color the pixels associated with benignity in green and the pixels associated with malignancy in red

#### Inputs

* class_list --> a vector of classes

* outcome_list --> a character vector, corresponding to the levels of the outcome variable, indicating if the corresponding outcome is bening or malignant

* class_vector --> a vector containing the predicted class for each observation (usually, the column of the classification matrix with all the predictions)

#### Returns

* a list containing:
  * class_outcome_matrix --> a matrix that associates the levels of the response variable ('class_list') witht he outcome ('outcome_list')
  * class_vector_as_numeric --> the same 'class_vector' but with the classes replaced by a number (resembling the color of the pixels) according to the association with the outcome
  * legend_text --> the text of the figure legend
  * legend_fill --> the color code to fill the squares in the figure legend







# EXAMPLE WORKFLOW
## Import spectra
'
# Import spectra from files, in different formats
spectra <- import_spectra(filepath = "/folder/where/the/spectra/are", spectra_format = "imzML", mass_range = c(3000, 15000), allow_parallelization = FALSE, spectral_names = "name", replace_sample_name_field = TRUE, remove_empty_spectra = TRUE)
'

## Preprocess spectra
'
# Preprocess the imported spectra, by defining the preprocessing parameters
spectra <- preprocess_spectra(spectra, tof_mode = "linear", preprocessing_parameters = list(mass_range = c(3000, 15000), transformation_algorithm = NULL, smoothing_algorithm = "SavitzkyGolay", smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_algorithm_parameter = 100, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_algorithm = NULL, spectral_alignment_reference = "average_spectrum"), allow_parallelization = FALSE, tolerance_ppm = 1000)
'

## Feature extraction (i.e. peak picking and alignment)
'
# Perform peak picking onto the spectra
peaks <- peak_picking(spectra, peak_picking_algorithm = "SuperSmoother", tof_mode = "linear", SNR = 3, allow_parallelization = FALSE, deisotope_peaklist = FALSE, envelope_peaklist = FALSE, signals_to_take = 50)

# Perform peak alignment and filtering
peaks <- align_and_filter_peaks(peaks, peak_picking_algorithm = "SuperSmoother", tof_mode = "linear", peak_filtering_frequency_threshold_percent = 5, class_vector_for_peak_filtering = NULL, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", reference_peaklist = NULL, spectra = spectra, alignment_iterations = 5, allow_parallelization = FALSE, tolerance_ppm = 1000)

# Generate the peaklist statistical matrix
peaklist_matrix <- intensityMatrix(peaks = peaks, spectra = spectra)
# Add the "Class" and the "Sample" columns (additional info to the statistical matrix for statistical analysis)
peaklist_matrix <- matrix_add_class_and_sample(peaklist_matrix, peaks = peaks, class_list = c("benign", "malignant"), spectra_format = "imzML", sample_output = TRUE, class_output = TRUE, row_labels = "Sample")
'

## Model ensemble training and tuning
'
# Train and tune the classifiers of the ensemble on the entire peaklist matrix
model_ensemble_tuner <- model_ensemble_embedded_fs(training_set = peaklist_matrix, features_to_select = 20, common_features_to_select = 0, model_tuning = "after", selection_metric = "Accuracy", discriminant_attribute = "Class", non_features = c("Sample", "Class"), seed = 12345, automatically_select_features = FALSE, generate_plots = TRUE, cv_repeats_control = 5, k_fold_cv_control = 3, preprocessing = c("center", "scale"), allow_parallelization = TRUE, feature_reranking = FALSE, try_combination_of_parameters = FALSE, outcome_list = c("benign", "malignant"), progress_bar = NULL, test_set = NULL)

# Extract the variables
model_list <- model_ensemble_tuner$model_list
feature_list <- model_ensemble_tuner$feature_list
model_performance_matrix <- model_ensemble_tuner$model_performance_matrix
common_features_matrix <- model_ensemble_tuner$common_features_matrix
model_performance_parameter_list <- model_ensemble_tuner$model_performance_parameter_list
'

## Ensemble pixel-by-pixel classification
'
# Perform the spectral classification
classification_of_patients <- spectral_classification(spectra_path = "/folder/where/the/spectra/to/be/classified/are", filepath_R = NULL, model_list = model_list, model_performance_parameter_list = model_performance_parameter_list, classification_mode = "pixel", peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = FALSE, preprocessing_parameters = list(mass_range = c(4000,15000), transformation_algorithm = NULL, smoothing_algorithm = "SavitzkyGolay", smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_algorithm_parameter = 100, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_algorithm = NULL), tof_mode = "linear", allow_parallelization = FALSE, decision_method_ensemble = c("unweighted majority", "bayesian probabilities"), pixel_grouping = "single", seed = 12345, plot_figures = TRUE, plot_legends = c("sample name", "legend", "plot name"), progress_bar = NULL, tolerance_ppm = 2000, alpha = 0.05)
'