################## FEATURE SELECTION FUNCTIONS

###INSTALL THE REQUIRED PACKAGES
required_packages <- c ("parallel", "caret", "tcltk", "e1071", "doMC", "kernlab", "pROC")
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




############ LOAD THE REQUIRED PACKAGES
library (parallel)
library (caret)
library (tcltk)
library (doMC)
library (e1071)
library (kernlab)
library (pROC)




############################## MULTICORE
# Detect the number of cores
cpu_thread_Number <- detectCores(logical=TRUE)
cpu_core_number <- cpu_thread_number/2
# Register the foreach backend
registerDoMC(cores = cpu_core_number)
##############################





################################################################################
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
	library(p)
}
}





####################################### READ (SOURCE) SCRIPT FROM URL
source_github <- function(url, ...) {
  # Load required package (install it if not present)
  if ("RCurl" %in% installed.packages()[,1] == FALSE) {
	  install.packages("RCurl")
  }
  require(RCurl)
  # Parse and evaluate each R script (in the global environement)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}





####################################### FEATURE SELECTION BASED UPON CORRELATION
correlation_feature_selection <- function (peaklist, correlation_method="pearson", non_features=c("Sample", "Class", "THY")) {
# Take only the part of the matrix without Class and Sample
peaklist_features <- peaklist [,!(names(peaklist) %in% non_features)]
# Compute the correlation
feature_correlation <- cor(peaklist_features, method=correlation_method)
# Output the highly correlated features
highly_correlated <- find_correlation(feature_correlation, correlation_threshold)
# Keep only the part of the matrix that is not highly correlated
peaklist_features_low_correlation <- peaklist_features [,-highly_correlated]
# Sort the dataframe column by a non decreasing order
peaklist_features_low_correlation <- peaklist_features_low_correlation [,order(as.numeric(names(peaklist_features_low_correlation)))]
low_correlated_features <- names(peaklist_features_low_correlation)
# Add back the non Features (Sample and Class)
for (i in 1:length(non_features)) {
	peaklist_features_low_correlation <- data.frame (peaklist_features_low_correlation, peaklist[,non_features[i]])
}
names(peaklist_features_low_correlation) <- c(low_correlated_features, non_features)
# Return the dataframe
cat (paste("The number of selected features is:", length(low_correlated_features)))
return (peaklist_features_low_correlation)
}





################################################################################





################################ DEFINED SPECTRA GROUPING (PEAKLIST MATRIX)
# Generate a peaklist dataframe with a certain number of rows per patient
group_peaklist <- function (peaklist, rows_per_patient=100, seed=0, balanced=TRUE, discard_poor_samples=FALSE, discard_if_lower_than=100, non_features=c("Sample","Class"), algorithm="random") {
# Plant the seed only if a specified value is entered
if (seed != 0) {
	# Make the randomness reproducible
	set.seed(seed)
}
result_data_frame <- data.frame()
# Create the file vector and the patient vector
file_vector <- peaklist$Sample
patient_vector <- levels(file_vector)
#### If the dataset has to be balanced...
if (balanced == TRUE) {
	########## All the patients should have the same number of rows, always
	### Determine the shortest patient's spectra number
	lowest_observation_number <- NULL
	# Scroll the other patients
	for (p in 1:length(patient_vector)) {
		# Determine the number of observations for this patient
		observation_number <- nrow(peaklist [peaklist$Sample==patient_vector[p],])
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
		patient_data_frame <- peaklist [peaklist$Sample==patient_vector[p],]
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
		patient_data_frame <- peaklist [peaklist$Sample==patient_vector[p],]
		patient_data_frame_rows <- nrow(patient_data_frame)
		# Continue only if nothing has to be discarded or the dimension of the patient dataset is higher than the threshold
		if (discard_poor_samples == FALSE || (discard_poor_samples == TRUE && patient_data_frame_rows >= discard_if_lower_than)) {
			########################### RANDOMNESS
			if (algorithm == "random") {
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
				index <- createFolds(patient_data_frame$Sample, k = rows_per_patient)
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
			if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
				# Compute the distance matrix
				distance_matrix <- dist(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], method="euclidean")
				hca <- hclust(distance_matrix)
				# Associate to each row/spectrum the subgroup of the tree to which it belongs
				# k must be between 1 and 56
				if (rows_per_patient > 56) {
					rows_per_patient <- 56
				}
				hca_groups <- cutree(hca, k=rows_per_patient)
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
			if (algorithm == "kMeans" | algorithm == "kmeans" | algorithm == "k-Means") {
				# Compute the k-Means clustering
				k_means <- kmeans(patient_data_frame[,!(names(patient_data_frame) %in% non_features)], centers=rows_per_patient)
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
return (result_data_frame)
}





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





##################### MATRIX/DATAFRAME SPLITTING FUNCTION (TRAINING AND TESTING)
# Output which patient (for each class) should be taken for training and for testing
matrix_splitting_training_test <- function (peaklist, class_list=list(), seed=0, percentage_of_observation_for_training=50) {
# Plant the seed only if a specified value is entered
if (seed != 0) {
	# Make the randomness reproducible
	set.seed(seed)
}
############################# If there are no classes or one class, just randomly select rows on the entire dataset
if (length(class_list) <= 1) {
index_training <- sample(nrow(peaklist), size=round(nrow(peaklist)*percentage_of_observation_for_training/100))
}
############################# If there are more classes, randomly select the rows for each class
if (length(class_list) > 1) {
class_list <- sort(class_list)
# Automatically generate the corresponding class if not specified in the input matrix
########### Create the dataframe with only patients and classes
# Create two dataframe columns with the patients
sample_and_class <- data.frame(unique(peaklist$Sample), unique(peaklist$Sample))
# Give the correct names to the columns
colnames(sample_and_class) <- c ("Sample", "Class")
# Define the sample column as character (not as factor)
sample_and_class$Sample <- as.character (sample_and_class$Sample)
sample_and_class$Class <- as.character (sample_and_class$Class)
# Replace the patients in the second column with the class
for (j in 1:length(class_list)) {
	for (i in 1:length(sample_and_class$Class)) {
		if (length(grep(class_list[j], sample_and_class$Class[i], ignore.case=TRUE)) != 0) {
			sample_and_class$Class[i] <- class_list [j]
		}
	}
}
# Create a list of dataframes, one for each class
class_data_frame_list <- list()
for (j in 1:length(class_list)) {
	class_data_frame_list[[j]] <- sample_and_class [sample_and_class$Class == class_list[j],]
}
############ Create a list containing the training and testing indexes of the class dataframes
index_training <- list()
for (i in 1:length(class_data_frame_list)) {
	index_training[[i]] <- createDataPartition(y = class_data_frame_list[[i]]$Class, p = (percentage_of_observation_for_training/100), list = FALSE)
}
names(index_training) <- class_list
}
# Now we have randomly selected some patients for the training and some others for the testing, for each class
####################
# The indexTraining is integer if there are no classes or one class, and rows are selected randomly
if (is.integer(index_training)) {
	training_dataset <- peaklist [index_training,]
	test_dataset <- peaklist [-index_training,]
}
# The indexTraining is a list if there are more than one classes (one element per class),
# and rows are selected randomly for each class
if (is.list(index_training)) {
##### Create the training and the testing patient dataframes
training_data_frame_list <- list()
test_data_frame_list <- list()
# For each class (and so for each element of the list with the indexes per class)
for (i in 1:length(index_training)) {
	# Create the two complementary dataframes: randomly select the patients for each class
	# for training and testing
	training_data_frame_list[[i]] <- class_data_frame_list[[i]][index_training[[i]],]
	test_data_frame_list[[i]] <- class_data_frame_list[[i]][-index_training[[i]],]
}
#### In the original matrix, select only the rows of the selected samples
# Turn the patient lists into a data frame and store all the sample values in a vector
patientsForTraining <- character ()
for (i in 1:length(training_data_frame_list)) {
	patients_for_training <- append(patients_for_training, training_data_frame_list[[i]]$Sample)
}
patients_for_test <- character ()
for (i in 1:length(test_data_frame_list)) {
	patients_for_test <- append(patients_for_test, test_data_frame_list[[i]]$Sample)
}
# Create the two datasets
training_dataset <- data.frame()
test_dataset <- data.frame()
for (p in 1:length(patients_for_training)) {
	training_dataset <- rbind(training_dataset, peaklist[(peaklist$Sample==patients_for_training[p]),])
}
for (p in 1:length(patients_for_test)) {
	test_dataset <- rbind(test_dataset, peaklist[(peaklist$Sample==patients_for_test[p]),])
}
}
return (list(training_dataset = training_dataset, testDataset = test_dataset))
}





################################################################################





############################## RECURSIVE FEATURE ELIMINATION - FEATURE SELECTION
recursive_feature_elimination <- function (peaklist, features_to_select=20, selection_method="pls", discriminant="Class", non_features=c("Sample", "Class", "THY"), seed=0) {
if (seed != 0) {
	set.seed(seed)
}
########################################################### RFE MODEL (with PLS)
rfe_ctrl <- rfeControl (functions = caretFuncs,
	method = "cv",
	repeats = 3,
	number = 10
	)
# The simulation will fit models with subset sizes: (the subset size is the number of predictors to use)
if (selection_method != "") {
	if (seed != 0) {
		set.seed(seed)
	}
	rfe_model <- rfe (x = peaklist [,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant],
		sizes = features_to_select,
		rfeControl = rfe_ctrl,
		method = selection_method,
		)
} else {
	if (seed != 0) {
		set.seed(seed)
	}
	rfe_model <- rfe (x = peaklist [,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant],
		sizes = features_to_select,
		rfeControl = rfe_ctrl,
		)
}
# Output the best predictors after the RFE
predictors_rfe <- predictors (rfe_model) [1:features_to_select]
# Take the selected features
peaklist_rfe <- peaklist [,predictors_rfe]
# Add the non features back
for (i in 1:length(non_features)) {
	peaklist_rfe <- cbind(peaklist_rfe , peaklist[,non_features[i]])
}
names(peaklist_rfe) <- c(as.character(predictors_rfe), non_features)
#
return (peaklist_rfe)
}





################################################################################





############################################### TRUNCATE PEAKLIST MATRIX
truncate_peaklist <- function (peaklist, range=c(4000, 15000), non_features=c("Sample", "Class")){
features_to_keep <- numeric()
# Do not do anything if the values are set to zero
if (range[1] == 0 & range[2] == 0) {
	# Keep all the features
	features_to_keep <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
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
peaklistTruncated <- data.frame (peaklist [,as.character(features_to_keep)])
# Add the non features back
for (i in 1:length(non_features)) {
	peaklistTruncated <- data.frame (peaklistTruncated, peaklist[,non_features[i]])
}
names(peaklistTruncated) <- c(features_to_keep, non_features)
#
return (peaklistTruncated)
}





################################################################################





##################### CROSS-VALIDATION FOR A SUPPORT VECTOR MACHINE MODEL
################ One row per patient (average)
cross_validation_svm <- function (training_dataset, seed=0, k_fold=10, repeats=5, svm_model, non_features=c("Sample","Class","THY"), positive_class=levels(training_dataset$Class)[1]) {
if (seed != 0) {
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
result_matrix <- matrix (NA, ncol=14, nrow=0)
colnames(result_matrix) <- c("Accuracy", "Kappa", "No information rate", "Accuracy p-value", "McNiemar p-value", "Sensitivity", "Specificity", "PPV", "NPV", "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy", "ROC AUC")
#### For each repetition...
for (i in 1:repeats) {
	if (seed != 0) {
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
		if (seed != 0) {
			set.seed(seed)
		}
		model <- svm (x=training_predictors, y=training_outcomes, scale=TRUE, kernel=svm_kernel, degree=svm_degree, gamma=svm_gamma, cost=svm_cost)
		# Use the model to predict the testing subset
		if (seed != 0) {
			set.seed(seed)
		}
		predicted_classes <- predict(model, newdata=test_predictors)
		# Compute the performances
		performances <- confusionMatrix(data=predicted_classes, test_outcomes, positive=positive_class)
		confusion_matrix_list [[k*i]] <- performances
		#### ROC analysis
		svm_roc <- roc(response=test_outcomes, predictor=as.numeric(predicted_classes))
		# Create a matrix row with the partial results
		result_matrix_row <- matrix (0, ncol=ncol(result_matrix), nrow=1)
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
result_matrix_average <- apply(result_matrix, 2, mean, na.rm=TRUE)
return (list(confusion_matrix=confusion_matrix_list, result_matrix=result_matrix, result_matrix_average=result_matrix_average, seed=seed, features=names(training_dataset)[!(names(training_dataset) %in% non_features)]))
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





############################# SVM TUNING AND VALIDATION
svm_tuning_and_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), gamma=10^(-5:5), cost=10^(-5:5), epsilon=seq(1,2,by=1), degree=1:5, kernel="radial", k_fold_cv=2, repeats_cv=5, positive_class_cv="HP", seed=0, pca=FALSE, numer_of_components=3) {
if (seed != 0) {
	set.seed(seed)
}
############################################################################ PCA
if (pca == TRUE) {
# Compute the PCs
pca_training <- prcomp (as.matrix(peaklist_training [,!(names(peaklist_training) %in% non_features)]))
pca_test <- prcomp (as.matrix(peaklist_test [,!(names(peaklist_test) %in% non_features)]))
#################### SVM tuning
if (seed != 0) {
	set.seed(seed)
}
svm_tuning <- tune.svm (as.data.frame(pca_training$x[,1:number_of_components]), factor(peaklist_training$Class), gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
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
parameters_output <- list (kernel=svm_kernel, cost=svm_tuning$best.model$cost, degree=svm_tuning$best.model$degree, epsilon=svm_tuning$best.model$epsilon, gamma=svm_tuning$best.model$gamma)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predictthe outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed(seed)
	}
	predicted_classes_svm <- predict(svm_model, newdata = as.data.frame(pca_test$x[,1:number_of_components]))
	# Create the outcomes dataframe
	classification_results_svm <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	svm_roc <- list()
	roc_curve <- roc(response=classification_results_svm$true, predictor=as.numeric(classification_results_svm$predicted))
	svm_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste ("ROC area under the curve:", roc_curve$auc)
	legend ("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	svm_roc[[2]] <- recordPlot()
	# Output the results
	return (list (model=svm_model, svm_parameters=parameters_output, classification_results=classification_results_svm, performances=test_performances_svm, roc=svm_roc))
}	else {return (list (model=svm_model, svm_parameters=parameters_output, cross_validation=cv_svm_model))}
} else {
####################################################################### FEATURES
################### Optimise according to the kernel
# Find the best tuning parameters for the SVM
if (kernel == "radial") {
	if (seed != 0) {
		set.seed(seed)
	}
	svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
}
if (kernel == "polynomial") {
	if (seed != 0) {
		set.seed(seed)
	}
	svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon, degree=degree)
}
if (kernel == "linear") {
	if (seed != 0) {
		set.seed(seed)
	}
	svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma=gamma, cost=cost, kernel=kernel)
}
if (kernel == "sigmoid") {
	if (seed != 0) {
		set.seed(seed)
	}
	svm_tuning <- tune.svm(peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class), gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
}
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
parameters_output <- list (kernel=svm_kernel, cost=svm_tuning$best.model$cost, degree=svm_tuning$best.model$degree, epsilon=svm_tuning$best.model$epsilon, gamma=svm_tuning$best.model$gamma)
#################### CROSS-VALIDATION
cv_svm_model <- cross_validation_svm (peaklist_training, k_fold=k_fold_cv, repeats=repeats_cv, svm_model=svm_model, non_features=non_features, positive_class=positive_class_cv, seed=seed)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed(seed)
	}
	predicted_classes_svm <- predict(svm_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
	# Create the outcomes dataframe
	classification_results_svm <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	svm_roc <- list()
	roc_curve <- roc(response=classification_results_svm$true, predictor=as.numeric(classification_results_svm$predicted))
	svm_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	svm_roc[[2]] <- recordPlot()
	# Output the results
	return (list(model=svm_model, cross_validation=cv_svm_model, classification_results=classification_results_svm, parameters=parameters_output, performances=test_performances_svm, roc=svm_roc))
}	else {return (list(model=svm_model, parameters=parameters_output, cross_validation=cv_svm_model))}
	}
}





################################################################################





############################# SVM VALIDATION
svm_validation <- function (peaklist_training, peaklist_test=NULL, non_features=c("Sample","Class","THY"), gamma=0.1, cost=10, epsilon=0.1, degree=3, kernel="radial", k_fold_cv=10, repeats_cv=2, positive_class_cv="HP", seed=0, pca=FALSE, numer_of_components=3) {
if (seed != 0) {
	set.seed(seed)
}
############################################################################ PCA
if (pca == TRUE) {
# Compute the PCs
pca_training <- prcomp(as.matrix(peaklist_training [,!(names(peaklist_training) %in% non_features)]))
pca_test <- prcomp(as.matrix(peaklist_test [,!(names(peaklist_test) %in% non_features)]))
# SVM tuning
if (seed != 0) {
	set.seed(seed)
}
# SVM with defined parameters
svm_model <- svm (pca_training, factor(peaklist_training$Class), kernel=kernel, cost=cost, epsilon=epsilon, gamma=gamma)
parameters_output <- list (kernel=kernel, cost=cost, degree=degree, epsilon=epsilon, gamma=gamma)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed(seed)
	}
	predicted_classes_svm <- predict(svm_model, newdata = as.data.frame(pca_test$x[,1:number_of_components]))
	# Create the outcomes dataframe
	classification_results_svm <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	svm_roc <- list()
	roc_curve <- roc (response=classification_results_svm$true, predictor=as.numeric(classification_results_svm$predicted))
	svm_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	svm_roc[[2]] <- recordPlot()
	# Output the results
	return (list (model=svm_model, svm_parameters=parameters_output, classification_results=classification_results_svm, performances=test_performances_svm, roc=svm_roc))
}	else {return (list (model=svm_model, svm_parameters=parameters_output, cross_validation=cvsvm_model))}
} else {
####################################################################### FEATURES
# SVM with defined parameters
svm_model <- svm (peaklist_training [,!(names(peaklist_training) %in% non_features)], factor(peaklist_training$Class),
		kernel=kernel, cost=cost, epsilon=epsilon, gamma=gamma)
parameters_output <- list (Kernel=kernel, Cost=cost, Degree=degree, Epsilon=epsilon, Gamma=gamma)
#################### CROSS-VALIDATION
cvsvm_model <- cross_validation_svm (peaklist_training, k_fold=k_fold_cv, repeats=repeats_cv, svm_model=svm_model, non_features=non_features, positive_class=positive_class_cv, seed=seed)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix(peaklist_test) || is.data.frame(peaklist_test)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed(seed)
	}
	predicted_classes_svm <- predict(svm_model, newdata = peaklist_test [,!(names(peaklist_test) %in% non_features)])
	# Create the outcomes dataframe
	classification_results_svm <- data.frame (sample=peaklist_test$Sample, predicted=predicted_classes_svm, true=peaklist_test$Class)
	### Generate the confusion matrix to evaluate the performances
	test_performances_svm <- confusionMatrix(data = predicted_classes_svm, peaklist_test$Class, positive=positive_class_cv)
	#### ROC analysis
	svm_roc <- list()
	roc_curve <- roc (response=classification_results_svm$true, predictor=as.numeric(classification_results_svm$predicted))
	svm_roc[[1]] <- roc_curve$auc
	plot (roc_curve)
	roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
	legend("bottomright", legend=roc_legend, xjust=0.5, yjust=0.5)
	svm_roc[[2]] <- recordPlot()
	# Output the results
	return (list (model=svm_model, cross_validation=cvsvm_model, classification_results=classification_results_svm, parameters=parameters_output, performances=test_performances_svm, roc=svm_roc))
}	else {return (list (model=svm_model, parameters=parameters_output, cross_validation=cvsvm_model))}
	}
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





################################## ROUND NUMERIC FEATURES PEAKLIST
roundFeaturesPeaklist <- function (peaklist, decimal_digits=3, non_features=c("Sample","Class","THY")) {
# Take the features
feature_list <- as.numeric(names(peaklist [,!(names(peaklist) %in% non_features)]))
# Round the features
for (i in 1:length(feature_list)) {
	feature_list[i] <- round(feature_list[i], digits=decimal_digits)
}
# Regather the features
features <- c(feature_list, non_features)
names(peaklist) <- features
#
return (peaklist)
}
