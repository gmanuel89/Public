################## FEATURE SELECTION FUNCTIONS

###INSTALL THE REQUIRED PACKAGES
requiredPackages <- c ("parallel", "caret", "tcltk", "e1071", "doMC", "kernlab", "pROC")
installedPackages <- installed.packages () [,1]
missingPackages <- character()
for (p in 1:length(requiredPackages)) {
	if ((requiredPackages[p] %in% installedPackages) == FALSE) {
		missingPackages <- append (missingPackages, requiredPackages[p])
	}
}
if (length(missingPackages) > 0) {
	install.packages (missingPackages)
}




############ LOAD THE REQUIRED PACKAGES
library ("parallel")
library ("caret")
library ("tcltk")
library ("doMC")
library ("e1071")
library ("kernlab")
library ("pROC")




############################## MULTICORE
# Detect the number of cores
CPUthreadNumber <- detectCores (logical=TRUE)
CPUcoreNumber <- CPUthreadNumber/2
# Register the foreach backend
registerDoMC (cores = CPUcoreNumber)
##############################





################################################################################





####################################### FEATURE SELECTION BASED UPON CORRELATION
correlationFeatureSelection <- function (peaklist, correlationMethod="pearson",
	nonFeatures=c("Sample", "Class", "THY")) {
# Take only the part of the matrix without Class and Sample
peaklistFeatures <- peaklist [,!(names(peaklist) %in% nonFeatures)]
# Compute the correlation
featureCorrelation <- cor (peaklistFeatures, method=correlationMethod)
# Output the highly correlated features
highlyCorrelated <- findCorrelation (featureCorrelation, correlationThreshold)
# Keep only the part of the matrix that is not highly correlated
peaklistFeaturesLowCorrelation <- peaklistFeatures [,-highlyCorrelated]
# Sort the dataframe column by a non decreasing order
peaklistFeaturesLowCorrelation <- peaklistFeaturesLowCorrelation [,order(as.numeric(names(peaklistFeaturesLowCorrelation)))]
lowCorrelatedFeatures <- names (peaklistFeaturesLowCorrelation)
# Add back the non Features (Sample and Class)
for (i in 1:length(nonFeatures)) {
	peaklistFeaturesLowCorrelation <- data.frame (peaklistFeaturesLowCorrelation, peaklist[,nonFeatures[i]])
}
names (peaklistFeaturesLowCorrelation) <- c(lowCorrelatedFeatures, nonFeatures)
# Return the dataframe
cat (paste("The number of selected features is:", length(lowCorrelatedFeatures)))
return (peaklistFeaturesLowCorrelation)
}





################################################################################





################################ DEFINED SPECTRA GROUPING (PEAKLIST MATRIX)
# Generate a peaklist dataframe with a certain number of rows per patient
groupPeaklist <- function (peaklist, rowsPerPatient=100, seed=0, balanced=TRUE,
	discardPoorSamples=FALSE, discardIfLowerThan=100, nonFeatures=c("Sample","Class"), algorithm="random") {
# Plant the seed only if a specified value is entered
if (seed != 0) {
	# Make the randomness reproducible
	set.seed (seed)
}
resultDataframe <- data.frame()
# Create the file vector and the patient vector
fileVector <- peaklist$Sample
patientVector <- levels (fileVector)
#### If the dataset has to be balanced...
if (balanced == TRUE) {
	########## All the patients should have the same number of rows, always
	### Determine the shortest patient's spectra number
	# Set the lowest as the row number of the first patient
	lowestObservationNumber <- nrow (peaklist [peaklist$Sample==patientVector[1],])
	# Scroll the other patients
	for (p in 2:length(patientVector)) {
		# Determine the number of observations for this patient
		observationNumber <- nrow (peaklist [peaklist$Sample==patientVector[p],])
		# If its lower than the reference
		if (observationNumber < lowestObservationNumber) {
			lowestObservationNumber <- observationNumber
		}
	}
	#### Rows per patient should not be more than the minimum number of rows
	if (rowsPerPatient >= lowestObservationNumber) {
		rowsPerPatient <- lowestObservationNumber
	}
	if (rowsPerPatient < lowestObservationNumber) {
		rowsPerPatient <- rowsPerPatient
	}
}
# If the spectra per patient is equal to one, do the normal averaging
if (rowsPerPatient == 0 | rowsPerPatient == 1) {
	# For each patient
	for (p in 1:length(patientVector)) {
		# Subset the entire dataframe, isolating only the part corresponding to that patient
		patientDataframe <- peaklist [peaklist$Sample==patientVector[p],]
		patientDataframeRows <- nrow (patientDataframe)
		# Continue only if nothing has to be discarded or the dimension of the patient dataset is higher than the threshold
		if (discardPoorSamples == FALSE || (discardPoorSamples == TRUE && patientDataframeRows >= discardIfLowerThan)) {
		# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
		patientDataframeMatrix <- as.matrix (patientDataframe [,!(names(patientDataframe) %in% nonFeatures)])
		# Generate a single matrix row with the means of each column (feature)
		patientDataframeAVG <- apply (patientDataframeMatrix, 2, mean)
		# Put this row together again with the Sample and the Class back into a dataframe
		patientDataframeFinal <- data.frame (rbind(patientDataframeAVG))
		for (i in 1:length(nonFeatures)) {
			patientDataframeFinal <- data.frame (patientDataframeFinal, patientDataframe[,nonFeatures[i]][1])
		}
		# Fix the column names
		names (patientDataframeFinal) <- names(patientDataframe)
		# Attach this row to the final dataframe
		resultDataframe <- rbind (resultDataframe, patientDataframeFinal)
		}
	}
}
if (rowsPerPatient > 1) {
	# For each patient
	for (p in 1:length(patientVector)) {
		# Subset the entire dataframe, isolating only the part corresponding to that patient
		patientDataframe <- peaklist [peaklist$Sample==patientVector[p],]
		patientDataframeRows <- nrow (patientDataframe)
		# Continue only if nothing has to be discarded or the dimension of the patient dataset is higher than the threshold
		if (discardPoorSamples == FALSE || (discardPoorSamples == TRUE && patientDataframeRows >= discardIfLowerThan)) {
			########################### RANDOMNESS
			if (algorithm == "random") {
				# Generate random folds (one for each representative Spectrum)
				# Index of the spectra to be allocated in the fold, randomly selected
				index <- createFolds (patientDataframe$Sample, k = rowsPerPatient)
				# For each fold
				for (k in 1:length(index)) {
					# Take the corresponding spectra subset (of the selected fold)
					foldIndex <- index[[k]]
					patientDataframeTemp <- patientDataframe [foldIndex,]
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patientDataframeMatrix <- as.matrix (patientDataframeTemp [,!(names(patientDataframeTemp) %in% nonFeatures)])
					# Generate a single matrix row with the means of each column (feature)
					patientDataframeAVG <- apply (patientDataframeMatrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patientDataframeFinal <- data.frame (rbind(patientDataframeAVG))
					for (i in 1:length(nonFeatures)) {
						patientDataframeFinal <- data.frame (patientDataframeFinal, patientDataframe[,nonFeatures[i]][1])
					}
					# Fix the column names
					names (patientDataframeFinal) <- names(patientDataframeTemp)
					# Attach this row to the final dataframe
					resultDataframe <- rbind (resultDataframe, patientDataframeFinal)
				}
			}
			############################ HCA
			if (algorithm == "hierarchicalClustering" | algorithm == "hca" | algorithm == "HCA") {
				# Compute the distance matrix
				distanceMatrix <- dist (patientDataframe[,!(names(patientDataframe) %in% nonFeatures)], method="euclidean")
				hca <- hclust (distanceMatrix)
				# Associate to each row/spectrum the subgroup of the tree to which it belongs
				# k must be between 1 and 56
				if (rowsPerPatient > 56) {
					rowsPerPatient <- 56
				}
				hcaGroups <- cutree (hca, k=rowsPerPatient)
				# For each subgroup to be isolated...
				for (s in 1:rowsPerPatient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which (hcaGroups == s)
					patientDataframeTemp <- patientDataframe [index,]
					# Average the spectra in this HCA subgroup
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patientDataframeMatrix <- as.matrix (patientDataframeTemp [,!(names(patientDataframeTemp) %in% nonFeatures)])
					# Generate a single matrix row with the means of each column (feature)
					patientDataframeAVG <- apply (patientDataframeMatrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patientDataframeFinal <- data.frame (rbind(patientDataframeAVG))
					for (i in 1:length(nonFeatures)) {
						patientDataframeFinal <- data.frame (patientDataframeFinal, patientDataframe[,nonFeatures[i]][1])
					}
					# Fix the column names
					names (patientDataframeFinal) <- names(patientDataframeTemp)
					# Attach this row to the final dataframe
					resultDataframe <- rbind (resultDataframe, patientDataframeFinal)
				}
			}
			######################### K-MEANS
			if (algorithm == "kMeans" | algorithm == "kmeans" | algorithm == "k-Means") {
				# Compute the k-Means clustering
				kMeans <- kmeans (patientDataframe[,!(names(patientDataframe) %in% nonFeatures)], centers=rowsPerPatient)
				# Associate to each row/spectrum the subgroup/cluster to which it belongs
				kMeansGroups <- kMeans$cluster
				# For each subgroup to be isolated...
				for (s in 1:rowsPerPatient) {
					# Index the spectra under in the selected subgroup of the HCA
					index <- which (kMeansGroups == s)
					patientDataframeTemp <- patientDataframe [index,]
					# Average the spectra in this kMeans subgroup
					# Transform it into a matrix for "apply" (only the features, excluding the Sample and the Class)
					patientDataframeMatrix <- as.matrix (patientDataframeTemp [,!(names(patientDataframeTemp) %in% nonFeatures)])
					# Generate a single matrix row with the means of each column (feature)
					patientDataframeAVG <- apply (patientDataframeMatrix, 2, mean)
					# Put this row together again with the Sample and the Class back into a dataframe
					patientDataframeFinal <- data.frame (rbind(patientDataframeAVG))
					for (i in 1:length(nonFeatures)) {
						patientDataframeFinal <- data.frame (patientDataframeFinal, patientDataframe[,nonFeatures[i]][1])
					}
					# Fix the column names
					names (patientDataframeFinal) <- names(patientDataframeTemp)
					# Attach this row to the final dataframe
					resultDataframe <- rbind (resultDataframe, patientDataframeFinal)
				}
			}
		}
	}
}
return (resultDataframe)
}





################################################################################





################################################## NO HEMOGLOBIN IN THE PEAKLIST
noHemoglobinPeaklist <- function (peaklist, hemoglobinMass=15400, tolerancePPM=2000,
	nonFeatures=c("Sample", "Class")) {
# Identify the monocharged and the bicharged ion
hemoglobin <- c(hemoglobinMass, hemoglobinMass/2)
# Extract the feature names (as numeric)
featureList <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
# Create the index for the hemoglobin mass(es)
index <- integer()
# Scroll the feature list
for (f in 1:length(featureList)) {
	for (h in 1:length(hemoglobin)) {
		# Identify the hemoglobin
		if (abs(hemoglobin[h] - featureList[f])*10^6/hemoglobinMass <= tolerancePPM) {
			index <- append (index, f)
		}
	}
}
# Keep only the part of the dataframe without the hemoglobin column(s)
if (length(index) > 0) {
	peaklistNoHemoglobin <- peaklist [,-index]
}	else {peaklistNoHemoglobin <- peaklist}
return (peaklistNoHemoglobin)
}





################################################################################





##################### MATRIX/DATAFRAME SPLITTING FUNCTION (TRAINING AND TESTING)
# Output which patient (for each class) should be taken for training and for testing
matrixSplittingTrainingTesting <- function (peaklist, classList=list(), seed=0,
	percentageOfObservationForTraining=50) {
# Plant the seed only if a specified value is entered
if (seed != 0) {
	# Make the randomness reproducible
	set.seed (seed)
}
############################# If there are no classes or one class, just randomly select rows on the entire dataset
if (length(classList) <= 1) {
indexTraining <- sample (nrow(peaklist), size=round(nrow(peaklist)*percentageOfObservationForTraining/100))
}
############################# If there are more classes, randomly select the rows for each class
if (length(classList) > 1) {
classList <- sort (classList)
# Automatically generate the corresponding class if not specified in the input matrix
########### Create the dataframe with only patients and classes
# Create two dataframe columns with the patients
sampleAndClass <- data.frame (unique(peaklist$Sample), unique(peaklist$Sample))
# Give the correct names to the columns
colnames (sampleAndClass) <- c ("Sample", "Class")
# Define the sample column as character (not as factor)
sampleAndClass$Sample <- as.character (sampleAndClass$Sample)
sampleAndClass$Class <- as.character (sampleAndClass$Class)
# Replace the patients in the second column with the class
for (j in 1:length(classList)) {
	for (i in 1:length(sampleAndClass$Class)) {
		if (length(grep(classList[j], sampleAndClass$Class[i], ignore.case=TRUE)) != 0) {
			sampleAndClass$Class[i] <- classList [j]
		}
	}
}
# Create a list of dataframes, one for each class
classDataFrameList <- list()
for (j in 1:length(classList)) {
	classDataFrameList[[j]] <- sampleAndClass [sampleAndClass$Class == classList[j],]
}
############ Create a list containing the training and testing indexes of the class dataframes
indexTraining <- list()
for (i in 1:length(classDataFrameList)) {
	indexTraining[[i]] <- createDataPartition (y = classDataFrameList[[i]]$Class, p = (percentageOfObservationForTraining/100), list = FALSE)
}
names (indexTraining) <- classList
}
# Now we have randomly selected some patients for the training and some others for the testing, for each class
####################
# The indexTraining is integer if there are no classes or one class, and rows are selected randomly
if (is.integer(indexTraining)) {
	trainingDataset <- peaklist [indexTraining,]
	testingDataset <- peaklist [-indexTraining,]
}
# The indexTraining is a list if there are more than one classes (one element per class),
# and rows are selected randomly for each class
if (is.list(indexTraining)) {
##### Create the training and the testing patient dataframes
trainingDataFrameList <- list()
testingDataFrameList <- list()
# For each class (and so for each element of the list with the indexes per class)
for (i in 1:length(indexTraining)) {
	# Create the two complementary dataframes: randomly select the patients for each class
	# for training and testing
	trainingDataFrameList[[i]] <- classDataFrameList[[i]][indexTraining[[i]],]
	testingDataFrameList[[i]] <- classDataFrameList[[i]][-indexTraining[[i]],]
}
#### In the original matrix, select only the rows of the selected samples
# Turn the patient lists into a data frame and store all the sample values in a vector
patientsForTraining <- character ()
for (i in 1:length(trainingDataFrameList)) {
	patientsForTraining <- append (patientsForTraining, trainingDataFrameList[[i]]$Sample)
}
patientsForTesting <- character ()
for (i in 1:length(testingDataFrameList)) {
	patientsForTesting <- append (patientsForTesting, testingDataFrameList[[i]]$Sample)
}
# Create the two datasets
trainingDataset <- data.frame()
testingDataset <- data.frame()
for (p in 1:length(patientsForTraining)) {
	trainingDataset <- rbind (trainingDataset, peaklist[(peaklist$Sample==patientsForTraining[p]),])
}
for (p in 1:length(patientsForTesting)) {
	testingDataset <- rbind (testingDataset, peaklist[(peaklist$Sample==patientsForTesting[p]),])
}
}
return (list(trainingDataset = trainingDataset, testingDataset = testingDataset))
}





################################################################################





############################## RECURSIVE FEATURE ELIMINATION - FEATURE SELECTION
recursiveFeatureElimination <- function (peaklist, featuresToSelect=20, selectionMethod="pls",
	discriminant="Class", nonFeatures=c("Sample", "Class", "THY"), seed=0) {
if (seed != 0) {
	set.seed (seed)
}
########################################################### RFE MODEL (with PLS)
rfeCtrl <- rfeControl (functions = caretFuncs,
	method = "cv",
	repeats = 3,
	number = 10
	)
# The simulation will fit models with subset sizes: (the subset size is the number of predictors to use)
if (selectionMethod != "") {
	if (seed != 0) {
		set.seed (seed)
	}
	RFEmodel <- rfe (x = peaklist [,!(names(peaklist) %in% nonFeatures)], y = peaklist[,discriminant],
		sizes = featuresToSelect,
		rfeControl = rfeCtrl,
		method = selectionMethod,
		)
} else {
	if (seed != 0) {
		set.seed (seed)
	}
	RFEmodel <- rfe (x = peaklist [,!(names(peaklist) %in% nonFeatures)], y = peaklist[,discriminant],
		sizes = featuresToSelect,
		rfeControl = rfeCtrl,
		)
}
# Output the best predictors after the RFE
predictorsRFE <- predictors (RFEmodel) [1:featuresToSelect]
# Take the selected features
peaklistRFE <- peaklist [,predictorsRFE]
# Add the non features back
for (i in 1:length(nonFeatures)) {
	peaklistRFE <- cbind (peaklistRFE , peaklist[,nonFeatures[i]])
}
names (peaklistRFE) <- c(as.character (predictorsRFE), nonFeatures)
#
return (peaklistRFE)
}





################################################################################





############################################### TRUNCATE PEAKLIST MATRIX
truncatePeaklist <- function (peaklist, range=c(4000, 15000), nonFeatures=c("Sample", "Class")){
featuresToKeep <- numeric()
# Do not do anything if the values are set to zero
if (range[1] == 0 & range[2] == 0) {
	# Keep all the features
	featuresToKeep <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
}
#
if (range[1] != 0 & range[2] == 0) {
	featureList <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
	# Keep only the features that are above the lower threshold
	for (f in 1:length(featureList)) {
		if (featureList[f] >= range[1]) {
			featuresToKeep <- append (featuresToKeep, featureList[f])
		}
	}
}
#
if (range[1] != 0 & range[2] != 0) {
	featureList <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
	# Keep only the features that are above the lower threshold and below the upper threshold
	for (f in 1:length(featureList)) {
		if (featureList[f] >= range[1] && featureList[f] <= range[2]) {
			featuresToKeep <- append (featuresToKeep, featureList[f])
		}
	}
}
# Keep only the features of interest
peaklistTruncated <- data.frame (peaklist [,as.character(featuresToKeep)])
# Add the non features back
for (i in 1:length(nonFeatures)) {
	peaklistTruncated <- data.frame (peaklistTruncated, peaklist[,nonFeatures[i]])
}
names (peaklistTruncated) <- c(featuresToKeep, nonFeatures)
#
return (peaklistTruncated)
}





################################################################################





##################### CROSS-VALIDATION FOR A SUPPORT VECTOR MACHINE MODEL
################ One row per patient (average)
crossValidationSVM <- function (trainingDataset, seed=0, kFold=10, repeats=5, svmModel, nonFeatures=c("Sample","Class","THY"), positiveClass=levels(trainingDataset$Class)[1]) {
if (seed != 0) {
	set.seed (seed)
}
### Extract the svm model parameters
# Kernel
svmKernel <- svmModel$kernel
if (svmKernel == 0) {
	svmKernel <- "linear"
}
if (svmKernel == 1) {
	svmKernel <- "polynomial"
}
if (svmKernel == 2) {
	svmKernel <- "radial"
}
if (svmKernel == 3) {
	svmKernel <- "sigmoid"
}
svmCost <- svmModel$cost
svmDegree <- svmModel$degree
svmGamma <- svmModel$gamma
### Output matrix
confusionMatrixList <- list()
resultMatrix <- matrix (NA, ncol=14, nrow=0)
colnames (resultMatrix) <- c("Accuracy","Kappa","No information rate", "Accuracy p-value",
	"McNiemar p-value","Sensitivity","Specificity","PPV","NPV","Prevalence","Detection Rate","Detection Prevalence","Balanced Accuracy","ROC AUC")
#### For each repetition...
for (i in 1:repeats) {
	if (seed != 0) {
		set.seed (seed)
	}
	# Index of the spectra to be allocated in the folds, randomly selected
	# k equal-size partitions are created
	index <- createFolds (trainingDataset$Class, k = kFold)
	# Each single fold has to be used as a training and the rest as testing
	for (k in 1:length(index)) {
		# Generate the training and the testing subsets (based upon the index of the kFold)
		# The kFold indexed with k is used each time as testing
		testingSubset <- trainingDataset [index[[k]],]
		testingPredictors <- testingSubset [,!(names(testingSubset) %in% nonFeatures)]
		testingOutcomes <- testingSubset$Class
		# All the other ones are used as training
		trainingSubset <- trainingDataset [-index[[k]],]
		trainingPredictors <- trainingSubset [,!(names(trainingSubset) %in% nonFeatures)]
		trainingOutcomes <- trainingSubset$Class
		# Now the model has to be built using the training and tested onto the testing dataset
		if (seed != 0) {
			set.seed (seed)
		}
		model <- svm (x=trainingPredictors, y=trainingOutcomes, scale=TRUE, kernel=svmKernel,
			degree=svmDegree, gamma=svmGamma, cost=svmCost)
		# Use the model to predict the testing subset
		if (seed != 0) {
			set.seed (seed)
		}
		predictedClasses <- predict (model, newdata=testingPredictors)
		# Compute the performances
		performances <- confusionMatrix (data=predictedClasses, testingOutcomes, positive=positiveClass)
		confusionMatrixList [[k*i]] <- performances
		#### ROC analysis
		SVMroc <- roc (response=testingOutcomes, predictor=as.numeric(predictedClasses))
		# Create a matrix row with the partial results
		resultMatrixRow <- matrix (0, ncol=ncol(resultMatrix), nrow=1)
		resultMatrixRow [1,1] <- as.numeric (performances$overall[1])
		resultMatrixRow [1,2] <- as.numeric (performances$overall[2])
		resultMatrixRow [1,3] <- as.numeric (performances$overall[5])
		resultMatrixRow [1,4] <- as.numeric (performances$overall[6])
		resultMatrixRow [1,5] <- as.numeric (performances$overall[7])
		resultMatrixRow [1,6] <- as.numeric (performances$byClass[1])
		resultMatrixRow [1,7] <- as.numeric (performances$byClass[2])
		resultMatrixRow [1,8] <- as.numeric (performances$byClass[3])
		resultMatrixRow [1,9] <- as.numeric (performances$byClass[4])
		resultMatrixRow [1,10] <- as.numeric (performances$byClass[5])
		resultMatrixRow [1,11] <- as.numeric (performances$byClass[6])
		resultMatrixRow [1,12] <- as.numeric (performances$byClass[7])
		resultMatrixRow [1,13] <- as.numeric (performances$byClass[8])
		resultMatrixRow [1,14] <- as.numeric (SVMroc$auc)
		# Attach the row to the result matrix
		resultMatrix <- rbind (resultMatrix, resultMatrixRow)
	}
}
# Compute the average of each column
resultMatrixAVG <- apply (resultMatrix, 2, mean, na.rm=TRUE)
return (list (confusionMatrix=confusionMatrixList, resultMatrix=resultMatrix, resultMatrixAVG=resultMatrixAVG, seed=seed, Features=names(trainingDataset)[!(names(trainingDataset) %in% nonFeatures)]))
}





################################################################################





################################################## DEISOTOPING IN PEAKLIST
deisotopePeaklist <- function (peaklist, nonFeatures=c("Sample", "Class")) {
ions <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
# Deisotoped Ions
deisotopedIons <- numeric()
# Take each ion and evaluate the consecutive ions
isotopeBin <- numeric ()
# Scroll the ion list
for (i in 1:(length(ions)-1)) {
	# Create the isotope bin, if the consecutive peaks belong to a isotope cluster
	if (abs(ions[i+1] - ions[i]) < 1.5) {
		isotopeBin <- append (isotopeBin, ions[i])
	}
	# When a major distance is found... We are out of the cluster
	if (abs(ions[i+1] - ions[i]) >= 1.5) {
		# Add this mass to the isotope bin (it will be the last belonging to a this isotope cluster
		isotopeBin <- append (isotopeBin, ions[i])
		# Store the first isotope of the isotope Bin created so far in a new vector
		deisotopedIons <- append (deisotopedIons, isotopeBin[1])
		# Empty the isotope Bin in order for it to become available for the future cycle
		isotopeBin <- numeric()

	}
}
# Keep only the part of the matrix with the selected ions
deisotopedPeaklist <- peaklist [,as.character(deisotopedIons)]
# Add the non features back
columnNames <- names (deisotopedPeaklist)
for (i in 1:length(nonFeatures)) {
	deisotopedPeaklist <- cbind (deisotopedPeaklist, peaklist[,nonFeatures[i]])
}
names (deisotopedPeaklist) <- c(columnNames, nonFeatures)
return (deisotopedPeaklist)
}





################################################################################





################################################## NO TRYPSIN IN THE PEAKLIST
noTrypsinPeaklist <- function (peaklist, trypsinMass=c(841.50, 905.50, 1005.48, 1044.56, 1468.72, 1735.84, 1767.79, 2157.02, 2210.10, 2282.17, 3012.32, 4488.11),
	tolerancePPM=200, nonFeatures=c("Sample", "Class")) {
# Extract the feature names (as numeric)
featureList <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
# Create the index for the hemoglobin mass(es)
index <- integer()
# Scroll the feature list
for (f in 1:length(featureList)) {
	for (t in 1:length(trypsinMass)) {
		# Identify the hemoglobin
		if (abs(trypsinMass[t] - featureList[f])*10^6/trypsinMass[t] <= tolerancePPM) {
			index <- append (index, f)
		}
	}
}
# Keep only the part of the dataframe without the hemoglobin column(s)
if (length(index) > 0) {
	peaklistNoTrypsin <- peaklist [,-index]
}	else {peaklistNoTrypsin <- peaklist}
return (peaklistNoTrypsin)
}





################################################################################





############################# SVM TUNING AND TESTING
svmTuningAndClassification <- function (peaklistTraining, peaklistTesting=NULL,
	nonFeatures=c("Sample","Class","THY"),	gamma=10^(-5:5), cost=10^(-5:5),
	epsilon=seq(1,2,by=1), degree=1:5, kernel="radial", kFoldCV=2, repeatsCV=5,
	positiveClassCV="HP", seed=0, PCA=FALSE, numerOfComponents=3) {
if (seed != 0) {
	set.seed (seed)
}
############################################################################ PCA
if (PCA == TRUE) {
# Compute the PCs
pcaTraining <- prcomp (as.matrix(peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)]))
pcaTesting <- prcomp (as.matrix(peaklistTesting [,!(names(peaklistTesting) %in% nonFeatures)]))
#################### SVM tuning
if (seed != 0) {
	set.seed (seed)
}
SVMtuning <- tune.svm (as.data.frame(pcaTraining$x[,1:numberOfComponents]), factor(peaklistTraining$Class),
	gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
# Select automatically the best model after the tuning
SVMModel <- SVMtuning$best.model
# Parameters output
svmKernel <- SVMtuning$best.model$kernel
if (svmKernel == 0) {
	svmKernel <- "linear"
}
if (svmKernel == 1) {
	svmKernel <- "polynomial"
}
if (svmKernel == 2) {
	svmKernel <- "radial"
}
if (svmKernel == 3) {
	svmKernel <- "sigmoid"
}
parametersOutput <- list (Kernel=svmKernel, Cost=SVMtuning$best.model$cost, Degree=SVMtuning$best.model$degree, Epsilon=SVMtuning$best.model$epsilon, Gamma=SVMtuning$best.model$gamma)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix (peaklistTesting) || is.data.frame (peaklistTesting)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed (seed)
	}
	predictedClassesSVM <- predict (SVMModel, newdata = as.data.frame(pcaTesting$x[,1:numberOfComponents]))
	# Create the outcomes dataframe
	classificationResultsSVM <- data.frame (Sample=peaklistTesting$Sample, Predicted=predictedClassesSVM, True=peaklistTesting$Class)
	### Generate the confusion matrix to evaluate the performances
	testPerformancesSVM <- confusionMatrix (data = predictedClassesSVM, peaklistTesting$Class, positive=positiveClassCV)
	#### ROC analysis
	SVMroc <- list()
	rocCurve <- roc (response=classificationResultsSVM$True, predictor=as.numeric(classificationResultsSVM$Predicted))
	SVMroc[[1]] <- rocCurve$auc
	plot (rocCurve)
	rocLegend <- paste ("ROC area under the curve:", rocCurve$auc)
	legend ("bottomright", legend=rocLegend, xjust=0.5, yjust=0.5)
	SVMroc[[2]] <- recordPlot()
	# Output the results
	return (list (model=SVMModel, svmParameters=parametersOutput, classificationResults=classificationResultsSVM, performances=testPerformancesSVM, ROC=SVMroc))
}	else {return (list (model=SVMModel, svmParameters=parametersOutput, crossValidation=cvSVMModel))}
} else {
####################################################################### FEATURES
################### Optimise according to the kernel
# Find the best tuning parameters for the SVM
if (kernel == "radial") {
	if (seed != 0) {
		set.seed (seed)
	}
	SVMtuning <- tune.svm (peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)], factor(peaklistTraining$Class),
		gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
}
if (kernel == "polynomial") {
	if (seed != 0) {
		set.seed (seed)
	}
	SVMtuning <- tune.svm (peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)], factor(peaklistTraining$Class),
		gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon, degree=degree)
}
if (kernel == "linear") {
	if (seed != 0) {
		set.seed (seed)
	}
	SVMtuning <- tune.svm (peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)], factor(peaklistTraining$Class),
		gamma=gamma, cost=cost, kernel=kernel)
}
if (kernel == "sigmoid") {
	if (seed != 0) {
		set.seed (seed)
	}
	SVMtuning <- tune.svm (peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)], factor(peaklistTraining$Class),
		gamma=gamma, cost=cost, kernel=kernel, epsilon=epsilon)
}
# Select automatically the best model after the tuning
SVMModel <- SVMtuning$best.model
# Parameters output
svmKernel <- SVMtuning$best.model$kernel
if (svmKernel == 0) {
	svmKernel <- "linear"
}
if (svmKernel == 1) {
	svmKernel <- "polynomial"
}
if (svmKernel == 2) {
	svmKernel <- "radial"
}
if (svmKernel == 3) {
	svmKernel <- "sigmoid"
}
parametersOutput <- list (Kernel=svmKernel, Cost=SVMtuning$best.model$cost, Degree=SVMtuning$best.model$degree, Epsilon=SVMtuning$best.model$epsilon, Gamma=SVMtuning$best.model$gamma)
#################### CROSS-VALIDATION
cvSVMModel <- crossValidationSVM (peaklistTraining, kFold=kFoldCV, repeats=repeatsCV, svmModel=SVMModel,
	nonFeatures=nonFeatures, positiveClass=positiveClassCV, seed=seed)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix (peaklistTesting) || is.data.frame (peaklistTesting)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed (seed)
	}
	predictedClassesSVM <- predict (SVMModel, newdata = peaklistTesting [,!(names(peaklistTesting) %in% nonFeatures)])
	# Create the outcomes dataframe
	classificationResultsSVM <- data.frame (Sample=peaklistTesting$Sample, Predicted=predictedClassesSVM, True=peaklistTesting$Class)
	### Generate the confusion matrix to evaluate the performances
	testPerformancesSVM <- confusionMatrix (data = predictedClassesSVM, peaklistTesting$Class, positive=positiveClassCV)
	#### ROC analysis
	SVMroc <- list()
	rocCurve <- roc (response=classificationResultsSVM$True, predictor=as.numeric(classificationResultsSVM$Predicted))
	SVMroc[[1]] <- rocCurve$auc
	plot (rocCurve)
	rocLegend <- paste ("ROC area under the curve:", rocCurve$auc)
	legend ("bottomright", legend=rocLegend, xjust=0.5, yjust=0.5)
	SVMroc[[2]] <- recordPlot()
	# Output the results
	return (list (model=SVMModel, crossValidation=cvSVMModel, classificationResults=classificationResultsSVM, parameters=parametersOutput, performances=testPerformancesSVM, ROC=SVMroc))
}	else {return (list (model=SVMModel, parameters=parametersOutput, crossValidation=cvSVMModel))}
	}
}





################################################################################





############################# SVM TUNING AND TESTING
svmClassification <- function (peaklistTraining, peaklistTesting=NULL,
	nonFeatures=c("Sample","Class","THY"),	gamma=0.1, cost=10,
	epsilon=0.1, degree=3, kernel="radial", kFoldCV=10, repeatsCV=2,
	positiveClassCV="HP", seed=0, PCA=FALSE, numerOfComponents=3) {
if (seed != 0) {
	set.seed (seed)
}
############################################################################ PCA
if (PCA == TRUE) {
# Compute the PCs
pcaTraining <- prcomp (as.matrix(peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)]))
pcaTesting <- prcomp (as.matrix(peaklistTesting [,!(names(peaklistTesting) %in% nonFeatures)]))
# SVM tuning
if (seed != 0) {
	set.seed (seed)
}
# SVM with defined parameters
SVMModel <- svm (pcaTraining, factor(peaklistTraining$Class),
		kernel=kernel, cost=cost, epsilon=epsilon, gamma=gamma)
parametersOutput <- list (Kernel=kernel, Cost=cost, Degree=degree, Epsilon=epsilon, Gamma=gamma)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix (peaklistTesting) || is.data.frame (peaklistTesting)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed (seed)
	}
	predictedClassesSVM <- predict (SVMModel, newdata = as.data.frame(pcaTesting$x[,1:numberOfComponents]))
	# Create the outcomes dataframe
	classificationResultsSVM <- data.frame (Sample=peaklistTesting$Sample, Predicted=predictedClassesSVM, True=peaklistTesting$Class)
	### Generate the confusion matrix to evaluate the performances
	testPerformancesSVM <- confusionMatrix (data = predictedClassesSVM, peaklistTesting$Class, positive=positiveClassCV)
	#### ROC analysis
	SVMroc <- list()
	rocCurve <- roc (response=classificationResultsSVM$True, predictor=as.numeric(classificationResultsSVM$Predicted))
	SVMroc[[1]] <- rocCurve$auc
	plot (rocCurve)
	rocLegend <- paste ("ROC area under the curve:", rocCurve$auc)
	legend ("bottomright", legend=rocLegend, xjust=0.5, yjust=0.5)
	SVMroc[[2]] <- recordPlot()
	# Output the results
	return (list (model=SVMModel, svmParameters=parametersOutput, classificationResults=classificationResultsSVM, performances=testPerformancesSVM, ROC=SVMroc))
}	else {return (list (model=SVMModel, svmParameters=parametersOutput, crossValidation=cvSVMModel))}
} else {
####################################################################### FEATURES
# SVM with defined parameters
SVMModel <- svm (peaklistTraining [,!(names(peaklistTraining) %in% nonFeatures)], factor(peaklistTraining$Class),
		kernel=kernel, cost=cost, epsilon=epsilon, gamma=gamma)
parametersOutput <- list (Kernel=kernel, Cost=cost, Degree=degree, Epsilon=epsilon, Gamma=gamma)
#################### CROSS-VALIDATION
cvSVMModel <- crossValidationSVM (peaklistTraining, kFold=kFoldCV, repeats=repeatsCV, svmModel=SVMModel,
	nonFeatures=nonFeatures, positiveClass=positiveClassCV, seed=seed)
#################### EXTERNAL VALIDATION (If a dataset is provided)
if (is.matrix (peaklistTesting) || is.data.frame (peaklistTesting)) {
	#### Use the model to predict the outcome of the testing set (the new data must have only the predictors)
	if (seed != 0) {
		set.seed (seed)
	}
	predictedClassesSVM <- predict (SVMModel, newdata = peaklistTesting [,!(names(peaklistTesting) %in% nonFeatures)])
	# Create the outcomes dataframe
	classificationResultsSVM <- data.frame (Sample=peaklistTesting$Sample, Predicted=predictedClassesSVM, True=peaklistTesting$Class)
	### Generate the confusion matrix to evaluate the performances
	testPerformancesSVM <- confusionMatrix (data = predictedClassesSVM, peaklistTesting$Class, positive=positiveClassCV)
	#### ROC analysis
	SVMroc <- list()
	rocCurve <- roc (response=classificationResultsSVM$True, predictor=as.numeric(classificationResultsSVM$Predicted))
	SVMroc[[1]] <- rocCurve$auc
	plot (rocCurve)
	rocLegend <- paste ("ROC area under the curve:", rocCurve$auc)
	legend ("bottomright", legend=rocLegend, xjust=0.5, yjust=0.5)
	SVMroc[[2]] <- recordPlot()
	# Output the results
	return (list (model=SVMModel, crossValidation=cvSVMModel, classificationResults=classificationResultsSVM, parameters=parametersOutput, performances=testPerformancesSVM, ROC=SVMroc))
}	else {return (list (model=SVMModel, parameters=parametersOutput, crossValidation=cvSVMModel))}
	}
}





################################################################################





############################################ OUTLIERS REMOVAL
outliersRemoval <- function (vector, replaceWith="") {
summaryVector <- summary (vector)
# Calculate the interquartile range
interQuartileRange <- summaryVector[5] - summaryVector[2]
# Calculate the fences, beyond which the spectrum is an outlier
IQRfences <- c((summaryVector[2] - 1.5*interQuartileRange), (summaryVector[5] + 1.5*interQuartileRange))
# Find the outliers based on the fences condition
outliersPosition <- which (vector < IQRfences[1] | vector > IQRfences[2])
# If the outliers have to be discarded...
if (replaceWith == "") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vector <- vector [-outliersPosition]
	}
}
if (replaceWith == "mean") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vectorNoOutliers <- vector [-outliersPosition]
		# Replace the outliers with the vector mean (without outliers)
		vector [outliersPosition] <- mean (vectorNoOutliers)
	}
}
if (replaceWith == "median") {
	if (length(outliersPosition) > 0) {
		# Remove the correspondent spectra from the dataset
		vectorNoOutliers <- vector [-outliersPosition]
		# Replace the outliers with the vector median
		vector [outliersPosition] <- median (vectorNoOutliers)
	}
}
return (list (vector=vector, outliersPosition=outliersPosition))
}





################################################################################





############################## PEAK STATISTICS (on the final peaklist)
peakStatisticsPeaklist <- function (peaklist, nonFeatures=c("Sample","Class","THY"), removeOutliers=TRUE, replaceOutliersWith="") {
classList <- levels (peaklist$Class)
# Determine the number of classes
if (length(classList) == 0 | length(classList) == 1) {
numberOfClasses <- 1
}
if (length(classList) > 1) {
	numberOfClasses <- length(classList)
}
# Peak vector
#peakVector <- as.numeric (names(peaklist))
############################################################## ONE CLASS
if (numberOfClasses == 1) {
	# Output matrix
	peakStatMatrix <- matrix (0, nrow=(ncol(peaklist)-length(nonFeatures)), ncol=8)
	rownames (peakStatMatrix) <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
	colnames (peakStatMatrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter")
	# For each peak
	for (p in 1:(ncol(peaklist)-length(nonFeatures))) {
		intensityVector <- as.numeric (peaklist[,p])
		if (removeOutliers == TRUE) {
			intensityVector <- outliersRemoval (intensityVector, replaceWith=replaceOutliersWith)$vector
		}
		# Calculate the statistical parameters on the intensity values in the vector
		# Normality
		if (length(intensityVector) >= 3 & length(intensityVector) <= 5000) {
		shapiroTest <- shapiro.test (intensityVector)
			if (shapiroTest$p.value < 0.05) {
			distributionType <- "Non-normal"
			}
			if (shapiroTest$p.value >= 0.05) {
			distributionType <- "Normal"
			}
		}
		if (length(intensityVector) < 3) {
		distributionType <- "Not determinable, number of samples too low"
		}
		if (length(intensityVector) > 5000) {
		distributionType <- "Number of samples too high, assume it is normal"
		}
		# Other parameters
		#spectraNames <- levels (peaklist$Sample)
		stDevIntensity <- sd (intensityVector)
		summaryIntensityVector <- summary (intensityVector)
		meanIntensity <- summaryIntensityVector [4]
		coeffVariation <- (stDevIntensity / meanIntensity) *100
		medianIntensity <- summaryIntensityVector [3]
		firstQuartile <- summaryIntensityVector [2]
		thirdQuartile <- summaryIntensityVector [5]
		interQuartileRange <- thirdQuartile - firstQuartile
		spectraCounter <- length (intensityVector)
		# Fill the matrix with the values
		peakStatMatrix [p,1] <- distributionType
		peakStatMatrix [p,2] <- meanIntensity
		peakStatMatrix [p,3] <- stDevIntensity
		peakStatMatrix [p,4] <- coeffVariation
		peakStatMatrix [p,5] <- medianIntensity
		peakStatMatrix [p,6] <- interQuartileRange
		peakStatMatrix [p,7] <- paste ("1st quartile", firstQuartile, "; 3rd quartile", thirdQuartile)
		peakStatMatrix [p,8] <- spectraCounter
		#peakStatMatrix [p,9] <- spectraNames
	}
}
############################################################ TWO OR MORE CLASSES
# Every variable now is a list, each element of which corresponds to a certain value from a class
# So every variable is a list with the same length of the class list (each element of the list
# is referred to a class
if (numberOfClasses > 1) {
	# Output matrix
	peakStatMatrix <- matrix (0, nrow=(ncol(peaklist)-length(nonFeatures)), ncol=14)
	rownames (peakStatMatrix) <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
	colnames (peakStatMatrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %",
		"Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)",
		"t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
	# For each peak
	for (p in 1:(ncol(peaklist)-length(nonFeatures))) {
		# Put the intensity of that peak into one vector per class (in a global list)
		intensityVector <- list()
		# Scroll the peaklists and Add the peak intensity to a vector (one for each class)
		for (l in 1:length(classList)) {
			# Allocate in the intensity vector the rows for that peak belonging to the certain class
			intensityVector [[l]] <- peaklist [peaklist$Class == classList[l],p]
		}
		if (removeOutliers == TRUE) {
			for (i in 1:length(intensityVector)) {
				intensityVector[[l]] <- outliersRemoval (intensityVector[[l]], replaceWith=replaceOutliersWith)$vector
			}
		}
		######################## STATISTICAL PARAMETERS
		############################################### Normality for each class
		shapiroTest <- list()
		distributionType <- list()
		for (l in 1:length(classList)) {
			if (length(intensityVector[[l]]) >= 3 && length(intensityVector[[l]]) <=5000) {
				shapiroTest[[l]] <- shapiro.test (intensityVector[[l]])
				if (shapiroTest[[l]]$p.value < 0.05) {
				distributionType[[l]] <- "Non-normal"
				}
				if (shapiroTest[[l]]$p.value >= 0.05) {
				distributionType[[l]] <- "Normal"
				}
			}
			if (length(intensityVector[[l]]) < 3) {
			distributionType[[l]] <- "Not determinable, number of samples too low"
			}
			if (length(intensityVector) > 5000) {
			distributionType[[l]] <- "Number of samples too high, assume it is normal"
			}
		}
		##################################################### Homoscedasticity
		if (length(classList) == 2) {
			varianceTestPar <- var.test (intensityVector[[1]], intensityVector[[2]])
		}
		if (length(classList) >= 2) {
			varianceTestNonPar <- bartlett.test (as.numeric(peaklist[,p]), g=peaklist$Class)
		}
		########################################### Other parameters (per class)
		stDevIntensity <- list()
		summaryIntensityVector <- list()
		meanIntensity <- list()
		coeffVariation <- list()
		medianIntensity <- list()
		firstQuartile <- list()
		thirdQuartile <- list()
		interQuartileRange <- list()
		spectraCounter <- list()
		variance <- list()
		for (l in 1:length(classList)) {
			stDevIntensity[[l]] <- sd (intensityVector[[l]])
			summaryIntensityVector [[l]] <- summary (intensityVector[[l]])
			meanIntensity[[l]] <- summaryIntensityVector[[l]] [4]
			coeffVariation[[l]] <- (stDevIntensity[[l]] / meanIntensity[[l]]) *100
			medianIntensity[[l]] <- summaryIntensityVector[[l]] [3]
			firstQuartile[[l]] <- summaryIntensityVector[[l]] [2]
			thirdQuartile[[l]] <- summaryIntensityVector[[l]] [5]
			interQuartileRange[[l]] <- thirdQuartile[[l]] - firstQuartile[[l]]
			spectraCounter[[l]] <- length (intensityVector[[l]])
			variance[[l]] <- var (intensityVector[[l]])
		}
		############################################# Parameters between classes
		# T-test (Two classes, parametric)
		if (length(classList) == 2) {
			tTest <- t.test (intensityVector[[1]], intensityVector[[2]])
		}
		# ANOVA TEST (More than two classes, parametric)
		if (length(classList) >= 2) {
		anovaTest <- aov (peaklist[,p] ~ peaklist$Class)
		}
		# WILCOXON - MANN-WHITNEY TEST (Two classes, non parametric)
		if (length(classList) == 2) {
			WilcoxonTest <- wilcox.test (intensityVector[[1]], intensityVector[[2]])
		}
		# KRUSKAL-WALLIS TEST (more than two classes, non parametric)
		if (length(classList) >= 2) {
			KruskalWallisTest <- kruskal.test (peaklist[,p], g=peaklist$Class)
		}
		######################################## Fill the matrix with the values
		# Distribution Type
		distributionTypeName <- character()
		for (l in length(classList)) {
			distributionTypeName <- paste (distributionTypeName, " ", distributionType[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,1] <- paste (distributionTypeName)
		# Mean
		meanIntensityName <- character()
		for (l in length(classList)) {
			meanIntensityName <- paste (meanIntensityName, " ", meanIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,2] <- meanIntensityName
		# Standard Deviation
		StDevName <- character()
		for (l in length(classList)) {
			StDevName <- paste (StDevName, " ", stDevIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,3] <- StDevName
		# Coefficient of Variation
		coeffVariationName <- character()
		for (l in length(classList)) {
			coeffVariationName <- paste (coeffVariationName, " ", coeffVariation[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,4] <- coeffVariationName
		# Median
		medianIntensityName <- character()
		for (l in length(classList)) {
			medianIntensityName <- paste (medianIntensityName, " ", medianIntensity[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,5] <- medianIntensityName
		# Interquartile Range (IQR)
		interQuartileRangeName <- character()
		for (l in length(classList)) {
			interQuartileRangeName <- paste (interQuartileRangeName, " ", interQuartileRange[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,6] <- interQuartileRangeName
		# Spectra counter
		spectraCounterName <- character()
		for (l in length(classList)) {
			spectraCounterName <- paste (interQuartileRangeName, " ", spectraCounter[[l]], " - ", classList[l], sep="")
		}
		peakStatMatrix [p,7] <- spectraCounterName
		# Class
		#className <- classList
		#peakStatMatrix [p,8] <- className
		# Homoscedasticity (Parametric)
		if (varianceTestPar$p.value < 0.05) {
		homoscedasticityP <- paste ("Non homoscedastic data", "(p-value:", varianceTestPar$p.value, ")")
		}
		if (varianceTestPar$p.value >= 0.05) {
		homoscedasticityP <- paste ("Homoscedastic data", "(p-value:", varianceTestPar$p.value, ")")
		}
		if (varianceTestNonPar$p.value < 0.05) {
		homoscedasticityNP <- paste ("Non homoscedastic data", "(p-value:", varianceTestNonPar$p.value, ")")
		}
		if (varianceTestNonPar$p.value >= 0.05) {
		homoscedasticityNP <- paste ("Homoscedastic data", "(p-value:", varianceTestNonPar$p.value, ")")
		}
		peakStatMatrix [p,9] <- homoscedasticityP
		peakStatMatrix [p,10] <- homoscedasticityNP
		# t-Test
		peakStatMatrix [p,11] <- tTest$p.value
		# ANOVA
		peakStatMatrix [p,12] <- summary(anovaTest)[[1]]$"Pr(>F)"[1]
		# Wilcoxon / Mann-Whitney test
		peakStatMatrix [p,13] <- WilcoxonTest$p.value
		# Kruskal-Wallis test
		peakStatMatrix [p,14] <- KruskalWallisTest$p.value
	}
}
return (peakStatMatrix)
}



################################################################################





################################## ROUND NUMERIC FEATURES PEAKLIST
roundFeaturesPeaklist <- function (peaklist, decimalDigits=3, nonFeatures=c("Sample","Class","THY")) {
# Take the features
featureList <- as.numeric (names (peaklist [,!(names(peaklist) %in% nonFeatures)]))
# Round the features
for (i in 1:length(featureList)) {
	featureList[i] <- round (featureList[i], digits=decimalDigits)
}
# Regather the features
features <- c (featureList, nonFeatures)
names (peaklist) <- features
#
return (peaklist)
}


