############################################## LC-MS URINE STATISTICS 2017.02.14



















#################### INSTALL AND LOAD THE REQUIRED PACKAGES
##### Update the packages
update.packages(repos = "http://cran.mirror.garr.it/mirrors/CRAN/", ask = FALSE)
##### Install the required packages if not already installed
# FUNCTION
install_and_load_required_packages <- function(required_packages, repository="http://cran.mirror.garr.it/mirrors/CRAN/") {
	# Retrieve the installed packages
	installed_packages <- installed.packages() [,1]
	# Determine the missing packages
	missing_packages <- character()
	for (p in 1:length(required_packages)) {
		if ((required_packages[p] %in% installed_packages) == FALSE) {
			missing_packages <- append(missing_packages, required_packages[p])
		}
	}
	# If a repository is specified
	if (repository != "" || !is.null(repository)) {
		if (length(missing_packages) > 0) {
			install.packages(missing_packages, repos=repository)
		}
	} else {
		if (length(missing_packages) > 0) {
			install.packages(missing_packages)
		}
	}
	# Load the packages
	for (i in 1:length(required_packages)) {
		library (required_packages[i], character.only=TRUE)
	}
}
# RUN THE FUNCTION
install_and_load_required_packages(c("XLConnect", "rattle", "phia", "MASS", "ggplot2", "lawstat", "coin", "multcomp", "agricolae", "tcltk", "Hmisc")) # "Rcmdr", "RcmdrPlugin.coin"
# The package lawstat is used for levene test for non parametric anova































###############################################################################

############################# USER PARAMETERS (GUI)

##### Input file
filepath_import_select <- tkmessageBox(title = "Input file", message = "Select the file containing all the mass spectrometric information for the statistics", icon = "info")
input_file <- tclvalue(tkgetOpenFile(filetypes = "{{Comma Separated Value files} {.csv}} {{Microsoft Excel files} {.xls .xlsx}}"))
if (!nchar(input_file)) {
	tkmessageBox(message = "No file selected")
} else {
	tkmessageBox(message = paste("The following file will be read:", input_file))
}

##### Output folder
output_folder_select <- tkmessageBox(title = "Output folder", message = "Select the folder where all the outputs are saved", icon = "info")
output_folder <- tclvalue(tkchooseDirectory())
if (!nchar(output_folder)) {
	tkmessageBox(message = "No folder selected")
}	else {
	tkmessageBox(message = paste("Every file will be saved in", output_folder))
}
setwd(output_folder)


##### Output Format
output_format_select <- tkmessageBox(title = "Output format", message = "Select the file format for the output files", icon = "info")
output_format <- select.list(c("Comma Separated Values (.csv)", "Microsoft Excel (.xls)", "Microsoft Excel (.xlsx)"), title = "Choose")
# Fix the file format
if (output_format == "Comma Separated Values (.csv)" || output_format == "") {
	file_format <- "csv"
} else if (output_format == "Microsoft Excel (.xlsx)") {
	file_format <- "xlsx"
} else if (output_format == "Microsoft Excel (.xls)") {
	file_format <- "xls"
}

##### Image Format
image_format_select <- tkmessageBox(title = "Image format", message = "Select the file format for the image files", icon = "info")
image_format <- select.list(c("JPG (.jpg)", "PNG (.png)", "TIFF (.tiff)"), title = "Choose")
# Fix the file format
if (image_format == "JPG (.jpg)") {
	image_format <- ".jpg"
} else if (image_format == "PNG (.png)" || image_format == "") {
	image_format <- ".png"
} else if (image_format == "TIFF (.tiff)") {
	image_format <- ".tiff"
}



































# PVALUE
pvalue_expression <- 0.05
pvalue_tests <- 0.05




# Sampling parameters
TestPer_Base <- 0.17	#fair for class
TestPer_Adv <- 0.19		#fair for class
minimum_number_of_patients <- 3



### Data Management
data_record <- TRUE
### Correlation Analysis
correlation_analysis <- TRUE
# Outlier estimation and removal
remove_outliers <- FALSE
### Effect Analysis (2-Level Effects)
# One class vs all the others (if more classes are present); one class vs the other class (if only two classes are present) 
two_level_effect_analysis  <- TRUE
# The sampling = TRUE performs the sampling of around 20% of patients/ctrls to be used in R... The other 80% is exported for RapidMiner analysis (with the filtered signals)
sampling <- FALSE
# Outlier estimation in the 2-level effect analysis
remove_outliers_two_level_effect_analysis <- FALSE
### Effect Analysis (Multi-Level Effects) Stima pT (pT4 e pT5 sono uniti) e Grade
multi_level_effect_analysis <- TRUE
remove_outliers_multi_level_effect_analysis <- FALSE
# Age binning
age_binning <- TRUE
# Non signals (column)
non_signals <- c("Class", "Age", "Age_BINNED", "pT", "pT_2009", "Grade", "Dim", "No", "Centro", "DIA_1_2")
# Discriminant column
discriminant_feature <- "Class" ## allow to select it from a list after reading the file (list of column names)
####################################################




























########## Exception handling
if (isTRUE(two_level_effect_analysis)) {
	ifelse(isTRUE(sampling), print("The Two-Level Effect analysis is performed with sampling for RM"), print("WARNING!!! The Two-Level Effect analysis is performed without sampling!"))
}



























################################################################################
################################# FUNCTIONS


########## General Statistic Function: definition ###############################
### skip this !!!!!!!!!!!!!!!!!!
signif <- function(object, pvalue_tests) {
	if (object$p.value <= pvalue_tests) {
		object$p.value
	} else {
		NULL
	}
}
################################

assumption_p <- function(objt) {
	ifelse(is.object(objt), return(objt$p.value), return(NA))
}
assumption_anova <- function(objt) {
	ifelse(is.object(objt), return(objt$Pr[1]), return(NA))
}
assumption_permutation <- function(objt) {
	ifelse(is.object(objt), return(pvalue(objt)[1]), return(NA))
}


#################  Post hoc 
write_posthoc_file <- function(file_name, data, file_format = "xlsx"){
	# Store the original filename
	original_file_name <- file_name
	# The sheet name must not contain more than 31 characters
	length_sheet_name <- length(unlist(strsplit(original_file_name, "")))
	# Fix the sheet name if it is too long...
	if (length_sheet_name > 31) {
		sheet_name <- unlist(strsplit(original_file_name, ""))[1:31]
		final_sheet_name <- ""
		for (ch in 1:length(sheet_name)) {
			final_sheet_name <- paste(final_sheet_name, sheet_name[ch], sep = "")
		}
	} else {
		final_sheet_name <- original_file_name
	}
	# Fix the extension
	file_name <- paste(file_name, ".", file_format, sep = "")
	### Excel
	if (file_format == "xls" || file_format == "xlsx") {
		wb = loadWorkbook(file = file_name, create = TRUE)
		# Create a new sheet
		createSheet(wb, name = final_sheet_name)
		# cumulative length (rows) of matrices
		# +2 = 1 for list names, 1 for header row
		cumlen = cumsum(c(1, head(sapply(data, nrow), n = - 1) + 2))
		# Write data rows (implicitly vectorized!)
		writeWorksheet(wb, data = data, sheet = final_sheet_name, startRow = cumlen + 1, header = TRUE)
		# Write list names
		writeWorksheet(wb, data = as.list(names(data)), sheet = final_sheet_name, startRow = cumlen, header = FALSE)
		saveWorkbook(wb)
	} else if (file_format == "csv") {
		output_matrix <- matrix(nrow = 1, ncol = length(data))
		# Fill the matrix
		for (l in 1:length(data)) {
			output_matrix[1, l] <- data[[l]]
		}
		colnames(output_matrix) <- names(sapply(data, nrow))
		write.csv(output_matrix, file = file_name)
	}
}
###################################################

sammy <- function(pn, minpos, maxpos){
	sam <- NULL
	while (length(sam) < pn){
		z <- round(runif(1,min=minpos,max=maxpos))
		if (!(z %in% sam))
			sam <- append(sam, z)
		}
	return(sam)
	}	
	
	
##### Split the dependent variable according to the factor variable
group_dependent_variable <- function(dependent_variable, factor_variable) {
	# Extract the levels of the factor variable
	factor_levels <- levels(factor_variable)
	# Determine the number of levels
	number_of_levels <- length(factor_levels)
	# Generate a list in which each element is relative to the observations of the dependent variable with a certain value of ther factor variable
	dependent_variable_split <- list()
	# For each level of the factor variable...
	for (i in 1:number_of_levels){
		# Extract the dependent variable observations with that factor variable value
		dependent_variable_split[[i]] <- dependent_variable[factor_variable == factor_levels[i]]
	}
	# Return
	return(dependent_variable_split)
}





###################################################
# Function to export the data files
write_file <- function(file_name, data, file_format = "xlsx"){
	# Store the original filename
	original_file_name <- file_name
	# The sheet name must not contain more than 31 characters
	length_sheet_name <- length(unlist(strsplit(original_file_name, "")))
	# Fix the sheet name if it is too long...
	if (length_sheet_name > 31) {
		sheet_name <- unlist(strsplit(original_file_name, ""))[1:31]
		final_sheet_name <- ""
		for (ch in 1:length(sheet_name)) {
			final_sheet_name <- paste(final_sheet_name, sheet_name[ch], sep = "")
		}
	} else {
		final_sheet_name <- original_file_name
	}
	# Fix the extension
	file_name <- paste(file_name, ".", file_format, sep = "")
	if (file_format == "xls" || file_format == "xlsx") {
		wb <- loadWorkbook(file_name, create = TRUE)
		# Create a new sheet
		createSheet(wb, name = final_sheet_name)
		# Write data rows (implicitly vectorized!)
		writeWorksheet(wb, data = data, sheet = final_sheet_name, header = TRUE)
		saveWorkbook(wb)
	} else if (file_format == "csv") {
		write.csv(data, file = file_name)
	}
}

SamplingForRM <- function(Per_Base , Per_Adv){
			
	##### sampling indexes from group 0
	p <- Per_Base
	minval <- 1
	maxval <- dim(BaseTFD)[1]
	print("numb. of base subj")
	print(maxval)
	pn <- round(p*maxval)
	print("numb. of base subj to sample")
	print(pn)
	set.seed(123)
	n1 <- sammy(pn,minval, maxval)
	print("Selected base-group indeces for Sig.an.")
	print(n1)
	##### sampling indexes from group 1
	p <- Per_Adv
	minval <- 1
	maxval <- dim(AdvTFD)[1]
	print("numb. of adv subj")
	print(maxval)
	pn <- round(p*maxval)
	print("numb. of adv subj to sample")
	print(pn)
	n2 <- sammy(pn,minval, maxval)	
	print("Selected adv-group indeces for Sig.an.")
	print(n2)
	return (list(N1=n1,N2=n2))
	}
		
		
################################################################################





























#################### IMPORT THE DATA FROM THE FILE
#XLSInputDataName <- input_file
#InputDataName <- sprintf("%s%s",output_folder,XLSInputDataName)
#excel.InputFile <- file.path(InputDataName)
##### XLS or XLSX
if (length(grep(".xls", input_file, fixed = TRUE)) > 0) {
	input_data <- readWorksheetFromFile(input_file, sheet = 1)
	input_data <- as.data.frame(input_data)
} else if (length(grep(".csv", input_file, fixed = TRUE)) > 0) {
##### CSV
	input_data <- read.csv(input_file)
}
## Rownames
rownames(input_data) <- input_data$No
## Data type
input_data$Age <- as.numeric(input_data$Age)
input_data$pT <- as.numeric(input_data$pT)
input_data$pT_2009 <- as.factor(input_data$pT_2009)
input_data$Dim <- as.numeric(input_data$Dim)
input_data$Grade <- as.numeric(input_data$Grade)
input_data$Class <- as.factor(input_data$Class)
input_data$Centro <- as.character(input_data$Centro)
## Age binning
if (isTRUE(age_binning)) {
	age_bins <- 6
	Age_BINNED <- binning(input_data$Age, age_bins, method = "quantile", ordered=TRUE)
	levels_age_binned <- 
	levels(Age_BINNED) <- c(1:age_bins)
	# Replace the Age with the Age_BINNED
	input_data$Age <- cbind(Age_BINNED)
}
## Class list
class_list <- levels(input_data[,discriminant_feature])


##### Separate the mass spectrometric data from the demographic data
non_signals_data <- input_data[,(names(input_data) %in% non_signals)]
signals_data <- input_data[,!(names(input_data) %in% non_signals)]
number_of_signals = ncol(signals_data)
# Matrix of signal names
signal_name <- data.matrix(colnames(signals_data))
colnames(signal_name) <- "Signal Name"




























################ Data from input_data #######################
if (isTRUE(data_record)) {
	print("########## Data Management ##########")
	# Generate a list of data frames and corresponding filenames to dump
	list_of_dataframes <- list()
	list_of_filenames <- character()
	# If there are only two classes...
	if (length(class_list <= 2)) {
		# Generate a dataframe with those two classes, otherwise dump one dataframe per each class
		list_of_dataframes[[1]] <- input_data[,!(names(input_data) %in% non_signals)]
		list_of_filenames <- "Control_vs_Disease"
	} else if (length(class_list > 2)) {
		# If there are more than two classes...
		for (cl in 1:length(class_list)) {
			# Put every dataframe in the list
			list_of_dataframes[[cl]] <- input_data[(input_data[,discriminant_feature] == class_list[cl]), !(names(input_data) %in% non_signals)]
			list_of_filenames <- append(list_of_filenames, paste("Class", class_list[cl]))
		}
	}

#TFD_CntrRcc <- as.data.frame(subset(input_data, input_data$Class == 1 | input_data$Class == 2 , select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade) ))
#BaseClass <- TFD_CntrRcc[["Class"]]
#BaseClass <- as.data.frame(BaseClass)
# Rename "2" to "1" and "1" to "0" 
#TFD_CntrRcc$Class[TFD_CntrRcc$Class==1] <- 0
#TFD_CntrRcc$Class[TFD_CntrRcc$Class==2] <- 1
#TFD_CntrRcc <- cbind(BaseClass,TFD_CntrRcc)

#TFD_DIA <- as.data.frame(subset(input_data, input_data$Class == 3, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade) ))	
#BaseClass <- TFD_DIA[["Class"]]
#BaseClass <- as.data.frame(BaseClass)
#Class <- rep(0,nrow(TFD_DIA))
#TFD_DIA$Class <- Class
#TFD_DIA <- cbind(BaseClass,TFD_DIA)

#TFD_DIAMA <- as.data.frame(subset(input_data, input_data$Class == 4, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade) ))	
#BaseClass <- TFD_DIAMA[["Class"]]
#BaseClass <- as.data.frame(BaseClass)
#Class <- rep(0,nrow(TFD_DIAMA))
#TFD_DIAMA$Class <- Class
#TFD_DIAMA <- cbind(BaseClass,TFD_DIAMA)

#TFD_M <- as.data.frame(subset(input_data, input_data$Class == 5, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade) ))							
#BaseClass <- TFD_M[["Class"]]
#BaseClass <- as.data.frame(BaseClass)
#Class <- rep(1,nrow(TFD_M))
#TFD_M$Class <- Class
#TFD_M <- cbind(BaseClass,TFD_M)

#TFD_B <- as.data.frame(subset(input_data, input_data$Class == 6, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade) ))			
#BaseClass <- TFD_B[["Class"]]
#BaseClass <- as.data.frame(BaseClass)
#Class <- rep(0,nrow(TFD_B))
#TFD_B$Class <- Class
#TFD_B <- cbind(BaseClass,TFD_B)

# Create the folder containing all the data frames (and go to it temporarily to dump the files)
dir.create(file.path(output_folder, "Data"))
setwd(file.path(output_folder, "Data"))
# Dump the dataframes
for (d in 1:length(list_of_dataframes)) {
	write_file(file_name = list_of_filenames[d], data = list_of_dataframes[[d]], file_format = file_format)
}
# Go back to the original output folder
setwd(output_folder)
}























################################################################################
#################### CORRELATION ANALYSIS (Done only with patients affected by RCC, because controls do not have pT, grade, dimension, etc...)
if (isTRUE(correlation_analysis)) {
	# Features for correlation analysis (sorted for reproducibility reasons)
	non_signals_for_correlation_analysis <- sort(c("Age", "pT", "Grade", "Dim"))
	features_for_correlation_analysis <- append(non_signals_for_correlation_analysis, names(signals_data))
	print("########## Correlation Analysis ##########")
	##### Create the folder for the correlation data
	dir.create(file.path(output_folder, "Correlation"))
	setwd(file.path(output_folder, "Correlation"))

	########## Task #1: -- RCC correlation analysis 
	#Rcc Age VS Signal intensity
	#Rcc Dim VS Signal intensity
	#Rcc pT (four levels) VS Signal intensity
	#Rcc Grade (two levels) VS Signal intensity
	# Temporary data filter for task #1   
	### For each class...
	for (cl in 1:length(class_list)) {
		# Message
		print(paste("Task 1 - RCC correlation analysis, class:", class_list[cl]))
		# Filter the dataframe with only the data for that selected class
		data_frame_correlation <- input_data[input_data[, discriminant_feature] == class_list[cl], features_for_correlation_analysis]
		# Store the original (after it will be modified by the outlier exclusion)
		data_frame_correlation_original <- data_frame_correlation
		### The correlation list is a list of correlation lists: one list for each correlation, each one of them is a list of the correlation with the signal intensity and the variable (correlation_list[[non_signal1]][[signal1]], correlation_list[[non_signal1]][[signal2]], correlation_list[[non_signal1]][[signal3]], correlation_list[[non_signal2]][[signal1]], correlation_list[[non_signal2]][[signal2]], correlation_list[[non_signal2]][[signal3]], ...)
		correlation_list <- list()
		# List of outliers (one dataframe for each mass)
		outlier_list <- list()
		# List of the correlation matrices
		correlation_matrix_list <- list()
		### For each non-signal...
		for (f in 1:length(non_signals_for_correlation_analysis)) {
			ns <- non_signals_for_correlation_analysis[f]
			# Generate the sublist...
			correlation_sublist <- list()
			# Isolate the column corresponding to the non-signal...
			non_signal_column <- data_frame_correlation[, ns]
			##### Scroll the masses...
			for (ms in 1:ncol(signals_data)) {
				m <- colnames(signals_data)[ms]
				# Isolate the column corresponding to the mass (as a vector)
				mass_x <- data_frame_correlation[, m]
				### Outlier detection
				if (isTRUE(remove_outliers)) {
					# Detect the outliers (fence)
					outliers <- boxplot(mass_x, plot = FALSE)$out
					# If there are some...
					if (length(outliers) > 0) {
						# Build the dataframe to be dumped
						outliers_dataframe <- data_frame_correlation[mass_x %in% outliers, c(non_signals_for_correlation_analysis, m)]
						# Signal name
						signal_name <- rep(m, dim(outliers_dataframe)[1])
						# Patients
						patients <- rownames(outliers_dataframe)
						# Append the two columns...
						outliers_dataframe <- cbind(outliers_dataframe, patients, signal_name)
						# Fix the column names
						colnames(outliers_dataframe) <- c(non_signals_for_correlation_analysis,"Intensity","Patient","Mass")
						# Store this in the final list of outlier dataframes
						outlier_list[[ms]] <- outliers_dataframe
						# Replace the intensity in the original dataframe with NA (so that they are excluded)
						data_frame_correlation[mass_x %in% outliers, m] <- NA
						mass_x[mass_x %in% outliers] <- NA
					}
				}
				### Compute the results
				correlation_result_vector <- c(cor.test(mass_x, non_signal_column, method = "spearman")$estimate, cor.test(mass_x, non_signal_column, method = "spearman")$p.value)
				### Correlation (signal m with the non-signal ns)
				# Append to the final correlation list
				correlation_list[[ns]][[m]] <- correlation_result_vector
			}
		}
		### Correlation matrices
		# For each element of the correlation list...
		for (l in 1:length(correlation_list)) {
			# Generate a matrix from the correlation list (one matrix for each non-signal)
			corr_matrix <- matrix(0, nrow = ncol(signals_data), ncol = 2)
			corr_matrix <- do.call(rbind, correlation_list[[l]])
			colnames(corr_matrix) <- c(paste("Corr", class_list[cl], non_signals_for_correlation_analysis[l]), "p_value")
			correlation_matrix_list[[l]] <- corr_matrix
		}
		### Build the final correlation matrix
		final_correlation_matrix <- as.matrix(cbind(signal_name))
		for (l in 1:length(correlation_matrix_list)) {
			final_correlation_matrix <- cbind(final_correlation_matrix, correlation_matrix_list[[l]])
		}
		### Outlier matrices
		if (isTRUE(remove_outliers)) {
			final_outlier_matrix <- matrix(ncol = (length(non_signals_for_correlation_analysis)) + 3)
			final_outlier_matrix <- do.call(rbind, outlier_list)
			colnames(final_outlier_matrix) <- c(non_signals_for_correlation_analysis , "Intensity", "Patient", "Mass")
		}
		#### Save the outputs
		write_file(file_name = paste(class_list[cl], "- Correlation matrix"), data = final_correlation_matrix, file_format = file_format)	
		if (isTRUE(remove_outliers)) {
			write_file(file_name = paste(class_list[cl], "- Outlier matrix"), data = final_outlier_matrix, file_format = file_format)
		}
	}
}



























#################################################################################
#################### 2-LEVEL EFFECT ANALYSIS OVER SIGNAL INTENSITY
#####  Simple effects (2 levels): 
##### (1 - CntrRcc) Cntrl / Rcc group ,
##### (2 - Rcc PT1) Rcc-pt 1 / Rcc-pt 2 + 3 + 4 + 5
##### (3 - Rcc PT2) Rcc-pt 2 / Rcc-pt 3 + 4 + 5
##### (4 - RccGrade1) Rcc-grade 1 / Rcc-pt 2 + 3 + 4
##### (5 - RccGrade2) Rcc-grade 2 / Rcc-pt 3 + 4


if (isTRUE(two_level_effect_analysis)) {
	# Establish the variables for the two-level effect analysis
	two_level_effect_analysis_non_features <- sort(c(discriminant_feature, "pT", "Grade"))
	#features_for_correlation_analysis <- append(non_signals_for_correlation_analysis, names(signals_data))
	print("########## Simple Effect Analysis (2 LEVELS) ##########")
	### Build up the combinations of possibilities
	# Empty vector
	combination_vector <- character()
	# For each non-signal variable to be considered...
	for (t in 1:length(sort(two_level_effect_analysis_non_features))) {
		# Isolate the dataframe column
		data_column <- input_data[, two_level_effect_analysis_non_features[t]]
		# Extract the (unique) values
		for (ch in sort(unique(as.character(data_column)))) {
			# Append the combination ID to the vector...
			combination <- ch
			names(combination) <- two_level_effect_analysis_non_features[t]
			combination_vector <- append(combination_vector, combination)
		}
	}
	# Generate the combination names (each element of the vector vs the others)
	combination_names <- character()
	for (l in 1:length(combination_vector)) {
		combination_name <- paste(names(combination_vector[l]), combination_vector[l], "vs OTHERS")
		combination_names <- append(combination_names, combination_name)
	}
	##### For each combination... (everytime fall on the 0-1 case, by replacing values)
	for (comb in 1:length(combination_vector)) {
		# Message
		print(paste("Processing ->", combination_names[comb]))
		########## Each time put one class with the value 0 and the others with 1
		# Non signal variable
		non_signal_variable <- names(combination_vector[comb])
		# Non signal variable level selected (it will be put = 0)
		l <- combination_vector[comb]
		# Extract the selected factor level (= variable value)
		#l <- levels(as.factor(input_data[, non_signal_variable]))[lv]
		# Isolate the dataframe with the signals + the selected non-signal variable
		temp_data_frame <- input_data[, c(non_signal_variable, names(signals_data))]
		# Store the original data frame
		temp_data_frame_original <- temp_data_frame
		# Convert the selected non-feature variable to character
		temp_data_frame[, non_signal_variable] <- as.character(temp_data_frame[, non_signal_variable])
		# Replace one level (the current level l) with 0 and the other ones with 1
		temp_data_frame[temp_data_frame[, non_signal_variable] == l, non_signal_variable] <- "0"
		temp_data_frame[temp_data_frame[, non_signal_variable] != "0", non_signal_variable] <- "1"
		##### Sampling
		if (isTRUE(sampling)){
			
			base_sampling_df <- temp_data_frame[temp_data_frame[, non_signal_variable] == 0, ]
			adv_sampling_df <- temp_data_frame[temp_data_frame[, non_signal_variable] == 1, ]
			
			Samlist <- SamplingForRM(TestPer_Base,TestPer_Adv)	
			n1 <- Samlist$N1
			n2 <- Samlist$N2
			
			TFD <- rbind(BaseTFD[n1,],AdvTFD[n2,])
			DIAGTFD <- rbind(BaseTFD[-n1,],AdvTFD[-n2,])
		}
		# List of outputs
		outlier_list <- list()
		# Shapiro tests for class 0 and class 1
		shapiro_list_0 <- list()
		shapiro_list_1 <- list()
		number_of_patients_list <- list()
		variance_test_list <- list()
		Leven_test_list <- list()
		EqTTest_list <- list()
		UneqTTest_list <- list()
		MWUTest_list <- list()
		KSTest_list <- list()
		# For each signal...
		for (ms in 1:ncol(signals_data)) {
			# Extract the name of the signal...
			m <- colnames(signals_data)[ms]
			## Extract the column of the signal and the response variable...
			mass_x <- temp_data_frame[, m]
			response_variable <- as.integer(temp_data_frame[, non_signal_variable])
			### Outlier analysis 
			if (isTRUE(remove_outliers_two_level_effect_analysis)) {
				## Outlier detection
				# Detect the outliers (fence)
				outliers <- boxplot(mass_x ~ response_variable, plot = FALSE)$out
				# If there are some...
				if (length(outliers) > 0) {
					# Build the dataframe to be dumped
					outliers_dataframe <- temp_data_frame[mass_x %in% outliers, c(non_signal_variable, m)]
					# Signal name
					signal_name <- rep(m, dim(outliers_dataframe)[1])
					# Patients
					patients <- rownames(outliers_dataframe)
					# Append the two columns...
					outliers_dataframe <- cbind(outliers_dataframe, patients, signal_name)
					# Fix the column names
					colnames(outliers_dataframe) <- c(non_signal_variable,"Intensity","Patient","Mass")
					# Store this in the final list of outlier dataframes
					outlier_list[[ms]] <- outliers_dataframe
					# Replace the intensity in the original dataframe with NA (so that they are excluded)
					temp_data_frame[mass_x %in% outliers, m] <- NA
					mass_x[mass_x %in% outliers] <- NA
				}
			}
			### Clean the NA values from the mass column and the response variable column
			response_variable <- response_variable[!is.na(mass_x)]
			mass_x <- mass_x[!is.na(mass_x)]
			## Group the dependent variable (the mass column) according to the response variable
			mass_x_split <- group_dependent_variable(mass_x, as.factor(response_variable))
			# Define the groups (response variable = 0, response variable = 1)
			group_0 <- mass_x_split[[1]]
			group_1 <- mass_x_split[[2]]
			# Define the names and the number of patients
			name_group_0 <- sprintf("%s%s", m, "_0")
			name_group_1 <- sprintf("%s%s", m, "_1")
			number_of_patients_list[[m]] <- c(length(group_0),length(group_1))
			### If the number of patients in the two lists is more than the minimum number of patients allowed and all the elements of the two groups are not the same...
			if (length(group_0) >= minimum_number_of_patients && length(group_1) >= minimum_number_of_patients && !all(group_0 == group_0[1]) && !all(group_1 == group_1[1]) ) {
				# Checking for normal distributed data in the two groups (class 0 and class 1)
				shapiro_list_0[[name_group_0]] <- shapiro.test(group_0)
				shapiro_list_1[[name_group_1]] <- shapiro.test(group_1)
				# Checking variances (Normal data)
				variance_test_list[[m]] <- var.test(group_0, group_1)
				# Checking variances (Non-Normal data)
				Leven_test_list[[m]] <- levene.test(response_variable, mass_x, location = "mean", kruskal.test = T)
				# Normal Data, Equal Variances --> Equal Variance t-test
				EqTTest_list[[m]] <- t.test(group_0, group_1)
				# Normal Data, Unequal Variances --> Unequal Variance t-test:
				UneqTTest_list[[m]] <- t.test(group_0, group_1,var.equal = FALSE)
				# Non-normal Data, Equal Variances --> Mann-Whitney U-test (Wilcoxon)
				MWUTest_list[[m]] <- wilcox.test(group_0, group_1,paired = FALSE)
				# Non-normal Data, Unequal Variances --> Kolmogorov-Smirnov test
				KSTest_list[[m]] <- ks.test(group_0, group_1)
			} else {
				### Otherwise all the results are NA...
				shapiro_list_0[[name_group_0]] <- NA
				shapiro_list_1[[name_group_1]] <- NA
				variance_test_list[[m]] <- NA
				Leven_test_list[[m]] <- NA	
				EqTTest_list[[m]] <- NA
				UneqTTest_list[[m]] <- NA
				MWUTest_list[[m]] <- NA
				KSTest_list[[m]] <- NA
			}
		}
		# Extract the pvalue list from the test lists...
		list_of_shapiro_list_0 <- lapply(shapiro_list_0, assumption_p)
		list_of_shapiro_list_1 <- lapply(shapiro_list_1, assumption_p)
		list_of_variance_test_list <- lapply(variance_test_list, assumption_p)
		list_of_Leven_test_list <- lapply(Leven_test_list, assumption_p)
		list_of_EqTTest_list <- lapply(EqTTest_list, assumption_p)
		list_of_UneqTTest_list <- lapply(UneqTTest_list, assumption_p)
		list_of_MWUTest_list <- lapply(MWUTest_list, assumption_p)
		list_of_KSTest_list <- lapply(KSTest_list, assumption_p)
		# Vectorize the lists...
		vector_of_shapiro_list_0 <- unlist(list_of_shapiro_list_0)
		vector_of_shapiro_list_1 <- unlist(list_of_shapiro_list_1)
		vector_of_variance_test_list <- unlist(list_of_variance_test_list)
		vector_of_Leven_test_list <- unlist(list_of_Leven_test_list)
		vector_of_EqTTest_list <- unlist(list_of_EqTTest_list)
		vector_of_UneqTTest_list <- unlist(list_of_UneqTTest_list)
		vector_of_MWUTest_list <- unlist(list_of_MWUTest_list)
		vector_of_KSTest_list <- unlist(list_of_KSTest_list)
		# Matrix of patient numbers...
		patient_number_matrix <- matrix(ncol = 2)	
		patient_number_matrix <- do.call(rbind, number_of_patients_list)
		patient_number_matrix <- cbind(rownames(patient_number_matrix),patient_number_matrix)
		colnames(patient_number_matrix) <- c("Mass","Group 0","Group 1")
		# Outlier matrix
		outlier_matrix <- do.call(rbind, outlier_list)
		# pvalue  matrix
		pvalue_matrix <- data.matrix(cbind(vector_of_shapiro_list_0, vector_of_shapiro_list_1, vector_of_variance_test_list, vector_of_Leven_test_list, vector_of_EqTTest_list, vector_of_UneqTTest_list, vector_of_MWUTest_list, vector_of_KSTest_list))
		rownames(pvalue_matrix) <- names(signals_data)
		### Extract the data according to the conditions
		# Extract the normal and non-normal data
		normal_data <- as.data.frame(subset(pvalue_matrix, vector_of_shapiro_list_0 > pvalue_tests & vector_of_shapiro_list_1 > pvalue_tests))
		non_normal_data <- as.data.frame(subset(pvalue_matrix, vector_of_shapiro_list_0 <= pvalue_tests | vector_of_shapiro_list_1 <= pvalue_tests))
		# Extract the normal homoschedastic data
		normal_homoschedastic_data <- as.data.frame(subset(normal_data, vector_of_variance_test_list > pvalue_tests))
		# Extract the normal heteroschedastic data
		normal_heteroschedastic_data <- as.data.frame(subset(normal_data, vector_of_variance_test_list <= pvalue_tests))
		# Extract the non-normal homoschedastic data
		non_normal_homoschedastic_data <-	as.data.frame(subset(non_normal_data, vector_of_Leven_test_list > pvalue_tests))
		# Extract the non-normal heteroschedastic data
		non_normal_heteroschedastic_data <- as.data.frame(subset(non_normal_data, vector_of_Leven_test_list <= pvalue_tests))
		# Differentially expressed signals: normal homoschedastic data
		diff_normal_homoschedastic_data <- as.data.frame(subset(normal_homoschedastic_data, vector_of_EqTTest_list <= pvalue_expression))
		method_diff_norm_homo <- rep("Equal Variance t-test", nrow(diff_normal_homoschedastic_data))
		matrix_diff_norm_homo <- as.data.frame(subset(diff_normal_homoschedastic_data, select = vector_of_EqTTest_list))
		matrix_diff_norm_homo <- cbind(rownames(matrix_diff_norm_homo), matrix_diff_norm_homo, method_diff_norm_homo)
		colnames(matrix_diff_norm_homo) <- c("Signal","pvalue","Test")
		# Differentially expressed signals: normal heteroschedastic data
		diff_normal_heteroschedastic_data <- as.data.frame(subset(normal_heteroschedastic_data, vector_of_UneqTTest_list <= pvalue_expression))
		method_diff_norm_hetero <- rep("Unequal Variance t-test", nrow(diff_normal_heteroschedastic_data))
		matrix_diff_norm_hetero <- as.data.frame(subset(diff_normal_heteroschedastic_data, select = vector_of_UneqTTest_list))
		matrix_diff_norm_hetero <- cbind(rownames(matrix_diff_norm_hetero), matrix_diff_norm_hetero, method_diff_norm_hetero)
		colnames(matrix_diff_norm_hetero) <- c("Signal","pvalue","Test")
		# Differentially expressed signals: non_normal homoschedastic data
		diff_non_normal_homoschedastic_data <- as.data.frame(subset(non_normal_homoschedastic_data, vector_of_MWUTest_list <= pvalue_expression))
		method_diff_non_norm_homo <- rep("Wilcoxon",nrow(diff_non_normal_homoschedastic_data))
		matrix_diff_non_norm_homo <- as.data.frame(subset(diff_non_normal_homoschedastic_data, select = vector_of_MWUTest_list))
		matrix_diff_non_norm_homo <- cbind(rownames(matrix_diff_non_norm_homo), matrix_diff_non_norm_homo, method_diff_non_norm_homo)
		colnames(matrix_diff_non_norm_homo) <- c("Signal","pvalue","Test")
		# Differentially expressed signals: non_normal heteroschedastic data
		diff_non_normal_heteroschedastic_data <- as.data.frame(subset(non_normal_heteroschedastic_data, vector_of_KSTest_list <= pvalue_expression))
		method_diff_non_norm_hetero <- rep("Kolmogorov-Smirnov test", nrow(diff_non_normal_heteroschedastic_data))
		matrix_diff_non_norm_hetero <- as.data.frame(subset(diff_non_normal_heteroschedastic_data, select = vector_of_KSTest_list))
		matrix_diff_non_norm_hetero <- cbind(rownames(matrix_diff_non_norm_hetero),matrix_diff_non_norm_hetero, method_diff_non_norm_hetero)
		colnames(matrix_diff_non_norm_hetero) <- c("Signal","pvalue","Test")
		## Selected signals for inference
		selected_signals_for_inference <- c(rownames(diff_normal_homoschedastic_data), rownames(diff_normal_heteroschedastic_data), rownames(diff_non_normal_homoschedastic_data), rownames(diff_non_normal_heteroschedastic_data))
		# Print the message...
		print("Differentially expressed signals for inference")
		print(selected_signals_for_inference)
		# Generate a matrix with the signals for inference and the method used to identify them...
		inference_signals_method_matrix <- rbind(matrix_diff_norm_homo, matrix_diff_norm_hetero, matrix_diff_non_norm_homo, matrix_diff_non_norm_hetero)
		# Select the signals to be used for inference and the non-signal variable
		selected_signals_for_inference_intensity_df <- subset(temp_data_frame, select = c(non_signal_variable, selected_signals_for_inference))
		
		### data for testing: for diagnostic analysis 
		if (isTRUE(sampling)) TEST <- DIAGTFD[c("No",BaseEffName,non_signal_variable,SelectedSignsForEffect)]
		
		############################### Dump the files
		# Create the folder where to dump the files and go to it...
		dir.create(file.path(output_folder, combination_names[comb]))
		setwd(file.path(output_folder, combination_names[comb]))
		# For each signals of inference...
		for (s in selected_signals_for_inference){
			# Extract the intensity
			signal_intensity <- selected_signals_for_inference_intensity_df[[s]]
			# Extract the non-signal variable as an ordered factor
			non_signal_as_ordered_factor <- ordered(temp_data_frame[, non_signal_variable])
			# Remove possible NA values
			non_signal_as_ordered_factor <- non_signal_as_ordered_factor[!is.na(signal_intensity)]
			signal_intensity <- signal_intensity[!is.na(signal_intensity)]
			# Number the observations
			IDs <- c(1:length(signal_intensity))
			# Generate a matrix with the ID, the intensities of the selected signal for inference and the non-signal variable
			signal_dataframe <- data.frame(IDs, signal_intensity, non_signal_as_ordered_factor)
			##### Jitter plot
			plot_name <- sprintf("%s%s", non_signal_variable,"_VS_Intensity")
			file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
			jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5))
			ggsave(jitter_plot, file = file_name, width = 4, height = 4)
			
			##### Box plot
			plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
			file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
			box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot")
			ggsave(box_plot, file = file_name, width = 4, height = 4)
			
			##### Scatter plot
			plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
			# Sort the dataframe rows according to the values of the non-signal variable
			ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
			ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
			colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
			file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
			graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
			scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors)
			ggsave(scatter_plot, file = file_name , width = 4, height = 4)
		}
		### Dump the files
		write_file(file_name = "Patient number matrix", data = patient_number_matrix, file_format = file_format)
		write_file(file_name = "Outliers", data = outlier_matrix, file_format = file_format)
		write_file(file_name = "Methods", data = inference_signals_method_matrix, file_format = file_format)
		write_file(file_name = "Differentially expressed signals", data = selected_signals_for_inference_intensity_df, file_format = file_format)
		write_file(file_name = "All signals", data = temp_data_frame_original, file_format = file_format)
		write_file(file_name = "Test p-values", data = pvalue_matrix, file_format = file_format)
		### Sampling
		if (isTRUE(sampling)) {
			filename <-  "RMdataTEST"
			dataname <- TEST
			WriteFile(filename,dataname, file_format)		
			
			TFD_DIA <- as.data.frame(subset(input_data, input_data$Class == 3, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade,Class) ))	
			TFD_DIA <- TFD_DIA[c("No",SelectedSignsForEffect)]
			RowNumb <- dim(TFD_DIA)[1]
			Class <- rep(0,RowNumb)
			BaseClass <- rep(3,RowNumb)
			TFD_DIA <- cbind(BaseClass,Class,TFD_DIA)
			filename <-  "RMTESTDIABS"
			dataname <- TFD_DIA
			WriteFile(filename,dataname, file_format)		
			
			TFD_DIAMA <- as.data.frame(subset(input_data, input_data$Class == 4, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade,Class) ))			
			TFD_DIAMA <- TFD_DIAMA[c("No",SelectedSignsForEffect)]
			RowNumb <- dim(TFD_DIAMA)[1]
			Class <- rep(0,RowNumb)
			BaseClass <- rep(4,RowNumb)		
			TFD_DIAMA <- cbind(BaseClass,Class,TFD_DIAMA)		
			filename <-  "RMTESTDIAMA.xlsx"
			dataname <- TFD_DIAMA
			WriteFile(filename,dataname)		
			
			TFD_M <- as.data.frame(subset(input_data, input_data$Class == 5, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade,Class) ))							
			TFD_M <- TFD_M[c("No",SelectedSignsForEffect)]
			RowNumb <- dim(TFD_M)[1]
			Class <- rep(1,RowNumb)
			BaseClass <- rep(5,RowNumb)		
			TFD_M <- cbind(BaseClass,Class,TFD_M)				
			filename <-  "RMTESTMAL"
			dataname <- TFD_M
			WriteFile(filename,dataname, file_format)		
			
			TFD_B <- as.data.frame(subset(input_data, input_data$Class == 6, select=-c(DIA_1_2,Centro,Age,Age_BINNED,pT_2009,pT,Dim,Grade,Class) ))			
			TFD_B <- TFD_B[c("No",SelectedSignsForEffect)]
			RowNumb <- dim(TFD_B)[1]
			Class <- rep(0,RowNumb)
			BaseClass <- rep(6,RowNumb)		
			TFD_B <- cbind(BaseClass,Class,TFD_B)				
			filename <-  "RMTESTBEN"
			dataname <- TFD_B
			WriteFile(filename,dataname, file_format)				
			
		}
		# Go back to the output folder
		setwd(output_folder)
	}
}









#################################################################################
#################### MULTI-LEVEL EFFECT ANALYSIS OVER SIGNAL INTENSITY
#####  Simple effects (4 levels): 
##### (2 - Rcc PT) Rcc-pt 1 2 3 4+5
##### (4 - Rcc Grade) Rcc-grade 1 2 3 4
########################################################################################################################
########################################################################################################################

if (isTRUE(multi_level_effect_analysis)) {
	# Establish the variables for the multi-level effect analysis
	multi_level_effect_analysis_non_features <- sort(c(discriminant_feature, "pT", "Grade"))
	print("########## Simple Effect Analysis (MULTI-LEVEL) ##########")
	### Build up the combinations of possibilities (the possibilities are simply the name of the non-signal variables, since in each test, the test automatically reads all the levels of the variable)
	combination_vector <- multi_level_effect_analysis_non_features
	##### For each combination... (everytime read all the levels)
	for (comb in 1:length(combination_vector)) {
		print(paste("Simple effects (multi-level) -> ", combination_vector[comb]))
		# Non signal variable
		non_signal_variable <- combination_vector[comb]
		# Isolate the dataframe with the signals + the selected non-signal variable
		temp_data_frame <- input_data[, c(non_signal_variable, names(signals_data))]
		# Store the original data frame
		temp_data_frame_original <- temp_data_frame
		# Convert the selected non-feature variable to character
		temp_data_frame[, non_signal_variable] <- as.character(temp_data_frame[, non_signal_variable])
		# Levels of the non-signal variable
		levels_non_signal_variable <- levels(as.factor(temp_data_frame[, non_signal_variable]))
		### Output initialization
		shapiro_list_multi <- list() #one element per level of the non-signal variable
		for (lv in 1:length(levels_non_signal_variable)) {
			shapiro_list_multi[[lv]] <- list()
		}
		bartlett_list_multi <- list()
		kruskal_list_multi <- list()
		leven_list_multi <- list()
		welch_list_multi <- list()
		anova_list_multi <- list()
		permutation_test_list_multi <- list()
		number_of_patients_list_multi <- list()
		outlier_list_multi <- list()
		### For each signal...
		for (ms in 1:ncol(signals_data)) {
			# Extract the name of the signal...
			m <- colnames(signals_data)[ms]
			## Extract the column of the signal and the response variable...
			mass_x <- temp_data_frame[, m]
			response_variable <- temp_data_frame[, non_signal_variable]
			if (isTRUE(remove_outliers_multi_level_effect_analysis)) {
				## Outlier detection
				# Detect the outliers (fence)
				outliers <- boxplot(mass_x ~ response_variable, plot = FALSE)$out
				# If there are some...
				if (length(outliers) > 0) {
					# Build the dataframe to be dumped
					outliers_dataframe <- temp_data_frame[mass_x %in% outliers, c(non_signal_variable, m)]
					# Signal name
					signal_name <- rep(m, dim(outliers_dataframe)[1])
					# Patients
					patients <- rownames(outliers_dataframe)
					# Append the two columns...
					outliers_dataframe <- cbind(outliers_dataframe, patients, signal_name)
					# Fix the column names
					colnames(outliers_dataframe) <- c(non_signal_variable,"Intensity","Patient","Mass")
					# Store this in the final list of outlier dataframes
					outlier_list_multi[[ms]] <- outliers_dataframe
					# Replace the intensity in the original dataframe with NA (so that they are excluded)
					temp_data_frame[mass_x %in% outliers, m] <- NA
					mass_x[mass_x %in% outliers] <- NA
				}
			}
			### Clean the NA values (outliers) from the mass column and the response variable column
			response_variable <- as.factor(response_variable[!is.na(mass_x)])
			mass_x <- mass_x[!is.na(mass_x)]
			## Group the dependent variable (the mass column) according to the response variable
			mass_x_split <- group_dependent_variable(mass_x, as.factor(response_variable))
			# Define the names and the number of patients
			number_of_patients_x_multi <- integer()
			for (x in 1:length(mass_x_split)) {
				number_of_patients_x_multi <- append(number_of_patients_x_multi, length(mass_x_split[[x]]))
			}
			number_of_patients_list_multi[[m]] <- number_of_patients_x_multi
			### If the number of patients in the two lists is more than the minimum number of patients allowed and all the elements of the two groups are not the same...
			minimum_number_of_patients_for_all <- FALSE
			if (min(number_of_patients_list_multi[[m]]) >= minimum_number_of_patients) {
				minimum_number_of_patients_for_all <- TRUE
			}
			all_alements_are_the_same <- FALSE
			for (x in 1:length(mass_x_split)) {
				if (length(unique(mass_x_split[[x]])) == 1) {
					all_alements_are_the_same <- TRUE
				}
			}
			# If the condition is satisfied (enough patients per group)...
			if (isTRUE(minimum_number_of_patients_for_all) && !isTRUE(all_alements_are_the_same)) {
				# Checking for normal distributed data in the groups
				for (x in 1:length(mass_x_split)) {
					shapiro_list_multi[[x]][[length(shapiro_list_multi[[x]]) + 1]] <- shapiro.test(mass_x_split[[x]])
				}
				# Checking variances
				bartlett_list_multi[[m]] <- bartlett.test(mass_x ~ response_variable)
				### Anova
				anova_list_multi[[m]] <- anova(lm(formula = mass_x ~ response_variable))
				### Leven
				leven_list_multi[[m]] <- levene.test(mass_x, response_variable, location = "mean", kruskal.test = T)
				### Kruskal
				kruskal_list_multi[[m]] <- kruskal.test(mass_x ~ response_variable)
				### Welch
				welch_list_multi[[m]] <- oneway.test(mass_x ~ response_variable, var.equal = F)
				### k-sets Permutation Test for non parametric heteroschedastic data
				permutation_test_list_multi[[m]] <- oneway_test(mass_x ~ response_variable, alternative = 'two.sided')
			} else {
				for (x in 1:length(mass_x_split)) {
					shapiro_list_multi[[x]][[length(shapiro_list_multi[[x]]) + 1]] <- NA
				}
				bartlett_list_multi[[m]] <- NA
				anova_list_multi[[m]] <- NA
				leven_list_multi[[m]] <- NA	
				kruskal_list_multi[[m]] <- NA
				welch_list_multi[[m]] <- NA
				permutation_test_list_multi[[m]] <- NA
			}
		}
		# Extract the pvalue list from the test lists...
		list_of_shapiro_pvalue <- list()
		for (s in 1:length(shapiro_list_multi)) {
			list_of_shapiro_pvalue[[s]] <- lapply(shapiro_list_multi[[s]], assumption_p)
		}
		list_of_bartlett_pvalue <- lapply(bartlett_list_multi, assumption_p)
		list_of_leven_pvalue <- lapply(leven_list_multi, assumption_p)
		list_of_anova_pvalue <- lapply(anova_list_multi, assumption_anova)
		list_of_kruskal_pvalue <- lapply(kruskal_list_multi, assumption_p)
		list_of_permutation_pvalue <- lapply(permutation_test_list_multi, assumption_permutation)
		list_of_welch_pvalue <- lapply(welch_list_multi, assumption_p)
		# Vectorize the lists...
		vector_of_shapiro_pvalue_list <- list()
		for (l in 1:length(list_of_shapiro_pvalue)) {
			vector_of_shapiro_pvalue_list[[l]] <- unlist(list_of_shapiro_pvalue[[l]])
		}
		for (l in 1:length(vector_of_shapiro_pvalue_list)) {
			names(vector_of_shapiro_pvalue_list[[l]]) <- colnames(signals_data)
		}
		vector_of_bartlett_pvalue <- unlist(list_of_bartlett_pvalue)
		vector_of_leven_pvalue <- unlist(list_of_leven_pvalue)
		vector_of_kruskal_pvalue <- unlist(list_of_kruskal_pvalue)
		vector_of_welch_pvalue <- unlist(list_of_welch_pvalue)
		vector_of_permutation_pvalue <- unlist(list_of_permutation_pvalue)
		vector_of_anova_pvalue <- unlist(list_of_anova_pvalue)
		# Matrix of patient numbers...
		patient_number_matrix_multi <- matrix(ncol = 4)
		patient_number_matrix_multi <- do.call(rbind, number_of_patients_list_multi)
		patient_number_matrix_multi <- cbind(rownames(patient_number_matrix_multi),patient_number_matrix_multi)
		patient_number_matrix_multi_colnames <- character()
		for (lv in 1:length(levels_non_signal_variable)) {
			patient_number_matrix_multi_colnames <- append(patient_number_matrix_multi_colnames, paste(non_signal_variable, levels_non_signal_variable[lv]))
		}
		colnames(patient_number_matrix_multi) <- c("Mass", patient_number_matrix_multi_colnames)
		# Outlier matrix
		outlier_matrix_multi <- do.call(rbind, outlier_list_multi)
		# pvalue  matrix
		p_value_shapiro <- matrix(unlist(vector_of_shapiro_pvalue_list), nrow = ncol(signals_data), byrow = F)
		rownames(p_value_shapiro) <- colnames(signals_data)
		p_value_shapiro_colnames <- character()
		for (i in 1:ncol(p_value_shapiro)) {
			p_value_shapiro_colnames <- append(p_value_shapiro_colnames, paste("Shapiro", i))
		}
		colnames(p_value_shapiro) <- p_value_shapiro_colnames
		pvalue_matrix_multi <- data.matrix(cbind(p_value_shapiro, vector_of_bartlett_pvalue, vector_of_leven_pvalue, vector_of_kruskal_pvalue, vector_of_welch_pvalue, vector_of_permutation_pvalue, vector_of_anova_pvalue))
		rownames(pvalue_matrix_multi) <- names(signals_data)
		# Generate a vector of features explaining if they are normally distributed
		vector_of_normally_distributied_signals <- character()
		vector_of_non_normally_distributied_signals <- character()
		# For each signal...
		for (p in 1:nrow(p_value_shapiro)) {
			# By default the signal is not normally distributed
			signal_is_normally_distributed <- TRUE
			# Scroll the column to change the distribution type according to the pvalue
			for (cl in 1:ncol(p_value_shapiro)) {
				if (!is.na(p_value_shapiro[p, cl])) {
					if (p_value_shapiro[p, cl] <= pvalue_tests) {
						signal_is_normally_distributed <- FALSE
					}
				} else if (is.na(p_value_shapiro[p, cl])) {
					signal_is_normally_distributed <- NULL
				}
			}
			# Add it to the right vector
			if (isTRUE(signal_is_normally_distributed)) {
				# Add this to the list of normally distributed signals
				vector_of_normally_distributied_signals <- append(vector_of_normally_distributied_signals, rownames(p_value_shapiro)[p])
			} else if (is.null(signal_is_normally_distributed)) {
				vector_of_normally_distributied_signals <- NULL
				vector_of_non_normally_distributied_signals <- NULL
			} else if (!isTRUE(signal_is_normally_distributed)) {
				# Add this to the list of non-normally distributed signals
				vector_of_non_normally_distributied_signals <- append(vector_of_non_normally_distributied_signals, rownames(p_value_shapiro)[p])
			} 
		}
		### Extract the data according to the conditions
		# Extract the normal and non-normal data
		normal_data <- as.data.frame(pvalue_matrix_multi[rownames(pvalue_matrix_multi) %in% vector_of_normally_distributied_signals, ])
		non_normal_data <- as.data.frame(pvalue_matrix_multi[rownames(pvalue_matrix_multi) %in% vector_of_non_normally_distributied_signals, ])
		# Extract the normal homoschedastic data
		normal_homoschedastic_data <- as.data.frame(subset(normal_data, vector_of_bartlett_pvalue > pvalue_tests))
		# Extract the normal heteroschedastic data
		normal_heteroschedastic_data <- as.data.frame(subset(normal_data, vector_of_bartlett_pvalue <= pvalue_tests))
		# Extract the non-normal homoschedastic data
		non_normal_homoschedastic_data <-	as.data.frame(subset(non_normal_data, vector_of_leven_pvalue > pvalue_tests))
		# Extract the non-normal heteroschedastic data
		non_normal_heteroschedastic_data <- as.data.frame(subset(non_normal_data, vector_of_leven_pvalue <= pvalue_tests))
		# Differentially expressed signals: normal homoschedastic data
		diff_normal_homoschedastic_data <- as.data.frame(subset(normal_homoschedastic_data, vector_of_anova_pvalue <= pvalue_expression))
		method_diff_norm_homo <- rep("ANOVA", nrow(diff_normal_homoschedastic_data))
		matrix_diff_norm_homo <- as.data.frame(subset(diff_normal_homoschedastic_data, select = vector_of_anova_pvalue))
		matrix_diff_norm_homo <- cbind(rownames(matrix_diff_norm_homo), matrix_diff_norm_homo, method_diff_norm_homo)
		colnames(matrix_diff_norm_homo) <- c("Signal","pvalue","Test")
		# Differentially expressed signals: normal heteroschedastic data
		diff_normal_heteroschedastic_data <- as.data.frame(subset(normal_heteroschedastic_data, vector_of_welch_pvalue <= pvalue_expression))
		method_diff_norm_hetero <- rep("Welch", nrow(diff_normal_heteroschedastic_data))
		matrix_diff_norm_hetero <- as.data.frame(subset(diff_normal_heteroschedastic_data, select = vector_of_welch_pvalue))
		matrix_diff_norm_hetero <- cbind(rownames(matrix_diff_norm_hetero), matrix_diff_norm_hetero, method_diff_norm_hetero)
		colnames(matrix_diff_norm_hetero) <- c("Signal","pvalue","Test")
		# Differentially expressed signals: non_normal homoschedastic data
		diff_non_normal_homoschedastic_data <- as.data.frame(subset(non_normal_homoschedastic_data, vector_of_kruskal_pvalue <= pvalue_expression))
		method_diff_non_norm_homo <- rep("Kruskal-Wallis", nrow(diff_non_normal_homoschedastic_data))
		matrix_diff_non_norm_homo <- as.data.frame(subset(diff_non_normal_homoschedastic_data, select = vector_of_kruskal_pvalue))
		matrix_diff_non_norm_homo <- cbind(rownames(matrix_diff_non_norm_homo), matrix_diff_non_norm_homo, method_diff_non_norm_homo)
		colnames(matrix_diff_non_norm_homo) <- c("Signal","pvalue","Test")	
		# Differentially expressed signals: non-normal heteroschedastic data
		diff_non_normal_heteroschedastic_data <- as.data.frame(subset(non_normal_heteroschedastic_data, vector_of_permutation_pvalue <= pvalue_expression))
		method_diff_non_norm_hetero <- rep("Permutation", nrow(diff_non_normal_heteroschedastic_data))
		matrix_diff_non_norm_hetero <- as.data.frame(subset(diff_non_normal_heteroschedastic_data, select = vector_of_permutation_pvalue))
		matrix_diff_non_norm_hetero <- cbind(rownames(matrix_diff_non_norm_hetero), matrix_diff_non_norm_hetero, method_diff_non_norm_hetero)
		colnames(matrix_diff_non_norm_hetero) <- c("Signal","pvalue","Test")
		########## POST-HOC TESTS
		### ANOVA POST-HOC	
		post_hoc_anova_list <- list()
		if (nrow(diff_normal_homoschedastic_data)){
			# For each signal...
			for (m in rownames(diff_normal_homoschedastic_data)){
				## Extract the column of the signal and the response variable...
				mass_x <- temp_data_frame[[m]]
				response_variable <- temp_data_frame[[non_signal_variable]]
				### Clean the NA values from the mass column and the response variable column
				response_variable <- response_variable[!is.na(mass_x)]
				mass_x <- mass_x[!is.na(mass_x)]
				### Anova post hoc analysis: Tukey HSD test
				anova_test <- aov(lm(formula = mass_x ~ response_variable))
				# Post Hoc Tukey list creation
				tukey_post <- TukeyHSD(anova_test)
				rn <- rownames(tukey_post$response_variable)
				anova_table <- cbind(rn, tukey_post$response_variable)
				post_hoc_anova_list[[m]] <- anova_table
			}
		}
		### KRUSKAL-WALLIS POST-HOC	
		post_hoc_kruskal_list <- list()
		if (nrow(diff_non_normal_homoschedastic_data)){
			# For each signal...
			for (m in rownames(diff_non_normal_homoschedastic_data)){
				## Extract the column of the signal and the response variable...
				mass_x <- temp_data_frame[[m]]
				response_variable <- temp_data_frame[[non_signal_variable]]
				### Clean the NA values from the mass column and the response variable column
				response_variable <- response_variable[!is.na(mass_x)]
				mass_x <- mass_x[!is.na(mass_x)]
				### Kruskal post hoc analysis: Nemenyi-Damico-Wolfe-Dunn test - coin required!
				signal_data_frame <- data.frame(mass_x, response_variable)
				# Post hoc test
				NDWD <- oneway_test(mass_x ~ response_variable, data = signal_data_frame, ytrafo = function(data) trafo(data, numeric_trafo = rank), xtrafo = function(data) trafo(data, factor_trafo = function(mass_x) model.matrix(~mass_x - 1) %*% t(contrMat(table(mass_x), "Tukey"))), teststat = "max", distribution = approximate(B = 90000))
				post_hoc_pvalue <- pvalue(NDWD, method = "single-step")
				rn <- rownames(post_hoc_pvalue)
				post_hoc_table <-  cbind(rn, post_hoc_pvalue)
				post_hoc_kruskal_list[[m]] <- post_hoc_table
			}
		}
		### WELCH POST-HOC	
		post_hoc_welch_list <- list()
		if (nrow(diff_normal_heteroschedastic_data)) {
			# For each signal...
			for (m in rownames(diff_normal_heteroschedastic_data)) {
				## Extract the column of the signal and the response variable...
				mass_x <- temp_data_frame[[m]]
				response_variable <- temp_data_frame[[non_signal_variable]]
				### Clean the NA values from the mass column and the response variable column
				response_variable <- response_variable[!is.na(mass_x)]
				mass_x <- mass_x[!is.na(mass_x)]
				### Welch post hoc analysis: Duncan-Waller post-hoc test - agricolae required!
				degrees_of_freedom <- df.residual(lm(mass_x ~ response_variable))
				MSerror <- deviance(lm(mass_x ~ response_variable))/degrees_of_freedom
				Fc <- anova(lm(mass_x ~ response_variable))[1,4]
				comparison <- waller.test(mass_x, response_variable, degrees_of_freedom, MSerror, Fc, group = FALSE)
				rn <- rownames(as.data.frame(comparison[4]))
				ndf <- cbind(as.data.frame(comparison[4]), rn)
				post_hoc_welch_list[[m]] <- ndf
			}
		}
		### PERMUTATION POST-HOC	
		post_hoc_permutation_list <- list()
		if (nrow(diff_non_normal_heteroschedastic_data)) {	
			# For each signal...
			for (m in rownames(diff_non_normal_heteroschedastic_data)) {
				## Extract the column of the signal and the response variable...
				mass_x <- temp_data_frame[[m]]
				response_variable <- temp_data_frame[[non_signal_variable]]
				### Clean the NA values from the mass column and the response variable column
				response_variable <- response_variable[!is.na(mass_x)]
				mass_x <- mass_x[!is.na(mass_x)]
				### Permutation post hoc analysis: still follow Nemenyi-Damico-Wolfe-Dunn test - coin required!
				signal_data_frame <- data.frame(mass_x, response_variable)
				
				NDWD <- oneway_test(mass_x ~ response_variable, data = signal_data_frame, ytrafo = function(data) trafo(data, numeric_trafo = rank), xtrafo = function(data) trafo(data, factor_trafo = function(mass_x) model.matrix(~mass_x - 1) %*% t(contrMat(table(mass_x), "Tukey"))), teststat = "max", distribution = approximate(B = 90000))
				
				post_hoc_pvalue  <- pvalue(NDWD, method = "single-step")
				rn <- rownames(post_hoc_pvalue)
				post_hoc_table <-  cbind(rn, post_hoc_pvalue)
				post_hoc_permutation_list[[m]] <- post_hoc_table
			}
		}
		##### Selected signals for inference
		selected_signals_for_inference <- c(rownames(diff_normal_homoschedastic_data), rownames(diff_normal_heteroschedastic_data), rownames(diff_non_normal_homoschedastic_data), rownames(diff_non_normal_heteroschedastic_data))
		# Print the message...
		print("Differentially expressed signals for inference")
		print(selected_signals_for_inference)
		# Generate a matrix with the signals for inference and the method used to identify them...
		inference_signals_method_matrix <- rbind(matrix_diff_norm_homo, matrix_diff_norm_hetero, matrix_diff_non_norm_homo, matrix_diff_non_norm_hetero)
		# Select the signals to be used for inference and the non-signal variable
		selected_signals_for_inference_intensity_df <- subset(temp_data_frame, select = c(non_signal_variable, selected_signals_for_inference))
		########### Dump the files
		# Create the folder where to dump the files and go to it...
		subfolder <- paste(combination_vector[comb], "multi-levels")
		dir.create(file.path(output_folder, subfolder))
		setwd(file.path(output_folder, subfolder))
		# For each signals of inference...
		for (s in selected_signals_for_inference) {
			# Extract the intensity
			signal_intensity <- selected_signals_for_inference_intensity_df[[s]]
			# Extract the non-signal variable as an ordered factor
			non_signal_as_ordered_factor <- ordered(temp_data_frame[, non_signal_variable])
			# Remove possible NA values
			non_signal_as_ordered_factor <- non_signal_as_ordered_factor[!is.na(signal_intensity)]
			signal_intensity <- signal_intensity[!is.na(signal_intensity)]
			# Number the observations
			IDs <- c(1:length(signal_intensity))
			# Generate a matrix with the ID, the intensities of the selected signal for inference and the non-signal variable
			signal_dataframe <- data.frame(IDs, signal_intensity, non_signal_as_ordered_factor)
			##### Jitter plot	
			plot_name <- sprintf("%s%s", non_signal_variable,"_VS_Intensity")
			file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
			jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5))
			ggsave(jitter_plot, file = file_name, width = 4, height = 4)
			##### Box plot
			plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
			file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
			box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot")
			ggsave(box_plot, file = file_name, width = 4, height = 4)
			##### Scatter plot
			plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
			# Sort the dataframe rows according to the values of the non-signal variable
			ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
			ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
			colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
			file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
			graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
			scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors)
			ggsave(scatter_plot, file = file_name , width = 4, height = 4)
		}
		### Dump the files
		write_file(file_name = "Patient number matrix", data = patient_number_matrix_multi, file_format = file_format)
		write_file(file_name = "Outliers", data = outlier_matrix_multi, file_format = file_format)
		write_file(file_name = "Methods", data = inference_signals_method_matrix, file_format = file_format)
		write_file(file_name = "Differentially expressed signals", data = selected_signals_for_inference_intensity_df, file_format = file_format)
		write_file(file_name = "All signals", data = temp_data_frame_original, file_format = file_format)
		write_file(file_name = "Test p-values", data = pvalue_matrix_multi, file_format = file_format)
		if (nrow(diff_normal_homoschedastic_data)){
			write_posthoc_file(file_name = "PostHoc ANOVA", data = post_hoc_anova_list, file_format = file_format)
		}
		if (nrow(diff_non_normal_homoschedastic_data)) {
			write_posthoc_file(file_name = "PostHoc Kruskal-Wallis", data = post_hoc_kruskal_list, file_format = file_format)
		}
		if (nrow(diff_normal_heteroschedastic_data)) {
			write_posthoc_file(file_name = "PostHoc Welch", data = post_hoc_welch_list, file_format = file_format)
		}
		if(nrow(diff_non_normal_heteroschedastic_data)) {
			write_posthoc_file(file_name = "PostHoc Permutation", data = post_hoc_permutation_list, file_format = file_format)
		}
		# Go back to the output folder
		setwd(output_folder)
		# Restore the original data frame
		temp_data_frame <- temp_data_frame_original
	}
}








##### FINISH
finish <- tkmessageBox(title = "Operation completed", message = "All the statistical operations have been performed and the files have been dumped", icon = "info")
