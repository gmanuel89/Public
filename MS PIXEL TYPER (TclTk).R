#################### FUNCTIONS - MASS SPECTROMETRY 2017.04.04 ####################

########################################################################## MISC

###################################################### CHECK INTERNET CONNECTION
# This function checks if there is internet connection, by pinging a website. It returns TRUE or FALSE.
# Two methods are available: 'ping' tries to ping the website, while 'getURL' connects to it directly. The default is 'getURL', since it is more reliable than ping.
check_internet_connection <- function(method = "getURL", website_to_ping = "www.google.it") {
    ##### Start with getURL...
    there_is_internet <- FALSE
    ##### GET URL
    if (method == "getURL") {
        try({
            # Install the RCurl package if not installed
            if ("RCurl" %in% installed.packages()[,1]) {
                library(RCurl)
            } else {
                install.packages("RCurl", repos = "http://cran.mirror.garr.it/mirrors/CRAN/", quiet = TRUE, verbose = FALSE)
                library(RCurl)
            }
        }, silent = TRUE)
        there_is_internet <- FALSE
        try({
            there_is_internet <- is.character(getURL(u = website_to_ping, followLocation = TRUE, .opts = list(timeout = 1, maxredirs = 2, verbose = FALSE)))
            }, silent = TRUE)
    }
    ##### If getURL failed... Go back to ping (which should never fail)
    ##### PING
    if (method == "ping" || there_is_internet == FALSE) {
        if (Sys.info()[1] == "Linux") {
            # -c: number of packets sent/received (attempts) ; -W timeout in seconds
            there_is_internet <- !as.logical(system(command = paste("ping -c 1 -W 2", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        } else if (Sys.info()[1] == "Windows") {
            # -n: number of packets sent/received (attempts) ; -w timeout in milliseconds
            there_is_internet <- !as.logical(system(command = paste("ping -n 1 -w 2000", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        } else {
            there_is_internet <- !as.logical(system(command = paste("ping", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        }
    }
    return(there_is_internet)
}





################################################################################





##################################################### INSTALL REQUIRED PACKAGES
# This function installs and loads the selected packages
install_and_load_required_packages <- function(required_packages, repository = "http://cran.mirror.garr.it/mirrors/CRAN/", update_packages = FALSE) {
    ### Check internet connection
    there_is_internet <- check_internet_connection(method = "getURL", website_to_ping = "www.google.it")
    ########## Update all the packages (if there is internet connection)
    if (update_packages == TRUE) {
        if (there_is_internet == TRUE) {
            ##### If a repository is specified
            if (repository != "" || !is.null(repository)) {
                update.packages(repos = repository, ask = FALSE, checkBuilt = TRUE, quiet = TRUE, verbose = FALSE)
            } else {
                update.packages(ask = FALSE, checkBuilt = TRUE, quiet = TRUE, verbose = FALSE)
            }
            print("Packages updated")
        } else {
            print("Packages cannot be updated due to internet connection problems")
        }
    }
    ##### Retrieve the installed packages
    installed_packages <- installed.packages()[,1]
    ##### Determine the missing packages
    missing_packages <- character()
    for (p in 1:length(required_packages)) {
        if ((required_packages[p] %in% installed_packages) == FALSE) {
            missing_packages <- append(missing_packages, required_packages[p])
        }
    }
    ##### If there are packages to install...
    if (length(missing_packages) > 0) {
        ### If there is internet...
        if (there_is_internet == TRUE) {
            ### If a repository is specified
            if (repository != "" || !is.null(repository)) {
                install.packages(missing_packages, repos = repository, quiet = TRUE, verbose = FALSE)
            } else {
                ### If NO repository is specified
                install.packages(missing_packages, quiet = TRUE, verbose = FALSE)
            }
            print("All the required packages have been installed")
        } else {
            ### If there is NO internet...
            print("Some packages cannot be installed due to internet connection problems")
        }
    } else {
        print("All the required packages are installed")
    }
    ##### Load the packages (if there are all the packages)
    if ((length(missing_packages) > 0 && there_is_internet == TRUE) || length(missing_packages) == 0) {
        for (i in 1:length(required_packages)) {
            library(required_packages[i], character.only = TRUE)
        }
    } else {
        print("Packages cannot be installed/loaded... Expect issues...")
    }
}





###############################################################################





##################################################### R WORKSPACE DATA RETRIEVER
# This function takes an RData file path as input and returns the list of variables in the R workspace and the matrix with the type (R class) for each variable.
R_workspace_data_retriever <- function(filepath_R) {
    ### Create a temporary environment
    input_R_workspace <- new.env()
    ### Load the workspace
    load(filepath_R, envir = input_R_workspace)
    ### Extract the variable list (sorted)
    variable_list <- sort(ls(name = input_R_workspace))
    ### Extract the class list (the type of each variable)
    class_list <- character()
    ### If there are variables...
    if (length(variable_list) > 0) {
        # For each variable in the input workspace...
        for (v in 1:length(variable_list)) {
            ### Extract each variable to determine the class
            # Load the variable in it
            variable_x <- get(variable_list[v], pos = input_R_workspace)
            # Retrieve the class
            class_x <- as.character(class(variable_x)[1])
            # Add this to the final vector
            class_list <- append(class_list, class_x)
        }
        ## Generate the output matrix
        # Fill the matrix
        output_matrix <- cbind(variable_list, class_list)
        colnames(output_matrix) <- c("Variable name", "Variable type")
    }
    ### Return
    return(list(variable_list = variable_list, variable_matrix = output_matrix))
}





###############################################################################





########################################################### ENSEMBLE VOTE MATRIX
# The function takes as input the result matrix of an ensemble classification: each row is an observation/spectrum (patient or pixel) and each column is the predicted class of that observation by one model.
# The function returns a single column matrix with the ensemble classification results computed according to the input parameters (such as vote weights and method).
ensemble_vote_classification <- function(classification_matrix, class_list = NULL, decision_method = "majority", vote_weights = "equal") {
    ### Class list
    # Retrieve the class list according to the present classes (if not specified):
    if (is.null(class_list) || length(class_list) == 0) {
        # Initialize the class vector
        class_vector <- character()
        # Fill the class vector with the classification matrix columns
        for (cl in 1:ncol(classification_matrix)) {
            class_vector <- append(class_vector, as.character(classification_matrix[, cl]))
        }
        # Convert it into a factor
        class_vector <- as.factor(class_vector)
        # Extract the levels
        class_list <- levels(class_vector)
    }

    ########## Vote
    ##### Majority vote
    if (decision_method == "majority" && vote_weights == "equal") {
        # Function for matrix apply (x = row)
        majority_vote_function <- function(x, class_list) {
            # Generate the vote vector (same length as the class list, with the number of the votes for each class, labeled)
            votes <- integer(length = length(class_list))
            names(votes) <- class_list
            # Count the votes for each class
            for (class in class_list) {
                votes[which(class_list == class)] <- length(which(x == class))
            }
            # Determine the final majority vote
            final_vote <- names(votes)[which(votes == max(votes))]
            # Even vote
            if (length(final_vote) != 1) {
                final_vote <- NA
            }
            # Return the vote
            return(final_vote)
        }
        # For each spectrum (matrix row), establish the final majority vote
        classification_ensemble_matrix <- cbind(apply(X = classification_matrix, MARGIN = 1, FUN = function(x) majority_vote_function(x, class_list)))
        colnames(classification_ensemble_matrix) <- "Ensemble classification"
    }
    return(classification_ensemble_matrix)
}





###############################################################################





############################### ADD THE CLASS AND THE SAMPLE NAME TO THE MATRIX
# This function adds two column to the peaklist matrix (rows: spectra/patients, columns: aligned peaks): Sample and Class, according to the file name.
### The name of the rows will be either the sample name or the class name (depending on the function parameter).
# If the rows are named according to the sample name, an additional column for the class is added
matrix_add_class_and_sample <- function(signal_matrix, peaks = list(), class_list = list(), spectra_format = "imzml", sample_output = TRUE, class_output = TRUE, row_labels = "Sample") {
    # Convert the input matrix/dataframe into a matrix
    if (!is.matrix(signal_matrix)) {
        signal_matrix <- as.matrix(signal_matrix)
    }
    # Determine the number of spectra/peaklists
    if (isMassPeaksList(peaks)) {
        number_of_spectra <- length(peaks)
    } else if (isMassPeaks(peaks)) {
        number_of_spectra <- 1
    }
    ####################################### PATH and FILE VECTORS
    ##### Path vector
    # Create the empty vector
    path_vector <- character()
    # Add the file names recursively, scrolling the whole spectral dataset
    if (isMassPeaksList(peaks)) {
        for (i in 1:length(peaks)) {
            if (spectra_format == "imzml" || spectra_format == "imzML") {
                path_vector <- append(path_vector, peaks[[i]]@metaData$file[1])
            } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
                path_vector <- append(path_vector, peaks[[i]]@metaData$sampleName[1])
            }
        }
    } else if (isMassPeaks(peaks)) {
        if (spectra_format == "imzml" || spectra_format == "imzML") {
            path_vector <- append(path_vector, peaks@metaData$file[1])
        } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
            path_vector <- append(path_vector, peaks@metaData$sampleName[1])
        }
    }
    ##### File vector
    # Replace the path with the sample name in the peaks list
    peaks <- replace_sample_name(peaks, spectra_format = spectra_format)
    # Create the empty vector
    file_vector <- character()
    # Add the file names recursively, scrolling the whole spectral dataset
    if (isMassPeaksList(peaks)) {
        for (i in 1:length(peaks)) {
            if (spectra_format == "imzml" || spectra_format == "imzML") {
                file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
            } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
                file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName[1])
            }
        }
    } else if (isMassPeaks(peaks)) {
        if (spectra_format == "imzml" || spectra_format == "imzML") {
            file_vector <- append(file_vector, peaks@metaData$file[1])
        } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
            file_vector <- append(file_vector, peaks@metaData$sampleName[1])
        }
    }
    ################################## Only the sample
    if ((class_output == FALSE && sample_output == TRUE) || (class_output == TRUE && length(class_list) == 0 && sample_output == TRUE)) {
        # Create the sample matrix column and append it to the global matrix
        # Sample as rownames
        if ("sample" %in% row_labels || "Sample" %in% row_labels) {
            rownames(signal_matrix) <- file_vector
        } else {
            sample_column <- matrix("", ncol = 1, nrow = number_of_spectra)
            colnames(sample_column) <- "Sample"
            sample_column[,1] <- cbind(file_vector)
            signal_matrix <- cbind(signal_matrix, sample_column)
        }
    }
    ################################## Both the class and the sample
    if (class_output == TRUE && length(class_list) >= 1 && sample_output == TRUE) {
        # Create the sample matrix column and append it to the global matrix
        # Sample as rownames
        if ("sample" %in% row_labels || "Sample" %in% row_labels) {
        rownames(signal_matrix) <- file_vector
        } else {
            sample_column <- matrix("", ncol = 1, nrow = number_of_spectra)
            colnames(sample_column) <- "Sample"
            sample_column[,1] <- cbind(file_vector)
            signal_matrix <- cbind(signal_matrix, sample_column)
        }
        ### Add the class column
        class_list <- sort(class_list)
        # Rename the classes according to the class_list vector (the match should be /class/ to avoid catching the name of the class in previous folders)
        class_vector <- path_vector
        for (p in 1:length(class_vector)) {
            for (w in 1:length(class_list)) {
                if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                    if (length(grep(paste("/", class_list[w], "/", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
                        class_vector[p] <- class_list[w]
                    }
                } else if (Sys.info()[1] == "Windows") {
                    if (length(grep(paste("\\\\", class_list[w], "\\\\", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
                        class_vector[p] <- class_list[w]
                    }
                }
            }
        }
        # Class as rownames
        if ("class" %in% row_labels || "Class" %in% row_labels) {
            rownames(signal_matrix) <- class_vector
        } else {
            class_column <- matrix("", ncol = 1, nrow = number_of_spectra)
            colnames(class_column) <- "Class"
            # Fill in the matrix column with the file_vector classes and samples
            class_column[,1] <- cbind(class_vector)
            signal_matrix <- cbind(signal_matrix, class_column)
        }
    }
    ################################## Only the class
    if (class_output == TRUE && length(class_list) >= 1 && sample_output == FALSE) {
        class_list <- sort(class_list)
        # Rename the classes according to the class_list vector (the match should be /class/ to avoid catching the name of the class in previous folders)
        class_vector <- path_vector
        for (p in 1:length(class_vector)) {
            for (w in 1:length(class_list)) {
                if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                    if (length(grep(paste("/", class_list[w], "/", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
                        class_vector[p] <- class_list[w]
                    }
                } else if (Sys.info()[1] == "Windows") {
                    if (length(grep(paste("\\\\", class_list[w], "\\\\", sep = ""), class_vector[p], ignore.case = TRUE)) > 0) {
                        class_vector[p] <- class_list[w]
                    }
                }
            }
        }
        if ("class" %in% row_labels || "Class" %in% row_labels) {
            rownames(signal_matrix) <- class_vector
        } else {
            ### Add the class column
            class_column <- matrix("", ncol = 1, nrow = number_of_spectra)
            colnames(class_column) <- "Class"
            # Fill in the matrix column with the file_vector classes and samples
            class_column[,1] <- cbind(class_vector)
            signal_matrix <- cbind(signal_matrix, class_column)
        }
    }
    ### Add these matrix columns to the peaklist matrix
    return(signal_matrix)
}





################################################################################





###################################### FROM CLASS AND OUTCOME TO NUMBERS FOR MSI
# This function takes a list (character vector) of classes and a list (character vector) of corresponding outcomes and returns a dataframe with the classes, the outcomes and a number corresponding to the outcomes (0.5 = benign, 1= malignant), to be used to replace the intensities in the spectra list for plotting purposes.
# The outcome list must correspond to the class list, otherwise the outcomes are matched to the first elements of the class list, while the remaining class list elements are linked to the outcome 'other'; or the other outcomes are turned into a class.
# The outcome list can be: benign/ben/b - malignant/mal/m - other/oth/o
# If a vector is provided, the function will return also the same vector with the classes converted in numbers according to the outcome (benign = 0.5, malignant = 1)
outcome_and_class_to_MS <- function(class_list = c("HP", "PTC"), outcome_list = c("b", "m"), class_vector = NULL) {
    # Extract the unique values in case of duplicates
    class_list <- unique(class_list)
    outcome_list <- unique(outcome_list)
    # Fix the outcome names to a universal name (benign, malignant, other)
    for (ou in 1:length(outcome_list)) {
        if (length(grep("ben", outcome_list[ou])) > 0 || outcome_list[ou] == "b") {
            outcome_list[ou] <- "benign"
        } else if (length(grep("mal", outcome_list[ou])) > 0 || outcome_list[ou] == "m") {
            outcome_list[ou] <- "malignant"
        } else {
            outcome_list[ou] <- "other"
        }
    }
    # Establish the greater number of elements (class or outcome)
    if (length(class_list) >= length(outcome_list)) {
        greater_number <- length(class_list)
    } else {
        greater_number <- length(outcome_list)
    }
    # Get all the vectors to have the same dimension
    if (length(class_list) > length(outcome_list)) {
        for (i in 1:abs(length(class_list) - length(outcome_list))) {
            outcome_list <- append(outcome_list, "other")
        }
    } else if (length(class_list) < length(outcome_list)) {
        for (i in 1:abs(length(outcome_list) - length(class_list))) {
            class_list <- append(class_list, outcome_list[length(class_list) + i])
        }
    }
    ### Generate the matrix (classes)
    class_outcome_matrix <- matrix("", nrow = greater_number, ncol = 4)
    # Define the colnames
    colnames(class_outcome_matrix) <- c("Class", "Outcome", "Number", "Color")
    # Fill the matrix
    class_outcome_matrix[,1] <- cbind(class_list)
    class_outcome_matrix[,2] <- cbind(outcome_list)
    ### Generate the numbers
    outcome_list_as_number <- outcome_list
    for (ou in 1:length(outcome_list)) {
        # Benign = 0.5 (green pixels)
        if (length(grep("ben", outcome_list[ou])) > 0 || outcome_list[ou] == "b") {
            outcome_list_as_number[ou] <- 0.5
        } else if (length(grep("mal", outcome_list[ou])) > 0 || outcome_list[ou] == "m") {
        # Malignant = 1 (red pixels)
            outcome_list_as_number[ou] <- 1
        } else if (is.na(outcome_list[ou])) {
        # Other cases = 0 (black pixels)
            outcome_list_as_number[ou] <- 0
        } else {
        # Other cases = 0 (black pixels)
            outcome_list_as_number[ou] <- 0
        }
    }
    # Fill in the matrix column (outcome as number + NA)
    class_outcome_matrix[,3] <- cbind(outcome_list_as_number)
    ### Generate the colors
    outcome_list_as_color <- outcome_list
    for (ou in 1:length(outcome_list)) {
        # Benign = 0.5 (green pixels)
        if (length(grep("ben", outcome_list[ou])) > 0 || outcome_list[ou] == "b") {
            outcome_list_as_color[ou] <- "green"
        } else if (length(grep("mal", outcome_list[ou])) > 0 || outcome_list[ou] == "m") {
            # Malignant = 1 (red pixels)
            outcome_list_as_color[ou] <- "red"
        } else if (is.na(outcome_list[ou])) {
            # Other cases = 0 (black pixels)
            outcome_list_as_color[ou] <- "black"
        } else {
            # Other cases = 0 (black pixels)
            outcome_list_as_color[ou] <- "black"
        }
    }
    # Fill in the matrix column (outcome as number + NA)
    class_outcome_matrix[,4] <- cbind(outcome_list_as_color)
    ### Convert the class vector (if not null)
    if (!is.null(class_vector)) {
        for (cv in 1:length(class_vector)) {
            for (ou in 1:length(class_outcome_matrix[, "Class"])) {
                if (is.na(class_vector[cv])) {
                    class_vector[cv] <- as.numeric(0)
                } else if (class_vector[cv] == class_outcome_matrix[, "Class"][ou]) {
                    class_vector[cv] <- as.numeric(class_outcome_matrix[, "Number"][ou])
                }
            }
        }
    }
    # Convert the final vector into numeric
    class_vector <- as.numeric(class_vector)
    ###
    ### Return
    return(list(class_outcome_matrix = class_outcome_matrix, class_vector_as_numeric = class_vector, legend_text = class_list, legend_fill = outcome_list_as_color))
}





###################################### ADD THE THY CLASS TO THE MATRIX (THYROID)
# This function adds the THY column to the peaklist matrix, by reading the THY value from the sample name (imzML)
matrix_add_thy <- function(signal_matrix, peaks, spectra_format = "imzml") {
    number_of_spectra <- length(peaks)
    # Create the empty vector
    file_vector <- character()
    # Add the file names recursively, scrolling the whole spectral dataset
    for (i in 1:length(peaks)) {
        if (spectra_format == "imzml" || spectra_format == "imzML") {
            file_vector <- append(file_vector, peaks[[i]]@metaData$file[1])
        }
        if (spectra_format == "brukerflex" || spectra_format == "xmass") {
            file_vector <- append(file_vector, peaks[[i]]@metaData$sampleName)
        }
    }
    # Create the thy vector, symmetrical to the file vector
    thy_vector <- file_vector
    # Find the "THY" in the file name and store its value
    for (t in 1:length(thy_vector)) {
        ############## THY
        # THY1
        if (length(grep("THY1", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 1
        }
        # THY2
        if (length(grep("THY2", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 2
        }
        # THY3
        if (length(grep("THY3", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 3
        }
        # THY4
        if (length(grep("THY4", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 4
        }
        # THY5
        if (length(grep("THY5", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 5
        }
        ############## TIR
        # TIR1
        if (length(grep("TIR1", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 1
        }
        # TIR2
        if (length(grep("TIR2", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 2
        }
        # TIR3
        if (length(grep("TIR3", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 3
        }
        # TIR4
        if (length(grep("TIR4", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 4
        }
        # TIR5
        if (length(grep("TIR5", thy_vector[t], ignore.case = TRUE)) == 1) {
            thy_vector[t] <- 5
        }
    }
    # Create the sample matrix column and appendit to the global matrix
    thy_column <- matrix (0, ncol = 1, nrow = number_of_spectra)
    colnames(thy_column) <- "THY"
    # Fill in the matrix thy column with the thy_vector and attach it to the matrix
    thy_column [,1] <- cbind(thy_vector)
    signal_matrix <- cbind(signal_matrix, thy_column)
    # Return
    return(signal_matrix)
}





################################################################################





##################################################### REMOVE LOW INTENSITY PEAKS
# This function removes low-intensity peaks (in terms of level of intensity compared with the most intense peak in the peaklist) from the list of provided peaks (MALDIquant).
# If the method is selected to be "element-wise", each element of the peaklist is evaluated, and the intensity threshold is calculated over the peaks of only that element. Otherwise, if "whole" is selected, the threshold is calculated on all the peaks in the dataset.
remove_low_intensity_peaks <- function(peaks, low_intensity_peak_removal_threshold_percent = 0.1, low_intensity_peak_removal_threshold_method = "element-wise", allow_parallelization = FALSE) {
    ### Load the required libraries
    install_and_load_required_packages(c("parallel", "MALDIquant"))
    ### Fix the percentage value
    if (low_intensity_peak_removal_threshold_percent < 0) {
        low_intensity_peak_removal_threshold_percent <- 0
    } else if (low_intensity_peak_removal_threshold_percent > 100) {
        low_intensity_peak_removal_threshold_percent <- 100
    }
    ### If there is only one peaklist, there is no point in doing the 'whole' method, but only the element-wise.
    if (!isMassPeaksList(peaks)) {
        low_intensity_peak_removal_threshold_method <- "element-wise"
    }
    ########## Do everything only if there is a reasonable value of the percentage
    if (low_intensity_peak_removal_threshold_percent > 0 && low_intensity_peak_removal_threshold_percent < 100) {
        ########## ELEMENT-WISE
        if (low_intensity_peak_removal_threshold_method == "element-wise") {
            ##### INTENSITY FILTERING FUNCTION (ELEMENT-WISE)
            intensity_filtering_subfunction_element <- function(peaks, low_intensity_peak_removal_threshold_percent) {
                # Filter out the peaks whose intensity is below a certain threshold
                # Store mass and intensity into vectors
                intensity_values <- peaks@intensity
                mass_values <- peaks@mass
                snr_values <- peaks@snr
                # Identify the positions of the values to be discarded
                values_to_be_discarded <- intensity_values[((intensity_values * 100 / max(intensity_values, na.rm = TRUE)) < low_intensity_peak_removal_threshold_percent)]
                # If there are values to be discarded...
                if (length(values_to_be_discarded) > 0) {
                    # Identify the positions
                    positions_to_be_discarded <- numeric()
                    for (i in 1:length(values_to_be_discarded)) {
                        value_position <- which(intensity_values == values_to_be_discarded[i])
                        positions_to_be_discarded <- append(positions_to_be_discarded, value_position)
                    }
                    # Discard the values from the vectors
                    intensity_values <- intensity_values [-positions_to_be_discarded]
                    mass_values <- mass_values[-positions_to_be_discarded]
                    snr_values <- snr_values[-positions_to_be_discarded]
                    # Put the values back into the MALDIquant list
                    peaks@mass <- mass_values
                    peaks@intensity <- intensity_values
                    peaks@snr <- snr_values
                } else {
                    # If there aren't any values to be discarded...
                    peaks@mass <- mass_values
                    peaks@intensity <- intensity_values
                    peaks@snr <- snr_values
                }
                return(peaks)
            }
            ##### Multiple peaks elements
            if (isMassPeaksList(peaks)) {
                ### MULTICORE
                if (allow_parallelization == TRUE) {
                    # Detect the number of cores
                    cpu_thread_number <- detectCores(logical = TRUE)
                    cpu_thread_number <- cpu_thread_number / 2
                    if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                        peaks_filtered <- mclapply(peaks, FUN = function (peaks) intensity_filtering_subfunction_element(peaks, low_intensity_peak_removal_threshold_percent), mc.cores = cpu_thread_number)
                    } else if (Sys.info()[1] == "Windows") {
                        # Make the CPU cluster for parallelisation
                        cl <- makeCluster(cpu_thread_number)
                        # Make the cluster use the custom functions and the package functions along with their parameters
                        clusterEvalQ(cl, {library(MALDIquant)})
                        # Pass the variables to the cluster for running the function
                        clusterExport(cl = cl, varlist = c("peaks", "low_intensity_peak_removal_threshold_percent", "intensity_filtering_subfunction_element"), envir = environment())
                        # Apply the multicore function
                        peaks_filtered <- parLapply(cl, peaks, fun = function (peaks) intensity_filtering_subfunction_element(peaks, low_intensity_peak_removal_threshold_percent))
                        stopCluster(cl)
                    } else {
                        peaks_filtered <- lapply(peaks, FUN = function (peaks) intensity_filtering_subfunction_element(peaks, low_intensity_peak_removal_threshold_percent))
                    }
                } else {
                    ### SINGLE CORE
                    peaks_filtered <- lapply(peaks, FUN = function (peaks) intensity_filtering_subfunction_element(peaks, low_intensity_peak_removal_threshold_percent))
                }
            } else {
                ##### Single peaks element
                peaks_filtered <- intensity_filtering_subfunction_element(peaks, low_intensity_peak_removal_threshold_percent)
            }
        }
        ########## WHOLE DATASET
        if (low_intensity_peak_removal_threshold_method == "whole") {
            ##### INTENSITY FILTERING FUNCTION (ELEMENT-WISE)
            intensity_filtering_subfunction_whole <- function(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity) {
                # Filter out the peaks whose intensity is below a certain threshold
                # Store mass and intensity into vectors
                intensity_values <- peaks@intensity
                mass_values <- peaks@mass
                snr_values <- peaks@snr
                # Identify the positions of the values to be discarded
                values_to_be_discarded <- intensity_values[((intensity_values * 100 / highest_intensity) < low_intensity_peak_removal_threshold_percent)]
                # If there are values to be discarded...
                if (length(values_to_be_discarded) > 0) {
                    # Identify the positions
                    positions_to_be_discarded <- numeric()
                    for (i in 1:length(values_to_be_discarded)) {
                        value_position <- which(intensity_values == values_to_be_discarded[i])
                        positions_to_be_discarded <- append(positions_to_be_discarded, value_position)
                    }
                    # Discard the values from the vectors
                    intensity_values <- intensity_values [-positions_to_be_discarded]
                    mass_values <- mass_values[-positions_to_be_discarded]
                    snr_values <- snr_values[-positions_to_be_discarded]
                    # Put the values back into the MALDIquant list
                    peaks@mass <- mass_values
                    peaks@intensity <- intensity_values
                    peaks@snr <- snr_values
                } else {
                    # If there aren't any values to be discarded...
                    peaks@mass <- mass_values
                    peaks@intensity <- intensity_values
                    peaks@snr <- snr_values
                }
                return(peaks)
            }
            ### Determine the highest peak in the dataset
            highest_peak <- NULL
            highest_intensity <- 0
            for (p in 1:length(peaks)) {
                if (length(peaks[[p]]@mass) > 0) {
                    for (m in 1:length(peaks[[p]]@mass)) {
                        if (highest_intensity == 0 || peaks[[p]]@intensity[m] > highest_intensity) {
                            highest_intensity <- peaks[[p]]@intensity[m]
                            highest_peak <- peaks[[p]]@mass[m]
                        }
                    }
                }
            }
            ### Filter the peaks
            ##### Multiple peaks elements
            if (isMassPeaksList(peaks)) {
                ### MULTICORE
                if (allow_parallelization == TRUE) {
                    # Detect the number of cores
                    cpu_thread_number <- detectCores(logical = TRUE)
                    cpu_thread_number <- cpu_thread_number / 2
                    if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                        peaks_filtered <- mclapply(peaks, FUN = function(peaks) intensity_filtering_subfunction_whole(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity), mc.cores = cpu_thread_number)
                    } else if (Sys.info()[1] == "Windows") {
                        # Make the CPU cluster for parallelisation
                        cl <- makeCluster(cpu_thread_number)
                        # Make the cluster use the custom functions and the package functions along with their parameters
                        clusterEvalQ(cl, {library(MALDIquant)})
                        # Pass the variables to the cluster for running the function
                        clusterExport(cl = cl, varlist = c("peaks", "low_intensity_peak_removal_threshold_percent", "intensity_filtering_subfunction_whole", "highest_intensity"), envir = environment())
                        # Apply the multicore function
                        peaks_filtered <- parLapply(cl, peaks, fun = function(peaks) intensity_filtering_subfunction_whole(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity))
                        stopCluster(cl)
                    } else {
                        peaks_filtered <- lapply(peaks, FUN = function(peaks) intensity_filtering_subfunction_whole(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity))
                    }
                } else {
                    ### SINGLE CORE
                    peaks_filtered <- lapply(peaks, FUN = function(peaks) intensity_filtering_subfunction_whole(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity))
                }
            } else {
                ########## Single peaks element
                peaks_filtered <- intensity_filtering_subfunction_whole(peaks, low_intensity_peak_removal_threshold_percent, highest_intensity)
            }
        }
        return(peaks_filtered)
    } else {
        return(peaks)
    }
}





################################################################################





########################################################### SPECTRA FILES READER
# This function reads all the files from a folder and returns only the imzML files or the fid files.
read_spectra_files <- function(folder, spectra_format = "imzml", full_path = TRUE) {
    if (spectra_format == "imzml" || spectra_format == "imzML") {
        # Read all the files
        folder_all_files <- list.files(folder, full.names = full_path, recursive = TRUE)
        # Create the empty vector in which only the imzML files will be listed
        spectra_files <- character()
        # Put the imzML files in the new vector discarding the ibd files
        for (l in 1:length(folder_all_files)) {
            if (length(grep(".imzML", folder_all_files[l], fixed = TRUE)) == 1) {
                spectra_files <- append(spectra_files, folder_all_files[l])
            }
        }
    } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # Read all the folders
        sample_directory_list_all_folders <- dir(folder, ignore.case = TRUE, full.names = full_path, recursive = TRUE, include.dirs = TRUE)
        # Take only the fid paths (spectra)
        spectra_files <- character()
        for (i in 1:length(sample_directory_list_all_folders)) {
            if ((length(grep("fid", sample_directory_list_all_folders[i], ignore.case = TRUE)) == 1)) {
                spectra_files <- append(spectra_files, sample_directory_list_all_folders[i])
            }
        }
    }
    return(spectra_files)
}





################################################################################





######################### REPLACE THE SNR WITH THE STDEV IN THE PEAKS
# This function computes the standard deviation of each peak of an average spectrum peaklist, by replacing the existing SNR slot with the SD or CV: all the peaks (average and dataset) are aligned and each peak of the average peaklist is searched across the dataset thanks to the intensity matrix.
replace_SNR_in_avg_peaklist <- function (spectra, peak_picking_algorithm = "SuperSmoother", SNR = 5, tof_mode = "linear", tolerance_ppm = 2000, spectra_format = "imzml", replace_snr_with = "std") {
    install_and_load_required_packages(c("MALDIquant", "stats"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Check if it is not a single spectrum
    # If the spectra are many...
    if (length(spectra) > 0 && isMassSpectrumList(spectra)) {
        # Average the spectra
        avg_spectrum <- averageMassSpectra(spectra, method = "mean")
        avg_spectrum <- removeBaseline(avg_spectrum, method = "TopHat")
        avg_spectrum <- calibrateIntensity(avg_spectrum, method = "TIC")
        # Peak picking
        peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
        peaks_avg <- peak_picking(avg_spectrum, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
        # Merge for the alignment
        global_peaks <- append(peaks_avg, peaks)
        global_peaks <- align_and_filter_peaks(global_peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
        peaks_avg <- global_peaks[[1]]
        peaks <- global_peaks [2:length(global_peaks)]
        # Compute the intensity matrix
        intensity_matrix <- intensityMatrix(peaks, spectra)
        # Scroll the peaks of the average peaklist...
        for (p in 1:length(peaks_avg@mass)) {
            # Search for the column in the intensity_matrix
            for (z in 1:ncol(intensity_matrix)) {
                # Match...
                if (peaks_avg@mass[p] == colnames(intensity_matrix)[z]) {
                    # Create an empty vector where to allocate the intensity of this peak into the dataset
                    intensity_vector <- intensity_matrix[,z]
                    # Compute the standard deviation
                    if (replace_snr_with == "std" || replace_snr_with == "sd") {
                        std_int <- sd(intensity_vector)
                        # Replacement
                        peaks_avg@snr[p] <- std_int
                    }
                    if (replace_snr_with == "cv" || replace_snr_with == "CV") {
                        std_int <- sd(intensity_vector)
                        mean_int <- mean(intensity_vector)
                        cv_int <- std_int / mean_int
                        # Replacement
                        peaks_avg@snr[p] <- cv_int
                    }
                    # Save time avoiding to scroll to the end when found
                    break
                } else {peaks_avg@snr[p] <- NaN}
            }
        }
        return (peaks_avg)
    } else {
        return(print("The spectra list is empty or the spectra list contains only one spectrum"))
    }
}





################################################################################





############################################################### OUTLIERS REMOVAL
# This function takes a vector on values as input and returns the same vector without the outliers, calculated based upon the interquartile range (fence rule).
# The outliers can be replaced with nothing (so they are removed from the vector) or with some value
outliers_removal <- function (v, replace_with = "") {
    summary_vector <- summary(v)
    # Calculate the interquartile range
    inter_quartile_range <- summary_vector[5] - summary_vector[2]
    # Calculate the fences, beyond which the spectrum is an outlier
    iqr_fences <- c((summary_vector[2] - 1.5*inter_quartile_range), (summary_vector[5] + 1.5*inter_quartile_range))
    # Find the outliers based on the fences condition
    outliers_position <- which(v < iqr_fences[1] | v > iqr_fences[2])
    # If the outliers have to be discarded...
    if (replace_with == "") {
        if (length(outliers_position) > 0) {
            # Remove the correspondent elements from the dataset
            v <- v [-outliers_position]
        }
    }
    # If the outliers have to be replaced...
    if (is.numeric(replace_with)) {
        if (length(outliers_position) > 0) {
            # Replace the outliers with the replacement
            v [outliers_position] <- replace_with
        }
    }
    if (replace_with == 0 || replace_with == "zero") {
        if (length(outliers_position) > 0) {
            # Replace the outliers with the replacement
            v [outliers_position] <- 0
        }
    }
    if (replace_with == "NA" || is.na(replace_with)) {
        if (length(outliers_position) > 0) {
            # Replace the outliers with the replacement
            v [outliers_position] <- NA
        }
    }
    if (replace_with == "mean") {
        if (length(outliers_position) > 0) {
            # Remove the correspondent elements from the dataset
            vector_no_outliers <- v [-outliers_position]
            # Replace the outliers with the vector mean (no outliers)
            v [outliers_position] <- mean(vector_no_outliers)
        }
    }
    if (replace_with == "median") {
        if (length(outliers_position) > 0) {
            # Remove the correspondent elements from the dataset
            vector_no_outliers <- v [-outliers_position]
            # Replace the outliers with the vector median (no outliers)
            v [outliers_position] <- median(vector_no_outliers)
        }
    }
    return (list (vector = v, outliers_position = outliers_position))
}





################################################################################





######################################### PEAK STATISTICS (on processed Spectra)
# This function computes the peak statistics onto a selected spectra dataset (or to the provided peaks), both when the spectra belong to no (or one) class and more classes.
# It returns a NULL value if the peak statistics cannot be performed.
peak_statistics <- function(spectra, peaks = NULL, SNR = 3, peak_picking_algorithm = "SuperSmoother", class_list = NULL, class_in_file_name = TRUE, tof_mode = "linear", spectra_format = "imzml", exclude_spectra_without_peak = FALSE, alignment_iterations = 5, peak_filtering_frequency_threshold_percent = 25, remove_outliers = TRUE, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element_wise") {
    ########## Load the required libraries
    install_and_load_required_packages(c("MALDIquant", "stats"))
    ########## Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ########## Define the tolerance in PPM
    if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
        tolerance_ppm <- 2000
    } else if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
        tolerance_ppm <- 200
    }
    ########## Determine the number of classes
    if (length(class_list) == 0 || length(class_list) == 1 || is.null(class_list)) {
    number_of_classes <- 1
    } else if (length(class_list) > 1) {
        number_of_classes <- length(class_list)
    }
    ########## Detect (if not already provided) and Align Peaks
    if (is.null(peaks)) {
        peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
    }
    peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, alignment_iterations = alignment_iterations, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, reference_peaklist = NULL, spectra = spectra)
    # Generate the matrix (and convert it into a data frame)
    if (exclude_spectra_without_peak == FALSE) {
        signal_matrix <- intensityMatrix(peaks, spectra)
    } else if (exclude_spectra_without_peak == TRUE) {
        signal_matrix <- intensityMatrix(peaks)
    }
    # Peak vector
    if (is.matrix(signal_matrix)) {
        peak_vector <- as.numeric(colnames(signal_matrix))
    } else if (is.data.frame(signal_matrix)) {
        peak_vector <- as.numeric(names(signal_matrix))
    }
    ############################################################## ONE CLASS
    if (number_of_classes == 1) {
        ################################# FUNCTION for matrix APPLY (it will applied for each matrix column, for each peak)
        peak_statistcs_function <- function (signal_matrix_column, signal_matrix, remove_outliers) {
            # Generate the output matrix row
            peak_stat_matrix_row <- matrix (0, nrow = 1, ncol = 7)
            rownames(peak_stat_matrix_row) <- as.numeric(colnames(signal_matrix_column))
            colnames(peak_stat_matrix_row) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles")
            # Start the calculation
            intensity_vector <- as.numeric(signal_matrix_column)
            if (remove_outliers == TRUE) {
                intensity_vector <- outliers_removal(intensity_vector)$vector
            }
            # Calculate the statistical parameters on the intensity values in the vector
            # Normality
            if (length(intensity_vector) >= 3 & length(intensity_vector) <= 5000) {
                shapiro_test <- shapiro.test(intensity_vector)
                if (shapiro_test$p.value < 0.05) {
                    distribution_type <- paste("Non-normal", "(Shapiro p-value:", round(shapiro_test$p.value,3), ")")
                }
                if (shapiro_test$p.value >= 0.05) {
                    distribution_type <- paste("Normal", "(Shapiro p-value:", round(shapiro_test$p.value,3), ")")
                }
            } else if (length(intensity_vector) < 3) {
                distribution_type <- "Not determinable, number of samples too low"
            } else if (length(intensity_vector) > 5000) {
                distribution_type <- "Number of samples too high, assume it is normal"
            }
            # Other parameters
            st_dev_intensity <- sd(intensity_vector, na.rm = TRUE)
            summary_intensity_vector <- summary(intensity_vector)
            mean_intensity <- summary_intensity_vector [4]
            coeff_variation <- (st_dev_intensity / mean_intensity) *100
            median_intensityensity <- summary_intensity_vector [3]
            first_quartile <- summary_intensity_vector [2]
            third_quartile <- summary_intensity_vector [5]
            inter_quartile_range <- third_quartile - first_quartile
            # Fill the matrix with the values
            peak_stat_matrix_row [,1] <- distribution_type
            peak_stat_matrix_row [,2] <- as.numeric(mean_intensity)
            peak_stat_matrix_row [,3] <- as.numeric(st_dev_intensity)
            peak_stat_matrix_row [,4] <- as.numeric(coeff_variation)
            peak_stat_matrix_row [,5] <- as.numeric(median_intensityensity)
            peak_stat_matrix_row [,6] <- as.numeric(inter_quartile_range)
            peak_stat_matrix_row [,7] <- paste("1st quartile", first_quartile, "; 3rd quartile", third_quartile)
            return (peak_stat_matrix_row)
        }
        ###############
        # Fix the signal_matrix (Add the sample column)
        signal_matrix <- matrix_add_class_and_sample(signal_matrix, peaks = peaks, spectra_format = spectra_format, sample_output = TRUE, class_output = FALSE)
        # Output matrix
        peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-1), ncol = 8)
        rownames(peak_stat_matrix) <- as.numeric(colnames(signal_matrix)[1:(ncol(signal_matrix)-1)])
        colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Sample")
        # Only peaks
        signal_matrix_peaks <- signal_matrix [,1:(ncol(signal_matrix)-1)]
        # Apply the function (transpose the result matrix)
        peak_stat_matrix <- t(apply(signal_matrix_peaks, MARGIN = 2, FUN = function(x) peak_statistcs_function(x, signal_matrix, remove_outliers = remove_outliers)))
        # Generate the intensity matrix with NA if the peak is not present in the spectra
        intensity_matrix_with_na <- intensityMatrix(peaks)
        spectra_counter_vector <- numeric()
        for (pk in 1:ncol(intensity_matrix_with_na)) {
            intensity_vector <- intensity_matrix_with_na[,pk]
            spectra_counter_vector <- append(spectra_counter_vector, length(intensity_vector[!is.na(intensity_vector)]))
        }
        peak_stat_matrix <- cbind(peak_stat_matrix, spectra_counter_vector)
        # Fix the column names
        colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median",  "Interquartile Range (IQR)", "Quartiles", "Spectra counter")
        ## Return
        return(peak_stat_matrix)
    } else if (number_of_classes > 1) {
        ############################################################ TWO OR MORE CLASSES
        # Every variable now is a list, each element of which corresponds to a certain value from a class
        # So every variable is a list with the same length of the class list (each element of the list
        # is referred to a class
        # Fix the signal_matrix (Add the sample column)
        signal_matrix <- matrix_add_class_and_sample(signal_matrix, peaks = peaks, class_list = class_list, spectra_format = spectra_format, sample_output = TRUE, class_output = TRUE)
        # Check if there is a sufficient number of observations per class
        observations_per_class <- numeric()
        for (i in 1:length(class_list)) {
            observations_per_class <- append(observations_per_class, length(which(signal_matrix[,ncol(signal_matrix)] == class_list[i])))
        }
        sufficient_number_of_observations_per_class <- TRUE
        for (i in 1:length(observations_per_class)) {
            if (observations_per_class[i] < 3) {
                sufficient_number_of_observations_per_class <- FALSE
            }
        }
        ##### Run only if there is a sufficient number of samples
        if (sufficient_number_of_observations_per_class == TRUE) {
            # Output matrix
            peak_stat_matrix <- matrix (0, nrow=(ncol(signal_matrix)-2), ncol = 14)
            rownames(peak_stat_matrix) <- as.numeric(colnames(signal_matrix)[1:(ncol(signal_matrix)-2)])
            colnames(peak_stat_matrix) <- c("Intensity distribution type", "Mean", "Standard deviation", "Coefficient of Variation %", "Median", "Interquartile Range (IQR)", "Spectra counter", "Class", "Homoscedasticity (parametric)", "Homoscedasticity (non-parametric)", "t-Test", "ANOVA", "Wilcoxon - Mann-Whitney test", "Kruskal-Wallis test")
            # For each peak
            for (p in 1:(ncol(signal_matrix)-2)) {
                # Put the intensity of that peak into one vector per class (in a global list)
                intensity_vector <- list()
                # Scroll the peaklists and Add the peak intensity to a vector(one for each class)
                for (l in 1:length(class_list)) {
                    # Allocate in the intensity vector the rows for that peak belonging to the certain class
                    intensity_vector[[l]] <- as.numeric(signal_matrix [signal_matrix[,ncol(signal_matrix)] == class_list[l],p])
                }
                if (remove_outliers == TRUE) {
                    for (i in 1:length(intensity_vector)) {
                        intensity_vector[[i]] <- outliers_removal(intensity_vector[[i]])
                        intensity_vector[[i]] <- intensity_vector[[i]]$vector
                    }
                }
                ######################## STATISTICAL PARAMETERS
                ############################################### Normality for each class
                shapiro_test <- list()
                distribution_type <- list()
                for (l in 1:length(class_list)) {
                    if (length(intensity_vector[[l]]) >= 3 && length(intensity_vector[[l]]) <= 5000) {
                        shapiro_test[[l]] <- shapiro.test(intensity_vector[[l]])
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
                    variance_test_parametric <- var.test(intensity_vector[[1]], intensity_vector[[2]])
                }
                if (length(class_list) >= 2) {
                    variance_test_non_parametric <- bartlett.test(as.numeric(signal_matrix[,p]), g = as.factor(signal_matrix[,ncol(signal_matrix)]))
                }
                ########################################### Other parameters (per class)
                st_dev_intensity <- list()
                summary_intensity_vector <- list()
                mean_intensity <- list()
                coeff_variation <- list()
                median_intensityensity <- list()
                first_quartile <- list()
                third_quartile <- list()
                inter_quartile_range <- list()
                spectra_counter <- list()
                variance <- list()
                for (l in 1:length(class_list)) {
                    st_dev_intensity[[l]] <- sd(intensity_vector[[l]])
                    summary_intensity_vector [[l]] <- summary(intensity_vector[[l]])
                    mean_intensity[[l]] <- summary_intensity_vector[[l]] [4]
                    coeff_variation[[l]] <- (st_dev_intensity[[l]] / mean_intensity[[l]]) *100
                    median_intensityensity[[l]] <- summary_intensity_vector[[l]] [3]
                    first_quartile[[l]] <- summary_intensity_vector[[l]] [2]
                    third_quartile[[l]] <- summary_intensity_vector[[l]] [5]
                    inter_quartile_range[[l]] <- third_quartile[[l]] - first_quartile[[l]]
                    spectra_counter[[l]] <- length(intensity_vector[[l]])
                    variance[[l]] <- var(intensity_vector[[l]])
                }
                ############################################# Parameters between classes
                # T-test
                if (length(class_list) == 2) {
                    t_test <- t.test(intensity_vector[[1]], intensity_vector[[2]])
                }
                # ANOVA TEST
                if (length(class_list) >= 2) {
                anova_test <- aov(signal_matrix[,p] ~ signal_matrix[,ncol(signal_matrix)])
                }
                # WILCOXON - MANN-WHITNEY TEST
                if (length(class_list) == 2) {
                    wilcoxon_test <- wilcox.test(intensity_vector[[1]], intensity_vector[[2]])
                }
                # KRUSKAL-WALLIS TEST
                if (length(class_list) >= 2) {
                    kruskal_wallis_test <- kruskal.test(signal_matrix[,p], g = as.factor(signal_matrix[,ncol(signal_matrix)]))
                }
                ######################################## Fill the matrix with the values
                # Distribution Type
                distribution_type_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(distribution_type_name)) {
                        distribution_type_name <- paste(distribution_type[[l]], " - ", class_list[l], sep = "")
                    } else {
                        distribution_type_name <- paste(distribution_type_name, " , ", distribution_type[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,1] <- paste(distribution_type_name)
                # Mean
                mean_intensity_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(mean_intensity_name)) {
                        mean_intensity_name <- paste(mean_intensity[[l]], " - ", class_list[l], sep = "")
                    } else {
                        mean_intensity_name <- paste(mean_intensity_name, " , ", mean_intensity[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,2] <- mean_intensity_name
                # Standard Deviation
                st_dev_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(st_dev_name)) {
                        st_dev_name <- paste(st_dev_intensity[[l]], " - ", class_list[l], sep = "")
                    } else {
                        st_dev_name <- paste(st_dev_name, " , ", st_dev_intensity[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,3] <- st_dev_name
                # Coefficient of Variation
                coeff_variation_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(coeff_variation_name)) {
                        coeff_variation_name <- paste(coeff_variation[[l]], " - ", class_list[l], sep = "")
                    } else {
                        coeff_variation_name <- paste(coeff_variation_name, " , ", coeff_variation[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,4] <- coeff_variation_name
                # Median
                median_intensityensity_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(median_intensityensity_name)) {
                        median_intensityensity_name <- paste(median_intensityensity[[l]], " - ", class_list[l], sep = "")
                    } else {
                        median_intensityensity_name <- paste(median_intensityensity_name, " , ", median_intensityensity[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,5] <- median_intensityensity_name
                # Interquartile Range (IQR)
                inter_quartile_range_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(inter_quartile_range_name)) {
                        inter_quartile_range_name <- paste(inter_quartile_range[[l]], " - ", class_list[l], sep = "")
                    } else {
                        inter_quartile_range_name <- paste(inter_quartile_range_name, " , ", inter_quartile_range[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,6] <- inter_quartile_range_name
                # Spectra counter
                spectra_counter_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(spectra_counter_name)) {
                        spectra_counter_name <- paste(spectra_counter[[l]], " - ", class_list[l], sep = "")
                    } else {
                        spectra_counter_name <- paste(spectra_counter_name, " , ", spectra_counter[[l]], " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,7] <- spectra_counter_name
                # Class
                class_name <- NULL
                for (l in 1:length(class_list)) {
                    if (is.null(class_name)) {
                        class_name <- class_list[l]
                    } else {
                        class_name <- paste(class_name, " - ", class_list[l], sep = "")
                    }
                }
                peak_stat_matrix [p,8] <- class_name
                # Homoscedasticity (Parametric)
                if (variance_test_parametric$p.value < 0.05) {
                homoscedasticity_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
                } else if (variance_test_parametric$p.value >= 0.05) {
                homoscedasticity_parametric <- paste("Homoscedastic data", "(p-value:", variance_test_parametric$p.value, ")")
                }
                if (variance_test_non_parametric$p.value < 0.05) {
                homoscedasticity_non_parametric <- paste("Non homoscedastic data", "(p-value:", variance_test_non_parametric$p.value, ")")
                } else if (variance_test_non_parametric$p.value >= 0.05) {
                homoscedasticity_non_parametric <- paste("Homoscedastic data", "(p-value:", variance_test_non_parametric$p.value, ")")
                }
                peak_stat_matrix [p,9] <- homoscedasticity_parametric
                peak_stat_matrix [p,10] <- homoscedasticity_non_parametric
                # t-Test
                peak_stat_matrix [p,11] <- t_test$p.value
                # ANOVA
                peak_stat_matrix [p,12] <- summary(anova_test)[[1]]$"Pr(>F)"[1]
                # Wilcoxon / Mann-Whitney test
                peak_stat_matrix [p,13] <- wilcoxon_test$p.value
                # Kruskal-Wallis test
                peak_stat_matrix [p,14] <- kruskal_wallis_test$p.value
            }
            ## Return
            return(peak_stat_matrix)
        } else {
            ## Return NULL
            return(NULL)
        }
    }
}





################################################################################


























































































############################################# SPECTRA







################################################################ SPECTRA BINNING
# The function performs the binning onto a selected spectra dataset (list of MALDIquant spectra objects)
resample_spectra <- function(spectra, final_data_points = lowest_data_points, binning_method = "sum", allow_parallelization = FALSE) {
    ####################################################### BINNING FUNCTION
    binning_subfunction <- function(spectra, final_data_points, binning_method) {
        # Create the new spectra_binned list
        spectra_binned <- spectra
        # Empty the mass and intensity values
        for (s in 1:length(spectra_binned)) {
            spectra_binned@mass <- numeric()
            spectra_binned@intensity <- numeric()
        }
        # Calculate the number of datapoints per bin
        data_points_per_bin <- floor(length(spectra@mass) / final_data_points)
        # Define the indexes
        index1 <- 1
        index2 <- data_points_per_bin
        # For each bin...
        for (d in 1:final_data_points) {
            # Create the (temporary) bin vectors (mass and intensity), where the data points will be stored for the binning
            bin_mass <- numeric()
            bin_intensity <- numeric()
            # Scroll the data points, grouping them by bins
            for (i in index1:index2) {
                bin_mass <- append(bin_mass, spectra@mass[i])
                bin_intensity <- append(bin_intensity, spectra@intensity[i])
            }
            # Calculate the value of the new data point
            data_point_mass <- mean(bin_mass)
            if (binning_method == "sum") {
                data_point_intensity <- sum(bin_intensity)
            } else if (binning_method == "mean") {
                data_point_intensity <- mean(bin_intensity)
            } else if (binning_method == "median") {
                data_point_intensity <- median(bin_intensity)
            } else if (binning_method == "RMS" | binning_method == "rms") {
                data_point_intensity <- sqrt(sum(bin_intensity^2))
            }
            # Store it in the new spectra_binned list
            spectra_binned@mass <- append(spectra_binned@mass, data_point_mass)
            spectra_binned@intensity <- append(spectra_binned@intensity, data_point_intensity)
            # Increase the indexes
            index1 <- index1 + data_points_per_bin
            index2 <- index2 + data_points_per_bin
        }
        return(spectra_binned)
    }
    ############# More spectra
    if (isMassSpectrumList(spectra)) {
        # Load the required libraries
        install_and_load_required_packages("parallel")
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        ########################
        # Calculate the lowest amount of data points, that corresponds to the maximum
        # number of data points that can be used for the binning
        # Use this value as default if final_data_points is not specified (Function)
        lowest_data_points <- NULL
        # Compare it with all the others
        for (s in 1:length(spectra)) {
            if (is.null(lowest_data_points) || length(spectra[[s]]@mass) < lowest_data_points) {
                lowest_data_points <- length(spectra[[s]]@mass)
            }
        }
        # Check the qualities of the spectral dataset
        data_points_table <- table(sapply(spectra, length))
        datapoints_dataset <- as.numeric(names(data_points_table))
        equality_data_points <- length(data_points_table)
        ######## Do not bin if all the spectra are of the same lengthand have the same number of datapoints as defined
        if ((equality_data_points == 1) && (datapoints_dataset == final_data_points)) {
            spectra_binned <- spectra
        }
        ######## Perform the binning if the spectra are not of the same lengthor
        # they are of the same lengthbut with a different number of datapoints than
        # the desired one
        if (!(equality_data_points == 1) || ((equality_data_points == 1) && !(datapoints_dataset == final_data_points))) {
            # Do the binning only if the number of final data points is lower than the
            # lowest number of original data points
            if (final_data_points > lowest_data_points) {
                final_data_points <- lowest_data_points
                print("Binning at this sample rate is not possible, the highest number of data points possible will be used")
                if (allow_parallelization == TRUE) {
                    if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                        spectra_binned <- mclapply(spectra, FUN = function (spectra) binning_subfunction(spectra, final_data_points, binning_method), mc.cores = cpu_thread_number)
                    } else if (Sys.info()[1] == "Windows") {
                        cl <- makeCluster(cpu_thread_number)
                        # Pass the variables to the cluster for running the function
                        clusterExport(cl = cl, varlist = c("final_data_points", "binning_method"), envir = environment())
                        spectra_binned <- parLapply(cl, spectra, fun = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                        stopCluster(cl)
                    } else {
                        spectra_binned <- lapply(spectra, FUN = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                    }
                } else {
                    spectra_binned <- lapply(spectra, FUN = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                }
            }
            if (final_data_points <= lowest_data_points) {
                if (allow_parallelization == TRUE) {
                    if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                        spectra_binned <- mclapply(spectra, fun = function (spectra) binning_subfunction(spectra, final_data_points, binning_method), mc.cores = cpu_thread_number)
                    } else if (Sys.info()[1] == "Windows") {
                        cl <- makeCluster(cpu_thread_number)
                        # Pass the variables to the cluster for running the function
                        clusterExport(cl = cl, varlist = c("final_data_points", "binning_method"), envir = environment())
                        spectra_binned <- parLapply(cl, spectra, fun = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                        stopCluster(cl)
                    } else {
                        spectra_binned <- lapply(spectra, fun = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                    }
                } else {
                    spectra_binned <- lapply(spectra, fun = function (spectra) binning_subfunction(spectra, final_data_points, binning_method))
                }
            }
        }
        print(table(sapply(spectra_binned, length)))
        print(paste("Equal distance between datapoints", (all(sapply(spectra_binned, isRegular)))))
    } else {
        # Retrieve the number of datapoints
        lowest_data_points <- length(spectra@mass)
        # Perform the binning only if necessary
        if (final_data_points < lowest_data_points) {
            spectra_binned <- binning_subfunction(spectra, final_data_points, binning_method)
        } else {spectra_binned <- spectra}
    }
    return(spectra_binned)
}





################################################################################





################################################# SAMPLE NAME REPLACING
# This function replaces the sample name field in the spectrum with the actual sample name (keeping only the last part of the file path and discarding the folder tree)
# The input can be both spectra or peaks (MALDIquant)
replace_sample_name <- function(spectra, spectra_format = "imzml", allow_parallelization = FALSE) {
    ##### Function for lapply
    name_replacing_subfunction <- function(spectra, spectra_format) {
        ### imzML
        if (spectra_format == "imzml" || spectra_format == "imzML") {
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                # Split the filepath at /
                sample_name <- unlist(strsplit(spectra@metaData$file[1],"/"))
                # The sample name is the last part of the path
                sample_name <- sample_name[length(sample_name)]
                # Detach the file extension
                sample_name <- unlist(strsplit(sample_name, ".imzML"))
                sample_name <- sample_name[1]
                # Put the name back into the spectra
                spectra@metaData$file <- sample_name
            } else if (Sys.info()[1] == "Windows") {
                # Split the filepath at \
                sample_name <- unlist(strsplit(spectra@metaData$file[1],"\\\\"))
                # The sample name is the last part of the path
                sample_name <- sample_name[length(sample_name)]
                # Detach the file extension
                sample_name <- unlist(strsplit(sample_name, ".imzML"))
                sample_name <- sample_name[1]
                # Put the name back into the spectra
                spectra@metaData$file <- sample_name
            } else {
                # Keep the name unaltered
                sample_name <- spectra@metaData$file[1]
            }
        } else if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        ### Xmass
            sample_name <- spectra@metaData$sampleName[1]
            # Put the name back into the spectra
            spectra@metaData$file <- sample_name
        }
        ### Return
        return(spectra)
    }
    ##### More elements
    if (isMassSpectrumList(spectra) || isMassPeaksList(spectra)) {
        ##### Apply the function
        if (allow_parallelization == TRUE) {
            # Detect the number of cores
            cpu_thread_number <- detectCores(logical = TRUE)
            cpu_thread_number <- cpu_thread_number / 2
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                spectra <- mclapply(spectra, FUN = function(spectra) name_replacing_subfunction(spectra, spectra_format = spectra_format), mc.cores = cpu_thread_number)
            } else if (Sys.info()[1] == "Windows") {
                # Make the CPU cluster for parallelisation
                cl <- makeCluster(cpu_thread_number)
                # Make the cluster use the custom functions and the package functions along with their parameters
                clusterEvalQ(cl, {library(MALDIquant)})
                # Pass the variables to the cluster for running the function
                clusterExport(cl = cl, varlist = c("spectra", "spectra_format"), envir = environment())
                # Apply the multicore function
                spectra <- parLapply(cl, spectra, fun = function(spectra) name_replacing_subfunction(spectra, spectra_format = spectra_format))
                stopCluster(cl)
            } else {
                spectra <- lapply(spectra, FUN = function(spectra) name_replacing_subfunction(spectra, spectra_format = spectra_format))
            }
        } else {
            spectra <- lapply(spectra, FUN = function(spectra) name_replacing_subfunction(spectra, spectra_format = spectra_format))
        }
    } else if (isMassSpectrum(spectra) || isMassPeaks(spectra)) {
        spectra <- name_replacing_subfunction(spectra, spectra_format = spectra_format)
    }
    ### Return
    return(spectra)
}




################################################################################





################################################# SAMPLE NAME REPLACING
# This function replaces the sample name field in the spectrum with the actual sample name (keeping only the last part of the file path and discarding the folder tree) and puts the sample name as the name of the list element, so that names(spectra) is not NULL and a better identification of the spectra is obtained (through a unique name rather than a simple number)
# The input can be both spectra or peaks (MALDIquant)
replace_sample_name_list <- function(spectra, spectra_format = "imzml") {
    ##### Replace the metaData$file with the name with the replace_sample_name() function
    spectra <- replace_sample_name(spectra, spectra_format = spectra_format, allow_parallelization = FALSE)
    ##### More spectra
    if (isMassSpectrumList(spectra) || isMassPeaksList(spectra)) {
        # Generate a vector with all the spectra names
        spectra_names_vector <- character()
        for (s in 1:length(spectra)) {
            spectra_names_vector <- append(spectra_names_vector, spectra[[s]]@metaData$file)
        }
        # Extract the unique values
        if (length(spectra_names_vector) != length(unique(spectra_names_vector))) {
            spectra_names_vector <- make.names(spectra_names_vector, unique = TRUE)
        }
        # Fix the spectral list names
        names(spectra) <- spectra_names_vector
    } else if (isMassSpectrum(spectra) || isMassPeaks(spectra)) {
        ##### One spectrum
        spectra <- spectra
    }
    ### Return
    return(spectra)
}




################################################################################





########################################################## CLASS NAME REPLACING
# This function replaces the sample name field in the spectrum with the class the spectrum belongs to. If the filename contains the class, it replaces the filename with the class name, otherwise a class list symmetrical to the spectra must be provided.
# The input can be both spectra or peaks (MALDIquant)
replace_class_name <- function (spectra, class_list = NULL, class_in_file_name = TRUE, spectra_format = "imzml") {
    # If a class list is provided...
    if (!is.null(class_list)) {
        # If the class name is in the name
        if (class_in_file_name == TRUE) {
            # Extract the file name from the path
            #spectra <- replace_sample_name(spectra, spectra_format = spectra_format)
            # Scroll the class list...
            for (w in 1:length(class_list)) {
                # Scroll the spectra...
                for (i in 1:length(spectra)) {
                    # If there is a match between the file name and the class
                    if (length(grep(class_list[w],spectra[[i]]@metaData$file[1], fixed = TRUE)) != 0) {
                        # Replace the file name with the class name
                        spectra[[i]]@metaData$file <- class_list [w]
                    }
                }
            }
        } else {
            # If the class name is not in the name
            # The class list and the spectra are symmetrical
            for (w in 1:length(class_list)) {
                spectra[[w]]@metaData$file <- class_list [w]
            }
        }
    }
    # If a class list is not provided...
    if (is.null(class_list)) {
        # It's guessed from the filepath
        if (class_in_file_name == TRUE) {
            #### imzML
            if (spectra_format == "imzml" || spectra_format == "imzML") {
                # Scroll the spectra...
                for (i in 1:length(spectra)) {
                    # Split the filepath at / (universal for different OSs)
                    class_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"/"))
                    # The sample name is the last part of the path
                    class_name <- class_name[(length(class_name)-1)]
                    # Put the name back into the spectra
                    spectra[[i]]@metaData$file <- class_name
                }
            }
            #### Xmass
            if (spectra_format == "brukerflex" || spectra_format == "xmass") {
                # Scroll the spectra...
                for (i in 1:length(spectra)) {
                    # Split the filepath at / (universal for different OSs)
                    class_name <- unlist(strsplit(spectra[[i]]@metaData$file[1],"/"))
                    # The sample name is the last part of the path
                    class_name <- class_name[length(class_name)-6]
                    # Put the name back into the spectra
                    spectra[[i]]@metaData$file <- class_name
                }
            }
        } else {
            # The worst case returns the list of spectra as it is.
            spectra <- spectra
        }
    }
    return (spectra)
}





###############################################################################





######################################################## SPECTRA PRE-PROCESSING
# The function runs the preprocessing on the selected spectra (smoothing, baseline subtraction and normalization). The function can be applied both to a spectra list or a single spectrum, allowing parallel computation.
# The function allows to select some additional parameters of the preprocessing.
# This version of the function whould be faster because each element of the spectral list is subjected to all the preprocessing step.
# If an algorithm is set to NULL, that preprocessing step will not be performed.
preprocess_spectra <- function(spectra, tof_mode = "linear", preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = NULL, smoothing_algorithm = "SavitzkyGolay", smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 100, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_algorithm = NULL, spectral_alignment_reference = "auto"), allow_parallelization = FALSE) {
    ##### Load the required libraries
    install_and_load_required_packages(c("MALDIquant", "parallel"))
    ##### Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ##### Extract the parameters from the input list
    mass_range <- preprocessing_parameters$mass_range
    transformation_algorithm <- preprocessing_parameters$transformation_algorithm
    smoothing_algorithm <- preprocessing_parameters$smoothing_algorithm
    smoothing_strength <- preprocessing_parameters$smoothing_strength
    baseline_subtraction_algorithm <- preprocessing_parameters$baseline_subtraction_algorithm
    baseline_subtraction_iterations <- preprocessing_parameters$baseline_subtraction_iterations
    normalization_algorithm <- preprocessing_parameters$normalization_algorithm
    normalization_mass_range <- preprocessing_parameters$normalization_mass_range
    preprocess_spectra_in_packages_of <- preprocessing_parameters$preprocess_spectra_in_packages_of
    spectral_alignment_algorithm <- preprocessing_parameters$spectral_alignment_algorithm
    spectral_alignment_reference <- preprocessing_parameters$spectral_alignment_reference
    ##### Fix the names
    if (!is.null(smoothing_algorithm) && (smoothing_algorithm == "SavitzkyGolay" || smoothing_algorithm == "Savitzky-Golay" || smoothing_algorithm == "SG")) {
        smoothing_algorithm <- "SavitzkyGolay"
    } else if (!is.null(smoothing_algorithm) && (smoothing_algorithm == "MovingAverage" || smoothing_algorithm == "Moving Average" || smoothing_algorithm == "MA")) {
        smoothing_algorithm <- "MovingAverage"
    }
    ##### Define the smoothing half wondow size
    smoothing_half_window_size <- NULL
    if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
        if (!is.null(smoothing_strength) && smoothing_strength == "medium") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 10
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 2
            }
        } else if (!is.null(smoothing_strength) && smoothing_strength == "strong") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 20
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 4
            }
        } else if (!is.null(smoothing_strength) && smoothing_strength == "stronger") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 30
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 6
            }
        }
    } else if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
        if (!is.null(smoothing_strength) && smoothing_strength == "medium") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 3
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 0.6
            }
        } else if (!is.null(smoothing_strength) && smoothing_strength == "strong") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 6
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 1.2
            }
        } else if (!is.null(smoothing_strength) && smoothing_strength == "stronger") {
            if (!is.null(smoothing_algorithm) && smoothing_algorithm == "SavitzkyGolay") {
                smoothing_half_window_size <- 9
            } else if (!is.null(smoothing_algorithm) && smoothing_algorithm == "MovingAverage") {
                smoothing_half_window_size <- 1.8
            }
        }
    }
    ##### Generate the preprocessing function to be applied to every element of the spectra_temp list (x = spectrum)
    preprocessing_subfunction <- function(x, mass_range, transformation_algorithm, smoothing_algorithm, smoothing_half_window_size, baseline_subtraction_algorithm, baseline_subtraction_iterations, normalization_algorithm, normalization_mass_range) {
        ### Trimming
        # Mass range specified
        if (!is.null(mass_range)) {
            x <- trim_spectra(x, range = mass_range)
        }
        ### Transformation
        if (!is.null(transformation_algorithm)) {
            x <- transformIntensity(x, method = transformation_algorithm)
        }
        ### Smoothing
        if (!is.null(smoothing_algorithm)) {
            x <- smoothIntensity(x, method = smoothing_algorithm, halfWindowSize = smoothing_half_window_size)
        }
        ### Baseline removal
        if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "TopHat") {
            x <- removeBaseline(x, method = baseline_subtraction_algorithm)
        } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "SNIP") {
            # Default value for the number of iterations
            if (is.null(baseline_subtraction_iterations) || baseline_subtraction_iterations <= 0) {
                baseline_subtraction_iterations <- 100
            }
            x <- removeBaseline(x, method = baseline_subtraction_algorithm, iterations = baseline_subtraction_iterations)
        }
        ### Normalization
        if (!is.null(normalization_algorithm)) {
            if (normalization_algorithm == "TIC") {
                if (!is.null(normalization_mass_range) && is.numeric(normalization_mass_range)) {
                    x <- calibrateIntensity(x, method = normalization_algorithm, range = normalization_mass_range)
                } else {
                    x <- calibrateIntensity(x, method = normalization_algorithm)
                }
            } else {
                x <- calibrateIntensity(x, method = normalization_algorithm)
            }
        }
        ### Return the preprocessed spectrum (x)
        return(x)
    }
    #################### Multiple spectra
    if (isMassSpectrumList(spectra)) {
        ##### Trimming (same mass range for all the dataset)
        if (is.null(mass_range)) {
            spectra <- trim_spectra(spectra)
        }
        ##### Preprocess in packages
        # Fix the packages variable
        if (is.null(preprocess_spectra_in_packages_of) || preprocess_spectra_in_packages_of <= 0 || preprocess_spectra_in_packages_of > length(spectra)) {
            preprocess_spectra_in_packages_of <- length(spectra)
        }
        # Create the list containing the processed spectra
        preprocessed_spectra <- list()
        # Define the number of folds
        number_of_folds <- as.integer(ceiling(length(spectra) / preprocess_spectra_in_packages_of))
        # Define to which fold the spectra belong to
        if (length(number_of_folds) > 1) {
            spectra_id_folds <- cut(seq(1, length(spectra)), breaks = number_of_folds, labels = FALSE)
        } else if (length(number_of_folds) == 1) {
            spectra_id_folds <- rep(1, length(spectra))
        }
        # For each fold...
        for (f in 1:number_of_folds) {
            spectra_temp <- spectra[which(spectra_id_folds == f)]
            # Apply the function to the list of spectra_temp
            if (allow_parallelization == TRUE) {
                # Detect the number of cores
                cpu_thread_number <- detectCores(logical = TRUE)
                cpu_thread_number <- cpu_thread_number / 2
                if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                    spectra_temp <- mclapply(spectra_temp, FUN = function(spectra_temp) preprocessing_subfunction(spectra_temp, mass_range = mass_range, transformation_algorithm = transformation_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_half_window_size = smoothing_half_window_size, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range), mc.cores = cpu_thread_number)
                } else if (Sys.info()[1] == "Windows") {
                    cl <- makeCluster(cpu_thread_number)
                    clusterEvalQ(cl, {library(MALDIquant)})
                    clusterExport(cl = cl, varlist = c("mass_range", "transformation_algorithm", "smoothing_algorithm", "smoothing_half_window_size", "baseline_subtraction_algorithm", "baseline_subtraction_iterations", "normalization_algorithm", "normalization_mass_range", "preprocessing_subfunction"), envir = environment())
                    spectra_temp <- parLapply(cl, spectra_temp, fun = function(spectra_temp) preprocessing_subfunction(spectra_temp, mass_range = mass_range, transformation_algorithm = transformation_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_half_window_size = smoothing_half_window_size, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range))
                    stopCluster(cl)
                } else {
                    spectra_temp <- lapply(spectra_temp, FUN = function(spectra_temp) preprocessing_subfunction(spectra_temp, mass_range = mass_range, transformation_algorithm = transformation_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_half_window_size = smoothing_half_window_size, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range))
                }
            } else {
                spectra_temp <- lapply(spectra_temp, FUN = function(spectra_temp) preprocessing_subfunction(spectra_temp, mass_range = mass_range, transformation_algorithm = transformation_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_half_window_size = smoothing_half_window_size, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range))
            }
            # Add to the final preprocessed spectral dataset
            preprocessed_spectra <- append(preprocessed_spectra, spectra_temp)
        }
    } else if (isMassSpectrum(spectra)) {
    #################### Single spectrum
        spectra <- preprocessing_subfunction(spectra, mass_range, transformation_algorithm, smoothing_algorithm, smoothing_half_window_size, baseline_subtraction_algorithm, baseline_subtraction_iterations, normalization_algorithm, normalization_mass_range)
        # Add to the final preprocessed spectral dataset
        preprocessed_spectra <- spectra
    }
    #################### SPECTRAL ALIGNMENT (Multiple spectra)
    if (!is.null(spectral_alignment_algorithm)) {
        if (isMassSpectrumList(preprocessed_spectra)) {
            spectral_alignment_performed <- FALSE
            try({
                # Perform the alignment
                preprocessed_spectra <- align_spectra(preprocessed_spectra, spectral_alignment_algorithm = spectral_alignment_algorithm, spectral_alignment_reference = spectral_alignment_reference, tof_mode = tof_mode)
                spectral_alignment_performed <- TRUE
            })
            # Return message
            if (spectral_alignment_performed == TRUE) {
                print("The spectral aligment has been performed successfully!")
            } else {
                print("The spectral aligment could not be performed!")
            }
        }
    }
    ########## Return preprocessed (and aligned) spectra
    return(preprocessed_spectra)
}





###############################################################################





############################################################# SPECTRAL ALIGNMENT
# The function allows to perform the alignment among spectra. By selecting 'auto', the reference peaks will be automatically determined; if 'average spectrum' is selected, the reference peaklist is the peaklist of the average spectrum; if 'skyline spectrum' is selected, the reference peaklist is the peaklist of the skyline spectrum.
align_spectra <- function(spectra, spectral_alignment_algorithm = "cubic", spectral_alignment_reference = "auto", tof_mode = "linear", deisotope_peaklist = FALSE) {
    # Perform only if a method is specified
    if (!is.null(spectral_alignment_algorithm)) {
        if (isMassSpectrumList(spectra)) {
            if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
                half_window_alignment <- 20
                tolerance_ppm <- 2000
            } else if (tof_mode == "reflector" || tof_mode == "reflectron" || tof_mode == "R") {
                half_window_alignment <- 5
                tolerance_ppm <- 200
            }
            ##### Perform the alignment: Automatic computation of the reference peaklist
            if (!is.null(spectral_alignment_reference) && spectral_alignment_reference == "auto") {
                aligned_spectra <- alignSpectra(spectra, noiseMethod = "SuperSmoother", halfWindowSize = half_window_alignment, SNR = 2, tolerance = (tolerance_ppm/10^6), warpingMethod = spectral_alignment_algorithm)
            } else if (!is.null(spectral_alignment_reference) && spectral_alignment_reference == "average spectrum") {
                ##### Perform the alignment: Average spectrum
                # Ganerate the average spectrum
                average_spectrum <- averageMassSpectra(spectra, method = "mean")
                # Preprocess the average spectrum
                if (tof_mode == "reflectron") {
                    smoothing_algorithm_avg <- NULL
                } else {
                    smoothing_algorithm_avg <- "SavitzkyGolay"
                }
                average_spectrum <- preprocess_spectra(average_spectrum, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = NULL, smoothing_algorithm = smoothing_algorithm_avg, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 200, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_algorithm = NULL, spectral_alignment_reference = "auto"))
                # Detect peaks onto the average spectrum: refference peaklist
                average_spectrum_peaks <- detectPeaks(average_spectrum, halfWindowSize = half_window_alignment, method = "SuperSmoother", SNR = 2)
                # Deisotope peaklist
                if (deisotope_peaklist == TRUE) {
                    average_spectrum_peaks <- deisotope_peaks(average_spectrum_peaks)
                }
                # Align the spectra
                aligned_spectra <- alignSpectra(spectra, noiseMethod = "SuperSmoother", halfWindowSize = half_window_alignment, SNR = 2, tolerance = (tolerance_ppm/10^6), warpingMethod = spectral_alignment_algorithm, reference = average_spectrum_peaks)
            } else if (!is.null(spectral_alignment_reference) && spectral_alignment_reference == "skyline spectrum") {
                ##### Perform the alignment: Skyline spectrum
                # Ganerate the skyline spectrum
                skyline_spectrum <- averageMassSpectra(spectra, method = "sum")
                # Preprocess the skyline spectrum
                if (tof_mode == "reflectron") {
                    smoothing_algorithm_skyline <- NULL
                } else {
                    smoothing_algorithm_skyline <- "SavitzkyGolay"
                }
                skyline_spectrum <- preprocess_spectra(skyline_spectrum, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = NULL, smoothing_algorithm = smoothing_algorithm_avg, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 200, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_algorithm = NULL, spectral_alignment_reference = "auto"))
                # Detect peaks onto the skyline spectrum: refference peaklist
                skyline_spectrum_peaks <- detectPeaks(skyline_spectrum, halfWindowSize = half_window_alignment, method = "SuperSmoother", SNR = 2)
                # Deisotope peaklist
                if (deisotope_peaklist == TRUE) {
                    skyline_spectrum_peaks <- deisotope_peaks(skyline_spectrum_peaks)
                }
                # Align the spectra
                aligned_spectra <- alignSpectra(spectra, noiseMethod = "SuperSmoother", halfWindowSize = half_window_alignment, SNR = 2, tolerance = (tolerance_ppm/10^6), warpingMethod = spectral_alignment_algorithm, reference = skyline_spectrum_peaks)
            }
            ### Return the aligned spectra
            return(aligned_spectra)
        } else {
            ### Return the original spectra if there is only one spectrum
            return(spectra)
        }
    } else {
        ### Return the original spectra if the alignment algorithm is not specified
        return(spectra)
    }
}





################################################################################





######################## PLOT THE MEAN SPECTRUM WITH THE SD BARS ON THE AVERAGE
# This function takes a list of spectra (MALDIquant) as input and returns the average spectrum of the provided dataset with bars onto the peaks, after calculating the standard deviation.
average_spectrum_bars <- function(spectra, SNR = 5, peak_picking_algorithm = "SuperSmoother", tolerance_ppm = 2000, mass_range_plot = c(4000,15000), graph_title = "Spectrum", average_spectrum_color = "black", peak_points = "yes", points_color = "red", bar_width = 40, bar_color = "blue") {
    # Load the required libraries
    install_and_load_required_packages("MALDIquant")
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Generate the average spectrum
    average_spectrum <- averageMassSpectra(spectra, method = "mean")
    average_spectrum <- removeBaseline(average_spectrum, method = "TopHat")
    # Peak picking on the average spectrum (for plotting)
    peaks_average <- peak_picking(average_spectrum, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
    # Peak picking on the dataset
    peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = SNR)
    # Alignment: merge the two lists and align them all
    peaks_all <- append(peaks, peaks_average)
    peaks_all <- align_and_filter_peaks(peaks_all, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 0)
    # Empty the lists
    peaks_average <- list()
    peaks <- list()
    # Re-fill the lists with the aligned peaklists
    for (i in 1:(length(peaks_all) - 1)) {
        peaks <- append(peaks, peaks_all[[i]])
    }
    peaks_average <- peaks_all[[length(peaks_all)]]
    ######## Generate a vector with the standard deviations of the peaks in the average peaklist
    st_dev_intensity <- vector(length = 0)
    # For each peak in the average peaklist
    for (a in 1:length(peaks_average@mass)) {
        intensity_vector <- vector(length = 0)
        # Search for it in the dataset peaklist
        for (p in 1:length(peaks)) {
            for (i in 1:length(peaks[[p]]@mass)) {
                if (abs(peaks[[p]]@mass[i] == peaks_average@mass[a])) {
                    intensity_vector <- append(intensity_vector, peaks[[p]]@intensity[i])
                }
            }
        }
        st_dev_intensity <- append(st_dev_intensity, sd(intensity_vector))
    }
    ####### Plot
    # Spectrum
    plot(average_spectrum, main = graph_title, col.main = average_spectrum_color, xlab = "m/z", ylab = "Intensity (a.i.)", xlim = mass_range_plot, col = average_spectrum_color)
    # Peaks
    if (peak_points == "yes") {
        points(peaks_average, pch = 4, col = points_color)
    }
    # Bars
    # epsilon: length of the horizontal segment
    epsilon = bar_width
    for (i in 1:length(peaks_average@mass)) {
        # Define the upper and the lower limit of the vertical segment
        up <- peaks_average@intensity[i] + st_dev_intensity [i]
        low <- peaks_average@intensity[i] - st_dev_intensity [i]
        # Vertical bar (x,y x,y)
        segments(peaks_average@mass[i], low, peaks_average@mass[i], up, col = bar_color)
        # Horizontal segments(x,y , x,y)
        segments(peaks_average@mass[i]-epsilon, low, peaks_average@mass[i]+epsilon, low, col = bar_color)
        segments(peaks_average@mass[i]-epsilon, up, peaks_average@mass[i]+epsilon, up, col = bar_color)
    }
    avg_spectrum_with_bars <- recordPlot()
    return(avg_spectrum_with_bars)
}





################################################################################





############################################# MOST INTENSE PEAKS IN PEAK PICKING
# This function returns a peak list containing only the most intense peaks per spectrum. If the input is a list of spectra, the function computes the peak picking and keeps only the most intense ones, if it's a list of peaklists, it applies the filtering function directly on the peaks.
most_intense_signals <- function(spectra, signals_to_take = 20, tof_mode = "linear", peak_picking_algorithm = "SuperSmoother", allow_parallelization = FALSE, deisotope_peaklist = FALSE) {
    # Load the required libraries
    install_and_load_required_packages(c("parallel", "MALDIquant"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ####################################################### PICKING FUNCTION
    picking_subfunction <- function (peaks, signals_to_take) {
        # Create a dataframe with mass and intensity
        peaks_data_frame <- data.frame (mass = peaks@mass, intensity = peaks@intensity, SNR = peaks@snr)
        # Check if the provided number does not exceed the number of available signals
        if (signals_to_take > nrow(peaks_data_frame) || signals_to_take <= 0) {
            signals_to_take <- nrow(peaks_data_frame)
        }
        # Sort the dataframe according to the SNR
        peaks_data_frame <- peaks_data_frame [order(-peaks_data_frame$SNR),]
        # Select only the first most intense signals
        selected_signals <- peaks_data_frame [1:signals_to_take,]
        # Sort the dataframe back according to mass
        selected_signals <- selected_signals [order(selected_signals$mass),]
        # Put these signals back into the peaklist
        peaks@mass <- selected_signals$mass
        peaks@intensity <- selected_signals$intensity
        peaks@snr <- selected_signals$SNR
        return (peaks)
    }
    ########################################################################
    # Peak picking
    if (isMassSpectrumList(spectra) || isMassSpectrum(spectra)) {
        peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = 3, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
    } else if (isMassPeaksList(spectra) || isMassPeaks(spectra)) {
        peaks <- spectra
    }
    # Most intense signals
    if (isMassPeaksList(peaks)) {
        if (allow_parallelization == TRUE) {
            # Detect the number of cores
            cpu_thread_number <- detectCores(logical = TRUE)
            cpu_thread_number <- cpu_thread_number / 2
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                most_intense_peaks <- mclapply(peaks, FUN = function(peaks) picking_subfunction(peaks, signals_to_take = signals_to_take), mc.cores = cpu_thread_number)
            } else if (Sys.info()[1] == "Windows") {
                # Make the CPU cluster for parallelisation
                cl <- makeCluster(cpu_thread_number)
                # Pass the variables to the cluster for running the function
                clusterExport(cl = cl, varlist = "signals_to_take", envir = environment())
                most_intense_peaks <- parLapply(cl, peaks, fun = function(peaks) picking_subfunction(peaks, signals_to_take = signals_to_take))
                stopCluster(cl)
            } else {
                most_intense_peaks <- lapply(peaks, FUN = function(peaks) picking_subfunction(peaks, signals_to_take = signals_to_take))
            }
        } else {
            most_intense_peaks <- lapply(peaks, FUN = function(peaks) picking_subfunction(peaks, signals_to_take = signals_to_take))
        }
    } else if (isMassPeaks(peaks)) {
        most_intense_peaks <- picking_subfunction(peaks, signals_to_take)
    }
    return(most_intense_peaks)
}





################################################################################





############################################# AVERAGE THE REPLICATES (BY FOLDER)
# This function averages the spectra contained in the same folder (more suitable for brukerflex format)
average_replicates_by_folder <- function(spectra, folder, spectra_format = "brukerflex") {
    # Load the required libraries
    install_and_load_required_packages("MALDIquant")
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Count the number of spectra
    if (isMassSpectrumList(spectra)) {
        number_of_spectra <- length(spectra)
    } else {
        number_of_spectra <- 1
    }
    # List the spectra files (list the imzML files or the FID files)
    folder_files <- read_spectra_files(folder, spectra_format = spectra_format, full_path = TRUE)
    if (spectra_format == "brukerflex") {
        # Split the path into individual folders (list, each element is a vector with the path splitted for that spectrum)
        folder_files_splitted <- list()
        # Split the paths into folders
        for (f in 1:number_of_spectra) {
            folder_files_splitted[f] <- strsplit(folder_files[f], "/")
        }
        #### SITUATION #1
        # The n-5 folder (where n is the length of the folder_splitted) is not unique, since it can represent the treatment folder; the n-4 folder is not unique, since one can call the replicates all in the same way (e.g. rep1, rep2), because they are placed in the sample folder; the combination n-5/n-4 should be almost unique, but to guarantee the uniqueness it is better to use the n-6/n-5/n-4 combination, so that it can be folder/sample/replicate or sample/treatment/replicate.
        #### SITUATION #2
        # The second case is when the sample folder is the same, but the coordinates change, so the replicates are the spectra under the same coordinate folder, while different coordinates mean different samples.
        # The folder in which the samples are should be the same, there shouldn't be a more complicate folder structure, in that case this method should be revised.
        # n-6 = sample name
        sample_name <- character()
        for (f in 1:number_of_spectra) {
            sample_name <- append(sample_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-6])
        }
        # n-5 = treatment / class
        treatment_name <- character()
        for (f in 1:number_of_spectra) {
            treatment_name <- append(treatment_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-5])
        }
        # n-4 = replicate_name
        #replicate_name <- character()
        #for (f in 1:number_of_spectra) {
            #replicate_name <- append(replicate_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-4])
        #}
        # n-3 = coordinates
        #coordinate_name <- character()
        #for (f in 1:number_of_spectra) {
            #coordinate_name <- append(coordinate_name, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-3])
        #}
        # n-2 = replicate_coordinates
        #replicate_coordinates <- character()
        #for (f in 1:number_of_spectra) {
            #replicate_coordinates <- append(replicate_coordinates, folder_files_splitted[[f]][length(folder_files_splitted[[f]])-2])
        #}
        # Generate a file_vector with only the folder to identify the replicates
        unique_sample_name <- character(length = number_of_spectra)
        # Check for replicates in the same folder
        for (f in 1:(number_of_spectra-1)) {
            ### Situation #1 or Situation #2
            if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] == treatment_name[f])) {
                unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], sep = "/")
                unique_sample_name[f+1] <- paste(sample_name[f], treatment_name[f], sep = "/")
            } else {
                unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], sep = "/")
                unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], sep = "/")
            }
            #if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] == treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
                ### Situation #1
                #replicate_name[f+1] <- replicate_name[f]
                #unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], replicate_name[f], sep = "/")
                #unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep = "/")
            #} else if ((sample_name[f+1] == sample_name[f]) && (treatment_name[f+1] != treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
                #unique_sample_name[f] <- paste(sample_name[f], treatment_name[f], replicate_name[f], sep = "/")
                #unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep = "/")
            #} else if ((treatment_name[f+1] == treatment_name[f]) && (replicate_name[f+1] != replicate_name[f])) {
                # Situation #1 + Situation #2
                #unique_sample_name[f] <- replicate_name[f]
                #unique_sample_name[f+1] <- paste(sample_name[f+1], treatment_name[f+1], replicate_name[f+1], sep = "/")
            #}
        }
    } else if (spectra_format == "imzML" || spectra_format == "imzml") {
        unique_sample_name <- character()
        for (f in 1:number_of_spectra) {
            unique_sample_name <- append(unique_sample_name, spectra[[f]]@metaData$file[1])
        }
    }
    # Average the mass spectra, grouping them according to the sample_vector
    spectra_replicates_averaged <- averageMassSpectra(spectra, labels = unique_sample_name)
    return(spectra_replicates_averaged)
}





################################################################################





################################################################### PEAK PICKING
# This function takes a list of spectra (MALDIquant) and computes the peak picking.
peak_picking <- function(spectra, peak_picking_algorithm = "SuperSmoother", tof_mode = "linear", SNR = 3, allow_parallelization = FALSE, deisotope_peaklist = FALSE) {
    ########## Load the required libraries
    install_and_load_required_packages(c("MALDIquant", "parallel"))
    ########## Multi-core
    if (allow_parallelization == TRUE) {
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
    }
    ##### TOF-MODE
    if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
        half_window_size <- 20
    } else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
        half_window_size <- 5
    }
    ########## SINGLE SPECTRUM
    if (isMassSpectrum(spectra)) {
        peaks <- detectPeaks(spectra, method = peak_picking_algorithm, halfWindowSize = half_window_size, SNR = SNR)
    }
    ########## MULTIPLE SPECTRA
    if (isMassSpectrumList(spectra)) {
        # Peak detection
        if (allow_parallelization == TRUE) {
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                peaks <- mclapply(spectra, FUN = function(spectra) detectPeaks(spectra, method = peak_picking_algorithm, halfWindowSize = half_window_size, SNR = SNR), mc.cores = cpu_thread_number)
            } else if (Sys.info()[1] == "Windows") {
                # Make the cluster (one for each core/thread)
                cl <- makeCluster(cpu_thread_number)
                clusterEvalQ(cl, {library(MALDIquant)})
                clusterExport(cl = cl, varlist = c("peak_picking_algorithm", "half_window_size", "SNR"), envir = environment())
                peaks <- parLapply(cl, spectra, fun = function(spectra) detectPeaks(spectra, method = peak_picking_algorithm, halfWindowSize = half_window_size, SNR = SNR))
                stopCluster(cl)
            } else {
                peaks <- lapply(spectra, FUN = function(spectra) detectPeaks(spectra, method = peak_picking_algorithm, halfWindowSize = half_window_size, SNR = SNR))
            }
        } else {
            peaks <- lapply(spectra, FUN = function(spectra) detectPeaks(spectra, method = peak_picking_algorithm, halfWindowSize = half_window_size, SNR = SNR))
        }
    }
    ##### Deisotope peaklist
    if (deisotope_peaklist == TRUE && (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R")) {
        peaks <- deisotope_peaks(peaks, pattern_model_correlation = 0.95, isotopic_tolerance = 10^(-4), isotope_pattern_distance = 1.00235, isotopic_pattern_size = 3L:10L, allow_parallelization = allow_parallelization)
    }
    ##### Return
    return(peaks)
}





################################################################################





################################################################### PEAK PICKING
# This function takes a list of peaks (MALDIquant) and returns the same peak list without isotopic clusters, only monoisotopic peaks.
deisotope_peaks <- function(peaks, pattern_model_correlation = 0.95, isotopic_tolerance = 10^(-4), isotope_pattern_distance = 1.00235, isotopic_pattern_size = 3L:5L, allow_parallelization = FALSE) {
    ##### Load the required packages
    install_and_load_required_packages(c("MALDIquant", "parallel"))
    ##### Multiple peaks
    if (isMassPeaksList(peaks)) {
        ##### Multiple cores
        if (allow_parallelization == TRUE) {
            # Detect the number of cores
            cpu_thread_number <- detectCores(logical = TRUE)
            cpu_thread_number <- cpu_thread_number / 2
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                peaks_deisotoped <- mclapply(peaks, FUN = function(peaks) monoisotopicPeaks(peaks, minCor = pattern_model_correlation, tolerance = isotopic_tolerance, distance = isotope_pattern_distance, size = isotopic_pattern_size), mc.cores = cpu_thread_number)
            } else if (Sys.info()[1] == "Windows") {
                # Make the CPU cluster for parallelisation
                cl <- makeCluster(cpu_thread_number)
                # Make the cluster use the custom functions and the package functions along with their parameters
                clusterEvalQ(cl, {library(MALDIquant)})
                # Pass the variables to the cluster for running the function
                clusterExport(cl = cl, varlist = c("peaks", "pattern_model_correlation", "isotopic_tolerance", "isotope_pattern_distance", "isotopic_pattern_size"), envir = environment())
                # Apply the multicore function
                peaks_deisotoped <- parLapply(cl, peaks, fun = function(peaks) monoisotopicPeaks(peaks, minCor = pattern_model_correlation, tolerance = isotopic_tolerance, distance = isotope_pattern_distance, size = isotopic_pattern_size))
                stopCluster(cl)
            } else {
                # Run the algorithm
                peaks_deisotoped <- lapply(peaks, FUN = function(peaks) monoisotopicPeaks(peaks, minCor = pattern_model_correlation, tolerance = isotopic_tolerance, distance = isotope_pattern_distance, size = isotopic_pattern_size))
            }
        } else {
            # Run the algorithm
            peaks_deisotoped <- lapply(peaks, FUN = function(peaks) monoisotopicPeaks(peaks, minCor = pattern_model_correlation, tolerance = isotopic_tolerance, distance = isotope_pattern_distance, size = isotopic_pattern_size))
        }
    } else {
        # Run the algorithm
        peaks_deisotoped <- monoisotopicPeaks(peaks, minCor = pattern_model_correlation, tolerance = isotopic_tolerance, distance = isotope_pattern_distance, size = isotopic_pattern_size)
    }
    # Return
    return(peaks_deisotoped)
}





################################################################# PEAK ALIGNMENT
# This function takes a list of peaks (MALDIquant) and computes the peak alignment, along with the false positive removal and the removal of low-intensity peaks.
align_and_filter_peaks <- function(peaks, peak_picking_algorithm = "SuperSmoother", tof_mode = "linear", peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", reference_peaklist = NULL, spectra = NULL, alignment_iterations = 5, allow_parallelization = FALSE) {
    ########## Determine the tolerance in PPM
    if (tof_mode == "linear" || tof_mode == "Linear" || tof_mode == "L") {
        tolerance_ppm <- 2000
    } else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
        tolerance_ppm <- 200
    }
    ########## Align only if there are many peaklists
    if (isMassPeaksList(peaks)) {
        ##### Load the required libraries
        install_and_load_required_packages(c("MALDIquant", "parallel"))
        ##### Rename the trim function
        trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
        ###### Peak alignment
        # Fix the iteration number
        if (alignment_iterations <= 0) {
            alignment_iterations <- 5
        }
        # Initialize the variable
        peaks_aligned <- NULL
        # For each iteration...
        for (iter in 1:alignment_iterations) {
            # If it is the first time that the alignment is run...
            if (is.null(peaks_aligned)) {
                # Run the alignment
                peaks_aligned <- binPeaks(peaks, method = "relaxed", tolerance = (tolerance_ppm/10^6))
            } else {
                # Otherwise run another alignment on the already aligned peaks
                peaks_aligned <- binPeaks(peaks_aligned, method = "relaxed", tolerance = (tolerance_ppm/10^6))
            }
        }
        ##### False positive removal
        if (!is.null(peak_filtering_frequency_threshold_percent) && peak_filtering_frequency_threshold_percent > 0) {
            peaks_aligned <- filterPeaks(peaks_aligned, minFrequency = (peak_filtering_frequency_threshold_percent/100))
        }
        ##### Low-intensity peaks removal
        if (!is.null(low_intensity_peak_removal_threshold_percent) && low_intensity_peak_removal_threshold_percent > 0 && low_intensity_peak_removal_threshold_percent < 100) {
            peaks_aligned <- remove_low_intensity_peaks(peaks_aligned, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization)
        }
        ##### Align to a reference peaklist: AVERAGE SPECTRUM (if a spectra list is provided)
        if (is.character(reference_peaklist) && reference_peaklist == "average" && !is.null(spectra)) {
            # Average the spectra
            average_spectrum <- averageMassSpectra(spectra, method = "mean")
            # Peak picking
            average_spectrum_peaks <- peak_picking(average_spectrum, peak_picking_algorithm = peak_picking_algorithm, SNR = 5, allow_parallelization = allow_parallelization)
            # The reference peaklist is the paklist of the average spectrum
            reference_peaklist <- createMassPeaks(mass = average_spectrum_peaks@mass, intensity = average_spectrum_peaks@intensity, snr = rep.int(5, length(average_spectrum_peaks@mass)), metaData = list(name = "Reference peaklist AVG"))
        } else if (is.character(reference_peaklist) && reference_peaklist == "average" && is.null(spectra)) {
            reference_peaklist <- NULL
        } else if (!is.null(reference_peaklist) && is.vector(reference_peaklist)) {
            reference_peaklist <- createMassPeaks(mass = as.numeric(reference_peaklist), intensity = rep.int(1, length(reference_peaklist)), snr = rep.int(5, length(reference_peaklist)), metaData = list(name = "Reference peaklist"))
        }
        ##### Align to the reference peaklist
        if (!is.null(reference_peaklist)) {
            try({
                ### Determine the warping function (peak alignment to a reference)
                warping_functions <- determineWarpingFunctions(peaks_aligned, reference = reference_peaklist, tolerance = (tolerance_ppm/10^6), method = "quadratic")
                ############ Function for lapply
                #align_peaks_subfunction <- function(peaks, reference_peaklist, tolerance_ppm) {
                #    mass_vector <- peaks@mass
                #    # For each reference peak
                #    for (ref in reference_peaklist) {
                #        # Replace the value in the vector with the reference value
                #        mass_vector[which(abs(mass_vector - ref)*10^6/ref <= tolerance_ppm)] <- ref
                #    }
                #    # Put the fixed mass vector back into the MALDIquant peaklist
                #    peaks@mass <- mass_vector
                #    return(peaks)
                #}
                ############# If there are many peaklists or one peaklist (use multicore)
                if (isMassPeaksList(peaks_aligned)) {
                    if (allow_parallelization == TRUE) {
                        # Detect the number of cores
                        cpu_thread_number <- detectCores(logical = TRUE)
                        cpu_thread_number <- cpu_thread_number / 2
                        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                            #peaks_aligned <- mclapply(peaks_aligned, FUN = function(peaks_aligned) align_peaks_subfunction(peaks_aligned, reference_peaklist, tolerance_ppm), mc.cores = cpu_thread_number)
                            peaks_aligned <- mclapply(peaks_aligned, FUN = function(peaks_aligned) warpMassPeaks(peaks_aligned, w = warping_functions), mc.cores = cpu_thread_number)
                        } else if (Sys.info()[1] == "Windows") {
                            # Make the CPU cluster for parallelisation
                            cl <- makeCluster(cpu_thread_number)
                            # Apply the multicore function
                            # Pass the variables to the cluster for running the function
                            clusterExport(cl = cl, varlist = c("reference_peaklist", "tolerance_ppm"), envir = environment())
                            peaks_aligned <- parLapply(cl, peaks_aligned, fun = function(peaks_aligned) warpMassPeaks(peaks_aligned, w = warping_functions))
                            stopCluster(cl)
                        }
                    } else {
                        peaks_aligned <- warpMassPeaks(peaks_aligned, w = warping_functions)
                    }
                } else {
                    peaks_aligned <- warpMassPeaks(peaks_aligned, w = warping_functions)
                }
            }, silent = TRUE)
        }
        # Low-intensity peaks removal
        if (!is.null(low_intensity_peak_removal_threshold_percent) && low_intensity_peak_removal_threshold_percent > 0 && low_intensity_peak_removal_threshold_percent < 100) {
            peaks <- remove_low_intensity_peaks(peaks_aligned, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization)
        }
        return(peaks_aligned)
    } else {
        # Low-intensity peaks removal
        if (!is.null(low_intensity_peak_removal_threshold_percent) && low_intensity_peak_removal_threshold_percent > 0 && low_intensity_peak_removal_threshold_percent < 100) {
            peaks <- remove_low_intensity_peaks(peaks, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization)
        }
        return(peaks)
    }
}





################################################################################






































































































################################################################ CLASSIFICATION

############################## SINGLE MODEL CLASSIFICATION (MULTICORE, ENSEMBLE)
# The function takes a folder in which there are imzML files (one for each patient) or an imzML file or a list of MALDIquant spectra files, the R workspace containing the models with the name of the model objects in the workspace, and allows the user to specify something regarding the preprocessing of the spectra to be classified.
# The function outputs a list containing: a matrix with the classification (pixel-by-pixel), MS images with the pixel-by-pixel classification, the model list, a matrix with the ensemble classification (pixel-by-pixel) and MS images with the pixel-by-pixel ensemble classification.
# Parallel computation implemented.
# It outputs NULL values if the classification cannot be performed due to incompatibilities between the model features and the spectral features.
# The pixel grouping cannot be 'graph', otherwise, when embedded in the pixel by pixel classification function, the graph segmentation is performed for each model before making the predictons.
single_model_classification_of_spectra <- function(spectra, model_x, model_name = "model", preprocessing_parameters = NULL, peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = FALSE, peak_picking_SNR = 5, peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 1, low_intensity_peak_removal_threshold_method = "element-wise", tof_mode = "linear", allow_parallelization = FALSE, pixel_grouping = c("single", "hca", "moving window average", "graph"), number_of_hca_nodes = 5, moving_window_size = 5, final_result_matrix = NULL, seed = 12345, correlation_method_for_adjacency_matrix = "pearson", correlation_threshold_for_adjacency_matrix = 0.95, pvalue_threshold_for_adjacency_matrix = 0.05, max_GA_generations = 10, iterations_with_no_change = 5, number_of_spectra_partitions = 1, partitioning_method = "space", plot_figures = TRUE, plot_graphs = TRUE, plot_legends = c("sample name", "legend", "plot name")) {
    # Class list (from the custom model entry)
    class_list <- model_x$class_list
    # Outcome list (from the custom model entry)
    outcome_list <- model_x$outcome_list
    # Retrieve the peaks used to create the model (from the custom model entry)
    features_model <- model_x$features_model
    # Model ID
    model_ID <- model_x$model_ID
    # Model object
    model_object <- model_x$model
    # Outcomes
    outcome_list <- model_x$outcome_list
    ########## MULTIPLE SPECTRA
    if (isMassSpectrumList(spectra)) {
        ########## SINGLE PIXEL CLASSIFICATION
        if (pixel_grouping == "single" || number_of_hca_nodes == 1 || moving_window_size >= length(spectra)) {
            ### Generate the intensity matrix with the features from the model
            final_sample_matrix <- generate_custom_intensity_matrix(spectra, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
            ### Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
            if (!is.null(final_sample_matrix)) {
                # Put the X at the beginning of the peak names
                for (n in 1:length(colnames(final_sample_matrix))) {
                    colnames(final_sample_matrix)[n] <- paste("X", colnames(final_sample_matrix)[n], sep = "")
                }
                ##### Predictions (spectra by spectra) (class, no probabilities)
                if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                    predicted_classes <- as.character(predict(model_object, newdata = final_sample_matrix, type = "raw"))
                } else {
                    predicted_classes <- as.character(predict(model_object, newdata = final_sample_matrix))
                    #predicted_classes <- as.character(apply(X = final_sample_matrix, MARGIN = 1, FUN = function(x) predict(model_object, as.matrix(rbind(x)))))
                }
                # Generate a matrix with the results
                result_matrix_model <- matrix(nrow = length(predicted_classes), ncol = 1)
                sample_name <- spectra[[1]]@metaData$file[1]
                rownames(result_matrix_model) <- cbind(rep(sample_name, length(spectra)))
                result_matrix_model[,1] <- cbind(predicted_classes)
                colnames(result_matrix_model) <- paste("Predicted Class '", model_name, "'", sep = "")
                ##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
                if (is.null(final_result_matrix)) {
                    final_result_matrix <- result_matrix_model
                } else {
                    final_result_matrix <- cbind(final_result_matrix, result_matrix_model)
                }
                ########## Generate a molecular image of the classification
                # Define the class as number depending on the outcome
                outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_classes)
                # Replace the spectra intensities with the class number for plotting purposes
                class_as_number <- outcome_and_class$class_vector_as_numeric
                spectra_for_plotting <- spectra
                for (s in 1:length(spectra_for_plotting)) {
                    spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s], length(spectra_for_plotting[[s]]@intensity))
                }
                # Generate the MS image
                slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                plotMsiSlice(slices, legend = FALSE, scale = F)
                if ("legend" %in% plot_legends) {
                    # Define the legend
                    legend_text <- outcome_and_class$legend_text
                    legend_fill <- outcome_and_class$legend_fill
                    legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                }
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = sample_name, xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                }
                # Store the plot into the list of images (for model)
                classification_msi_model <- recordPlot()
            } else {
                ### It is NULL if there are incompatibilities between the model features and the spectral features
                result_matrix_model <- NULL
                classification_msi_model <- NULL
            }
        } else if (pixel_grouping == "moving window average") {
            ########## MOVING WINDOW AVERAGE
            ##### Rearrange the spectra according to the space coordinates
            spectra <- rearrange_spectral_dataset(spectra, rearranging_method = "space")
            ##### Initialize the model result matrix
            result_matrix_model <- NULL
            ##### For each spectrum...
            for (s in 1:length(spectra)) {
                ### Define the indices
                index1 <- s - floor(moving_window_size/2)
                index2 <- s + floor(moving_window_size/2)
                ### Check the indices
                if (index1 <= 0) {
                    index1 <- 1
                } else if (index1 > length(spectra)) {
                    index1 <- length(spectra)
                }
                if (index2 <= 0) {
                    index2 <- 1
                } else if (index2 > length(spectra)) {
                    index2 <- length(spectra)
                }
                ### Isolate the spectra from the bin
                spectra_bin <- spectra[index1:index2]
                ### Generate the average spectrum for the bin
                average_spectrum_bin <- averageMassSpectra(spectra_bin)
                ### Preprocessing the AVG spectrum
                average_spectrum_bin <- preprocess_spectra(average_spectrum_bin, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                ### Peak picking on the AVG spectrum
                sample_bin_matrix <- generate_custom_intensity_matrix(average_spectrum_bin, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                ### Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
                if (!is.null(sample_bin_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (n in 1:length(colnames(sample_bin_matrix))) {
                        colnames(sample_bin_matrix)[n] <- paste("X", colnames(sample_bin_matrix)[n], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_classes <- as.character(predict(model_object, newdata = final_sample_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_bin_matrix))
                    }
                    # Generate a matrix with the results
                    result_matrix_model_bin <- matrix(nrow = 1, ncol = 1)
                    sample_name <- average_spectrum_bin@metaData$file[1]
                    rownames(result_matrix_model_bin) <- sample_name
                    result_matrix_model_bin[1,1] <- as.character(predicted_class_avg)
                    colnames(result_matrix_model_bin) <- paste("Predicted Class '", model_name, "'", sep = "")
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    result_matrix_model_bin <- NULL
                }
                ### Add the result matrix to the result matrix of the model
                if (is.null(result_matrix_model)) {
                    result_matrix_model <- result_matrix_model_bin
                } else {
                    result_matrix_model <- rbind(result_matrix_model, result_matrix_model_bin)
                }
            }
            ### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
            if (is.null(final_result_matrix)) {
                final_result_matrix <- result_matrix_model
            } else {
                final_result_matrix <- cbind(final_result_matrix, result_matrix_model)
            }
            ########## Generate a molecular image of the classification
            ### Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
            if (!is.null(result_matrix_model)) {
                # Replace the spectra intensities with the class number for plotting purposes
                predicted_classes <- as.character(result_matrix_model[,1])
                # Define the class as number depending on the outcome
                outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_classes)
                # Replace the spectra intensities with the class number for plotting purposes
                class_as_number <- outcome_and_class$class_vector_as_numeric
                spectra_for_plotting <- spectra
                for (s in 1:length(spectra_for_plotting)) {
                    spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s], length(spectra_for_plotting[[s]]@intensity))
                }
                # Generate the MS images
                slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                plotMsiSlice(slices, legend = FALSE, scale = F)
                if ("legend" %in% plot_legends) {
                    # Define the legend
                    legend_text <- outcome_and_class$legend_text
                    legend_fill <- outcome_and_class$legend_fill
                    legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                }
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = sample_name, xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                }
                # Store the plot into the list of images (for model)
                classification_msi_model <- recordPlot()
            } else {
                ### It is NULL if there are incompatibilities between the model features and the spectral features
                classification_msi_model <- NULL
            }
        } else if (pixel_grouping == "hca") {
            ########## HCA
            ### Initialize the model result matrix
            result_matrix_model <- matrix("class", nrow = length(spectra), ncol = 1)
            colnames(result_matrix_model) <- paste("Predicted Class '", model_name, "'", sep = "")
            sample_name <- spectra[[1]]@metaData$file[1]
            rownames(result_matrix_model) <- rep(sample_name, nrow(result_matrix_model))
            ### Detect and align peaks
            peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = peak_picking_SNR, allow_parallelization = allow_parallelization)
            peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent)
            ### Generate the peaklist matrix
            peaklist <- intensityMatrix(peaks, spectra)
            ### Compute the distance matrix
            distance_matrix <- dist(peaklist, method = "euclidean")
            # Generate the dendrogram
            hca <- hclust(distance_matrix)
            ### Cut the tree to generate K number of sub-clusters
            hca_groups <- cutree(hca, k = number_of_hca_nodes)
            ## Generate the spectra for plotting list
            spectra_for_plotting <- list()
            ## For each HCA node...
            for (n in 1:number_of_hca_nodes) {
                # Index the spectra under in the selected subgroup of the HCA
                index <- which(hca_groups == n)
                spectra_hca <- spectra[index]
                # Generate the average spectrum for these spectra under the node
                average_spectrum_hca <- averageMassSpectra(spectra_hca)
                ### Preprocessing the AVG spectrum
                average_spectrum_hca <- preprocess_spectra(average_spectrum_hca, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                ### Peak picking on the AVG spectrum
                sample_hca_matrix <- generate_custom_intensity_matrix(average_spectrum_hca, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                ### Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
                if (!is.null(sample_hca_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (x in 1:length(colnames(sample_hca_matrix))) {
                        colnames(sample_hca_matrix)[x] <- paste("X", colnames(sample_hca_matrix)[x], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_hca_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_hca_matrix))
                    }
                    # Fill the model matrix with the results
                    result_matrix_model[index,1] <- as.character(predicted_class_avg)
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    result_matrix_model <- NULL
                }
            }
            ### Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
            if (!is.null(result_matrix_model)) {
                # Replace the spectra intensities with the class number for plotting purposes
                predicted_classes <- as.character(result_matrix_model[,1])
                # Define the class as number depending on the outcome
                outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_classes)
                # Replace the spectra intensities with the class number for plotting purposes
                class_as_number <- outcome_and_class$class_vector_as_numeric
                ### Fix the spectra for plotting list (edit the intensity values)
                spectra_for_plotting <- spectra
                for (s in 1:length(spectra_for_plotting)) {
                    spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s], length(spectra_for_plotting[[s]]@intensity))
                }
                # Generate the MS images
                slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                plotMsiSlice(slices, legend = FALSE, scale = F)
                if ("legend" %in% plot_legends) {
                    # Define the legend
                    legend_text <- outcome_and_class$legend_text
                    legend_fill <- outcome_and_class$legend_fill
                    legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                }
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = sample_name, xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                }
                # Store the plot into the list of images (for SVM)
                classification_msi_model <- recordPlot()
            } else {
                ### It is NULL if there are incompatibilities between the model features and the spectral features
                classification_msi_model <- NULL
            }
            ##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
            if (is.null(final_result_matrix)) {
                final_result_matrix <- result_matrix_model
            } else {
                final_result_matrix <- cbind(final_result_matrix, result_matrix_model)
            }
        } else if (pixel_grouping == "graph") {
            # Initialize output
            spectra_for_plotting <- list()
            ### Initialize the model result matrix
            result_matrix_model <- matrix("class", nrow = length(spectra), ncol = 1)
            colnames(result_matrix_model) <- paste("Predicted Class '", model_name, "'", sep = "")
            sample_name <- spectra[[1]]@metaData$file[1]
            rownames(result_matrix_model) <- rep(sample_name, nrow(result_matrix_model))
            ########## GRAPH SEGMENTATION
            graph_segmentation <- graph_MSI_segmentation(filepath_imzml = spectra, preprocessing_parameters = NULL, allow_parallelization = allow_parallelization, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = deisotope_peaklist, SNR = peak_picking_SNR, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, custom_feature_vector = features_model, correlation_method_for_adjacency_matrix = correlation_method_for_adjacency_matrix, correlation_threshold_for_adjacency_matrix = correlation_threshold_for_adjacency_matrix, pvalue_threshold_for_adjacency_matrix = pvalue_threshold_for_adjacency_matrix, number_of_high_degree_vertices_for_subgraph = 0, vertices_not_in_induced_subgraph = "independent", max_GA_generations = max_GA_generations, iterations_with_no_change = iterations_with_no_change, plot_figures = plot_figures, plot_graphs = plot_graphs, number_of_spectra_partitions = number_of_spectra_partitions, partitioning_method = partitioning_method, seed = seed, plot_legends = plot_legends)
            ### Extract the spectra from the clique and from the independent set (mind that they are two lists if the spectra are taken all at the same time or they are two lists of lists if the spectra are partitioned)
            if (number_of_spectra_partitions > 1) {
                spectra_clique <- graph_segmentation$spectra_clique
                spectra_independent <- graph_segmentation$spectra_independent
                # Generate the average of the clique and the independent set
                spectra_clique_avg <- list()
                spectra_independent_avg <- list()
                for (l in 1:length(spectra_independent)) {
                    spectra_independent_avg <- append(spectra_independent_avg, averageMassSpectra(spectra_independent[[l]]))
                }
                for (l in 1:length(spectra_clique)) {
                    spectra_clique_avg <- append(spectra_clique_avg, averageMassSpectra(spectra_clique[[l]]))
                }
                # Preprocess the average spectra
                spectra_independent_avg <- preprocess_spectra(spectra_independent_avg, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                spectra_clique_avg <- preprocess_spectra(spectra_clique_avg, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                # Peak picking on the AVG spectrum
                sample_independent_matrix <- generate_custom_intensity_matrix(spectra_independent_avg, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                sample_clique_matrix <- generate_custom_intensity_matrix(spectra_clique_avg, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                ## Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
                if (!is.null(sample_clique_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (x in 1:length(colnames(sample_clique_matrix))) {
                        colnames(sample_clique_matrix)[x] <- paste("X", colnames(sample_clique_matrix)[x], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_clique_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_clique_matrix))
                    }
                    # Fill the model matrix with the results
                    #result_matrix_svm_clique[index,1] <- as.character(predicted_class_avg)
                    ### Fix the spectra for plotting list (edit the intensity values)
                    # Define the class as number depending on the outcome
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_class_avg)
                    # Replace the spectra intensities with the class number for plotting purposes
                    class_as_number <- outcome_and_class$class_vector_as_numeric
                    for (l in 1:length(spectra_clique)) {
                        spectra_for_plotting_clique <- spectra_clique[[l]]
                        for (s in 1:length(spectra_for_plotting_clique)) {
                            spectra_for_plotting_clique[[s]]@intensity <- rep(class_as_number[l], length(spectra_for_plotting_clique[[s]]@intensity))
                        }
                        # Append to the final list of spectra for plotting
                        spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_clique)
                    }
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    spectra_for_plotting_clique <- list()
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_clique)
                    #result_matrix_svm <- NULL
                }
                if (!is.null(sample_independent_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (x in 1:length(colnames(sample_independent_matrix))) {
                        colnames(sample_independent_matrix)[x] <- paste("X", colnames(sample_independent_matrix)[x], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_independent_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_independent_matrix))
                    }
                    # Fill the model matrix with the results
                    #result_matrix_svm_clique[index,1] <- as.character(predicted_class_avg)
                    ### Fix the spectra for plotting list (edit the intensity values)
                    # Define the class as number depending on the outcome
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_class_avg)
                    # Replace the spectra intensities with the class number for plotting purposes
                    class_as_number <- outcome_and_class$class_vector_as_numeric
                    for (l in 1:length(spectra_independent)) {
                        spectra_for_plotting_independent <- spectra_independent[[l]]
                        for (s in 1:length(spectra_for_plotting_independent)) {
                            spectra_for_plotting_independent[[s]]@intensity <- rep(class_as_number[l], length(spectra_for_plotting_independent[[s]]@intensity))
                        }
                        # Append to the final list of spectra for plotting
                        spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_independent)
                    }
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    spectra_for_plotting_independent <- list()
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_independent)
                    #result_matrix_svm <- NULL
                }
                ### Run only if there are some spectra classified
                if (length(spectra_for_plotting) > 0) {
                    # Generate the MS images
                    slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                    plotMsiSlice(slices, legend = FALSE, scale = F)
                    if ("legend" %in% plot_legends) {
                        # Define the legend
                        legend_text <- outcome_and_class$legend_text
                        legend_fill <- outcome_and_class$legend_fill
                        legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                    }
                    if ("sample name" %in% plot_legends) {
                        legend(x = "topright", legend = spectra_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                    }
                    if ("plot name" %in% plot_legends) {
                        legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                    }
                    # Store the plot into the list of images (for SVM)
                    classification_msi_model <- recordPlot()
                    # Generate the result matrix from the rearranged spectra (so that it is reproducible, because the spectra of the clique mught not be always the same) according to the intensity value
                    spectra_for_plotting <- rearrange_spectral_dataset(spectra_for_plotting, rearranging_method = "space")
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list)
                    predicted_classes <- character()
                    for (s in 1:length(spectra_for_plotting)) {
                        for (ou in 1:length(outcome_and_class$class_outcome_matrix$Number)) {
                            if (spectra_for_plotting[[s]]@intensity[1] == outcome_and_class$class_outcome_matrix$Number[ou]) {
                                predicted_classes <- append(predicted_classes, outcome_and_class$class_outcome_matrix$Class[ou])
                            }
                        }
                    }
                    result_matrix_model[,1] <- cbind(as.character(predicted_classes))
                    ##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
                    if (is.null(final_result_matrix)) {
                        final_result_matrix <- result_matrix_model
                    } else {
                        final_result_matrix <- cbind(final_result_matrix, result_matrix_model)
                    }
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    result_matrix_model <- NULL
                    classification_msi_model <- NULL
                }
            } else if (number_of_spectra_partitions <= 1) {
                ######### NO PARTITION SPECTRA
                spectra_clique <- graph_segmentation$spectra_clique
                spectra_independent <- graph_segmentation$spectra_independent
                # Generate the average of the clique and the independent set
                spectra_clique_avg <- averageMassSpectra(spectra_clique)
                spectra_independent_avg <- averageMassSpectra(spectra_independent)
                # Preprocess the average spectrum
                spectra_clique_avg <- preprocess_spectra(spectra_clique_avg, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                spectra_independent_avg <- preprocess_spectra(spectra_independent_avg, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                # Peak picking on the AVG spectrum
                sample_clique_matrix <- generate_custom_intensity_matrix(spectra_clique_avg, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                sample_independent_matrix <- generate_custom_intensity_matrix(spectra_independent_avg, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = peak_picking_SNR, peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
                ## Run only if the sample matrix is not NULL: it is NULL if there are incompatibilities between the model features and the spectral features
                if (!is.null(sample_clique_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (x in 1:length(colnames(sample_clique_matrix))) {
                        colnames(sample_clique_matrix)[x] <- paste("X", colnames(sample_clique_matrix)[x], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_clique_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_clique_matrix))
                    }
                    ### Fix the spectra for plotting list (edit the intensity values)
                    # Define the class as number depending on the outcome
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_class_avg)
                    # Replace the spectra intensities with the class number for plotting purposes
                    class_as_number <- outcome_and_class$class_vector_as_numeric
                    spectra_for_plotting_clique <- spectra_clique
                    for (s in 1:length(spectra_for_plotting_clique)) {
                        spectra_for_plotting_clique[[s]]@intensity <- rep(class_as_number, length(spectra_for_plotting_clique[[s]]@intensity))
                    }
                    # Append to the final list of spectra for plotting
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_clique)
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    spectra_for_plotting_clique <- list()
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_clique)
                }
                if (!is.null(sample_independent_matrix)) {
                    ### Put the X at the beginning of the peak names
                    for (x in 1:length(colnames(sample_independent_matrix))) {
                        colnames(sample_independent_matrix)[x] <- paste("X", colnames(sample_independent_matrix)[x], sep = "")
                    }
                    ##### Predictions (AVG spectrum) (class, no probabilities)
                    if (model_ID == "rf" || model_ID == "nbc" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda") {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_independent_matrix, type = "raw"))
                    } else {
                        predicted_class_avg <- as.character(predict(model_object, newdata = sample_independent_matrix))
                    }
                    ### Fix the spectra for plotting list (edit the intensity values)
                    # Define the class as number depending on the outcome
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list, class_vector = predicted_class_avg)
                    # Replace the spectra intensities with the class number for plotting purposes
                    class_as_number <- outcome_and_class$class_vector_as_numeric
                    spectra_for_plotting_independent <- spectra_independent
                    for (s in 1:length(spectra_for_plotting_independent)) {
                        spectra_for_plotting_independent[[s]]@intensity <- rep(class_as_number, length(spectra_for_plotting_independent[[s]]@intensity))
                    }
                    # Append to the final list of spectra for plotting
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_independent)
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    spectra_for_plotting_independent <- list()
                    spectra_for_plotting <- append(spectra_for_plotting, spectra_for_plotting_independent)
                }
                ### Run only if there are some spectra classified
                if (length(spectra_for_plotting) > 0) {
                    # Generate the MS images
                    slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                    plotMsiSlice(slices, legend = FALSE, scale = F)
                    if ("legend" %in% plot_legends) {
                        # Define the legend
                        legend_text <- outcome_and_class$legend_text
                        legend_fill <- outcome_and_class$legend_fill
                        legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                    }
                    if ("sample name" %in% plot_legends) {
                        legend(x = "topright", legend = spectra_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                    }
                    if ("plot name" %in% plot_legends) {
                        legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                    }
                    # Store the plot into the list of images
                    classification_msi_model <- recordPlot()
                    # Generate the result matrix from the rearranged spectra (so that it is reproducible, because the spectra of the clique mught not be always the same) according to the intensity value
                    spectra_for_plotting <- rearrange_spectral_dataset(spectra_for_plotting, rearranging_method = "space")
                    outcome_and_class <- outcome_and_class_to_MS(class_list = class_list, outcome_list = outcome_list)
                    predicted_classes <- character()
                    for (s in 1:length(spectra_for_plotting)) {
                        for (ou in 1:length(outcome_and_class$class_outcome_matrix$Number)) {
                            if (spectra_for_plotting[[s]]@intensity[1] == outcome_and_class$class_outcome_matrix$Number[ou]) {
                                predicted_classes <- append(predicted_classes, outcome_and_class$class_outcome_matrix$Class[ou])
                            }
                        }
                    }
                    result_matrix_model[,1] <- cbind(as.character(predicted_classes))
                    ##### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
                    if (is.null(final_result_matrix)) {
                        final_result_matrix <- result_matrix_model
                    } else {
                        final_result_matrix <- cbind(final_result_matrix, result_matrix_model)
                    }
                } else {
                    ### It is NULL if there are incompatibilities between the model features and the spectral features
                    result_matrix_model <- NULL
                    classification_msi_model <- NULL
                }
            }
        }
        ##### Return
        return(list(result_matrix_model = result_matrix_model, classification_msi_model = classification_msi_model, final_result_matrix = final_result_matrix))
    } else if (isMassSpectrum(spectra)) {
        ########## SINGLE (AVG) SPECTRUM
        ### Inizialize output
        result_matrix <- NULL
        average_spectrum_with_bars <- NULL
        ### Generate the intensity matrix with the features from the model
        final_sample_matrix <- generate_custom_intensity_matrix(spectra, custom_feature_vector = features_model, tof_mode = tof_mode, preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, peak_picking_SNR = 3, peak_filtering_frequency_threshold_percent = NULL, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
        ### Run only if there is compatibility between the spectral features and the model features (it is NULL if there is no compatibility)
        if (!is.null(final_sample_matrix)) {
            # Put the X at the beginning of the peak names
            for (n in 1:length(colnames(final_sample_matrix))) {
                name <- paste("X", colnames(final_sample_matrix)[n], sep = "")
                colnames(final_sample_matrix)[n] <- name
            }
            # Predictions
            if (model_ID == "rf" || model_ID == "knn" || model_ID == "nnet" || model_ID == "lda" || model_ID == "nbc") {
                predicted_classes <- as.character(predict(model_object, newdata = final_sample_matrix, type = "raw"))
            } else {
                predicted_classes <- as.character(predict(model_object, newdata = final_sample_matrix))
            }
            # Generate a matrix with the results
            result_matrix <- matrix(nrow = 1, ncol = 1)
            rownames(result_matrix) <- spectra@metaData$file[[1]]
            result_matrix[, 1] <- as.character(predicted_classes)
            colnames(result_matrix) <- paste("Predicted Class", model_name)
            #### Add the result matrix to a global final matrix (the final result matrix is the patient matrix if it is still non existent)
            if (is.null(final_result_matrix)) {
                final_result_matrix <- result_matrix
            } else {
                final_result_matrix <- cbind(final_result_matrix, result_matrix)
            }
            ################# Average spectrum with bars onto the signals used by the model
            for (f in 1:length(features_model)) {
                name_splitted <- unlist(strsplit(features_model[f],""))
                feature_def <- name_splitted[2]
                for (i in 3:length(name_splitted)) {
                    feature_def <- paste(feature_def, name_splitted[i], sep = "")
                }
                features_model[f] <- feature_def
            }
            # Average spectrum: sample_spectra_avg; model features: features_model; peaks: sample_peaks_avg
            # Detect peaks in the avg (SNR = 1)
            sample_peaks_avg_for_bars <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = 3)
            # Determine the coordinates of the bars
            coordinates_of_bars <- list(x = numeric(), y = numeric())
            # Check if the features used for the model are in the spectrum
            for (f in features_model) {
                presence_in_the_avg <- FALSE
                # Scroll the peaks
                for (z in 1:length(sample_peaks_avg_for_bars@mass)) {
                    # If there is a match...
                    if ((abs(sample_peaks_avg_for_bars@mass[z] - as.numeric(f))*10^6/as.numeric(f)) <= ifelse(tof_mode == "linear", 2000, 200)) {
                        # Add the intensity of this peak to the y coordinates of the bars
                        coordinates_of_bars$x = append(coordinates_of_bars$x,sample_peaks_avg_for_bars@mass[z])
                        coordinates_of_bars$y = append(coordinates_of_bars$y,sample_peaks_avg_for_bars@intensity[z])
                        # Set the presence in the peaklist to true
                        presence_in_the_avg <- TRUE
                        # Break the for cycle to avoid duplicates and to continue
                        break
                    }
                }
                if (presence_in_the_avg == FALSE) {
                    # If the feature is not in the peaklist, scroll the datapoints in the spectrum
                    for (j in 1:length(spectra@mass)) {
                        # If there is a match...
                        if ((abs(spectra@mass[j]-as.numeric(f))*10^6/as.numeric(f)) <= ifelse(tof_mode == "linear", 2000, 200)) {
                            # Add the intensity of this peak to the y coordinates of the bars
                            coordinates_of_bars$x = append(coordinates_of_bars$x, spectra@mass[j])
                            coordinates_of_bars$y = append(coordinates_of_bars$y, spectra@intensity[j])
                            # Break the for cycle to avoid duplicates and to continue
                            break
                        }
                    }
                }
            }
            if (plot_figures == TRUE) {
                plot(spectra, xlab = "m/z", ylab = "Intensity (a.i.)")
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = spectra@metaData$file[[1]], xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = model_name, xjust = 0.5, yjust = 0.5)
                }
                # Draw the bars
                for (s in 1:length(coordinates_of_bars$x)) {
                    # Vertical bars (x,y x,y)
                    segments(coordinates_of_bars$x[s], 0, coordinates_of_bars$x[s], coordinates_of_bars$y[s], col = "red", lwd = 2)
                    # Horizontal segments(x,y , x,y)
                    segments(coordinates_of_bars$x[s] - 20, 0, coordinates_of_bars$x[s] + 20, 0, col = "red", lwd = 2)
                    segments(coordinates_of_bars$x[s] - 20, coordinates_of_bars$y[s], coordinates_of_bars$x[s] + 20, coordinates_of_bars$y[s], col = "red", lwd = 2)
                }
                average_spectrum_with_bars <- recordPlot()
            } else {
                average_spectrum_with_bars <- NULL
            }
        } else {
            ### Return NULL in case of incompatibility
            result_matrix <- NULL
            average_spectrum_with_bars <- NULL
        }
        ### Return
        return(list(result_matrix_model = result_matrix, average_spectrum_with_bars = average_spectrum_with_bars, final_result_matrix = final_result_matrix))
    }
}





################################################################################





############ CLASSIFICATION: PIXEL-BY-PIXEL AND/OR PROFILE (MULTICORE, ENSEMBLE)
# The function takes a folder in which there are imzML files (one for each patient) or an imzML file or a list of MALDIquant spectra files, the R workspace containing the models with the name of the model objects in the workspace, and allows the user to specify something regarding the preprocessing of the spectra to be classified.
# The function outputs a list containing: a matrix with the classification (pixel-by-pixel and/or profile), MS images with the pixel-by-pixel classification, a matrix with the ensemble classification (pixel-by-pixel and/or profile), MS images with the pixel-by-pixel ensemble classification and the plot of the average spectrum with red bars to indicate the signals used for classification.
# Parallel computation implemented.
# It outputs NULL values if the classification cannot be performed due to incompatibilities between the model features and the spectral features.
spectral_classification <- function(spectra_path, filepath_R, model_list_object = "model_list", classification_mode = c("pixel", "profile"), peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = FALSE, preprocessing_parameters = list(mass_range = c(4000,15000), transformation_algorithm = NULL, smoothing_algorithm = "SavitzkyGolay", smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 100, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_method = NULL), tof_mode = "linear", allow_parallelization = FALSE, decision_method_ensemble = "majority", vote_weights_ensemble = "equal", pixel_grouping = c("single", "moving window average", "graph", "hca"), moving_window_size = 5, number_of_hca_nodes = 10, number_of_spectra_partitions_graph = 1, partitioning_method_graph = "space", correlation_method_for_adjacency_matrix = "pearson", correlation_threshold_for_adjacency_matrix = 0.95, pvalue_threshold_for_adjacency_matrix = 0.05, max_GA_generations = 10, iterations_with_no_change_GA = 5, seed = 12345, plot_figures = TRUE, plot_graphs = TRUE, plot_legends = c("sample name", "legend", "plot name")) {
    ### Install and load the required packages
    install_and_load_required_packages(c("MALDIquant", "MALDIquantForeign","stats", "parallel", "kernlab", "MASS", "klaR", "pls", "randomForest", "lda", "caret", "nnet"))
    ### Defaults
    if (pixel_grouping == "") {
        pixel_grouping <- "single"
    }
    if (length(classification_mode) == 1 && classification_mode == "") {
        classification_mode <- "profile"
    }
    ### Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ### TOF-MODE
    if (tof_mode == "linear") {
        tolerance_ppm <- 2000
    } else if (tof_mode == "reflector" || tof_mode == "reflectron") {
        tolerance_ppm <- 200
    }
    ### Mass range
    mass_range <- preprocessing_parameters$mass_range
    ########## List the imzML files in the selected folder (if the path provided is a folder): check if its is folder, imzML file or spectra list
    ## Multiple imzML file path provided (folder)
    if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed = TRUE)) == 0) {
        filepath_test_imzml <- read_spectra_files(spectra_path, spectra_format = "imzml", full_path = TRUE)
    } else if (!is.list(spectra_path) && length(grep(".imzML", spectra_path, fixed = TRUE)) > 0) {
        ## Single imzML filepath
        filepath_test_imzml <- spectra_path
    } else if (isMassSpectrumList(spectra_path) || isMassSpectrum(spectra_path)) {
        ## Spectra list
        # Assign a string to the filepath for the future IF cycles
        filepath_test_imzml <- "List of spectra"
        sample_spectra <- spectra_path
    }
    ########## Global OUTPUT Initialization
    ### Pixel by pixel
    final_result_matrix_msi_list <- list()
    classification_ms_images_list <- list()
    classification_ensemble_matrix_msi_all <- list()
    classification_ensemble_ms_image_list <- list()
    ### Profile
    final_result_matrix_profile_list <- list()
    classification_ensemble_matrix_profile_all <- list()
    average_spectrum_with_bars_profile_list <- list()
    ########### SPECTRA: Process the sample to be classified
    ##### For each imzML file...
    for (p in 1:length(filepath_test_imzml)) {
        ### Output initialization (patient)
        # Pixel by pixel
        final_result_matrix_msi_patient <- NULL
        classification_ms_images_patient <- list()
        # Profile
        final_result_matrix_profile_patient <- NULL
        average_spectrum_with_bars_patient <- list()
        ### Import the spectra
        if (!isMassSpectrumList(spectra_path) || !isMassSpectrum(spectra_path)) {
            # Import the spectra (one imzML at a time)
            if (!is.null(mass_range)) {
                sample_spectra <- importImzMl(filepath_test_imzml[p], massRange = mass_range)
            } else {
                sample_spectra <- importImzMl(filepath_test_imzml[p])
            }
        }
        ### Rearrange the spectra according to the space coordinates (for reproducibility purposes)
        #sample_spectra <- rearrange_spectral_dataset(sample_spectra, rearranging_method = "space")
        ### Replace the sample name (path) with the actual sample name
        sample_spectra <- replace_sample_name(sample_spectra, spectra_format = "imzml", allow_parallelization = allow_parallelization)
        ### Retrieve the sample name
        sample_name <- sample_spectra[[1]]@metaData$file[1]
        ### Preprocess spectra
        if (!is.null(preprocessing_parameters) && is.list(preprocessing_parameters) && length(preprocessing_parameters) > 0) {
            sample_spectra <- preprocess_spectra(sample_spectra, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
        }
        ### Average spectrum (protein profile)
        if ("profile" %in% classification_mode) {
            ## Generate the average spectrum
            if (isMassSpectrumList(sample_spectra)) {
                sample_spectra_avg <- averageMassSpectra(sample_spectra, method = "mean")
                ## Preprocess average spectrum
                if (!is.null(preprocessing_parameters) && is.list(preprocessing_parameters) && length(preprocessing_parameters) > 0) {
                    sample_spectra_avg <- preprocess_spectra(sample_spectra_avg, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
                }
            } else {
                sample_spectra_avg <- sample_spectra
            }
        }
        ### LOAD THE R WORKSPACE WITH THE MODEL LIST
        # Create a temporary environment
        temporary_environment <- new.env()
        # Load the workspace
        load(filepath_R, envir = temporary_environment)
        # Get the models (R objects) from the workspace
        model_list <- get(model_list_object, pos = temporary_environment)
        # Get the list of models
        list_of_models <- names(model_list)
        # For each model...
        for (md in 1:length(list_of_models)) {
            ### Pixel by pixel
            if ("pixel" %in% classification_mode) {
                # Perform the classification
                model_classification <- single_model_classification_of_spectra(spectra = sample_spectra, model_x = model_list[[md]], model_name = list_of_models[md], preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = deisotope_peaklist, peak_picking_SNR = 3, peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", tof_mode = tof_mode, allow_parallelization = allow_parallelization, pixel_grouping = pixel_grouping, number_of_hca_nodes = number_of_hca_nodes, moving_window_size = moving_window_size, final_result_matrix = final_result_matrix_msi_patient, seed = seed, correlation_method_for_adjacency_matrix = correlation_method_for_adjacency_matrix, correlation_threshold_for_adjacency_matrix = correlation_threshold_for_adjacency_matrix, pvalue_threshold_for_adjacency_matrix = pvalue_threshold_for_adjacency_matrix, max_GA_generations = max_GA_generations, iterations_with_no_change = iterations_with_no_change_GA, number_of_spectra_partitions = number_of_spectra_partitions_graph, partitioning_method = partitioning_method_graph, plot_figures = plot_figures, plot_graphs = plot_graphs, plot_legends = plot_legends)
                # MSI classification
                if (plot_figures == TRUE) {
                    classification_ms_images_model <- model_classification$classification_msi_model
                } else {
                    classification_ms_images_model <- NULL
                }
                # Add the model result matrix to the final matrix (the classification function automatically attach the result to a result matrix if provided)
                final_result_matrix_msi_patient <- model_classification$final_result_matrix
                # Append the classification image list to the final list
                classification_ms_images_patient[[(list_of_models[md])]] <- classification_ms_images_model
            }
            if ("profile" %in% classification_mode) {
                # Perform the classification
                model_classification_profile <- single_model_classification_of_spectra(spectra = sample_spectra_avg, model_x = model_list[[md]], model_name = list_of_models[md], preprocessing_parameters = NULL, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = deisotope_peaklist, peak_picking_SNR = 3, peak_filtering_frequency_threshold_percent = NULL, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", tof_mode = tof_mode, allow_parallelization = allow_parallelization, pixel_grouping = pixel_grouping, number_of_hca_nodes = number_of_hca_nodes, moving_window_size = moving_window_size, final_result_matrix = final_result_matrix_profile_patient, seed = seed, correlation_method_for_adjacency_matrix = correlation_method_for_adjacency_matrix, correlation_threshold_for_adjacency_matrix = correlation_threshold_for_adjacency_matrix, pvalue_threshold_for_adjacency_matrix = pvalue_threshold_for_adjacency_matrix, max_GA_generations = max_GA_generations, iterations_with_no_change = iterations_with_no_change_GA, number_of_spectra_partitions = number_of_spectra_partitions_graph, partitioning_method = partitioning_method_graph, plot_figures = plot_figures, plot_graphs = plot_graphs, plot_legends = plot_legends)
                # Plot AVG spectrum
                if (plot_figures == TRUE) {
                    average_spectrum_profile_with_bars_model <- model_classification_profile$average_spectrum_with_bars
                } else {
                    average_spectrum_profile_with_bars_model <- NULL
                }
                # Add the model result matrix to the final matrix (the classification function automatically attach the result to a result matrix if provided)
                final_result_matrix_profile_patient <- model_classification_profile$final_result_matrix
                # Append the classification image list to the final list
                average_spectrum_with_bars_patient[[(list_of_models[md])]] <- average_spectrum_profile_with_bars_model
            }
        }
        ### Add to the global final patient list
        if ("pixel" %in% classification_mode) {
            ##### Store the matrix with the classification from all the models in the final list of matrices
            final_result_matrix_msi_list[[sample_name]] <- final_result_matrix_msi_patient
            ##### Store all the MS images in the element of the final list of MS images
            classification_ms_images_list[[sample_name]] <- classification_ms_images_patient
        }
        if ("profile" %in% classification_mode) {
            ##### Store the matrix with the classification
            final_result_matrix_profile_list[[sample_name]] <- final_result_matrix_profile_patient
            ##### Store the plots with the AVG spectrum with bars
            average_spectrum_with_bars_profile_list[[sample_name]] <- average_spectrum_with_bars_patient
        }
        ######################################## ENSEMBLE VOTE
        ### The ensemble classification can be possible only if: there is the classifiation matrix, there are at least 3 models and if the classes/outcomes are the same for each model
        classes_are_the_same_for_each_model <- TRUE
        for (md in 1:(length(model_list) - 1)) {
            if (isTRUE(model_list[[md]]$class_list != model_list[[md + 1]]$class_list)) {
                classes_are_the_same_for_each_model <- FALSE
                break
            }
        }
        outcomes_are_the_same_for_each_model <- TRUE
        for (md in 1:(length(model_list) - 1)) {
            if (isTRUE(model_list[[md]]$outcome_list != model_list[[md + 1]]$outcome_list)) {
                outcomes_are_the_same_for_each_model <- FALSE
                break
            }
        }
        ########## Ensemble results
        if ("pixel" %in% classification_mode) {
            if (length(list_of_models) > 2 && !is.null(final_result_matrix_msi_patient) && classes_are_the_same_for_each_model == TRUE && outcomes_are_the_same_for_each_model == TRUE) {
                ### Classification matrix
                classification_ensemble_matrix_msi <- ensemble_vote_classification(classification_matrix = final_result_matrix_msi_patient, class_list = model_list[[1]]$class_list, decision_method = decision_method_ensemble, vote_weights = vote_weights_ensemble)
                # Store the ensemble classification matrix in the final output list
                classification_ensemble_matrix_msi_all[[sample_name]] <- classification_ensemble_matrix_msi
                ### Molecular image of the classification
                # Generate the "predicted classes" vector from the ensemble classification matrix
                predicted_classes <- as.character(classification_ensemble_matrix_msi)
                # Define the class as number depending on the outcome
                outcome_and_class <- outcome_and_class_to_MS(class_list = model_list[[1]]$class_list, outcome_list = model_list[[1]]$outcome_list, class_vector = predicted_classes)
                # Replace the spectra intensities with the class number for plotting purposes
                class_as_number <- outcome_and_class$class_vector_as_numeric
                spectra_for_plotting <- sample_spectra
                for (s in 1:length(spectra_for_plotting)) {
                    spectra_for_plotting[[s]]@intensity <- rep(class_as_number[s], length(spectra_for_plotting[[s]]@intensity))
                }
                # Generate the MS images
                slices <- msiSlices(spectra_for_plotting, center = spectra_for_plotting[[1]]@mass[(length(spectra_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
                plotMsiSlice(slices, legend = FALSE, scale = F)
                # Define the legend
                if ("legend" %in% plot_legends) {
                    legend_text <- outcome_and_class$legend_text
                    legend_fill <- outcome_and_class$legend_fill
                    legend(x = "bottomright", legend = legend_text, fill = legend_fill, xjust = 0.5, yjust = 0.5)
                }
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = spectra_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = "Ensemble classifier", xjust = 0.5, yjust = 0.5)
                }
                # Store the plot into the list of images
                classification_ensemble_ms_image_list[[sample_name]] <- recordPlot()
            } else {
                classification_ensemble_matrix_msi_all[[sample_name]] <- NULL
                classification_ensemble_ms_image_list[[sample_name]] <- NULL
            }
        }
        if ("profile" %in% classification_mode) {
            if (length(list_of_models) > 2 && !is.null(final_result_matrix_profile_patient) && classes_are_the_same_for_each_model == TRUE && outcomes_are_the_same_for_each_model == TRUE) {
                ########## Ensemble results
                classification_ensemble_matrix_profile <- ensemble_vote_classification(classification_matrix = final_result_matrix_profile_patient, class_list = model_list[[1]]$class_list, decision_method = decision_method_ensemble, vote_weights = vote_weights_ensemble)
                # Store the ensemble classification matrix in the final output list
                if (is.null(classification_ensemble_matrix_profile_all)) {
                    classification_ensemble_matrix_profile_all <- classification_ensemble_matrix_profile
                } else {
                    classification_ensemble_matrix_profile_all <- rbind(classification_ensemble_matrix_profile_all, classification_ensemble_matrix_profile)
                }
            } else {
                classification_ensemble_matrix_profile_all <- NULL
            }
        }
    }
    # Return the results
    return(list(final_result_matrix_msi_list = final_result_matrix_msi_list, final_result_matrix_profile_list = final_result_matrix_profile_list, classification_ms_images_list = classification_ms_images_list, classification_ensemble_matrix_msi_all = classification_ensemble_matrix_msi_all, classification_ensemble_matrix_profile_all = classification_ensemble_matrix_profile_all, classification_ensemble_ms_image_list = classification_ensemble_ms_image_list, average_spectrum_with_bars_profile_list = average_spectrum_with_bars_profile_list))
}











################################################################################











































##################################################################### STATISTICS

##################### MATRIX/DATAFRAME SPLITTING FUNCTION (TRAINING AND TESTING)
# This function takes a peaklist matrix as input and outputs a list containing the part of the dataset used for the training and the part to be used as test dataset. The split happens onto the discriminant (class) variable.
# The split is made according to the user desire and according to a selected column (usually the class column).
# All the entries (rows) are considered independent observations.
matrix_splitting_training_test <- function(peaklist, discriminant_feature = "Class", seed = 12345, percentage_of_observations_for_training = 66) {
    # Install the required packages
    install_and_load_required_packages("caret")
    if (!is.null(discriminant_feature)) {
        ##### Determine the class list
        class_list <- levels(as.factor(peaklist[,discriminant_feature]))
    } else {
        class_list <- character()
    }
    ### If there are no classes or one class, just randomly select rows on the entire dataset
    if (length(class_list) <= 1) {
        # Plant the seed only if a specified value is entered
        if (!is.null(seed)) {
            # Make the randomness reproducible
            set.seed(seed)
        }
        index_training <- sample(nrow(peaklist), size = round(nrow(peaklist)*percentage_of_observations_for_training/100))
    } else if (length(class_list) > 1) {
        ### If there are more classes, randomly select the rows for each class
        # Create a list of dataframes, one for each class
        class_data_frame_list <- list()
        for (j in 1:length(class_list)) {
            class_data_frame_list[[j]] <- peaklist[peaklist[,discriminant_feature] == class_list[j],]
        }
        ### Create a list containing the training and testing indexes of the class dataframes
        index_training <- list()
        for (i in 1:length(class_data_frame_list)) {
            # Plant the seed only if a specified value is entered
            if (!is.null(seed)) {
                # Make the randomness reproducible
                set.seed(seed)
            }
            index_training[[i]] <- createDataPartition(y = class_data_frame_list[[i]][,discriminant_feature], p = (percentage_of_observations_for_training/100), list = FALSE)
        }
        names(index_training) <- class_list
    }
    ### Now we have randomly selected some patients for the training and some others for the testing, for each class
    ### The indexTraining is integer if there are no classes or one class, and rows are selected randomly
    if (is.integer(index_training)) {
        training_dataset <- peaklist[index_training,]
        test_dataset <- peaklist[-index_training,]
    } else if (is.list(index_training)) {
        # The indexTraining is a list if there are more than one classes (one element per class), and rows are selected randomly for each class
        ### Create the training and the testing patient dataframes
        training_data_frame_list <- list()
        test_data_frame_list <- list()
        # For each class (and so for each element of the list with the indexes per class)
        for (i in 1:length(index_training)) {
            # Create the two complementary dataframes: randomly select the patients for each class
            # for training and testing
            training_data_frame_list[[i]] <- class_data_frame_list[[i]][index_training[[i]],]
            test_data_frame_list[[i]] <- class_data_frame_list[[i]][-index_training[[i]],]
        }
        ### Merge the two traininf and test sets
        training_dataset <- NULL
        test_dataset <- NULL
        for (p in 1:length(training_data_frame_list)) {
            if (is.null(training_dataset)) {
                training_dataset <- training_data_frame_list[[p]]
            } else {
                training_dataset <- rbind(training_dataset, training_data_frame_list[[p]])
            }
        }
        for (p in 1:length(test_data_frame_list)) {
            if (is.null(test_dataset)) {
                test_dataset <- test_data_frame_list[[p]]
            } else {
                test_dataset <- rbind(test_dataset, test_data_frame_list[[p]])
            }
        }
    }
    return(list(training_dataset = training_dataset, test_dataset = test_dataset))
}





################################################################################





############################################################## FEATURE SELECTION
# This function runs the feature selection algorithm onto the peaklist matrix, returning the peaklist without the redundant/non-informative features, the original peaklist and the list of selected features.
# The function allows for the use of several feature selection algorithms.
feature_selection <- function(peaklist, feature_selection_method = "ANOVA", features_to_select = 20, selection_method = "pls", selection_metric = "Kappa", correlation_method = "pearson", correlation_threshold = 0.75, auc_threshold = 0.7, cv_repeats_control = 5, k_fold_cv_control = 10, discriminant_attribute = "Class", non_features = c("Sample", "Class", "THY"), seed = NULL, automatically_select_features = FALSE, generate_plots = TRUE, preprocessing = c("center","scale"), allow_parallelization = FALSE, feature_reranking = FALSE) {
    # Load the required libraries
    install_and_load_required_packages(c("stats", "pROC"))
    if (allow_parallelization == TRUE) {
        ### PARALLEL BACKEND
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            install_and_load_required_packages("doMC")
            # Register the foreach backend
            registerDoMC(cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            install_and_load_required_packages("doParallel")
            # Register the foreach backend
            cl <- makeCluster(cpu_thread_number, type='PSOCK')
            registerDoParallel(cl)
        }
    }
    # Initialization
    feature_weights <- NULL
    variable_importance <- NULL
    fs_model_performance <- NULL
    fs_model <- NULL
    ######################################################################### ANOVA
    if (feature_selection_method == "ANOVA") {
        feature_list <- names(peaklist[,!(names(peaklist) %in% non_features)])
        features_ANOVA <- character()
        # For each feature, calculate the impact on the classification capability by fitting an ANOVA
        for (feature in feature_list) {
            # Isolate the peak intensities
            intensity_vector <- peaklist[, feature]
            # Compute an ANOVA based upon the discriminant attribute
            feature_ANOVA <- aov(intensity_vector ~ peaklist[, discriminant_attribute])
            # Extract the p-value
            feature_ANOVA_pvalue <- summary(feature_ANOVA)[[1]]$"Pr(>F)"[1]
            # Keep the feature if with a great impact
            if (feature_ANOVA_pvalue <= 0.05) {
                features_ANOVA <- append(features_ANOVA, feature)
            }
        }
        # Predictors
        predictors_feature_selection <- features_ANOVA
    }
    ################################################################ KRUSKAL-WALLIS
    if (feature_selection_method == "kruskal" || feature_selection_method == "kruskal-wallis" || feature_selection_method == "Kruskal-Wallis") {
        feature_list <- names(peaklist[,!(names(peaklist) %in% non_features)])
        features_kruskal <- character()
        # For each feature, calculate the impact on the classification capability by fitting an ANOVA
        for (feature in feature_list){
            # Isolate the peak intensities
            intensity_vector <- peaklist [,feature]
            # Compute an ANOVA based upon the discriminant attribute
            feature_kruskal <- kruskal.test(intensity_vector ~ peaklist[,discriminant_attribute])
            # Extract the p-value
            feature_kruskal_pvalue <- feature_kruskal$p.value
            # Keep the feature if with a great impact
            if (feature_kruskal_pvalue <= 0.05) {
                features_kruskal <- append(features_kruskal, feature)
            }
        }
        # Predictors
        predictors_feature_selection <- features_kruskal
    }
    ################################################# CORRELATION FEATURE SELECTION
    if (feature_selection_method == "correlation") {
        # Take only the part of the matrix without Class and Sample
        peaklist_features <- peaklist [,!(names(peaklist) %in% non_features)]
        # Compute the correlation
        feature_correlation <- cor(peaklist_features, method = correlation_method)
        # Output the highly correlated features
        highly_correlated <- findCorrelation(feature_correlation, correlation_threshold)
        # List the highly correlated features
        highly_correlated_features <- names(peaklist_features[,highly_correlated])
        # Features to keep
        low_correlation_features <- names(peaklist_features[,-highly_correlated])
        print(paste("The number of selected features is", length(low_correlation_features), "out of", length(peaklist_features)))
        # Predictors
        predictors_feature_selection <- low_correlation_features
    }
    #################################################################### IMPORTANCE
    if (feature_selection_method == "importance") {
        feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
        # For each feature, calculate the impact on the classification capability
        model_control <- trainControl(method = "repeatedcv", number = k_fold_cv_control, repeats = cv_repeats_control)
        # Compute a model based upon the discriminant attribute
        feature_model <- train(x = peaklist [,!(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], method = "pls", preProcess = "scale", trControl = model_control)
        # Estimate variable importance
        feature_importance <- varImp(feature_model, scale = FALSE)
        # Isolate the most important features (rank them according to their importance first!)
        feature_importance_df <- feature_importance$importance
        feature_importance_df$Features <- rownames(feature_importance_df)
        feature_importance_df <- feature_importance_df [order(-feature_importance_df$Overall),]
        # Predictors
        predictors_feature_selection <- feature_importance_df$Features [1:features_to_select]
    }
    ########################################################################### ROC
    if (feature_selection_method == "ROC" || feature_selection_method == "roc") {
        feature_list <- names(peaklist [,!(names(peaklist) %in% non_features)])
        # List of important features
        features_ROC <- character()
        ##### Automatically select features
        if (automatically_select_features == TRUE) {
            # For each feature, calculate the impact on the classification capability by computing a ROC
            for (feature in feature_list){
                # Compute the ROC of the feature
                feature_ROC <- roc(response = peaklist[,discriminant_attribute], predictor = peaklist[,feature])
                # Extract the AUC
                feature_ROC_AUC <- feature_ROC$auc
                # Keep the feature if with a great impact
                if (feature_ROC_AUC >= auc_threshold) {
                    features_ROC <- append(features_ROC, feature)
                }
            }
        } else {
            ##### Select the most N important features
            # For each feature, calculate the impact on the classification capability by computing a ROC
            feature_ROC_vector <- numeric(length = length(feature_list))
            names(feature_ROC_vector) <- feature_list
            for (f in 1:length(feature_list)) {
                # Compute the ROC of the feature
                feature_ROC <- roc(response = peaklist[,discriminant_attribute], predictor = peaklist[,feature_list[f]])
                # Append the AUC to a vector
                feature_ROC_vector[f] <- feature_ROC$auc
            }
            # Sort the vector
            feature_ROC_vector_sorted <- sort(feature_ROC_vector, decreasing = TRUE)
            # Take the first N features
            features_ROC <- names(feature_ROC_vector_sorted)[1:features_to_select]
        }
        # Predictors
        predictors_feature_selection <- features_ROC
    }
    ########################################################################## Plot
    if (generate_plots == TRUE) {
        # Initialize
        feature_selection_graphics <- NULL
        if (feature_selection_method == "correlation") {
        }
        if (feature_selection_method == "importance") {
            feature_selection_graphics <- plot(feature_importance)
        }
    } else {feature_selection_graphics <- NULL}
    #################################################### Take the selected features
    peaklist_feature_selection <- peaklist [,predictors_feature_selection]
    # Add the non features back
    for (i in 1:length(non_features)) {
        peaklist_feature_selection <- cbind(peaklist_feature_selection, peaklist[,non_features[i]])
    }
    names(peaklist_feature_selection) <- c(as.character(predictors_feature_selection), non_features)
    # Return the values
    return(list(peaklist_feature_selection = peaklist_feature_selection, predictors_feature_selection = predictors_feature_selection, feature_selection_graphics = feature_selection_graphics, feature_weights = feature_weights, variable_importance = variable_importance, fs_model_performance = fs_model_performance, fs_model = fs_model))
}





################################################################################





##################### EMBEDDED FEATURE SELECTION (RECURSIVE FEATURE ELIMINATION)
# This function runs the feature selection algorithm onto the peaklist matrix, returning the peaklist without the redundant/non-informative features, the original peaklist and the list of selected features, along with the model used for selecting the features.
# The function allows for the use of several feature selection algorithms.
embedded_rfe <- function(peaklist, features_to_select = 20, selection_method = "pls", model_tuning = TRUE, model_tuning_mode = c("embedded", "after"), model_tune_grid = list(), selection_metric = "Accuracy", cv_repeats_control = 5, k_fold_cv_control = 10, discriminant_attribute = "Class", non_features = c("Sample", "Class"), seed = NULL, automatically_select_features = TRUE, generate_plots = TRUE, preprocessing = c("center","scale"), allow_parallelization = FALSE, feature_reranking = TRUE, external_peaklist = NULL, positive_class_cv = "HP") {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "stats", "pROC", "nnet", "e1071", "kernlab", "randomForest", "klaR", "MASS", "pls", "iterators", "nnet", "SparseM", "stringi"))
    ### Define better the (helper) functions to be used in rfe control
    rfe_helper_functions <- list(summary = defaultSummary,
                                 fit = function(x, y, ...) {
                                     library("e1071")
                                     svm(x, y, ...)
                                 },
                                 pred = function(object, x) predict.svm(object, nwdata = x),
                                 rank = function(object, x, y) {
                                     vimp <- varImp(object)
                                     vimp <- vimp[order(vimp$Overall, decreasing = TRUE),, drop = FALSE]
                                     vimp$var <- rownames(vimp)
                                     vimp
                                 },
                                 selectSize = pickSizeBest,
                                 selectVar = pickVars
    )
    # Store th helper functions in the control function
    #rfe_ctrl$functions <- rfe_helper_functions
    if (allow_parallelization == TRUE) {
        ### PARALLEL BACKEND
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            install_and_load_required_packages("doMC")
            # Register the foreach backend
            registerDoMC(cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            install_and_load_required_packages("doParallel")
            # Register the foreach backend
            cl <- makeCluster(cpu_thread_number, type='PSOCK')
            registerDoParallel(cl)
        }
    }
    ### Initialization
    feature_weights <- NULL
    variable_importance <- NULL
    fs_model_performance <- NULL
    pie_chart_classification <- NULL
    model_roc <- NULL
    ### The simulation will fit models with subset sizes: (the subset size is the number of predictors to use)
    if (automatically_select_features == TRUE) {
        # Define the feature subset sizes
        subset_sizes <- seq(2, features_to_select, by = 1)
    } else {
        # Define the feature subset sizes
        subset_sizes <- features_to_select
    }
    # Plant the seed only if a specified value is entered
    if (!is.null(seed)) {
        # Make the randomness reproducible
        set.seed(seed)
    }
    ### Define the control function of the RFE
    rfe_ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, saveDetails = TRUE, allowParallel = allow_parallelization, rerank = feature_reranking, seeds = NULL)
    ### Run the RFE (model tuning during feature selection or after)
    if (model_tuning == TRUE && model_tuning_mode == "embedded") {
        rfe_model <- rfe(x = peaklist[, !(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, method = selection_method, metric = selection_metric, preProcess = preprocessing, tuneGrid = expand.grid(model_tune_grid))
    } else if (model_tuning == FALSE || (model_tuning == TRUE && model_tuning_mode == "after")) {
        rfe_model <- rfe(x = peaklist[, !(names(peaklist) %in% non_features)], y = peaklist[,discriminant_attribute], sizes = subset_sizes, rfeControl = rfe_ctrl, method = selection_method, metric = selection_metric, preProcess = preprocessing)
    }
    # Extract the model
    fs_model <- rfe_model$fit$finalModel
    # Variable importance
    variable_importance <- varImp(rfe_model$fit)
    # Feature weights
    feature_weights <- rfe_model$fit
    # Model performances
    if (selection_metric == "kappa" || selection_metric == "Kappa") {
        fs_model_performance <- as.numeric(max(rfe_model$fit$results$Kappa, na.rm = TRUE))
        names(fs_model_performance) <- "Kappa"
    } else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
        fs_model_performance <- as.numeric(max(rfe_model$fit$results$Accuracy, na.rm = TRUE))
        names(fs_model_performance) <- "Accuracy"
    }
    if (automatically_select_features == TRUE) {
        # Output the best predictors after the RFE
        predictors_rfe <- predictors(rfe_model)
    } else {
        # Output the best predictors after the RFE
        predictors_rfe <- predictors(rfe_model)[1:features_to_select]
    }
    # Predictors
    predictors_feature_selection <- predictors_rfe
    ##### Plots
    if (generate_plots == TRUE) {
        feature_selection_graphics <- plot(rfe_model, type = c("g","o"))
    } else {
        feature_selection_graphics <- NULL
    }
    #### Take the selected features
    peaklist_feature_selection <- peaklist[, predictors_feature_selection]
    ### Add the non features back
    for (i in 1:length(non_features)) {
        peaklist_feature_selection <- cbind(peaklist_feature_selection, peaklist[, non_features[i]])
    }
    names(peaklist_feature_selection) <- c(as.character(predictors_feature_selection), non_features)
    #### Tune the model after the features have been selected
    if (model_tuning == TRUE && model_tuning_mode == "after" && !is.null(model_tune_grid)) {
        # Make the randomness reproducible
        if (!is.null(seed)) {
            set.seed(seed)
        }
        # Define the control function
        train_ctrl <- trainControl(method = "repeatedcv", repeats = cv_repeats_control, number = k_fold_cv_control, allowParallel = allow_parallelization, seeds = NULL)
        # Define the model tuned
        fs_model <- train(x = peaklist_feature_selection[, !(names(peaklist_feature_selection) %in% non_features)], y = peaklist_feature_selection[, discriminant_attribute], method = selection_method, preProcess = preprocessing, tuneGrid = expand.grid(model_tune_grid), trControl = train_ctrl, metric = selection_metric)
        # Model performances
        if (selection_metric == "kappa" || selection_metric == "Kappa") {
            fs_model_performance <- as.numeric(max(fs_model$results$Kappa, na.rm = TRUE))
            names(fs_model_performance) <- "Kappa"
        } else if (selection_metric == "accuracy" || selection_metric == "Accuracy") {
            fs_model_performance <- as.numeric(max(fs_model$results$Accuracy, na.rm = TRUE))
            names(fs_model_performance) <- "Accuracy"
        }
    }
    if (allow_parallelization == TRUE) {
        # Close the parallelization cluster (on Windows)
        if (Sys.info()[1] == "Windows") {
            stopCluster(cl)
        }
    }
    ### External validation
    if (!is.null(external_peaklist) && (is.matrix(external_peaklist) || is.data.frame(external_peaklist))) {
        ## Use the model to predict the outcome of the testing set (the new data must have only the predictors)
        # Plant the seed only if a specified value is entered
        if (!is.null(seed)) {
            # Make the randomness reproducible
            set.seed(seed)
        }
        # Class prediction (with the same features!)
        external_peaklist <- external_peaklist[, names(peaklist_feature_selection)]
        predicted_classes_model <- predict(fs_model, newdata = external_peaklist[,!(names(external_peaklist) %in% non_features)])
        # Create the outcomes dataframe
        classification_results_model <- data.frame(Sample = external_peaklist$Sample, Predicted = predicted_classes_model, True = external_peaklist$Class)
        # Generate the confusion matrix to evaluate the performances
        if (is.null(positive_class_cv) || !(positive_class_cv %in% levels(as.factor(peaklist[, discriminant_attribute])))) {
            positive_class_cv <- levels(as.factor(peaklist[, discriminant_attribute]))[1]
        }
        test_performances_model <- confusionMatrix(data = predicted_classes_model, external_peaklist$Class, positive = positive_class_cv)
        ## ROC analysis
        model_roc <- list()
        roc_curve <- roc(response = as.numeric(classification_results_model$True), predictor = as.numeric(classification_results_model$Predicted))
        model_roc[[1]] <- roc_curve$auc
        if (generate_plots == TRUE) {
            plot(roc_curve)
            roc_legend <- paste("ROC area under the curve:", roc_curve$auc)
            legend("bottomright", legend = roc_legend, xjust = 0.5, yjust = 0.5)
            model_roc[[2]] <- recordPlot()
        } else {
            model_roc[[2]] <- NULL
        }
        ### Pie chart classification
        correctly_classified <- 0
        misclassified <- 0
        for (i in 1:nrow(classification_results_model)) {
            if (as.character(classification_results_model$Predicted[i]) == as.character(classification_results_model$True[i])) {
                correctly_classified <- correctly_classified + 1
            } else {
                misclassified <- misclassified + 1
            }
        }
        classification_pie <- c(correctly_classified, misclassified)
        if (generate_plots == TRUE) {
            pie(x = classification_pie, labels = c("Correctly classified", "Misclassified"), col = c("green","blue"))
            pie_chart_classification <- recordPlot()
        } else {
            pie_chart_classification <- NULL
        }
    }
    ### Return the values
    return(list(peaklist_feature_selection = peaklist_feature_selection, predictors_feature_selection = predictors_feature_selection, feature_selection_graphics = feature_selection_graphics, feature_weights = feature_weights, variable_importance = variable_importance, fs_model_performance = fs_model_performance, feature_selection_model = fs_model, class_list = levels(as.factor(peaklist[, discriminant_attribute])), pie_chart_classification = pie_chart_classification, model_roc = model_roc))
}





################################################################################





########### AUTOMATED EMBEDDED FEATURE SELECTION (RECURSIVE FEATURE ELIMINATION)
# This function iteratively runs the embedded-rfe feature selection function onto the same input objects as that function, in order to find the best combination of parameters (preprocessing, feature reranking) for the feature selection. It returns the same elements of the embedded_rfe function, but the best chosen after trying all of the parameter combinations.
# The function allows for the use of several feature selection algorithms.
automated_embedded_rfe <- function(peaklist, features_to_select = 20, selection_method = "pls", model_tuning = TRUE, model_tuning_mode = c("embedded", "after"), model_tune_grid = data.frame(ncomp = 1:5), selection_metric = "Accuracy", cv_repeats_control = 5, k_fold_cv_control = 10, discriminant_attribute = "Class", non_features = c("Sample", "Class"), seed = NULL, automatically_select_features = TRUE, generate_plots = TRUE, preprocessing = c("center","scale"), allow_parallelization = FALSE, feature_reranking = FALSE, external_peaklist = NULL, positive_class_cv = "HP", try_combination_of_parameters = TRUE) {
    # Load the required libraries
    install_and_load_required_packages(c("caret", "stats", "pROC", "nnet", "e1071", "kernlab", "randomForest", "klaR", "MASS", "pls", "iterators", "nnet", "SparseM", "stringi"))
    if (try_combination_of_parameters == TRUE) {
        ### Establish the combination of parameters (preprocessing + feature reranking) to establish the best model
        # Inizialize the output of combination of parameters
        parameter_combination <- list()
        # Define the parameters to be tested
        preprocessing_values <- list(NULL, c("center", "scale"))
        feature_reranking_values <- list(TRUE, FALSE)
        # Generate the combination list (each list element is a combination of values)
        for (p in 1:length(preprocessing_values)) {
            for (f in 1:length(feature_reranking_values)) {
                parameter_combination[[(length(parameter_combination) + 1)]] <- list(preprocessing = preprocessing_values[[p]], feature_reranking = feature_reranking_values[[f]])
            }
        }
        ### Test every combination...
        # Store the best performance value and the best model
        best_model_performance <- NULL
        best_rfe_model <- NULL
        # Run every combination, storing the result if good
        for (comb in 1:length(parameter_combination)) {
            single_rfe_model <- embedded_rfe(peaklist, features_to_select = features_to_select, selection_method = selection_method, model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = model_tune_grid, selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = parameter_combination[[comb]]$preprocessing, allow_parallelization = allow_parallelization, feature_reranking = parameter_combination[[comb]]$feature_reranking, external_peaklist = external_peaklist, positive_class_cv = positive_class_cv)
            ### Check (and store) the performance values and the model
            if (is.null(best_model_performance) || single_rfe_model$fs_model_performance > best_model_performance) {
                best_model_performance <- single_rfe_model$fs_model_performance
                best_rfe_model <- single_rfe_model
            }
        }
    } else {
    ### NO combinations, run the single function...
        best_rfe_model <- embedded_rfe(peaklist, features_to_select = features_to_select, selection_method = selection_method, model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = model_tune_grid, selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, external_peaklist = external_peaklist, positive_class_cv = positive_class_cv)
    }
    ### Return the values
    return(best_rfe_model)
}





################################################################################





############### MODEL ENSEMBLE TRAINING/TUNING (WITH EMBEDDED FEATURE SELECTION)
# The function takes a peaklist as input (samples/spectra as rows and aligned peaks as columns) and performs feature selection coupled with model training and tuning according to the input parameters. It can automatically select the number of features according to the model performances in the tuning phase or select a fix number of features in the tuning phase: corss-validation is used to adjust the tuning parameters. The best model is selected and returned.
# Several models undergo training/tuning with a selected number of features (each model selects its own).
# The function returns a list of models (with the model object, the model ID, the performances, the outcome list, the class list and the features), a matrix listing the cross-validation performances for each model and the feature list (common and for each model).
model_ensemble_embedded_fs <- function(peaklist, features_to_select = 20, model_tuning = TRUE, model_tuning_mode = "after", selection_metric = "Accuracy", discriminant_attribute = "Class", non_features = c("Sample", "Class"), seed = 12345, automatically_select_features = FALSE, generate_plots = TRUE, cv_repeats_control = 5, k_fold_cv_control = 3, preprocessing = c("center", "scale"), allow_parallelization = TRUE, feature_reranking = FALSE, try_combination_of_parameters = FALSE, outcome_list = c("benign", "malignant")) {
    ##### Feature selection (embedded method)
    # Partial Least Squares
    pls_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "pls", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = data.frame(ncomp = 1:5), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    pls_model <- pls_model_rfe$feature_selection_model
    pls_model_features <- pls_model_rfe$predictors_feature_selection
    pls_model_class_list <- pls_model_rfe$class_list
    pls_model_ID <- "pls"
    pls_model_performance <- pls_model_rfe$fs_model_performance
    print("Partial Least Squares")
    print(pls_model_features)
    print(pls_model_performance)
    # Support Vector Machine (with Radial Basis Kernel function)
    svmRadial_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "svmRadial", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = list(sigma = 10^(-5:5), C = 10^(-5:5)), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    svmRadial_model <- svmRadial_model_rfe$feature_selection_model
    svmRadial_model_features <- svmRadial_model_rfe$predictors_feature_selection
    svmRadial_model_class_list <- svmRadial_model_rfe$class_list
    svmRadial_model_ID <- "svm"
    svmRadial_model_performance <- svmRadial_model_rfe$fs_model_performance
    print("Support Vector Machine (with Radial Basis Kernel function)")
    print(svmRadial_model_features)
    print(svmRadial_model_performance)
    # Support Vector Machine (with Polynomial Kernel function)
    svmPoly_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "svmPoly", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = list(C = 10^(-5:5), degree = 1:5, scale = 1), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    svmPoly_model <- svmPoly_model_rfe$feature_selection_model
    svmPoly_model_features <- svmPoly_model_rfe$predictors_feature_selection
    svmPoly_model_class_list <- svmPoly_model_rfe$class_list
    svmPoly_model_ID <- "svm"
    svmPoly_model_performance <- svmPoly_model_rfe$fs_model_performance
    print("Support Vector Machine (with Polynomial Kernel function)")
    print(svmPoly_model_features)
    print(svmPoly_model_performance)
    # Random Forest
    rf_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "rf", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = list(mtry = seq(1,5, by = 1)), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    rf_model <- rf_model_rfe$feature_selection_model
    rf_model_features <- rf_model_rfe$predictors_feature_selection
    rf_model_class_list <- rf_model_rfe$class_list
    rf_model_ID <- "rf"
    rf_model_performance <- rf_model_rfe$fs_model_performance
    print("Random Forest")
    print(rf_model_features)
    print(rf_model_performance)
    # Naive Bayes Classifier
    nbc_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "nb", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = data.frame(fL = seq(0, 1, by = 0.2), usekernel = c(TRUE, FALSE), adjust = c(TRUE, FALSE)), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    nbc_model <- nbc_model_rfe$feature_selection_model
    nbc_model_features <- nbc_model_rfe$predictors_feature_selection
    nbc_model_class_list <- nbc_model_rfe$class_list
    nbc_model_ID <- "nbc"
    nbc_model_performance <- nbc_model_rfe$fs_model_performance
    print("Naive Bayes Classifier")
    print(nbc_model_features)
    print(nbc_model_performance)
    # K-Nearest Neighbor
    knn_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "knn", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = list(k = seq(1,15, by = 1)), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    knn_model <- knn_model_rfe$feature_selection_model
    knn_model_features <- knn_model_rfe$predictors_feature_selection
    knn_model_class_list <- knn_model_rfe$class_list
    knn_model_ID <- "knn"
    knn_model_performance <- knn_model_rfe$fs_model_performance
    print("K-Nearest Neighbor")
    print(knn_model_features)
    print(knn_model_performance)
    # Neural Network
    nnet_model_rfe <- automated_embedded_rfe(peaklist = peaklist, features_to_select = features_to_select, selection_method = "nnet", model_tuning = model_tuning, model_tuning_mode = model_tuning_mode, model_tune_grid = list(size = seq(1,5, by = 1), decay = seq(0,2,by = 1)), selection_metric = selection_metric, cv_repeats_control = cv_repeats_control, k_fold_cv_control = k_fold_cv_control, discriminant_attribute = discriminant_attribute, non_features = non_features, seed = seed, automatically_select_features = automatically_select_features, generate_plots = generate_plots, preprocessing = preprocessing, allow_parallelization = allow_parallelization, feature_reranking = feature_reranking, try_combination_of_parameters = try_combination_of_parameters)
    nnet_model <- nnet_model_rfe$feature_selection_model
    nnet_model_features <- nnet_model_rfe$predictors_feature_selection
    nnet_model_class_list <- nnet_model_rfe$class_list
    nnet_model_ID <- "nnet"
    nnet_model_performance <- nnet_model_rfe$fs_model_performance
    print("Neural Network")
    print(nnet_model_features)
    print(nnet_model_performance)
    ##### Elements for the RData
    ### Each model contains: the model object, the class list, the outcome list, the list of features, the model ID and the model performances
    # Partial Least Squares
    PLS_model_list <- list(model = pls_model, class_list = pls_model_class_list, outcome_list = outcome_list, features_model = pls_model_features, model_ID = pls_model_ID, model_performance = pls_model_performance)
    # Support Vector Machine (with Radial Basis Kernel function)
    RSVM_model_list <- list(model = svmRadial_model, class_list = svmRadial_model_class_list, outcome_list = outcome_list, features_model = svmRadial_model_features, model_ID = svmRadial_model_ID, model_performance = svmRadial_model_performance)
    # Support Vector Machine (with Polynomial Kernel function)
    PSVM_model_list <- list(model = svmPoly_model, class_list = svmPoly_model_class_list, outcome_list = outcome_list, features_model = svmPoly_model_features, model_ID = svmPoly_model_ID, model_performance = svmPoly_model_performance)
    # Random Forest
    RF_model_list <- list(model = rf_model, class_list = rf_model_class_list, outcome_list = outcome_list, features_model = rf_model_features, model_ID = rf_model_ID, model_performance = rf_model_performance)
    # Naive Bayes Classifier
    NBC_model_list <- list(model = nbc_model, class_list = nbc_model_class_list, outcome_list = outcome_list, features_model = nbc_model_features, model_ID = nbc_model_ID, model_performance = nbc_model_performance)
    # K-Nearest Neighbor
    KNN_model_list <- list(model = knn_model, class_list = knn_model_class_list, outcome_list = outcome_list, features_model = knn_model_features, model_ID = knn_model_ID, model_performance = knn_model_performance)
    # Neural Network
    NNET_model_list <- list(model = nnet_model, class_list = nnet_model_class_list, outcome_list = outcome_list, features_model = nnet_model_features, model_ID = nnet_model_ID, model_performance = nnet_model_performance)
    ### Build the final model list (to be exported) (each element has the proper name of the model)
    model_list <- list("SVM Radial Basis" = RSVM_model_list, "SVM Polynomial" = PSVM_model_list, "Partial Least Squares" = PLS_model_list, "Random Forest" = RF_model_list, "Naive Bayes Classifier" = NBC_model_list, "K-Nearest Neighbor" = KNN_model_list, "Neural Network" = NNET_model_list)
    ### Build the final feature vector
    feature_list <- extract_feature_list_from_model_list(filepath_R = model_list, model_list_object = "model_list", features_to_return = features_to_select)
    ### Build the matrix with the model performances
    model_performance_matrix <- extract_performance_matrix_from_model_list(filepath_R = model_list, model_list_object = "model_list")
    ##### Return
    return(list(model_list = model_list, feature_list = feature_list, model_performance_matrix = model_performance_matrix))
}





################################################################################





#################### SPECTRAL TYPER SCORE ACCORDING TO THE HIERARCHICAL DISTANCE
# This function computes the Spectral Typer score by comparing the test spectra with the library spectra, determining the similarity (through the euclidean distance) and assigning a category according to the distance.
# Each sample gets compared with all the entries in the database, simultaneously.
spectral_typer_score_hierarchical_distance <- function (spectra_database, spectra_test, peaks_database, peaks_test, class_list_library = NULL, peaks_filtering_percentage_threshold = 5, low_intensity_percentage_threshold = 0, low_intensity_threshold_method = "element-wise", tolerance_ppm = 2000, spectra_path_output = TRUE, score_only = TRUE, spectra_format = "brukerflex", normalise_distances = TRUE, normalization_method = "sum") {
    # Load the required libraries
    install_and_load_required_packages(c("MALDIquant", "stats", "ggplot2", "ggdendro"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Sample and Library size
    if (isMassPeaksList(peaks_test)) {
        number_of_samples <- length(peaks_test)
    } else if (isMassPeaks(peaks_test)) {
        number_of_samples <- 1
    }
    if (isMassPeaksList(peaks_database)) {
        database_size <- length(peaks_database)
    } else if (isMassPeaks(peaks_database)) {
        database_size <- 1
    }
    ####### Peak alignment
    # Merge the peaklists and the spectra
    peaks_all <- append(peaks_database, peaks_test)
    spectra_all <- append(spectra_database, spectra_test)
    # Align
    peaks_all <- align_and_filter_peaks(peaks_all, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peaks_filtering_percentage_threshold, low_intensity_peak_removal_threshold_percent = low_intensity_percentage_threshold, low_intensity_peak_removal_threshold_method = low_intensity_threshold_method)
    # Restore the lists
    peaks_database <- peaks_all [1:database_size]
    peaks_test <- peaks_all [(database_size+1):length(peaks_all)]
    #### Replace the sample name, both in the library and in the test set
    peaks_test <- replace_sample_name(peaks_test, spectra_format = spectra_format)
    peaks_database <- replace_class_name(peaks_database,  class_list = class_list_library, spectra_format = spectra_format)
    ####### Create the sample vector
    sample_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    ####### Create the library vector
    database_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    # Generate the path vector
    spectra_path_vector <- character()
    for (spectrum in spectra_test) {
        spectra_path_vector <- append(spectra_path_vector, spectrum@metaData$file[1])
    }
    # Generate the matrix (for hca)
    peaklist_matrix <- intensityMatrix(peaks_all, spectra_all)
    # Add additional info to the matrix
    peaklist_matrix <- matrix_add_class_and_sample(peaklist_matrix, peaks = peaks_all, class_list = list(), spectra_format = spectra_format, sample_output = TRUE, class_output = FALSE)
    rownames(peaklist_matrix) <- make.names(peaklist_matrix[,"Sample"], unique = TRUE)
    # Compute the hca
    distance_matrix <- dist(peaklist_matrix[,1:(ncol(peaklist_matrix)-1)], method = "euclidean")
    hierarchical_clustering <- hclust(distance_matrix)
    #plot(hierarchical_clustering, main = "Hierarchical clustering analysis - Spectral Typer, xlab = "Samples", ylab = "Tree height")
    hca_dendrogram <- ggdendrogram(hierarchical_clustering, segments = TRUE, labels = TRUE, leaf_labels = TRUE, rotate = TRUE, theme_dendro = TRUE)#, main = "Hierarchical clustering analysis - Spectral Typer", xlab = "Samples", ylab = "Tree height")
    #hca_dendrogram <- recordPlot()
    #
    distance_matrix <- as.matrix(distance_matrix)
    # The distance matrix displays the distance between the spectra
    colnames(distance_matrix) <- peaklist_matrix[,"Sample"]
    rownames(distance_matrix) <- peaklist_matrix[,"Sample"]
    # Remove the first rows (the spectra from the database)
    distance_matrix <- distance_matrix[(database_size+1):nrow(distance_matrix),]
    # Keep only the first columns (the spectra from the database)
    distance_matrix <- distance_matrix[,1:database_size]
    ### Normalise the euclidean distances
    if (normalise_distances == TRUE) {
        # TIC (SUM)
        if (normalization_method == "sum") {
            # Compute the sum of the rows
            row_sums <- apply(distance_matrix, MARGIN = 1, FUN = sum)
            # Divide each element of the matrix by the sum of the row
            for (r in 1:nrow(distance_matrix)) {
                distance_matrix [r,] <- distance_matrix[r,] / row_sums[r]
            }
            # Multiply everything by 10, to have more readable results
            distance_matrix <- distance_matrix * 10
            # The classification is made by comparing the single sample spectrum with the spectrum of the database class (the distance is displayed in the distance matrix): the closer the better
            # Scroll the rows, assign the class based upon the distance, create the output matrix for results (create a function to apply to each matrix row)
            scoring_function <- function (x) {
                if (x < 1) {
                    x <- paste("YES (", round(as.numeric(x),3), ")", sep = "")
                } else if (x >= 1 && x < 1.2) {
                    x <- paste("NI (", round(as.numeric(x),3), ")", sep = "")
                } else if (x >= 1.2) {
                    x <- paste("NO (", round(as.numeric(x),3), ")", sep = "")
                }
                return (x)
            }
            result_matrix <- apply(distance_matrix, MARGIN = c(1,2), FUN = function(x) scoring_function(x))
        }
        # MAX
        if (normalization_method == "max") {
            # Divide each element of the matrix by the maximum of the row
            for (r in 1:nrow(distance_matrix)) {
                distance_matrix [r,] <- distance_matrix[r,] / max(distance_matrix[r,])
            }
            # Multiply everything by 100, to have more readable results (percentage of the max)
            distance_matrix <- distance_matrix * 100
            # The classification is made by comparing the single sample spectrum with the spectrum of the database class (the distance is displayed in the distance matrix): the closer the better
            # Scroll the rows, assign the class based upon the distance, create the output matrix for results (create a function to apply to each matrix row)
            scoring_function <- function (x) {
                if (x < 50) {
                    x <- paste("YES (", round(as.numeric(x),3), ")", sep = "")
                } else if (x >= 50 && x < 75) {
                    x <- paste("NI (", round(as.numeric(x),3), ")", sep = "")
                } else if (x >= 75) {
                    x <- paste("NO (", round(as.numeric(x),3), ")", sep = "")
                }
                return (x)
            }
            result_matrix <- apply(distance_matrix, MARGIN = c(1,2), FUN = function(x) scoring_function(x))
        }
    }
    # Spectra path
    if (spectra_path_output == TRUE) {
        result_matrix <- cbind(result_matrix, sample_vector)
    }
    return(list(result_matrix = result_matrix, hca_dendrogram = hca_dendrogram))
}





################################################################################





####################################### SPECTRAL TYPER SCORE: CORRELATION MATRIX
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity symmetry via the correlation matrix.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_correlation_matrix <- function(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library = NULL, peaks_filtering_percentage_threshold = 5, low_intensity_percentage_threshold = 0, low_intensity_threshold_method = "element-wise", tolerance_ppm = 2000, intensity_correction_coefficient = 1, spectra_format = "brukerflex", spectra_path_output = TRUE, score_only = FALSE, allow_parallelization = FALSE) {
    install_and_load_required_packages(c("MALDIquant", "corrplot", "weights", "stats", "parallel"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    # Rename the trim function to avoid conflicts
    trim_weights <- get(x = "trim", pos = "package:weights")
    ### Folder lists
    #database_folder_list <- dir(filepath_database, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    #test_folder_list <- dir(filepath_test, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    # Sample and Library size
    if (isMassPeaksList(peaks_test)) {
        number_of_samples <- length(peaks_test)
    } else if (isMassPeaks(peaks_test)) {
        number_of_samples <- 1
    }
    if (isMassPeaksList(peaks_database)) {
        database_size <- length(peaks_database)
    } else if (isMassPeaks(peaks_database)) {
        database_size <- 1
    }
    #### Replace the sample name, both in the library and in the test set
    peaks_test <- replace_sample_name(peaks_test, spectra_format = spectra_format)
    peaks_database <- replace_class_name(peaks_database, class_list = class_list_library, class_in_file_name = TRUE, spectra_format = spectra_format)
    ####### Create the sample vector
    sample_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    ####### Create the library vector
    database_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    # Generate the path vector
    spectra_path_vector <- character()
    for (sp in 1:length(spectra_test)) {
        spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
    }
    ############################################################ SCORE (FRI)
    number_of_signals_database <- numeric(length = database_size)
    for (d in 1:database_size) {
        number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
    }
    ################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
    # Create a list to be used for lapply
    global_list <- list()
    for (spl in 1:number_of_samples) {
        # Extract the peaklist and the spectrum
        peaks_database_temp <- peaks_database
        spectra_database_temp <- spectra_database
        peaks_sample <- peaks_test[[spl]]
        spectrum_sample <- spectra_test[[spl]]
        # Generate the entry of the global list
        global_list_entry <- append(peaks_sample, peaks_database_temp)
        global_list_entry <- append(global_list_entry, spectrum_sample)
        global_list_entry <- append(global_list_entry, spectra_database_temp)
        global_list[[spl]] <- global_list_entry
    }
    ############################################## Define the function for parLapply
    # x = each element of the global list
    comparison_sample_db_subfunction_correlation <- function(x) {
        # Generate the matrix rows for the output
        matching_signals_matrix <- matrix(0, nrow = 1, ncol = database_size)
        number_of_signals_database_matrix <- matrix(0, nrow = 1, ncol = database_size)
        fit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        retrofit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        intensity_correlation_matrix <- matrix(0, ncol = database_size, nrow = 1)
        pvalue_matrix <- matrix(0, ncol = database_size, nrow = 1)
        slope_matrix <- matrix(0, ncol = database_size, nrow = 1)
        ###### Compare with all the elements in the library
        ### For each entry in the library...
        for (db in 1:database_size) {
            # Extract the peaklist and the spectrum
            peaks_sample <- x[[1]]
            spectrum_sample <- x[[(database_size+2)]]
            peaks_database_temp_all <- x[2:(database_size+1)]
            spectra_database_temp_all <- x[(database_size+3):length(x)]
            peaks_database_temp <- peaks_database_temp_all[[db]]
            spectra_database_temp <- spectra_database_temp_all[[db]]
            ####### Peak alignment
            # Merge the peaklists
            peaks_all <- append(peaks_database_temp, peaks_sample)
            spectra_all <- append(spectra_database_temp, spectrum_sample)
            # Align the peaks
            peaks_all <- align_and_filter_peaks(peaks_all, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peaks_filtering_percentage_threshold, low_intensity_peak_removal_threshold_percent = low_intensity_percentage_threshold, low_intensity_peak_removal_threshold_method = low_intensity_threshold_method)
            # Restore the lists
            peaks_database_temp <- peaks_all[[1]]
            peaks_sample <- peaks_all[[2]]
            #################### Number of signals
            number_of_signals_samples <- length(peaks_sample@mass)
            number_of_signals_database <- length(peaks_database_temp@mass)
            number_of_signals_database_matrix[1,db] <- number_of_signals_database
            ###### COUNTER 0 - MATCHING SIGNALS
            # Create a counter, symmetrical to the database Peaklist
            matching_signals_number <- 0
            matching_signals <- numeric()
            # For each peaklist in the Library
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                # Scroll the peaks in the sample
                for (i in 1:length(peaks_sample@mass)) {
                    # Scroll the peaklist in the library
                    for (j in 1:length(peaks_database_temp@mass)) {
                        # Count the number of closely matching signals with the reference in the sample peaklist
                        if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
                            matching_signals_number <- matching_signals_number + 1
                            matching_signals <- append(matching_signals, peaks_sample@mass[i])
                            break
                        }
                    }
                }
            } else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
                    matching_signals_number <- 0
            }
            # Append this row to the global matrix
            matching_signals_matrix [1,db] <- matching_signals_number
            ###### COUNTER 1 - FIT
            fit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                fit_sample <- matching_signals_number / length(peaks_sample@mass)
            }
            # Append this row to the global matrix
            fit_matrix[1,db] <- fit_sample
            ###### COUNTER 2 - RETRO FIT
            retrofit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                retrofit_sample <- matching_signals_number / length(peaks_database_temp@mass)
            }
            # Append this row to the global matrix
            retrofit_matrix[1,db] <- retrofit_sample
            ###### COUNTER 3
            # Symmetry -> comparison between intensities
            # Compute the correlation matrix with the library
            # Intensity matrix
            intensity_matrix_global <- intensityMatrix(peaks_all, spectra_all)
            # Keep only the matching signals
            #columns_to_keep <- as.character(matching_signals)
            #intensity_matrix_global <- intensity_matrix_global[,columns_to_keep]
            # Weighted correlation between samples (library + test samples) (samples must be as columns and features as test) - With weights
            if (intensity_correction_coefficient != 0 && intensity_correction_coefficient != 1) {
                # Compute the vector of weights
                weights_vector <- c(rep(1, length(database_vector)), rep(intensity_correction_coefficient, nrow(t(intensity_matrix_global))))
                correlation_sample <- wtd.cors(x = t(intensity_matrix_global), weight = weights_vector)
                intensity_correlation_sample <- as.matrix(correlation_matrix [(database_size+1):nrow(correlation_matrix), 1:database_size])
            } else if (intensity_correction_coefficient == 1) {
                t_intensity_matrix_global <- t(intensity_matrix_global)
                correlation_sample <- cor.test(t_intensity_matrix_global[,1], t_intensity_matrix_global[,2], method = "pearson")
                intensity_correlation_sample <- correlation_sample$estimate
                # pvalue
                pvalue <- correlation_sample$p.value
                pvalue_replacement_function <- function(x, number_of_digits) {
                    if (is.na(x)) {
                        x <- "Not available"
                    } else if (x < 0.00001) {
                        x <- "< 0.00001"
                    } else {
                        x <- as.character(round(x, digits = number_of_digits))
                    }
                    return (x)
                }
                pvalue <- pvalue_replacement_function(pvalue, number_of_digits = 6)
                # Append this row to the global matrix
                pvalue_matrix[1,db] <- pvalue
            } else if (intensity_correction_coefficient == 0) {
                intensity_correlation_sample <- 1
            }
            # Extract the absolute values and fix the NAs
            intensity_correlation_sample <- abs(intensity_correlation_sample)
            if (is.na(intensity_correlation_sample)) {
                intensity_correlation_sample <- 0
            }
            # Append this row to the global matrix
            intensity_correlation_matrix[1,db] <- intensity_correlation_sample
            ###### COUNTER 4 - REGRESSION CURVE
            t_intensity_matrix_global <- t(intensity_matrix_global)
            t_intensity_matrix_database <- rbind(as.matrix(t_intensity_matrix_global[,1]))
            t_intensity_matrix_test <- rbind(as.matrix(t_intensity_matrix_global[,2]))
            linear_regression <- lm(t_intensity_matrix_database[,1] ~ t_intensity_matrix_test[,1])
            regression_slope <- linear_regression$coefficients[2]
            regression_intercept <- linear_regression$coefficients[1]
            slope_sample <- round(regression_slope, digits = 3)
            # Append this row to the global matrix
            slope_matrix[1,db] <- slope_sample
        }
        # Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
        return(list(number_of_signals_samples = number_of_signals_samples, number_of_signals_database_matrix = number_of_signals_database_matrix, matching_signals_matrix = matching_signals_matrix, fit_matrix = fit_matrix, retrofit_matrix = retrofit_matrix, intensity_correlation_matrix = intensity_correlation_matrix, pvalue_matrix = pvalue_matrix, slope_matrix = slope_matrix))
    }
    if (allow_parallelization == TRUE) {
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            output_list <- mclapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_correlation(global_list), mc.cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            # Make the CPU cluster for parallelisation
            cls <- makeCluster(cpu_thread_number)
            # Make the cluster use the custom functions and the package functions along with their parameters
            clusterEvalQ(cls, {library(MALDIquant)})
            # Pass the variables to the cluster for running the function
            clusterExport(cl = cls, varlist = c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_correlation", "database_size", "tof_mode", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_percentage_threshold", "low_intensity_threshold_method", "intensity_correction_coefficient", "correlation_pvalue"), envir = environment())
            output_list <- parLapply(cls, global_list, fun = function(global_list) comparison_sample_db_subfunction_correlation(global_list))
            stopCluster(cls)
        }
    } else {
        output_list <- lapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_correlation(global_list))
    }
    ############################ Merge the matrix pieces together
    matching_signals_matrix_all <- NULL
    number_of_signals_database_matrix_all <- NULL
    fit_matrix_all <- NULL
    retrofit_matrix_all <- NULL
    intensity_correlation_matrix_all <- NULL
    pvalue_matrix_all <- NULL
    slope_matrix_all <- NULL
    for (ns in 1:number_of_samples) {
        # Matching signals
        if (is.null(matching_signals_matrix_all)) {
            matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
        } else {
            matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
        }
        # Number of signals database
        if (is.null(number_of_signals_database_matrix_all)) {
            number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
        } else {
            number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
        }
        # Fit
        if (is.null(fit_matrix_all)) {
            fit_matrix_all <- output_list[[ns]]$fit_matrix
        } else {
            fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
        }
        # Retrofit
        if (is.null(retrofit_matrix_all)) {
            retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
        } else {
            retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
        }
        # Intensity correlation
        if (is.null(intensity_correlation_matrix_all)) {
            intensity_correlation_matrix_all <- output_list[[ns]]$intensity_correlation_matrix
        } else {
            intensity_correlation_matrix_all <- rbind(intensity_correlation_matrix_all, output_list[[ns]]$intensity_correlation_matrix)
        }
        # pvalue
        if (is.null(pvalue_matrix_all)) {
            pvalue_matrix_all <- output_list[[ns]]$pvalue_matrix
        } else {
            pvalue_matrix_all <- rbind(pvalue_matrix_all, output_list[[ns]]$pvalue_matrix)
        }
        # Slope
        if (is.null(slope_matrix_all)) {
            slope_matrix_all <- output_list[[ns]]$slope_matrix
        } else {
            slope_matrix_all <- rbind(slope_matrix_all, output_list[[ns]]$slope_matrix)
        }
    }
    ######################################
    # Fix the rownames and colnames
    rownames(intensity_correlation_matrix_all) <- sample_vector
    colnames(intensity_correlation_matrix_all) <- database_vector
    ################### Score calculation
    if (intensity_correction_coefficient != 0) {
    score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_correlation_matrix_all*1000)
    } else {
        score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_correlation_matrix_all*100)
    }
    colnames(score) <- database_vector
    rownames(score) <- sample_vector
    #### Output the classification
    output <- matrix ("", nrow = number_of_samples, ncol = database_size)
    colnames(output) <- database_vector
    rownames(output) <- sample_vector
    if (spectra_path_output == TRUE) {
        output <- cbind(output, spectra_path_vector)
        colnames(output) <- c(database_vector, "Spectrum path")
    }
    if (score_only == TRUE) {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w]>=2) {
                    output[r,w] <- paste("YES","(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]<1.5) {
                    output[r,w] <- paste("NO", "(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]>=1.5 && score[r,w]<2) {
                    output[r,w] <- paste("NI","(", round(score[r,w], digits = 3), ")")
                }
            }
        }
    } else {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w] >= 2) {
                    output[r,w] <- paste("YES","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits = 3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
                if (score[r,w] < 1.5) {
                    output[r,w] <- paste("NO","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits = 3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
                if (score[r,w] >= 1.5 && score[r,w] < 2) {
                    output[r,w] <- paste("NI","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "Corr:", round(intensity_correlation_matrix_all[r,w], digits = 3), ",", "p:", pvalue_matrix_all[r,w], "sl:", slope_matrix_all[r,w], ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
            }
        }
    }
return(output)
}





################################################################################





######################################### SPECTRAL TYPER SCORE: SIGNAL INTENSITY
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity comparison.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_signal_intensity <- function(spectra_database, spectra_test, peaks_database, peaks_test, class_list_library = NULL, comparison = c("intensity percentage", "standard deviation"), peaks_filtering_percentage_threshold = 5, low_intensity_percentage_threshold = 0, low_intensity_threshold_method = "element-wise", tolerance_ppm = 2000, intensity_tolerance_percent_threshold = 50, spectra_format = "brukerflex", spectra_path_output = TRUE, score_only = TRUE, number_of_st_dev = 1, allow_parallelization = FALSE) {
    # Load the required libraries
    install_and_load_required_packages("MALDIquant","parallel")
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ### Folder lists
    #database_folder_list <- dir(filepath_database, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    #test_folder_list <- dir(filepath_test, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    # Sample and Library size
    if (isMassPeaksList(peaks_test)) {
        number_of_samples <- length(peaks_test)
    } else if (isMassPeaks(peaks_test)) {
        number_of_samples <- 1
    }
    if (isMassPeaksList(peaks_database)) {
        database_size <- length(peaks_database)
    } else if (isMassPeaks(peaks_database)) {
        database_size <- 1
    }
    #### Replace the sample name, both in the library and in the test set
    peaks_test <- replace_sample_name(peaks_test, spectra_format = spectra_format)
    peaks_database <- replace_class_name(peaks_database, class_list = class_list_library, class_in_file_name = TRUE, spectra_format = spectra_format)
    ####### Create the sample vector
    sample_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    ####### Create the library vector
    database_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    # Generate the path vector
    spectra_path_vector <- character()
    for (sp in 1:length(spectra_test)) {
        spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
    }
    ############################################################ SCORE (FRI)
    number_of_signals_database <- numeric(length = database_size)
    for (d in 1:database_size) {
        number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
    }
    ################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
    # Create a list to be used for lapply
    global_list <- list()
    for (spl in 1:number_of_samples) {
        # Extract the peaklist and the spectrum
        peaks_database_temp <- peaks_database
        spectra_database_temp <- spectra_database
        peaks_sample <- peaks_test[[spl]]
        spectrum_sample <- spectra_test[[spl]]
        # Generate the entry of the global list
        global_list_entry <- append(peaks_sample, peaks_database_temp)
        global_list_entry <- append(global_list_entry, spectrum_sample)
        global_list_entry <- append(global_list_entry, spectra_database_temp)
        global_list[[spl]] <- global_list_entry
    }
    ############################################## Define the function for parLapply
    # x = each element of the global list
    comparison_sample_db_subfunction_intensity <- function(x) {
        # Generate the matrix rows for the output
        matching_signals_matrix <- matrix(0, nrow = 1, ncol = database_size)
        number_of_signals_database_matrix <- matrix(0, nrow = 1, ncol = database_size)
        fit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        retrofit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        intensity_matching_matrix <- matrix(0, ncol = database_size, nrow = 1)
        ###### Compare with all the elements in the library
        ### For each entry in the database...
        for (db in 1:database_size) {
            # Extract the peaklist and the spectrum
            peaks_sample <- x[[1]]
            spectrum_sample <- x[[(database_size+2)]]
            peaks_database_temp_all <- x[2:(database_size+1)]
            spectra_database_temp_all <- x[(database_size+3):length(x)]
            peaks_database_temp <- peaks_database_temp_all[[db]]
            spectra_database_temp <- spectra_database_temp_all[[db]]
            ####### Peak alignment
            # Merge the peaklists
            peaks_all <- append(peaks_database_temp, peaks_sample)
            spectra_all <- append(spectra_database_temp, spectrum_sample)
            # Align the peaks
            peaks_all <- align_and_filter_peaks(peaks_all, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peaks_filtering_percentage_threshold, low_intensity_peak_removal_threshold_percent = low_intensity_percentage_threshold, low_intensity_peak_removal_threshold_method = low_intensity_threshold_method)
            # Restore the lists
            peaks_database_temp <- peaks_all[[1]]
            peaks_sample <- peaks_all[[2]]
            #################### Number of signals
            number_of_signals_samples <- length(peaks_sample@mass)
            number_of_signals_database <- length(peaks_database_temp@mass)
            number_of_signals_database_matrix[1,db] <- number_of_signals_database
            ###### COUNTER 0 - MATCHING SIGNALS
            # Create a counter, symmetrical to the database Peaklist
            matching_signals_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                # Scroll the peaks in the sample
                for (i in 1:length(peaks_sample@mass)) {
                    # Scroll the peaklist in the library
                    for (j in 1:length(peaks_database_temp@mass)) {
                        # Count the number of closely matching signals with the reference in the sample peaklist
                        if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
                            matching_signals_sample <- matching_signals_sample + 1
                            break
                        }
                    }
                }
            } else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
                matching_signals_sample <- 0
            }
            # Append this row to the global matrix
            matching_signals_matrix[1,db] <- matching_signals_sample
            ###### COUNTER 1 - FIT
            fit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                fit_sample <- matching_signals_sample / length(peaks_sample@mass)
            }
            # Append this row to the global matrix
            fit_matrix[1,db] <- fit_sample
            ###### COUNTER 2 - RETRO FIT
            retrofit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                retrofit_sample <- matching_signals_sample / length(peaks_database_temp@mass)
            }
            # Append this row to the global matrix
            retrofit_matrix[1,db] <- retrofit_sample
            ###### COUNTER 3 (INTENSITY MATCHING)
            if ("intensity percentage" %in% comparison && !("standard deviation" %in% comparison)) {
                # Create a counter, symmetrical to the database Peaklist
                intensity_matching_sample <- 0
                if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                    # Scroll the peaks in the sample
                    for (i in 1:length(peaks_sample@mass)) {
                        # Scroll the peaklist in the library
                        for (j in 1:length(peaks_database_temp@mass)) {
                            # Count the number of closely matching signals with the reference in the sample peaklist
                            if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
                                # Evaluate the difference in intensity
                                if ((abs(peaks_sample@intensity[i] - peaks_database_temp@intensity[j])*100/peaks_database_temp@intensity[j]) <= intensity_tolerance_percent_threshold) {
                                    intensity_matching_sample <- intensity_matching_sample + 1
                                }
                                break
                            }
                        }
                    }
                } else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
                    intensity_matching_sample <- 0
                }
                if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                    intensity_matching_sample <- intensity_matching_sample / matching_signals_sample
                    if (is.na(intensity_matching_sample)) {
                        intensity_matching_sample <- 0
                    }
                }
                # Append this row to the global matrix
                intensity_matching_matrix[1,db] <- intensity_matching_sample
            } else if (!("intensity percentage" %in% comparison) && "standard deviation" %in% comparison) {
                ########### To be implemented
                intensity_matching_sample <- NULL
            }
        }
        # Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
        return(list(number_of_signals_samples = number_of_signals_samples, number_of_signals_database_matrix = number_of_signals_database_matrix, matching_signals_matrix = matching_signals_matrix, fit_matrix = fit_matrix, retrofit_matrix = retrofit_matrix, intensity_matching_matrix = intensity_matching_matrix))
    }
    if (allow_parallelization == TRUE) {
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            output_list <- mclapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_intensity(global_list), mc.cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            # Make the CPU cluster for parallelisation
            cls <- makeCluster(cpu_thread_number)
            # Make the cluster use the custom functions and the package functions along with their parameters
            clusterEvalQ(cls, {library(MALDIquant)})
            clusterExport(cl = cls, varlist = c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_intensity", "database_size", "tof_mode", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_percentage_threshold", "low_intensity_threshold_method", "intensity_correction_coefficient", "correlation_pvalue", "intensity_tolerance_percent_threshold", "comparison"), envir = environment())
            output_list <- parLapply(cls, global_list, fun = function(global_list) comparison_sample_db_subfunction_intensity(global_list))
            stopCluster(cls)
        }
    } else {
        output_list <- lapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_intensity(global_list))
    }
    ############################ Merge the matrix pieces together
    matching_signals_matrix_all <- NULL
    number_of_signals_database_matrix_all <- NULL
    fit_matrix_all <- NULL
    retrofit_matrix_all <- NULL
    intensity_matching_matrix_all <- NULL
    pvalue_matrix_all <- NULL
    slope_matrix_all <- NULL
    for (ns in 1:number_of_samples) {
        # Matching signals
        if (is.null(matching_signals_matrix_all)) {
            matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
        } else {
            matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
        }
        # Number of signals database
        if (is.null(number_of_signals_database_matrix_all)) {
            number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
        } else {
            number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
        }
        # Fit
        if (is.null(fit_matrix_all)) {
            fit_matrix_all <- output_list[[ns]]$fit_matrix
        } else {
            fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
        }
        # Retrofit
        if (is.null(retrofit_matrix_all)) {
            retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
        } else {
            retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
        }
        # Intensity correlation
        if (is.null(intensity_matching_matrix_all)) {
            intensity_matching_matrix_all <- output_list[[ns]]$intensity_matching_matrix
        } else {
            intensity_matching_matrix_all <- rbind(intensity_matching_matrix_all, output_list[[ns]]$intensity_matching_matrix)
        }
    }
    ######################################
    # Fix the rownames and colnames
    colnames(intensity_matching_matrix_all) <- database_vector
    rownames(intensity_matching_matrix_all) <- sample_vector
    ### Score calculation
    score <- log10(fit_matrix_all*retrofit_matrix_all*intensity_matching_matrix_all*1000)
    colnames(score) <- database_vector
    rownames(score) <- sample_vector
    #### Output the classification
    output <- matrix ("NO", nrow = number_of_samples, ncol = database_size)
    colnames(output) <- database_vector
    rownames(output) <- sample_vector
    if (spectra_path_output == TRUE) {
        output <- cbind(output, spectra_path_vector)
        colnames(output) <- c(database_vector, "Spectrum path")
    }
    if (score_only == TRUE) {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w]>=2) {
                    output[r,w] <- paste("YES","(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]<1.5) {
                    output[r,w] <- paste("NO", "(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]>=1.5 && score[r,w]<2) {
                    output[r,w] <- paste("NI","(", round(score[r,w], digits = 3), ")")
                }
            }
        }
    } else {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w]>=2) {
                    output[r,w] <- paste("YES","(Score:", round(score[r,w], digits = 3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits = 3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
                if (score[r,w]<1.5) {
                    output[r,w] <- paste("NO","(Score:", round(score[r,w], digits = 3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits = 3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
                if (score[r,w]>=1.5 && score[r,w]<2) {
                    output[r,w] <- paste("NI","(Score:", round(score[r,w], digits = 3), "), ", "(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "IntMtch:", round(intensity_matching_matrix_all[r,w], digits = 3), ",", "ns:", matching_signals_matrix_all[r,w], ")")
                }
            }
        }
    }
    return(output)
}





######################################### SPECTRAL TYPER SCORE: SIMILARITY INDEX
# The function calculates the score for the Spectral Typer program, by comparing the test peaklist with the database peaklist, in terms of peak matching and intensity symmetry via the similarity index computation.
# Each sample gets compared with each entry in the database, separately.
# Parallel implemented.
spectral_typer_score_similarity_index <- function(spectra_database, spectra_test, peaks_database, peaks_test, filepath_database, filepath_test, class_list_library = NULL, peaks_filtering_percentage_threshold = 5, low_intensity_percentage_threshold = 0.1, low_intensity_threshold_method = "element-wise", tolerance_ppm = 2000, intensity_correction_coefficient = 1, spectra_format = "brukerflex", spectra_path_output = TRUE, score_only = FALSE, allow_parallelization = FALSE) {
    install_and_load_required_packages(c("MALDIquant", "stats", "parallel"))
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ### Folder lists
    #database_folder_list <- dir(filepath_database, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    #test_folder_list <- dir(filepath_test, ignore.case = TRUE, full.names = FALSE, recursive = FALSE, include.dirs = TRUE)
    # Sample and Library size
    if (isMassPeaksList(peaks_test)) {
        number_of_samples <- length(peaks_test)
    } else if (isMassPeaks(peaks_test)) {
        number_of_samples <- 1
    }
    if (isMassPeaksList(peaks_database)) {
        database_size <- length(peaks_database)
    } else if (isMassPeaks(peaks_database)) {
        database_size <- 1
    }
    #### Replace the sample name, both in the library and in the test set
    peaks_test <- replace_sample_name(peaks_test, spectra_format = spectra_format)
    peaks_database <- replace_class_name(peaks_database, class_list = class_list_library, class_in_file_name = TRUE, spectra_format = spectra_format)
    ####### Create the sample vector
    sample_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:number_of_samples) {
            sample_vector <- append(sample_vector, peaks_test[[s]]@metaData$file[1])
        }
    }
    ####### Create the library vector
    database_vector <- character()
    if (spectra_format == "brukerflex" || spectra_format == "xmass") {
        # If a spectrum is the result of the averaging of several spectra, take only the first name in the file name (the file name is a vector with all the names of the original spectra)
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    if (spectra_format == "imzml" | spectra_format == "imzML") {
        for (s in 1:database_size) {
            database_vector <- append(database_vector, peaks_database[[s]]@metaData$file[1])
        }
    }
    # Generate the path vector
    spectra_path_vector <- character()
    for (sp in 1:length(spectra_test)) {
        spectra_path_vector <- append(spectra_path_vector, spectra_test[[sp]]@metaData$file[1])
    }
    ############################################################ SCORE (FRI)
    number_of_signals_database <- numeric(length = database_size)
    for (d in 1:database_size) {
        number_of_signals_database[d] <- length(peaks_database[[d]]@mass)
    }
    ################### Each sample gets compared with the database (create a copy of the original database each time, otherwise it gets modified when processed together with the sample)
    # Create a list to be used for lapply
    global_list <- list()
    for (spl in 1:number_of_samples) {
        # Extract the peaklist and the spectrum
        peaks_database_temp <- peaks_database
        spectra_database_temp <- spectra_database
        peaks_sample <- peaks_test[[spl]]
        spectrum_sample <- spectra_test[[spl]]
        # Generate the entry of the global list
        global_list_entry <- append(peaks_sample, peaks_database_temp)
        global_list_entry <- append(global_list_entry, spectrum_sample)
        global_list_entry <- append(global_list_entry, spectra_database_temp)
        global_list[[spl]] <- global_list_entry
    }
    ############################################## Define the function for parLapply
    # x = each element of the global list
    comparison_sample_db_subfunction_similarity_index <- function(x) {
        # Generate the matrix rows for the output
        matching_signals_matrix <- matrix(0, nrow = 1, ncol = database_size)
        number_of_signals_database_matrix <- matrix(0, nrow = 1, ncol = database_size)
        fit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        retrofit_matrix <- matrix(0, ncol = database_size, nrow = 1)
        similarity_index_matrix <- matrix(0, ncol = database_size, nrow = 1)
        ###### Compare with all the elements in the library
        ### For each entry in the library...
        for (db in 1:database_size) {
            # Extract the peaklist and the spectrum
            peaks_sample <- x[[1]]
            spectrum_sample <- x[[(database_size+2)]]
            peaks_database_temp_all <- x[2:(database_size+1)]
            spectra_database_temp_all <- x[(database_size+3):length(x)]
            peaks_database_temp <- peaks_database_temp_all[[db]]
            spectra_database_temp <- spectra_database_temp_all[[db]]
            ####### Peak alignment
            # Merge the peaklists
            peaks_all <- append(peaks_database_temp, peaks_sample)
            spectra_all <- append(spectra_database_temp, spectrum_sample)
            # Align the peaks
            peaks_all <- align_and_filter_peaks(peaks_all, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peaks_filtering_percentage_threshold, low_intensity_peak_removal_threshold_percent = low_intensity_percentage_threshold, low_intensity_peak_removal_threshold_method = low_intensity_threshold_method)
            # Restore the lists
            peaks_database_temp <- peaks_all[[1]]
            peaks_sample <- peaks_all[[2]]
            #################### Number of signals
            number_of_signals_samples <- length(peaks_sample@mass)
            number_of_signals_database <- length(peaks_database_temp@mass)
            number_of_signals_database_matrix[1,db] <- number_of_signals_database
            ###### COUNTER 0 - MATCHING SIGNALS
            # Create a counter, symmetrical to the database Peaklist
            matching_signals_number <- 0
            matching_signals <- numeric()
            # For each peaklist in the Library
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                # Scroll the peaks in the sample
                for (i in 1:length(peaks_sample@mass)) {
                    # Scroll the peaklist in the library
                    for (j in 1:length(peaks_database_temp@mass)) {
                        # Count the number of closely matching signals with the reference in the sample peaklist
                        if (peaks_sample@mass[i] == peaks_database_temp@mass[j]) {
                            matching_signals_number <- matching_signals_number + 1
                            matching_signals <- append(matching_signals, peaks_sample@mass[i])
                            break
                        }
                    }
                }
            } else if (length(peaks_sample@mass) == 0 || length(peaks_database_temp@mass) == 0) {
                    matching_signals_number <- 0
            }
            # Append this row to the global matrix
            matching_signals_matrix [1,db] <- matching_signals_number
            ###### COUNTER 1 - FIT
            fit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                fit_sample <- matching_signals_number / length(peaks_sample@mass)
            }
            # Append this row to the global matrix
            fit_matrix[1,db] <- fit_sample
            ###### COUNTER 2 - RETRO FIT
            retrofit_sample <- 0
            if (length(peaks_sample@mass) > 0 && length(peaks_database_temp@mass) > 0) {
                retrofit_sample <- matching_signals_number / length(peaks_database_temp@mass)
            }
            # Append this row to the global matrix
            retrofit_matrix[1,db] <- retrofit_sample
            ###### COUNTER 3
            # Symmetry -> comparison between intensities
            # Compute the similaroty index with the library
            similarity_index_matrix_global <- intensityMatrix(peaks_all, spectra_all)
            # Similarity index (E Id * Ix / sqrt(E Id^2 * E Ix^2)) = A / sqrt (B * E)
            A <- 0
            for (z in 1:ncol(similarity_index_matrix_global)) {
                A <- A + (similarity_index_matrix_global[1,z]*similarity_index_matrix_global[2,z])
            }
            B <- 0
            for (z in 1:ncol(similarity_index_matrix_global)) {
                B <- B + (similarity_index_matrix_global[1,z])^2
            }
            E <- 0
            for (z in 1:ncol(similarity_index_matrix_global)) {
                E <- E + (similarity_index_matrix_global[2,z])^2
            }
            similarity_index_sample <- A / sqrt (B * E)
            # Append this row to the global matrix
            similarity_index_matrix[1,db] <- similarity_index_sample
        }
        # Return a list, each element of which is a matrix row. Finally, all the matrix rows will be rbind together.
        return(list(number_of_signals_samples = number_of_signals_samples, number_of_signals_database_matrix = number_of_signals_database_matrix, matching_signals_matrix = matching_signals_matrix, fit_matrix = fit_matrix, retrofit_matrix = retrofit_matrix, similarity_index_matrix = similarity_index_matrix))
    }
    if (allow_parallelization == TRUE) {
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            output_list <- mclapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_similarity_index(global_list), mc.cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            # Make the CPU cluster for parallelisation
            cls <- makeCluster(cpu_thread_number)
            # Make the cluster use the custom functions and the package functions along with their parameters
            clusterEvalQ(cls, {library(MALDIquant)})
            clusterExport(cl = cls, varlist = c("align_and_filter_peaks", "install_and_load_required_packages", "comparison_sample_db_subfunction_correlation", "database_size", "tof_mode", "peaks_filtering_percentage_threshold", "remove_low_intensity_peaks", "low_intensity_percentage_threshold", "low_intensity_threshold_method"), envir = environment())
            output_list <- parLapply(cls, global_list, fun = function(global_list) comparison_sample_db_subfunction_similarity_index(global_list))
            stopCluster(cls)
        }
    } else {
        output_list <- lapply(global_list, FUN = function(global_list) comparison_sample_db_subfunction_similarity_index(global_list))
    }
    ############################ Merge the matrix pieces together
    matching_signals_matrix_all <- NULL
    number_of_signals_database_matrix_all <- NULL
    fit_matrix_all <- NULL
    retrofit_matrix_all <- NULL
    similarity_index_matrix_all <- NULL
    for (ns in 1:number_of_samples) {
        # Matching signals
        if (is.null(matching_signals_matrix_all)) {
            matching_signals_matrix_all <- output_list[[ns]]$matching_signals_matrix
        } else {
            matching_signals_matrix_all <- rbind(matching_signals_matrix_all, output_list[[ns]]$matching_signals_matrix)
        }
        # Number of signals database
        if (is.null(number_of_signals_database_matrix_all)) {
            number_of_signals_database_matrix_all <- output_list[[ns]]$number_of_signals_database_matrix
        } else {
            number_of_signals_database_matrix_all <- rbind(number_of_signals_database_matrix_all, output_list[[ns]]$number_of_signals_database_matrix)
        }
        # Fit
        if (is.null(fit_matrix_all)) {
            fit_matrix_all <- output_list[[ns]]$fit_matrix
        } else {
            fit_matrix_all <- rbind(fit_matrix_all, output_list[[ns]]$fit_matrix)
        }
        # Retrofit
        if (is.null(retrofit_matrix_all)) {
            retrofit_matrix_all <- output_list[[ns]]$retrofit_matrix
        } else {
            retrofit_matrix_all <- rbind(retrofit_matrix_all, output_list[[ns]]$retrofit_matrix)
        }
        # Similarity index
        if (is.null(similarity_index_matrix_all)) {
            similarity_index_matrix_all <- output_list[[ns]]$similarity_index_matrix
        } else {
            similarity_index_matrix_all <- rbind(similarity_index_matrix_all, output_list[[ns]]$similarity_index_matrix)
        }
    }
    ######################################
    # Fix the rownames and colnames
    rownames(similarity_index_matrix_all) <- sample_vector
    colnames(similarity_index_matrix_all) <- database_vector
    ################### Score calculation
    score <- log10(fit_matrix_all*retrofit_matrix_all*similarity_index_matrix_all*1000)
    colnames(score) <- database_vector
    rownames(score) <- sample_vector
    #### Output the classification
    output <- matrix ("", nrow = number_of_samples, ncol = database_size)
    colnames(output) <- database_vector
    rownames(output) <- sample_vector
    if (spectra_path_output == TRUE) {
        output <- cbind(output, spectra_path_vector)
        colnames(output) <- c(database_vector, "Spectrum path")
    }
    if (score_only == TRUE) {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w]>=2) {
                    output[r,w] <- paste("YES","(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]<1.5) {
                    output[r,w] <- paste("NO", "(", round(score[r,w], digits = 3), ")")
                }
                if (score[r,w]>=1.5 && score[r,w]<2) {
                    output[r,w] <- paste("NI","(", round(score[r,w], digits = 3), ")")
                }
            }
        }
    } else {
        for (r in 1:number_of_samples) {
            for (w in 1:database_size) {
                if (score[r,w] >= 2) {
                    output[r,w] <- paste("YES","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits = 3), ")")
                }
                if (score[r,w] < 1.5) {
                    output[r,w] <- paste("NO","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits = 3), ")")
                }
                if (score[r,w] >= 1.5 && score[r,w] < 2) {
                    output[r,w] <- paste("NI","(Score:", round(score[r,w], digits = 3), "), ","(F:", matching_signals_matrix_all[r,w], "/", output_list[[r]]$number_of_signals_samples, " = ", round(fit_matrix_all[r,w], digits = 3), ",", "RF:", matching_signals_matrix_all[r,w], "/", number_of_signals_database_matrix_all[r,w], " = ", round(retrofit_matrix_all[r,w], digits = 3), ",", "SI:", round(similarity_index_matrix_all[r,w], digits = 3), ")")
                }
            }
        }
    }
    return(output)
}





################################################################################





#################################################### ADJACENCY MATRIX GENERATION
# The function generates an adjacency matrix from a peak list matrix, in order to generate a graph. The matrix is computed first by generating a correlation matrix and then replacing the correlation coefficients with 0 or 1 according to a threshold value.
generate_adjacency_matrix <- function(peaklist_matrix, correlation_method = "spearman", correlation_threshold = 0.90, pvalue_threshold = 0.05) {
    ##### Install the required packages
    install_and_load_required_packages("stats")
    ##### Transpose the peaklist matrix to compute the correlation between observations and not features
    peaklist_matrix_t <- t(peaklist_matrix)
    ##### Generate the correlation matrix
    correlation_matrix <- cor(peaklist_matrix_t, method = correlation_method)
    ### If the p-value is not considered...
    if (pvalue_threshold == 0 || is.null(pvalue_threshold)) {
        ##### Generate the function to apply to the matrix (x = matrix entry)
        matrix_replacement_subfunction <- function(x, threshold) {
            if (abs(x) >= threshold) {
                x <- 1
            } else {
                x <- 0
            }
        }
        ##### Generate the final adjacency matrix
        adjacency_matrix <- apply(correlation_matrix, MARGIN = c(1,2), FUN = function(x) matrix_replacement_subfunction(x, threshold = correlation_threshold))
        ##### Fix the column and row names
        rownames(adjacency_matrix) <- rownames(peaklist_matrix)
        colnames(adjacency_matrix) <- rownames(peaklist_matrix)
    } else {
        ##### Install the required packages
        install_and_load_required_packages("psych")
        ##### Generate the correlation matrix
        # Generate unique row and column names for the corr.test function
        rownames(peaklist_matrix_t) <- seq(1:nrow(peaklist_matrix_t))
        colnames(peaklist_matrix_t) <- seq(1:ncol(peaklist_matrix_t))
        correlation_matrix_pvalue <- corr.test(peaklist_matrix_t, method = correlation_method, adjust = "holm", ci = FALSE, alpha = pvalue_threshold)$p
        ##### Generate a global matrix (each entry displays "coefficient, pvalue")
        global_matrix <- matrix("", nrow = nrow(correlation_matrix_pvalue), ncol = ncol(correlation_matrix_pvalue))
        rownames(global_matrix) <- seq(1:nrow(global_matrix))
        colnames(global_matrix) <- seq(1:ncol(global_matrix))
        for (rw in 1:nrow(correlation_matrix_pvalue)) {
            for (cl in 1:ncol(correlation_matrix_pvalue)) {
                global_matrix[rw,cl] <- paste(correlation_matrix[rw,cl], correlation_matrix_pvalue[rw,cl], sep = ",")
            }
        }
        ##### Generate the function to apply to the matrix (x = matrix entry)
        matrix_replacement_subfunction <- function(x, coeff_threshold, p_threshold) {
            # Split the entry
            splitted_x <- as.numeric(unlist(strsplit(x, ",")))
            if (abs(splitted_x[1]) >= coeff_threshold && abs(splitted_x[2]) <= p_threshold) {
                x <- 1
            } else {
                x <- 0
            }
        }
        ##### Generate the final adjacency matrix
        adjacency_matrix <- apply(global_matrix, MARGIN = c(1,2), FUN = function(x) matrix_replacement_subfunction(x, coeff_threshold = correlation_threshold, p_threshold = pvalue_threshold))
        ##### Fix the column and row names
        rownames(adjacency_matrix) <- rownames(peaklist_matrix)
        colnames(adjacency_matrix) <- rownames(peaklist_matrix)
    }
    ##### Return the matrix
    return(adjacency_matrix)
}





################################################################################





############################################# GENERATE A CUSTOM INTENSITY MATRIX
# The function takes a list of spectra and a vector of custom features to be included in the generation of the final peaklist intensity matrix. The functions takes the spectra, preprocesses the spectra according to the specified parameters, performs the peak picking and outputs the intensity matrix only for the peaks specified as input (not all of those custom peaks if they are outside of the spectral mass range).
# If the range provided is too large, the function will return a NULL value, since some custom features cannot be found.
# This function is suited for aligning the spectral features (of an unknown dataset) with the model features.
generate_custom_intensity_matrix <- function(spectra, custom_feature_vector = NULL, tof_mode = "linear", preprocessing_parameters = list(mass_range = c(800,3000), transformation_algorithm = NULL, smoothing_algorithm = NULL, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 100, normalization_algorithm = "TIC", normalization_mass_range = NULL, preprocess_spectra_in_packages_of = 0, spectral_alignment_method = NULL), peak_picking_algorithm = "SuperSmoother", peak_picking_SNR = 5, peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", allow_parallelization = FALSE, deisotope_peaklist = FALSE) {
    ### Install the required packages
    install_and_load_required_packages("MALDIquant")
    # Rename the trim function
    trim_spectra <- get(x = "trim", pos = "package:MALDIquant")
    ### Define the tolerance
    if (tof_mode == "linear" || tof_mode == "L") {
        tolerance_ppm <- 2000
    } else if (tof_mode == "reflectron" || tof_mode == "reflector" || tof_mode == "R") {
        tolerance_ppm <- 200
    }
    ### Preprocessing
    if (!is.null(preprocessing_parameters) && is.list(preprocessing_parameters) && length(preprocessing_parameters) > 0) {
        spectra <- preprocess_spectra(spectra, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, allow_parallelization = allow_parallelization)
    }
    ### Peak picking and alignment (with the custom features)
    peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, SNR = peak_picking_SNR, allow_parallelization = allow_parallelization, deisotope_peaklist = deisotope_peaklist)
    ### Run the alignment only if the vector of custom features is not null
    if (!is.null(custom_feature_vector)) {
        # Check if there are X at the beginning of the feature numbers before converting into numbers
        if (unlist(strsplit(as.character(custom_feature_vector[1]),""))[1] == "X") {
            # Remove the X
            for (f in 1:length(custom_feature_vector)) {
                name_splitted <- unlist(strsplit(custom_feature_vector[f],""))
                feature_def <- name_splitted [2]
                for (i in 3:length(name_splitted)) {
                    feature_def <- paste(feature_def, name_splitted[i], sep = "")
                }
                custom_feature_vector[f] <- feature_def
            }
        }
        # Convert the custom feature vector in numeric
        custom_feature_vector <- sort(as.numeric(unique(custom_feature_vector)))
        # Retrieve the maximum and minimum data points in the spectra
        spectra_peaks <- numeric()
        if (isMassSpectrumList(spectra)) {
            for (sp in 1:length(spectra)) {
                spectra_peaks <- append(spectra_peaks, c(spectra[[sp]]@mass[1], spectra[[sp]]@mass[length(spectra[[sp]]@mass)]))
            }
        } else {
            spectra_peaks <- append(spectra_peaks, c(spectra@mass[1], spectra@mass[length(spectra@mass)]))
        }
        spectra_peaks <- sort(as.numeric(unique(spectra_peaks)))
        ### Check the compatibility between the spectra and the provided mass list
        if (spectra_peaks[1] <= custom_feature_vector[1] && spectra_peaks[length(spectra_peaks)] >= custom_feature_vector[length(custom_feature_vector)]) {
            feature_compatibility <- TRUE
        } else {
            feature_compatibility <- FALSE
        }
        ### If there is feature compatibility...
        if (feature_compatibility == TRUE) {
            ### Function of peak alignment for lapply
            peak_alignment_subfunction <- function(peaks, reference_masses, tolerance_ppm) {
                # Scroll the peaks
                for (l in 1:length(reference_masses)) {
                    # Identify the peak to be aligned to the reference
                    peak_id <- which(abs((as.numeric(peaks@mass) - as.numeric(reference_masses[l])) / as.numeric(reference_masses[l]) * 10^6) <= tolerance_ppm)
                    # If there is one peak, align it with the reference
                    if (length(peak_id == 1)) {
                        peaks@mass[peak_id] <- as.numeric(reference_masses[l])
                    } else if (length(peak_id) > 1) {
                        # If there is more than one peak, align only the first with the reference
                        peaks@mass[peak_id[1]] <- as.numeric(reference_masses[l])
                    }
                }
                # Return
                return(peaks)
            }
            # Fix the peak values (with LAPPLY)
            if (isMassPeaksList(peaks)) {
                if (allow_parallelization == TRUE) {
                    # Detect the number of cores
                    cpu_thread_number <- detectCores(logical = TRUE)
                    cpu_thread_number <- cpu_thread_number / 2
                    if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                        peaks <- mclapply(peaks, FUN = function(peaks) peak_alignment_subfunction(peaks = peaks, reference_masses = custom_feature_vector, tolerance_ppm = tolerance_ppm), mc.cores = cpu_thread_number)
                    } else if (Sys.info()[1] == "Windows") {
                        # Make the CPU cluster for parallelisation
                        cl <- makeCluster(cpu_thread_number)
                        # Make the cluster use the custom functions and the package functions along with their parameters
                        clusterEvalQ(cl, {library(MALDIquant)})
                        # Pass the variables to the cluster for running the function
                        clusterExport(cl = cl, varlist = c("peaks", "custom_feature_vector", "tolerance_ppm", "peak_alignment_subfunction"), envir = environment())
                        # Apply the multicore function
                        peaks <- parLapply(cl, peaks, fun = function(peaks) peak_alignment_subfunction(peaks = peaks, reference_masses = custom_feature_vector, tolerance_ppm = tolerance_ppm))
                        stopCluster(cl)
                    } else {
                        peaks <- lapply(peaks, FUN = function(peaks) peak_alignment_subfunction(peaks = peaks, reference_masses = custom_feature_vector, tolerance_ppm = tolerance_ppm))
                    }
                } else {
                    peaks <- lapply(peaks, FUN = function(peaks) peak_alignment_subfunction(peaks = peaks, reference_masses = custom_feature_vector, tolerance_ppm = tolerance_ppm))
                }
            } else if (isMassPeaks(peaks)) {
                peaks <- peak_alignment_subfunction(peaks = peaks, reference_masses = custom_feature_vector, tolerance_ppm = tolerance_ppm)
            }
            # Generate a fake spectrum and a fake peaklist with the custom features (to exploit the intensityMatrix function afterwards)
            if (isMassSpectrumList(spectra)) {
                fake_spectrum <- createMassSpectrum(mass = spectra[[1]]@mass, intensity = spectra[[1]]@intensity, metaData = list(name = "Fake spectrum"))
            } else if (isMassSpectrum(spectra)) {
                fake_spectrum <- createMassSpectrum(mass = spectra@mass, intensity = spectra@intensity, metaData = list(name = "Fake spectrum"))
            }
            fake_peaks <- createMassPeaks(mass = as.numeric(custom_feature_vector), intensity = rep(1, length(custom_feature_vector)), snr = rep(peak_picking_SNR, length(custom_feature_vector)), metaData = list(name = "Fake peaklist"))
            # Append the fake spectrum and the fake peaklist to the original lists (the fake will be the first element of the list)
            spectra_all <- append(fake_spectrum, spectra)
            peaks_all <- append(fake_peaks, peaks)
            # Generate the intensity matrix (with the custom features, which will be a little misaligned due to the alignment with the other peaks)
            intensity_matrix_all <- intensityMatrix(peaks_all, spectra_all)
            # Remove the first row (corresponding to the fake spectrum)
            intensity_matrix_all <- intensity_matrix_all[2:nrow(intensity_matrix_all), ]
            # Keep only the columns of interest
            if (is.matrix(intensity_matrix_all)) {
                final_peaklist_matrix <- intensity_matrix_all[, as.character(custom_feature_vector)]
            } else if (is.vector(intensity_matrix_all)) {
                final_peaklist_matrix <- as.matrix(rbind(intensity_matrix_all[names(intensity_matrix_all) %in% as.character(custom_feature_vector)]))
            }
            # Rownames
            if (is.null(names(spectra))) {
                rownames(final_peaklist_matrix) <- seq(1, nrow(final_peaklist_matrix), by = 1)
            } else {
                rownames(final_peaklist_matrix) <- names(spectra)
            }
            ### Return the final matrix with the custom features
            return(final_peaklist_matrix)
        } else {
            ### Return a NULL value if features are not compatible
            return(NULL)
        }
    } else {
        ### Return the simple peaklist matrix if no custom vector is provided
        peaks <- align_and_filter_peaks(peaks, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, reference_peaklist = NULL, spectra = spectra, alignment_iterations = 5, allow_parallelization = allow_parallelization)
        final_peaklist_matrix <- intensityMatrix(peaks, spectra)
        # Rownames
        if (is.null(names(spectra))) {
            rownames(final_peaklist_matrix) <- seq(1, nrow(final_peaklist_matrix), by = 1)
        } else {
            rownames(final_peaklist_matrix) <- names(spectra)
        }
        return(final_peaklist_matrix)
    }
}





################################################################################





##################################################### REARRANGE THE SPECTRAL DATASET
# The function takes a list of spectra and rearranges it in a certain way (defined by the user) on spatial, random or hierarchical basis. It returns a list of spectra, appropriately sorted.
# If space is selected, the spectral list is rearranged along the widest dimension: for example, if X is the widest dimension, the spectra are rearranged as X1,Y1 X2,Y1 X3,Y1 and so on... If the spectral list does not contain spatial coordinates (it is not an imaging dataset), nothing is performedin terms of space.
rearrange_spectral_dataset <- function(spectra, rearranging_method = c("space","random"), seed = NULL) {
    ##### Install and load the required packages
    install_and_load_required_packages("MALDIquant")
    ########## Do everything only if it is a list of spectra
    if (isMassSpectrumList(spectra)) {
        #################### SPACE
        if (rearranging_method == "space") {
            ########## Do everything only if there are spatial coordinates
            if (!is.null(spectra[[1]]@metaData$imaging$pos)) {
                ##### Generate the matrix of spectral coordinates
                spectral_coordinates <- matrix(0, ncol = 2, nrow = length(spectra))
                ### Fill in the matrix
                for (s in 1:length(spectra)) {
                    spectral_coordinates[s,] <- spectra[[s]]@metaData$imaging$pos
                }
                ### Fix the matrix column and row names
                rownames(spectral_coordinates) <- seq(1:nrow(spectral_coordinates))
                colnames(spectral_coordinates) <- c("x", "y")
                ### Convert it into a dataframe for sorting
                spectral_coordinates <- as.data.frame(spectral_coordinates)
                ### Sort along the greater dimension
                ## Define the rectangle containing the tissue section
                rectangle_coordinates <- matrix(0, nrow = 2, ncol = 2)
                rectangle_coordinates[1,1] <- min(spectral_coordinates[,1])
                rectangle_coordinates[2,1] <- max(spectral_coordinates[,1])
                rectangle_coordinates[1,2] <- min(spectral_coordinates[,2])
                rectangle_coordinates[2,2] <- max(spectral_coordinates[,2])
                rownames(rectangle_coordinates) <- c("min", "max")
                colnames(rectangle_coordinates) <- c("x", "y")
                ## Determine the WIDEST coordinate
                if (rectangle_coordinates[2,1] >= rectangle_coordinates[2,2]) {
                    widest_coordinate <- "x"
                } else if (rectangle_coordinates[2,2] > rectangle_coordinates[2,1]) {
                    widest_coordinate <- "y"
                }
                ## X coordinate is the widest
                if (widest_coordinate == "x") {
                    # Sort according the X
                    spectral_coordinates_sorted <- spectral_coordinates[order(spectral_coordinates[,2], spectral_coordinates[,1], decreasing = FALSE), ]
                    # Extract the row names, which are the ID numbers of the spectral list
                    spectra_ID <- as.integer(rownames(spectral_coordinates_sorted))
                    # Rearrange the spectral list
                    spectra_riarranged <- spectra[spectra_ID]
                } else if (widest_coordinate == "y") {
                    ### Y coordinate is the widest
                    # Sort according the Y
                    spectral_coordinates_sorted <- spectral_coordinates[order(spectral_coordinates[,1], spectral_coordinates[,2], decreasing = FALSE), ]
                    # Extract the row names, which are the ID numbers of the spectral list
                    rearranged_spectra_IDs <- as.integer(rownames(spectral_coordinates_sorted))
                    # Rearrange the spectral list
                    spectra_riarranged <- spectra[rearranged_spectra_IDs]
                }
            } else {
                spectra_riarranged <- spectra
            }
        } else if (rearranging_method == "random") {
            #################### RANDOM
            ##### Rearrange the spectra randomly in N parts
            ##### Install and load the required packages
            install_and_load_required_packages("caret")
            # Set the seed (make randomness reproducible)
            if (!is.null(seed)) {
                set.seed(seed)
            }
            # Generate a random list of numbers for random rearrangement
            if (!is.null(seed)) {
                set.seed(seed)
            }
            rearranged_spectra_IDs <- sample(1:length(spectra), size = length(spectra))
            # Rearrange the spectral list
            spectra_riarranged <- spectra[rearranged_spectra_IDs]
        } else if (rearranging_method == "hca") {
            #################### HCA
            NULL
        }
        #################### Output
        return(spectra_riarranged)
    } else if (isMassSpectrum(spectra)) {
        ########## Single spectrum
        return(spectra)
    }
}





################################################################################





##################################################### PARTITION THE SPECTRAL DATASET IN SUBSETS
# The function takes a list of spectra and partitions it in a certain number of subsets (defined by the user) on spatial, random or hierarchical basis. It returns a list of sub-lists of spectra.
partition_spectral_dataset <- function(spectra, partitioning_method = c("space","random", "hca"), number_of_partitions = 3, seed = NULL, tof_mode = "reflectron") {
    ##### Install and load the required packages
    install_and_load_required_packages("MALDIquant")
    ########## Do everything only if it is a list of spectra
    if (isMassSpectrumList(spectra)) {
        #################### SPACE
        if (partitioning_method == "space") {
            ########## Do everything only if there are spatial coordinates
            if (!is.null(spectra[[1]]@metaData$imaging$pos)) {
                ##### Generate the matrix of spectral coordinates
                spectral_coordinates <- matrix(0, ncol = 2, nrow = length(spectra))
                # Fill in the matrix
                for (s in 1:length(spectra)) {
                    spectral_coordinates[s,] <- spectra[[s]]@metaData$imaging$pos
                }
                # Fix the matrix column and row names
                rownames(spectral_coordinates) <- seq(1:nrow(spectral_coordinates))
                colnames(spectral_coordinates) <- c("x", "y")
                # Define the rectangle containing the tissue section
                rectangle_coordinates <- matrix(0, nrow = 2, ncol = 2)
                rectangle_coordinates[1,1] <- min(spectral_coordinates[,1])
                rectangle_coordinates[2,1] <- max(spectral_coordinates[,1])
                rectangle_coordinates[1,2] <- min(spectral_coordinates[,2])
                rectangle_coordinates[2,2] <- max(spectral_coordinates[,2])
                rownames(rectangle_coordinates) <- c("min", "max")
                colnames(rectangle_coordinates) <- c("x", "y")
                ##### Split the dataset based upon the the WIDEST coordinate
                if (abs(rectangle_coordinates[2,1] - rectangle_coordinates[1,1]) >= abs(rectangle_coordinates[2,2] - rectangle_coordinates[1,2])) {
                    widest_coordinate <- "x"
                } else {
                    widest_coordinate <- "y"
                }
                ### X coordinate is the widest
                if (widest_coordinate == "x") {
                    # Define the intervals
                    interval_width <- ceiling(abs(rectangle_coordinates[2,1] - rectangle_coordinates[1,1])/number_of_partitions)
                    # Define the separating points
                    interval_points <- seq.int(from = rectangle_coordinates[1,1], to = rectangle_coordinates[2,1], by = interval_width)
                    # Add the maximum coordinate to the interval points (if not already present)
                    if (interval_points[length(interval_points)] != rectangle_coordinates[2,1]) {
                        interval_points <- c(interval_points, rectangle_coordinates[2,1])
                    }
                    ## Generate the final list of partitioned spectra
                    # Initialize the output
                    spectra_partitioned <- list()
                    # Populate the final list
                    for (p in 1:(length(interval_points) - 1)) {
                        # Generate the sublist of the spectra for that interval
                        spectra_interval <- list()
                        # Scroll the spectra...
                        for (s in 1:length(spectra)) {
                            # Populate the list...
                            if (spectra[[s]]@metaData$imaging$pos[1] >= interval_points[p] && spectra[[s]]@metaData$imaging$pos[1] < interval_points[p + 1]) {
                                spectra_interval <- append(spectra_interval, spectra[[s]])
                            }
                        }
                        # Add the spectra with the max coordinate
                        if (p == (length(interval_points) - 1)) {
                            # Scroll the spectra...
                            for (s in 1:length(spectra)) {
                                # Populate the list...
                                if (spectra[[s]]@metaData$imaging$pos[1] == interval_points[p + 1]) {
                                    spectra_interval <- append(spectra_interval, spectra[[s]])
                                }
                            }
                        }
                        spectra_partitioned[[p]] <- spectra_interval
                    }
                } else if (widest_coordinate == "y") {
                    ### Y coordinate is the widest
                    # Define the intervals
                    interval_width <- ceiling((rectangle_coordinates[2,2] - rectangle_coordinates[1,2])/number_of_partitions)
                    # Define the separating points
                    interval_points <- seq.int(rectangle_coordinates[1,2], rectangle_coordinates[2,2], by = interval_width)
                    # Add the maximum coordinate to the interval points (if not already present)
                    if (interval_points[length(interval_points)] != rectangle_coordinates[2,2]) {
                        interval_points <- c(interval_points, rectangle_coordinates[2,2])
                    }
                    ## Generate the final list of partitioned spectra
                    # Initialize the output
                    spectra_partitioned <- list()
                    # Populate the final list
                    for (p in 1:(length(interval_points) - 1)) {
                        # Generate the sublist of the spectra for that interval
                        spectra_interval <- list()
                        # Scroll the spectra...
                        for (s in 1:length(spectra)) {
                            # Populate the list...
                            if (spectra[[s]]@metaData$imaging$pos[2] >= interval_points[p] && spectra[[s]]@metaData$imaging$pos[2] < interval_points[p + 1]) {
                                spectra_interval <- append(spectra_interval, spectra[[s]])
                            }
                        }
                        # Add the spectra with the max coordinate
                        if (p == (length(interval_points) - 1)) {
                            # Scroll the spectra...
                            for (s in 1:length(spectra)) {
                                # Populate the list...
                                if (spectra[[s]]@metaData$imaging$pos[2] == interval_points[p + 1]) {
                                    spectra_interval <- append(spectra_interval, spectra[[s]])
                                }
                            }
                        }
                        spectra_partitioned[[p]] <- spectra_interval
                    }
                }
            } else {
                spectra_partitioned <- spectra
            }
        } else if (partitioning_method == "random") {
            #################### RANDOM
            ##### Split the spectra randomly in N parts
            ##### Install and load the required packages
            install_and_load_required_packages("caret")
            # Set the seed (make randomness reproducible)
            if (!is.null(seed)) {
                set.seed(seed)
            }
            # Partition the spectral list
            spectra_partitioned_IDs <- createFolds(y = rep("spectra", length(spectra)), k = number_of_partitions, returnTrain = FALSE)
            # Generate the split sublist of spectra
            spectra_partitioned <- list()
            # For each element of the list containing the IDs of the fold...
            for (s in 1:length(spectra_partitioned_IDs)) {
                # Extract the corresponding spectra and put them in a final list
                spectra_partitioned[[s]] <- spectra[spectra_partitioned_IDs[[s]]]
            }
        } else if (partitioning_method == "hca") {
            #################### HCA
            ### Detect and align peaks
            peaks <- peak_picking(spectra, peak_picking_algorithm = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", "SuperSmoother", "MAD"), tof_mode = tof_mode, SNR = ifelse(tof_mode == "reflectron" || tof_mode == "reflector", 5, 3), allow_parallelization = FALSE)
            peaks <- align_and_filter_peaks(peaks, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = 5)
            ### Generate the peaklist matrix
            peaklist <- intensityMatrix(peaks, spectra)
            ### Compute the distance matrix
            distance_matrix <- dist(peaklist, method = "euclidean")
            # Generate the dendrogram
            hca <- hclust(distance_matrix)
            ### Cut the tree to generate K number of sub-clusters
            hca_groups <- cutree(hca, k = number_of_partitions)
            ### Initialize the output
            spectra_partitioned <- list()
            ### For each subgroup to be isolated...
            for (p in 1:number_of_partitions) {
                # Index the spectra under in the selected subgroup of the HCA
                index <- which(hca_groups == p)
                spectra_hca <- spectra[index]
                # Store them in the final list
                spectra_partitioned[[p]] <- spectra_hca
            }
        }
        #################### Output
        return(spectra_partitioned)
    } else if (isMassSpectrum(spectra)) {
        ########## Single spectrum
        return(spectra)
    }
}





################################################################################





############################################### EXTRACT FEATURES FROM MODEL LIST
# The function takes an RData workspace as input, in which there are all the models. The model_list is a list in which each element corresponds to a list containing the model object, the feature list, the class list and the outcome list and its tuning/cross-validation performance. If the input is the model list, the same operations are performed (only the import of the model_list from the RData is skipped).
# The function returns a vector of features, sorted according to the importance in all the models, and a dataframe listing the features for each model.
extract_feature_list_from_model_list <- function(filepath_R, model_list_object = "model_list", features_to_return = 20) {
    if (!is.list(filepath_R)) {
        ### LOAD THE R WORKSPACE WITH THE MODEL LIST
        # Create a temporary environment
        temporary_environment <- new.env()
        # Load the workspace
        load(filepath_R, envir = temporary_environment)
        # Get the models (R objects) from the workspace
        model_list <- get(model_list_object, pos = temporary_environment)
    } else if (is.list(filepath_R)) {
        model_list <- filepath_R
    }
    # Get the list of models
    list_of_models <- names(model_list)
    ### Initialize the matrix for output (the number of rows must be equal to the highest amount of features in the models)
    highest_number_of_features <- 0
    for (md in 1:length(list_of_models)) {
        if (length(model_list[[md]]$features_model) > highest_number_of_features) {
            highest_number_of_features <- length(model_list[[md]]$features_model)
        }
    }
    ### Initialize the feature matrix
    model_features_matrix <- matrix(0, nrow = highest_number_of_features, ncol = length(list_of_models))
    colnames(model_features_matrix) <- list_of_models
    ### Initialize the feature dataframe (list of feature with feature ranking)
    feature_dataframe <- NULL
    ### Extract the features for each model
    for (md in 1:length(list_of_models)) {
        # Extract the features from the model list
        feature_vector <- model_list[[md]]$features_model
        # Remove the X from the features
        if (unlist(strsplit(as.character(feature_vector[1]),""))[1] == "X") {
            for (f in 1:length(feature_vector)) {
                name_splitted <- unlist(strsplit(feature_vector[f],""))
                feature_def <- name_splitted [2]
                for (i in 3:length(name_splitted)) {
                    feature_def <- paste(feature_def, name_splitted[i], sep = "")
                }
                feature_vector[f] <- feature_def
            }
        }
        # Fill the matrix
        model_features_matrix[, md] <- cbind(as.numeric(feature_vector))
        # Fill the dataframe
        feature_df <- data.frame(features = cbind(as.numeric(feature_vector)), rank = cbind(seq(1, length(feature_vector), by = 1)))
        if (is.null(feature_dataframe)) {
            feature_dataframe <- feature_df
        } else {
            feature_dataframe <- rbind(feature_dataframe, feature_df)
        }
    }
    # Sort the feature dataframe according to the rank
    feature_dataframe <- feature_dataframe[order(feature_dataframe$rank), ]
    # Return to a vector and extract the unique values
    model_features <- unique(as.numeric(feature_dataframe$features))
    # Output a certain number of features
    if (features_to_return > 0 && features_to_return < length(model_features)) {
        model_features <- model_features[1:features_to_return]
    } else {
        model_features <- model_features
    }
    ### Return
    return(list(model_features = model_features, model_features_matrix = model_features_matrix))
}





################################################################################





############################################ EXTRACT PERFORMANCE FROM MODEL LIST
# The function takes an RData workspace as input, in which there are all the models. The model_list is a list in which each element corresponds to a list containing the model object, the feature list, the class list and the outcome list and its tuning/cross-validation performance. If the input is the model list, the same operations are performed (only the import of the model_list from the RData is skipped).
# The function returns a matrix, listing the performance for each model.
extract_performance_matrix_from_model_list <- function(filepath_R, model_list_object = "model_list") {
    if (!is.list(filepath_R)) {
        ### LOAD THE R WORKSPACE WITH THE MODEL LIST
        # Create a temporary environment
        temporary_environment <- new.env()
        # Load the workspace
        load(filepath_R, envir = temporary_environment)
        # Get the models (R objects) from the workspace
        model_list <- get(model_list_object, pos = temporary_environment)
    } else if (is.list(filepath_R)) {
        model_list <- filepath_R
    }
    # Get the list of models
    list_of_models <- names(model_list)
    ### Initialize the matrix for output
    model_performance_matrix <- NULL
    ### Extract the performance for each model
    for (md in 1:length(list_of_models)) {
        # Extract the features from the model list
        model_performance <- model_list[[md]]$model_performance
        # Fill the matrix
        if (is.null(model_performance_matrix)) {
            single_model_performance <- matrix(0, nrow = 1, ncol = 2)
            colnames(single_model_performance) <- c("Model", names(model_performance))
            rownames(single_model_performance) <- list_of_models[md]
            single_model_performance[1,1] <- list_of_models[md]
            single_model_performance[1,2] <- model_performance
            model_performance_matrix <- single_model_performance
        } else {
            single_model_performance <- matrix(0, nrow = 1, ncol = 2)
            colnames(single_model_performance) <- c("Model", names(model_performance))
            rownames(single_model_performance) <- list_of_models[md]
            single_model_performance[1,1] <- list_of_models[md]
            single_model_performance[1,2] <- model_performance
            model_performance_matrix <- rbind(model_performance_matrix, single_model_performance)
        }
    }
    ### Return
    return(model_performance_matrix)
}




































































































########################################################################## GRAPH THEORY

################################################### GENETIC ALGORITHM FOR GRAPHS
# The function takes an adjacency matrix as input, converts it into a graph (using the 'igraph' package) and runs the genetic algorithm on the chromosome population generated by the graph nodes. The aim of the genetic algorithm is to identify and isolate a set of highly correlated observations (clique) from the rest of the dataset (independent set) by performing vertex and triangle mutations (free low-grade vertices and bind high-grade vertices) and selecting (fitness function) the population with a minimal change compared to the original, in such a way that the maximum amount of information is preserved.
genetic_algorithm_graph <- function(input_adjacency_matrix, graph_type = "Preferential", vertex_mutation = TRUE, triangle_mutation = TRUE, allow_parallelization = FALSE, number_of_high_degree_vertices_for_subgraph = 0, vertices_not_in_induced_subgraph = c("independent", "reassigned"), vertex_independency_threshold = 200, iterations_with_no_change = 5, max_GA_generations = 10, seed = 12345) {
    ##### Install and load the required packages
    install_and_load_required_packages(c("MALDIquant", "parallel", "doParallel", "foreach", "iterators", "igraph", "GA"))
    ### Generate the graph
    initial_graph <- graph_from_adjacency_matrix(input_adjacency_matrix)
    ### Generate a subgraph with only the high-degree vertices
    if (!is.null(number_of_high_degree_vertices_for_subgraph) && number_of_high_degree_vertices_for_subgraph > 0 && number_of_high_degree_vertices_for_subgraph < length(V(initial_graph))) {
        # Generate the dataframe with the vertex IDs and their degree
        degree_dataframe <- data.frame(ID = seq(1:length(V(initial_graph))), degree = degree(initial_graph))
        # Sort it according to the degree
        degree_dataframe_sorted <- degree_dataframe[order(degree_dataframe$degree, decreasing = TRUE), ]
        input_graph <- induced_subgraph(initial_graph, vids = degree_dataframe_sorted$ID[1:number_of_high_degree_vertices_for_subgraph], impl = "auto")
        # Regenerate the adjacency matrix according to the new graph
        input_adjacency_matrix <- as.matrix(as_adjacency_matrix(input_graph))
    } else {
        input_graph <- initial_graph
        input_adjacency_matrix <- input_adjacency_matrix
    }
    ########## Vertex mutation parameters
    # The vertex number is the number of spectra in the graph
    vertex_number <- length(V(input_graph))
    # Define the number of vertices to be sampled
    independent_vertex_number_sampling <- round(0.28 * vertex_number)
    # Define the threshold (number of bridges) below which the vertex is considered quasi-independent and almost belonging to the independent set.
    degree_threshold_independency <- vertex_independency_threshold
    # Output the sequence of the degree (number of connections) for each vertex
    degree_sequence <- degree(input_graph)
    # Generate the matrix listing the degrees (for each vertex, each vertex being a row)
    degree_matrix <- cbind(ids = 1:vertex_number, degree_sequence)
    # Pick the quasi-independent vertices according to the degree threshold
    quasi_independent_vertices <- which(degree_matrix[, "degree_sequence"] <= degree_threshold_independency)
    # Compute the number of quasi-independent vertices
    independent_vertex_number_graph <- length(quasi_independent_vertices)
    ########## Triangle mutation parameters
    # Define the number of triangles to be sampled
    triangle_number_sampling <- round(0.75 * vertex_number)
    # Isolate the triangles in the graph
    graph_triangles <- triangles(input_graph)
    # Extract the array displaying the triangles (vector of vertices, grouped by 3 for triangles) (ex x1,x2,x3,y1,y2,y3,z1,z2,z3,...)
    triangle_array <- array(graph_triangles)
    # Arrange the vertices in a matrix (each row is a triangle, three columns listing the vertices)
    triangle_matrix <- matrix(triangle_array, nrow = length(triangle_array)/3, byrow = TRUE)
    # Compute the number of triangles
    triangle_number_graph <- nrow(triangle_matrix)
    # Population size
    population_size <- vertex_number * 10
    if (allow_parallelization == TRUE) {
        ### PARALLEL BACKEND
        # Detect the number of cores
        cpu_thread_number <- detectCores(logical = TRUE)
        cpu_thread_number <- cpu_thread_number / 2
        if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
            install_and_load_required_packages("doMC")
            # Register the foreach backend
            registerDoMC(cores = cpu_thread_number)
        } else if (Sys.info()[1] == "Windows") {
            install_and_load_required_packages("doParallel")
            # Register the foreach backend
            cls <- makeCluster(cpu_thread_number)
            registerDoParallel(cls)
        }
    }
    # Establish if the parallel computation can be performed or not in the GA (otherwise it returns an error)
    if (allow_parallelization == TRUE) {
        if (cpu_thread_number > 1) {
            GA_parallel <- TRUE
        } else {
            GA_parallel <- FALSE
        }
    } else if (allow_parallelization == FALSE) {
        GA_parallel <- FALSE
    }
    ########## Run the genetic algorithm
    ############################################################### FUNCTIONS ZOPPIS
    ###################### SPLIT GRAPH
    # The function takes a chromosome (likely the final, after all the optimizations performed by the genetic algorithm, the chromosome with the best fitness) and generates the graph, through the adjacency matrix.
    FITN_SplitGraph <- function(chromosome) {
        # Calculate the size of the chromosome
        chromosome_size <- length(chromosome)
        # Calculate the number of 1 (simply the sum of the numbers in the chromosome, which is composed only of 0 and 1)
        number_of_ones <- sum(chromosome)
        # print(chromosome)
        # print(number_of_ones)
        # Generate the adjacency matrix (n x n)
        chromosome_adjacency_matrix <- matrix(0, nrow = chromosome_size, ncol = chromosome_size)
        # Fill in the adjacency matrix
        chromosome_adjacency_matrix[which(chromosome == 1),] <- rep(chromosome, each = number_of_ones)
        # Set the diagonal to zero
        diag(chromosome_adjacency_matrix) <- 0
        # Generate the graph from the adjacency matrix
        chromosome_graph <- graph_from_adjacency_matrix(chromosome_adjacency_matrix, mode = "undirect")
        # Return the fitness value (the changes made, the difference between the original matrix and the matrix after the mutation: since the fitness has to be maximized, it is calculated as the negative sum of the differences between the matrices). Since the matrix is symmetrical, the difference has to be divided by 2, to avoid it being considered twice.
        fitness_GA <- -sum(abs(input_adjacency_matrix - chromosome_adjacency_matrix))/2
        return(fitness_GA)
    }
    ###################### SPLIT MUTATION
    # This mutation functions puts some vertices in a triangle and some vertices in the independent set.
    SplitMutation <- function(object, parent_ind) {
        parent_chromosome <- as.vector(object@population[parent_ind, ])
        #parentCR <- c(1,0,0,1,0,1,0,1,0,1)
        # Calculate the size of the chromosome
        chromosome_length <- length(parent_chromosome)
        # Create a copy of the parent chromosome, which will be subjected to mutation
        mutated_chromosome <- parent_chromosome
        ### Triangle Mutation
        if (triangle_mutation == TRUE) {
            # Compute the uniform distribution
            rsamp_RowIX <- runif(triangle_number_sampling, 1, triangle_number_graph)
            RoundedSampl_RowIX <- round(rsamp_RowIX)
            # For each triangle to be sampled, extract the triangle to be sampled from the triangle matrix (so this extracts the vertices of the triangles) and set them to 1 (be part of a triangle) in the mutated chromosome.
            for (i in 1:triangle_number_sampling) {
                triangle_sample <- triangle_matrix[RoundedSampl_RowIX[i],]
                mutated_chromosome[triangle_sample] <- 1
            }
        }
        ### Vertex Mutation
        if (vertex_mutation == TRUE) {
            # Compute the uniform distribution
            rsamp_IndVert <- runif(independent_vertex_number_sampling, 1, independent_vertex_number_graph)
            RoundedSampl_IndVertIX <- round(rsamp_IndVert)
            # Extracts the quasi-independent vertices to be extracted and be set to 0 (part of the independent set) in the mutatet chromosome
            vertex_sample <- quasi_independent_vertices[RoundedSampl_IndVertIX]
            mutated_chromosome[vertex_sample] <- 0
        }
        # output for user defined selection
        # mutated_chromosome <- list(population = object@population[sel, , drop = FALSE], fitness = f[sel])
        # Return the mutated chromosome
        return(as.vector(mutated_chromosome))
    }
    ###################### SPLIT GRAPH 2
    FinalSplitGraph <- function(chromosome) {
        # Calculate the size of the chromosome
        chromosome_size <- length(chromosome)
        # Calculate the number of 1 (simply the sum of the numbers in the chromosome, which is composed only of 0 and 1)
        number_of_ones <- sum(chromosome)
        # print(chromosome)
        # print(number_of_ones)
        # Generate the adjacency matrix (n x n)
        chromosome_adjacency_matrix <- matrix(0, nrow = chromosome_size, ncol = chromosome_size)
        # Fill in the adjacency matrix
        chromosome_adjacency_matrix[which(chromosome == 1),] <- rep(chromosome, each = number_of_ones)
        # Set the diagonal to zero
        diag(chromosome_adjacency_matrix) <- 0
        # Generate the graph from the adjacency matrix
        chromosome_graph <- graph_from_adjacency_matrix(chromosome_adjacency_matrix, mode = "undirect")
        # Return the graph
        return(chromosome_graph)
    }
    ################################################################################
    #ptm <- proc.time()
    # Pass the variables to the parallel backened
    #clusterExport(cls, varlist = c("input_adjacency_matrix", "FITN_SplitGraph", "vertex_number", "SplitMutation", "max_GA_generations", "population_size", "GA_parallel", "seed", "triangle_mutation", "vertex_mutation"))
    GA_model <- ga(type = "binary",
                   fitness = FITN_SplitGraph,
                   nBits = vertex_number,
                   selection = gabin_rwSelection,
                   mutation = SplitMutation,
                   maxiter = max_GA_generations,
                   run = iterations_with_no_change,
                   popSize = population_size,
                   pcrossover = 0.8,
                   pmutation = 0.1,
                   elitism = base::max(1, round(population_size*0.1)),
                   keepBest = TRUE,
                   parallel = GA_parallel,
                   seed = seed
    )
    #finaltime <- proc.time() - ptm
    #out <- summary(GA_model)
    #print(out)
    # Stop the cluster for parallel computing (only on Windows)
    if (Sys.info()[1] == "Windows" && allow_parallelization == TRUE) {
        stopCluster(cls)
    }
    # Extract the final graph
    best_solution <- GA_model@solution[1,]
    final_graph <- FinalSplitGraph(best_solution)
    # Generate the adjacency matrix again from the final graph
    output_adjacency_matrix <- as.matrix(as_adjacency_matrix(final_graph, type = "both"))
    rownames(output_adjacency_matrix) <- rownames(input_adjacency_matrix)
    colnames(output_adjacency_matrix) <- colnames(input_adjacency_matrix)
    ########## Extract the final chromosome (with the vertex IDs) and Reassign the leftover vertices
    ##### Subgraph generated
    if (!is.null(number_of_high_degree_vertices_for_subgraph) && number_of_high_degree_vertices_for_subgraph > 0 && number_of_high_degree_vertices_for_subgraph < length(V(initial_graph))) {
        # Extract the final chromosome (with the vertex IDs)
        if (is.matrix(GA_model@solution)) {
            final_chromosome <- as.vector(GA_model@solution[1,])
            names(final_chromosome) <- degree_dataframe_sorted$ID[1:number_of_high_degree_vertices_for_subgraph]
        } else {
            final_chromosome <- as.vector(GA_model@solution)
            names(final_chromosome) <- degree_dataframe_sorted$ID[1:number_of_high_degree_vertices_for_subgraph]
        }
        ##### Independent set
        if (vertices_not_in_induced_subgraph == "independent") {
            # Extract the vertex IDs in the chromosome
            vertex_IDs_chromosome <- as.numeric(names(final_chromosome))
            # Extract the vertex IDs not in the chromosome
            vertex_IDs_outside_chromosome <- as.numeric(V(initial_graph)[!(V(initial_graph) %in% vertex_IDs_chromosome)])
            # Generate the vector of 0 (all the vertices go into the independent set)
            vertices_outside_chromosome <- numeric(length=length(vertex_IDs_outside_chromosome))
            names(vertices_outside_chromosome) <- vertex_IDs_outside_chromosome
            # Add this vector to the final chromosome
            final_chromosome <- append(final_chromosome, vertices_outside_chromosome)
            # Convert it into a dataframe to sort it according to the names
            final_chromosome_df <- as.data.frame(cbind(ID = as.numeric(names(final_chromosome)), Value = as.numeric(final_chromosome)))
            final_chromosome_df <- final_chromosome_df[order(final_chromosome_df$ID), ]
            # Convert it back to a vector
            final_chromosome <- as.vector(final_chromosome_df$Value)
            names(final_chromosome) <- final_chromosome_df$ID
        }
    } else {
    ##### NO subgraph generated
        # Extract the final chromosome (with the vertex IDs)
        if (nrow(as.matrix(GA_model@solution)) > 1) {
            final_chromosome <- as.vector(GA_model@solution[1,])
            names(final_chromosome) <- colnames(output_adjacency_matrix)
        } else {
            final_chromosome <- as.vector(GA_model@solution)
            names(final_chromosome) <- colnames(output_adjacency_matrix)
        }
    }
    ### Return
    return(list(GA_model = GA_model, GA_model_solution = best_solution, input_graph = input_graph, final_graph = final_graph, output_adjacency_matrix = output_adjacency_matrix, final_chromosome = final_chromosome))
}





################################################################################





####################### FROM GENETIC ALGORITHM ON GRAPH TO SPECTRA AND MS IMAGES
# The function takes the result of the genetic algorithm optimization (the genetic model object) and the list of spectra used to generate the adjacency matrix for the function 'genetic_algorithm_graph', it extracts the final chromosome (generated after the optimization) and it mirrors it onto the list of spectra, so that the output is a list of spectra belonging to the clique and a list of spectra belonging to the independent set. In addition, a list of spectra for plotting purposes is returned.
from_GA_to_MS <- function(final_chromosome_GA, spectra, spectra_format = "imzml") {
    ##### Install and load the required packages
    install_and_load_required_packages(c("MALDIquant", "parallel", "doParallel", "igraph", "GA"))
    ########## Isolate the final chromosome (after mutations and fitness optimization)
    ## The spectra with identified with a 0 are the spectra outside the clique, while the spectra identified by a 1 belong to the clique (take only the first row if it is a matrix)
    ## Identify the spectra belonging to the clique and to the independent set
    spectra_clique_id <- which(final_chromosome_GA == 1)
    spectra_independent_id <- which(final_chromosome_GA == 0)
    ### Generate the lists (clique and independent)
    if (is.null(names(spectra))) {
        spectra_clique <- spectra[spectra_clique_id]
        spectra_independent <- spectra[spectra_independent_id]
    } else {
        spectra_clique <- spectra[names(spectra_clique_id)]
        spectra_independent <- spectra[names(spectra_independent_id)]
    }
    ### Replace intensities with a number resembling their belonging to the clique (red) or independent set (green)
    spectra_clique_for_plotting <- spectra_clique
    spectra_independent_for_plotting <- spectra_independent
    # If there are any spectra in the clique, replace the intensity values (for plotting) and add them to the final list
    if (isMassSpectrumList(spectra_clique_for_plotting)) {
        for (spp in 1:length(spectra_clique_for_plotting)) {
            spectra_clique_for_plotting[[spp]]@intensity <- rep(1, length(spectra_clique_for_plotting[[spp]]@intensity))
        }
    } else if (isMassSpectrum(spectra_clique_for_plotting)) {
        spectra_clique_for_plotting@intensity <- rep(1, length(spectra_clique_for_plotting@intensity))
    } else if (length(spectra_clique_for_plotting) == 0) {
        spectra_clique_for_plotting <- spectra_clique_for_plotting
    }
    # If there are any spectra in the independent set, replace the intensity values (for plotting) and add them to the final list
    if (isMassSpectrumList(spectra_independent_for_plotting)) {
        for (spp in 1:length(spectra_independent_for_plotting)) {
            spectra_independent_for_plotting[[spp]]@intensity <- rep(0.5, length(spectra_independent_for_plotting[[spp]]@intensity))
        }
    } else if (isMassSpectrum(spectra_independent_for_plotting)) {
        spectra_independent_for_plotting@intensity <- rep(0.5, length(spectra_independent_for_plotting@intensity))
    } else if (length(spectra_independent_for_plotting) == 0) {
        spectra_independent_for_plotting <- spectra_independent_for_plotting
    }
    ## Append to the final spectral lists
    spectra_all_for_plotting <- append(spectra_clique_for_plotting, spectra_independent_for_plotting)
    # Fix the spectra name (for output)
    #spectra_all_for_plotting <- replace_sample_name(spectra_all_for_plotting, spectra_format = spectra_format)
    ### Rearrange spectra for reproducibility
    #spectra_all_for_plotting <- rearrange_spectral_dataset(spectra = spectra_all_for_plotting, rearranging_method = "space")
    #file_explan <- "Performance"
    #FILE_NAME <- sprintf("%s%s%s%s%s",vertex_number,"_",graph_type,"_",file_explan)
    #sizeGrWindow(10,5)
    ##### Print outputs
    #BestFitv <- abs(GA_model@fitnessValue)
    #x <- (BestFitv/choose(vertex_number,2))*100
    #print("difference %")
    #print(x)
    #print("cpu time")
    #print(finaltime)
    #file_explan <- "SolPlot"
    #FILE_NAME <- sprintf("%s%s%s%s%s",vertex_number,"_",graph_type,"_",file_explan)
    ### Return
    return(list(spectra_clique = spectra_clique, spectra_independent = spectra_independent, spectra_for_plotting = spectra_all_for_plotting))
}





################################################################################





#################################################### GRAPH SEGMENTATION FUNCTION
# The function returns (for the imzML MSI dataset provided as the variable spectra) the list of spectra in the clique, the list of spectra in the independent set and the MS images (with pixels related to spectra in the clique in red and pixels related to spectra in the independent set in green), the initial and the final graph.
# It returns a NULL value if the segmentation is not possible due to incompatibilities between the features in the dataset and the ones provided by the model.
graph_MSI_segmentation <- function(filepath_imzml, preprocessing_parameters = list(mass_range = c(800,3000), transformation_algorithm = NULL, smoothing_algorithm = NULL, smoothing_strength = "medium", baseline_subtraction_algorithm = "SNIP", baseline_subtraction_iterations = 200, normalization_algorithm = "TIC", normalization_mass_range = NULL, process_spectra_in_packages_of = 0, spectral_alignment_method = NULL), allow_parallelization = FALSE, peak_picking_algorithm = "SuperSmoother", deisotope_peaklist = FALSE, SNR = 5, tof_mode = "reflectron", peak_filtering_frequency_threshold_percent = 5, low_intensity_peak_removal_threshold_percent = 0, low_intensity_peak_removal_threshold_method = "element-wise", custom_feature_vector = NULL, correlation_method_for_adjacency_matrix = "pearson", correlation_threshold_for_adjacency_matrix = 0.90, pvalue_threshold_for_adjacency_matrix = 0.05, number_of_high_degree_vertices_for_subgraph = 0, vertices_not_in_induced_subgraph = c("independent", "reassigned"), max_GA_generations = 10, iterations_with_no_change = 5, plot_figures = TRUE, plot_graphs = TRUE, plot_legends = c("sample name", "legend", "plot name"), number_of_spectra_partitions = 1, partitioning_method = "space", seed = 12345, spectra_format = "imzml") {
    # Install and load the required packages
    install_and_load_required_packages(c("MALDIquantForeign", "MALDIquant", "parallel", "caret", "pls", "tcltk", "kernlab", "pROC", "e1071", "igraph", "GA"))
    ### Import the dataset (if filepath_imzml is not already a list of spectra)
    if (!isMassSpectrumList(filepath_imzml) && !isMassSpectrum(filepath_imzml)) {
        if (!is.null(preprocessing_parameters$mass_range)) {
            spectra <- importImzMl(filepath_imzml, massRange = preprocessing_parameters$mass_range)
        } else {
            spectra <- importImzMl(filepath_imzml)
        }
    } else if (isMassSpectrumList(filepath_imzml) || isMassSpectrum(filepath_imzml)) {
        spectra <- filepath_imzml
    }
    ### Identify the original spectra by setting the names (better unique identification of the spectra)
    spectra <- replace_sample_name_list(spectra, spectra_format = spectra_format)
    ### Rearrange spectra for reproducibility
    #spectra <- rearrange_spectral_dataset(spectra = spectra, rearranging_method = "space")
    #################### Partition the spectra
    if (number_of_spectra_partitions > 1) {
        # Initialize the output
        final_spectra_clique <- list()
        final_spectra_independent <- list()
        graph_spectra_plot_list <- list()
        final_graph_plot_list <- list()
        ga_model_plot_list <- list()
        spectra_all_for_plotting <- list()
        ##### Generate the partitioned list
        spectra_partitioned <- partition_spectral_dataset(spectra, partitioning_method = partitioning_method, number_of_partitions = number_of_spectra_partitions, seed = seed, tof_mode = tof_mode)
        ##### For each fold...
        for (prt in 1:length(spectra_partitioned)) {
            spectra_partition <- unlist(spectra_partitioned[[prt]])
            ##### Do everything only of there are more than one spectra in the partition
            if (isMassSpectrumList(spectra_partition)) {
                ### Use only the features from the model
                peaklist_matrix <- generate_custom_intensity_matrix(spectra_partition, custom_feature_vector = custom_feature_vector, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = deisotope_peaklist, peak_picking_SNR = SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization)
                # Compute the adjacency matrix
                input_adjacency_matrix <- generate_adjacency_matrix(peaklist_matrix, correlation_method = correlation_method_for_adjacency_matrix, correlation_threshold = correlation_threshold_for_adjacency_matrix, pvalue_threshold = pvalue_threshold_for_adjacency_matrix)
                # Run the genetic algorithm
                GA_output <- genetic_algorithm_graph(input_adjacency_matrix, max_GA_generations = max_GA_generations, seed = seed, number_of_high_degree_vertices_for_subgraph = number_of_high_degree_vertices_for_subgraph, vertices_not_in_induced_subgraph = vertices_not_in_induced_subgraph, allow_parallelization = allow_parallelization, iterations_with_no_change = iterations_with_no_change)
                # Record the plots
                if (plot_graphs == TRUE) {
                    plot(GA_output$input_graph)
                    graph_spectra_plot_list[[prt]] <- recordPlot()
                    plot(GA_output$final_graph)
                    final_graph_plot_list[[prt]] <- recordPlot()
                } else {
                    graph_spectra_plot_list[[prt]] <- NULL
                    final_graph_plot_list[[prt]] <- NULL
                }
                if (plot_figures == TRUE) {
                    plot(GA_output$GA_model)
                    ga_model_plot_list[[prt]] <- recordPlot()
                } else {
                    ga_model_plot_list[[prt]] <- NULL
                }
                # From the optimised graph, extract the spectra for plotting, the adjacency matrix and the figures
                MS_from_GA <- from_GA_to_MS(GA_output$final_chromosome, spectra = spectra_partition)
                # Record the plot
                final_spectra_clique[[prt]] <- MS_from_GA$spectra_clique
                final_spectra_independent[[prt]] <- MS_from_GA$spectra_independent
                # Extract the spectra for plotting and the spectra in the clique and in the independent set
                spectra_all_for_plotting <- append(spectra_all_for_plotting, MS_from_GA$spectra_for_plotting)
            } else if (isMassSpectrum(spectra_partition)) {
                ##### The spectra from the partition contain actually one spectrum
                # Place it in the independent set
                final_spectra_independent[[prt]] <- spectra_partition
                # Replace its intensity with 0.5 (it will be plotted as independent)
                single_spectrum_partition_for_plotting <- spectra_partition
                single_spectrum_partition_for_plotting@intensity <- rep(0.5, length(single_spectrum_partition_for_plotting@intensity))
                # Add it to the final list of spectra_for_plotting
                spectra_all_for_plotting <- append(spectra_all_for_plotting, single_spectrum_partition_for_plotting)
            }
        }
        ### Generate the MS images (slices)
        slices <- msiSlices(spectra_all_for_plotting, center = spectra_all_for_plotting[[1]]@mass[(length(spectra_all_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
        ### Plot the MS image
        if (plot_figures == TRUE) {
            plotMsiSlice(slices, legend = FALSE, scale = F)
            if ("legend" %in% plot_legends) {
                ### Add the legend
                legend(x = "bottomright", legend = c("clique","independent set"), fill = c("red","green"), xjust = 0.5, yjust = 0.5)
            }
            if ("sample name" %in% plot_legends) {
                legend(x = "topright", legend = spectra_all_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
            }
            if ("plot name" %in% plot_legends) {
                legend(x = "topleft", legend = "Graph segmentation", xjust = 0.5, yjust = 0.5)
            }
            # Record the plot
            msi_segmentation <- recordPlot()
        }
        ########## Return
        return(list(spectra_for_plotting = spectra_all_for_plotting, spectra_clique = final_spectra_clique, spectra_independent = final_spectra_independent, msi_segmentation = msi_segmentation, initial_graph_spectra_plot_list = graph_spectra_plot_list, final_graph_plot_list = final_graph_plot_list, ga_model_plot = ga_model_plot_list))
    } else if (number_of_spectra_partitions <= 1) {
        #################### NO Partitioning of the spectra
        ### Initialize the outputs
        graph_spectra_plot <- NULL
        msi_segmentation <- NULL
        ga_model_plot <- NULL
        final_graph_plot <- NULL
        ### Use only the features from the model
        peaklist_matrix <- generate_custom_intensity_matrix(spectra, custom_feature_vector = custom_feature_vector, tof_mode = tof_mode, preprocessing_parameters = preprocessing_parameters, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = deisotope_peaklist, peak_picking_SNR = SNR, peak_filtering_frequency_threshold_percent = peak_filtering_frequency_threshold_percent, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, allow_parallelization = allow_parallelization)
        ### If a NULL value is returned, it means that the model is incompatible with the features in the dataset
        if (!is.null(peaklist_matrix)) {
            # Compute the adjacency matrix
            input_adjacency_matrix <- generate_adjacency_matrix(peaklist_matrix, correlation_method = correlation_method_for_adjacency_matrix, correlation_threshold = correlation_threshold_for_adjacency_matrix, pvalue_threshold = pvalue_threshold_for_adjacency_matrix)
            # Run the genetic algorithm
            GA_output <- genetic_algorithm_graph(input_adjacency_matrix, max_GA_generations = max_GA_generations, seed = seed, number_of_high_degree_vertices_for_subgraph = number_of_high_degree_vertices_for_subgraph, vertices_not_in_induced_subgraph = vertices_not_in_induced_subgraph, allow_parallelization = allow_parallelization, iterations_with_no_change = iterations_with_no_change)
            # Record the plots
            if (plot_graphs == TRUE) {
                plot(GA_output$input_graph)
                graph_spectra_plot <- recordPlot()
                plot(GA_output$final_graph)
                final_graph_plot <- recordPlot()
            } else {
                graph_spectra_plot <- NULL
                final_graph_plot <- NULL
            }
            if (plot_figures == TRUE) {
                plot(GA_output$GA_model)
                ga_model_plot <- recordPlot()
            } else {
                ga_model_plot <- NULL
            }
            # From the optimised graph, extract the spectra for plotting, the adjacency matrix and the figures
            MS_from_GA <- from_GA_to_MS(final_chromosome_GA = GA_output$final_chromosome, spectra = spectra)
            # Extract the spectra for plotting
            spectra_all_for_plotting <- MS_from_GA$spectra_for_plotting
            ### Generate the MS images (slices)
            slices <- msiSlices(spectra_all_for_plotting, center = spectra_all_for_plotting[[1]]@mass[(length(spectra_all_for_plotting[[1]]@mass)/2)], tolerance = 1, adjust = TRUE, method = "median")
            ### Plot the MS image
            if (plot_figures == TRUE) {
                plotMsiSlice(slices, legend = FALSE, scale = F)
                if ("legend" %in% plot_legends) {
                    ### Add the legend
                    legend(x = "bottomright", legend = c("clique","independent set"), fill = c("red","green"), xjust = 0.5, yjust = 0.5)
                }
                if ("sample name" %in% plot_legends) {
                    legend(x = "topright", legend = spectra_all_for_plotting[[1]]@metaData$file[1], xjust = 0.5, yjust = 0.5)
                }
                if ("plot name" %in% plot_legends) {
                    legend(x = "topleft", legend = "Graph segmentation", xjust = 0.5, yjust = 0.5)
                }
                # Record the plot
                msi_segmentation <- recordPlot()
            } else {
                msi_segmentation <- NULL
            }
            ### Return
            return(list(spectra_for_plotting = spectra_all_for_plotting, spectra_clique = MS_from_GA$spectra_clique, spectra_independent = MS_from_GA$spectra_independent, msi_segmentation = msi_segmentation, graph_spectra_plot = graph_spectra_plot, final_graph_plot = final_graph_plot, ga_model_plot = ga_model_plot))
        } else {
            return(NULL)
        }
    }
}





################################################################################
























































####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################































#################### MS PIXEL TYPER ####################





### Program version (Specified by the program writer!!!!)
R_script_version <- "2017.04.13.1"
### GitHub URL where the R file is
github_R_url <- "https://raw.githubusercontent.com/gmanuel89/Public-R-UNIMIB/master/MS%20PIXEL%20TYPER.R"
### Name of the file when downloaded
script_file_name <- "MS PIXEL TYPER"
# Change log
change_log <- "1. Bugfix\n2. Fixed GUI"






############## INSTALL AND LOAD THE REQUIRED PACKAGES
install_and_load_required_packages(c("tcltk", "parallel"), repository = "http://cran.mirror.garr.it/mirrors/CRAN/")






###################################### Initialize the variables (default values)
filepath_import <- NULL
tof_mode <- "linear"
output_folder <- getwd()
spectra_format <- "imzml"
peaks_filtering <- TRUE
low_intensity_peak_removal_threshold_method <- "element-wise"
peak_picking_algorithm <- "SuperSmoother"
file_type_export_matrix <- "csv"
file_type_export_images <- "png"
spectra <- NULL
peaks <- NULL
allow_parallelization <- FALSE
transform_data <- FALSE
transform_data_algorithm <- NULL
smoothing <- TRUE
smoothing_algorithm <- "SavitzkyGolay"
smoothing_strength <- "medium"
baseline_subtraction <- TRUE
baseline_subtraction_algorithm <- "SNIP"
baseline_subtraction_iterations <- 200
normalization <- TRUE
normalization_algorithm <- "TIC"
normalization_mass_range <- NULL
preprocess_spectra_in_packages_of <- 0
mass_range <- c(3000, 15000)
spectral_alignment <- FALSE
spectral_alignment_algorithm <- NULL
spectral_alignment_reference <- NULL
peaks_deisotoping <- FALSE
RData_file_integrity <- FALSE
classification_mode <- "pixel"
pixel_grouping <- "single"
decision_method_ensemble <- "majority"
vote_weights_ensemble <- "equal"
number_of_hca_nodes <- 4
moving_window_size <- 5




################## Values of the variables (for displaying and dumping purposes)
tof_mode_value <- "Linear"
filepath_import_value <- NULL
output_folder_value <- output_folder
spectra_format_value <- "imzML"
peaks_filtering_value <- "YES"
peak_picking_algorithm_value <- "Super Smoother"
low_intensity_peak_removal_threshold_method_value <- "element-wise"
allow_parallelization_value <- "NO"
transform_data_value <- "NO"
smoothing_value <- paste("YES", "\n( ", "SavitzkyGolay", ",\n" , "medium", " )", sep = "")
baseline_subtraction_value <- paste("YES", "\n( ", "SNIP", ",\niterations:", "200", " )", sep = "")
normalization_value <- paste("YES", "\n( ", "TIC", " )", sep = "")
spectral_alignment_value <- "NO"
spectral_alignment_algorithm_value <- ""
spectral_alignment_reference_value <- ""
peaks_deisotoping_value <- "NO"
RData_file_integrity_value <- "NO RDATA FILE\nSELECTED"
classification_mode_value <- "pixel"
decision_method_ensemble_value <- "majority"
vote_weights_ensemble_value <- "equal"






##################################################### DEFINE WHAT THE BUTTONS DO

##### Check for updates (from my GitHub page) (it just updates the label telling the user if there are updates) (it updates the check for updates value that is called by the label)
check_for_updates_function <- function() {
    ### Initialize the version number
    online_version_number <- NULL
    ### Initialize the variable that says if there are updates
    update_available <- FALSE
    ### Initialize the change log
    online_change_log <- "Bug fixes"
    # Check if there is internet connection by pinging a website
    there_is_internet <- check_internet_connection(method = "getURL", website_to_ping = "www.google.it")
    # Check for updates only in case of working internet connection
    if (there_is_internet == TRUE) {
        try({
            ### Read the file from the web (first 10 lines)
            online_file <- readLines(con = github_R_url)
            ### Retrieve the version number
            for (l in online_file) {
                if (length(grep("R_script_version", l, fixed = TRUE)) > 0) {
                    # Isolate the "variable" value
                    online_version_number <- unlist(strsplit(l, "R_script_version <- ", fixed = TRUE))[2]
                    # Remove the quotes
                    online_version_number <- unlist(strsplit(online_version_number, "\""))[2]
                    break
                }
            }
            ### Retrieve the change log
            for (l in online_file) {
                if (length(grep("change_log", l, fixed = TRUE)) > 0) {
                    # Isolate the "variable" value
                    online_change_log <- unlist(strsplit(l, "change_log <- ", fixed = TRUE))[2]
                    # Remove the quotes
                    online_change_log_split <- unlist(strsplit(online_change_log, "\""))[2]
                    # Split at the \n
                    online_change_log_split <- unlist(strsplit(online_change_log_split, "\\\\n"))
                    # Put it back to the character
                    online_change_log <- ""
                    for (o in online_change_log_split) {
                        online_change_log <- paste(online_change_log, o, sep = "\n")
                    }
                    break
                }
            }
            ### Split the version number in YYYY.MM.DD
            online_version_YYYYMMDDVV <- unlist(strsplit(online_version_number, ".", fixed = TRUE))
            ### Compare with the local version
            local_version_YYYYMMDDVV = unlist(strsplit(R_script_version, ".", fixed = TRUE))
            ### Check the versions (from the Year to the Day)
            # Check the year
            if (as.numeric(local_version_YYYYMMDDVV[1]) < as.numeric(online_version_YYYYMMDDVV[1])) {
                update_available <- TRUE
            }
            # If the year is the same (update is FALSE), check the month
            if (update_available == FALSE) {
                if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) < as.numeric(online_version_YYYYMMDDVV[2]))) {
                    update_available <- TRUE
                }
            }
            # If the month and the year are the same (update is FALSE), check the day
            if (update_available == FALSE) {
                if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) == as.numeric(online_version_YYYYMMDDVV[2])) && (as.numeric(local_version_YYYYMMDDVV[3]) < as.numeric(online_version_YYYYMMDDVV[3]))) {
                    update_available <- TRUE
                }
            }
            # If the day and the month and the year are the same (update is FALSE), check the daily version
            if (update_available == FALSE) {
                if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) == as.numeric(online_version_YYYYMMDDVV[2])) && (as.numeric(local_version_YYYYMMDDVV[3]) == as.numeric(online_version_YYYYMMDDVV[3])) && (as.numeric(local_version_YYYYMMDDVV[4]) < as.numeric(online_version_YYYYMMDDVV[4]))) {
                    update_available <- TRUE
                }
            }
            ### Return messages
            if (is.null(online_version_number)) {
                # The version number could not be ckecked due to internet problems
                # Update the label
                check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked: connection problems", sep = "")
            } else {
                if (update_available == TRUE) {
                    # Update the label
                    check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdate available: ", online_version_number, sep = "")
                } else {
                    # Update the label
                    check_for_updates_value <- paste("Version: ", R_script_version, "\nNo updates available", sep = "")
                }
            }
        }, silent = TRUE)
    }
    ### Something went wrong: library not installed, retrieving failed, errors in parsing the version number
    if (is.null(online_version_number)) {
        # Update the label
        check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked: connection problems", sep = "")
    }
    # Escape the function
    .GlobalEnv$update_available <- update_available
    .GlobalEnv$online_change_log <- online_change_log
    .GlobalEnv$check_for_updates_value <- check_for_updates_value
    .GlobalEnv$online_version_number <- online_version_number
}

##### Download the updated file (from my GitHub page)
download_updates_function <- function() {
    # Download updates only if there are updates available
    if (update_available == TRUE) {
        # Initialize the variable which says if the file has been downloaded successfully
        file_downloaded <- FALSE
        # Choose where to save the updated script
        tkmessageBox(title = "Download folder", message = "Select where to save the updated script file", icon = "info")
        download_folder <- tclvalue(tkchooseDirectory())
        if (!nchar(download_folder)) {
            # Get the output folder from the default working directory
            download_folder <- getwd()
        }
        # Go to the working directory
        setwd(download_folder)
        tkmessageBox(message = paste("The updated script file will be downloaded in:\n\n", download_folder, sep = ""))
        # Download the file
        try({
            download.file(url = github_R_url, destfile = paste(script_file_name, " (", online_version_number, ").R", sep = ""), method = "auto")
            file_downloaded <- TRUE
        }, silent = TRUE)
        if (file_downloaded == TRUE) {
            tkmessageBox(title = "Updated file downloaded!", message = paste("The updated script, named:\n\n", script_file_name, "\n\nhas been downloaded to:\n\n", download_folder, "\n\nClose everything, delete this file and run the script from the new file!", sep = ""), icon = "info")
            tkmessageBox(title = "Changelog", message = paste("The updated script contains the following changes:\n", online_change_log, sep = ""), icon = "info")
        } else {
            tkmessageBox(title = "Connection problem", message = paste("The updated script file could not be downloaded due to internet connection problems!\n\nManually download the updated script file at:\n\n", github_R_url, sep = ""), icon = "warning")
        }
    } else {
        tkmessageBox(title = "No update available", message = "NO UPDATES AVAILABLE!\n\nThe latest version is running!", icon = "info")
    }
}

##### Preprocessing window
preprocessing_window_function <- function() {
    ##### Functions
    # Transform the data
    transform_data_choice <- function() {
        # Catch the value from the menu
        transform_data <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "NO")
        # Default
        if (transform_data == "YES") {
            transform_data <- TRUE
            # Ask for the algorithm
            transform_data_algorithm <- select.list(c("sqrt","log","log2","log10"), title = "Choose", multiple = FALSE, preselect = "log10")
            # Default
            if (transform_data_algorithm == "") {
                transform_data_algorithm <- "sqrt"
            }
        } else if (transform_data == "NO" || transform_data == "") {
            transform_data <- FALSE
            transform_data_algorithm <- NULL
        }
        # Set the value of the displaying label
        if (transform_data == TRUE) {
            transform_data_value <- paste("YES", "\n( ", transform_data_algorithm, " )", sep = "")
        } else {
            transform_data_value <- "NO"
        }
        transform_data_value_label <- tklabel(preproc_window, text = transform_data_value, font = label_font, bg = "white", width = 20, height = 2)
        tkgrid(transform_data_value_label, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$transform_data <- transform_data
        .GlobalEnv$transform_data_algorithm <- transform_data_algorithm
        .GlobalEnv$transform_data_value <- transform_data_value
    }
    # Smoothing
    smoothing_choice <- function() {
        # Catch the value from the menu
        smoothing <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "YES")
        # Default
        if (smoothing == "YES" || smoothing == "") {
            smoothing <- TRUE
            # Ask for the algorithm
            smoothing_algorithm <- select.list(c("SavitzkyGolay","MovingAverage"), title = "Choose", multiple = FALSE, preselect = "SavitzkyGolay")
            # Default
            if (smoothing_algorithm == "") {
                smoothing_algorithm <- "SavitzkyGolay"
            }
            # Strength
            smoothing_strength <- select.list(c("medium","strong","stronger"), title = "Choose", multiple = FALSE, preselect = "medium")
            if (smoothing_strength == "") {
                smoothing_strength <- "medium"
            }
        } else if (smoothing == "NO") {
            smoothing <- FALSE
            smoothing_algorithm <- NULL
        }
        # Set the value of the displaying label
        if (smoothing == TRUE) {
            smoothing_value <- paste("YES", "\n( ", smoothing_algorithm, ",\n" , smoothing_strength, " )", sep = "")
        } else {
            smoothing_value <- "NO"
        }
        smoothing_value_label <- tklabel(preproc_window, text = smoothing_value, font = label_font, bg = "white", width = 20, height = 3)
        tkgrid(smoothing_value_label, row = 4, column = 2, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$smoothing <- smoothing
        .GlobalEnv$smoothing_strength <- smoothing_strength
        .GlobalEnv$smoothing_algorithm <- smoothing_algorithm
        .GlobalEnv$smoothing_value <- smoothing_value
    }
    # Baseline subtraction
    baseline_subtraction_choice <- function() {
        # Catch the value from the menu
        baseline_subtraction <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "YES")
        # Default
        if (baseline_subtraction == "YES" || baseline_subtraction == "") {
            baseline_subtraction <- TRUE
            # Ask for the algorithm
            baseline_subtraction_algorithm <- select.list(c("SNIP","TopHat","ConvexHull","median"), title = "Choose", multiple = FALSE, preselect = "SNIP")
            # SNIP
            if (baseline_subtraction_algorithm == "SNIP") {
                baseline_subtraction_iterations <- tclvalue(baseline_subtraction_iterations2)
                baseline_subtraction_iterations_value <- as.character(baseline_subtraction_iterations)
                baseline_subtraction_iterations <- as.integer(baseline_subtraction_iterations)
            }
            # Default
            if (baseline_subtraction_algorithm == "") {
                baseline_subtraction_algorithm <- "SNIP"
                baseline_subtraction_iterations <- 200
            }
        } else if (baseline_subtraction == "NO") {
            baseline_subtraction <- FALSE
            baseline_subtraction_algorithm <- NULL
        }
        # Set the value of the displaying label
        if (baseline_subtraction == TRUE && baseline_subtraction_algorithm != "SNIP") {
            baseline_subtraction_value <- paste("YES", "\n( ", baseline_subtraction_algorithm, " )\n", sep = "")
        } else if (baseline_subtraction == TRUE && baseline_subtraction_algorithm == "SNIP") {
            baseline_subtraction_value <- paste("YES", "\n( ", baseline_subtraction_algorithm, ",\niterations:", baseline_subtraction_iterations, " )", sep = "")
        } else {
            baseline_subtraction_value <- "NO"
        }
        baseline_subtraction_value_label <- tklabel(preproc_window, text = baseline_subtraction_value, font = label_font, bg = "white", width = 20, height = 3)
        tkgrid(baseline_subtraction_value_label, row = 5, column = 3, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$baseline_subtraction <- baseline_subtraction
        .GlobalEnv$baseline_subtraction_iterations <- baseline_subtraction_iterations
        .GlobalEnv$baseline_subtraction_algorithm <- baseline_subtraction_algorithm
        .GlobalEnv$baseline_subtraction_value <- baseline_subtraction_value
    }
    # Normalization
    normalization_choice <- function() {
        # Catch the value from the menu
        normalization <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "YES")
        # Default
        if (normalization == "YES" || normalization == "") {
            normalization <- TRUE
            # Ask for the algorithm
            normalization_algorithm <- select.list(c("TIC","PQN","median"), title = "Choose", multiple = FALSE, preselect = "TIC")
            # TIC
            if (normalization_algorithm == "TIC") {
                normalization_mass_range <- tclvalue(normalization_mass_range2)
                normalization_mass_range_value <- as.character(normalization_mass_range)
                if (normalization_mass_range != 0 && normalization_mass_range != "") {
                normalization_mass_range <- unlist(strsplit(normalization_mass_range, ","))
                normalization_mass_range <- as.numeric(normalization_mass_range)
                } else if (normalization_mass_range == 0 || normalization_mass_range == "") {
                    normalization_mass_range <- NULL
                }
            }
            # Default
            if (normalization_algorithm == "") {
                normalization_algorithm <- "TIC"
            }
        } else if (normalization == "NO") {
            normalization <- FALSE
            normalization_algorithm <- NULL
        }
        # Set the value of the displaying label
        if (normalization == TRUE && normalization_algorithm != "TIC") {
            normalization_value <- paste("YES", "\n( ", normalization_algorithm, " )\n", sep = "")
        } else if (normalization == TRUE && normalization_algorithm == "TIC") {
            if (!is.null(normalization_mass_range)) {
                normalization_value <- paste("YES", "\n( ", normalization_algorithm, ",\nrange:\n", normalization_mass_range_value, " )", sep = "")
            } else {
                normalization_value <- paste("YES", "\n( ", normalization_algorithm, " )", sep = "")
            }
        } else {
            normalization_value <- "NO"
        }
        normalization_value_label <- tklabel(preproc_window, text = normalization_value, font = label_font, bg = "white", width = 20, height = 4)
        tkgrid(normalization_value_label, row = 6, column = 3, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$normalization <- normalization
        .GlobalEnv$normalization_mass_range <- normalization_mass_range
        .GlobalEnv$normalization_algorithm <- normalization_algorithm
        .GlobalEnv$normalization_value <- normalization_value
    }
    # Spectral alignment
    spectral_alignment_choice <- function() {
        # Catch the value from the menu
        spectral_alignment <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "NO")
        # YES and Default
        if (spectral_alignment == "YES") {
            spectral_alignment <- TRUE
            # Ask for the algorithm
            spectral_alignment_algorithm <- select.list(c("cubic","quadratic","linear","lowess"), title = "Choose", multiple = FALSE, preselect = "cubic")
            # Default
            if (spectral_alignment_algorithm == "") {
                spectral_alignment_algorithm <- "cubic"
            }
            ## Ask for the reference peaklist
            spectral_alignment_reference <- select.list(c("auto","average spectrum", "skyline spectrum"), title = "Choose", multiple = FALSE, preselect = "auto")
            if (spectral_alignment_reference == "") {
                spectral_alignment_reference <- "auto"
            }
        } else if (spectral_alignment == "NO" || spectral_alignment == "") {
            spectral_alignment <- FALSE
            spectral_alignment_algorithm <- NULL
            spectral_alignment_reference <- NULL
        }
        # Set the value of the displaying label
        if (spectral_alignment == TRUE) {
            spectral_alignment_value <- paste("YES", "\n( ", spectral_alignment_algorithm, ",\n", spectral_alignment_reference, " )", sep = "")
        } else {
            spectral_alignment_value <- "NO"
        }
        spectral_alignment_value_label <- tklabel(preproc_window, text = spectral_alignment_value, font = label_font, bg = "white", width = 20, height = 3)
        tkgrid(spectral_alignment_value_label, row = 8, column = 2, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$spectral_alignment <- spectral_alignment
        .GlobalEnv$spectral_alignment_algorithm <- spectral_alignment_algorithm
        .GlobalEnv$spectral_alignment_reference <- spectral_alignment_reference
        .GlobalEnv$spectral_alignment_value <- spectral_alignment_value
    }
    # TOF mode
    tof_mode_choice <- function() {
        # Catch the value from the menu
        tof_mode <- select.list(c("Linear", "Reflector"), title = "Choose")
        # Default
        if (tof_mode == "" || tof_mode == "Linear") {
            tof_mode <- "linear"
        }
        if (tof_mode == "Reflector") {
            tof_mode <- "reflector"
        }
        # Set the value of the displaying label
        if (tof_mode == "linear") {
            tof_mode_value <- "Linear"
        } else if (tof_mode == "reflector") {
            tof_mode_value <- "Reflector"
        }
        tof_mode_value_label <- tklabel(preproc_window, text = tof_mode_value, font = label_font, bg = "white", width = 20)
        tkgrid(tof_mode_value_label, row = 2, column = 3, padx = c(10, 10), pady = c(10, 10))
        # Escape the function
        .GlobalEnv$tof_mode <- tof_mode
        .GlobalEnv$tof_mode_value <- tof_mode_value
    }
    # Commit preprocessing
    commit_preprocessing_function <- function() {
        # Get the values (they are filled with the default anyway)
        # Mass range
        mass_range <- tclvalue(mass_range2)
        mass_range_value <- as.character(mass_range)
        mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
        # Preprocessing
        preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of2)
        preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)
        preprocess_spectra_in_packages_of_value <- as.character(preprocess_spectra_in_packages_of)
        # Escape the function
        .GlobalEnv$mass_range <- mass_range
        .GlobalEnv$preprocess_spectra_in_packages_of <- preprocess_spectra_in_packages_of
        .GlobalEnv$preprocessing_parameters <- list(mass_range = mass_range, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of, spectral_alignment_algorithm = spectral_alignment_algorithm, spectral_alignment_reference = spectral_alignment_reference)
        # Destroy the window upon committing
        tkdestroy(preproc_window)
    }
    ##### List of variables, whose values are taken from the entries in the GUI (create new variables for the sub window, that will replace the ones in the global environment, only if the default are changed)
    mass_range2 <- tclVar("")
    preprocess_spectra_in_packages_of2 <- tclVar("")
    baseline_subtraction_iterations2 <- tclVar("")
    normalization_mass_range2 <- tclVar("")
    ##### Window
    preproc_window <- tktoplevel(bg = "white")
    tktitle(preproc_window) <- "Spectral preprocessing parameters"
    tkpack.propagate(preproc_window, FALSE)
    # Mass range
    mass_range_label <- tklabel(preproc_window, text = "Mass range", font = label_font, bg = "white", width = 20)
    mass_range_entry <- tkentry(preproc_window, width = 15, textvariable = mass_range2, font = entry_font, bg = "white", width = 20)
    tkinsert(mass_range_entry, "end", as.character(paste(mass_range[1],",",mass_range[2])))
    # Preprocessing (in packages of)
    preprocess_spectra_in_packages_of_label <- tklabel(preproc_window, text="Preprocess spectra\nin packages of", font = label_font, bg = "white", width = 20)
    preprocess_spectra_in_packages_of_entry <- tkentry(preproc_window, width = 10, textvariable = preprocess_spectra_in_packages_of2, font = entry_font, bg = "white", width = 20)
    tkinsert(preprocess_spectra_in_packages_of_entry, "end", as.character(preprocess_spectra_in_packages_of))
    # Tof mode
    tof_mode_label <- tklabel(preproc_window, text="Select the TOF mode", font = label_font, bg = "white", width = 20)
    tof_mode_entry <- tkbutton(preproc_window, text="Choose the TOF mode", command = tof_mode_choice, font = button_font, bg = "white", width = 20)
    # Transform the data
    transform_data_button <- tkbutton(preproc_window, text="Transform the data", command = transform_data_choice, font = button_font, bg = "white", width = 20)
    # Smoothing
    smoothing_button <- tkbutton(preproc_window, text="Smoothing", command = smoothing_choice, font = button_font, bg = "white", width = 20)
    # Baseline subtraction
    baseline_subtraction_button <- tkbutton(preproc_window, text="Baseline subtraction", command = baseline_subtraction_choice, font = button_font, bg = "white", width = 20)
    baseline_subtraction_iterations_entry <- tkentry(preproc_window, width = 15, textvariable = baseline_subtraction_iterations2, font = entry_font, bg = "white", width = 20)
    tkinsert(baseline_subtraction_iterations_entry, "end", as.character(baseline_subtraction_iterations))
    # Normalization
    normalization_button <- tkbutton(preproc_window, text="Normalization", command = normalization_choice, font = button_font, bg = "white", width = 20)
    normalization_mass_range_entry <- tkentry(preproc_window, width = 15, textvariable = normalization_mass_range2, font = entry_font, bg = "white", width = 20)
    tkinsert(normalization_mass_range_entry, "end", as.character(normalization_mass_range))
    # Spectral alignment
    spectral_alignment_button <- tkbutton(preproc_window, text="Align spectra", command = spectral_alignment_choice, font = button_font, bg = "white", width = 20)
    # Commit preprocessing
    commit_preprocessing_button <- tkbutton(preproc_window, text="Commit preprocessing", command = commit_preprocessing_function, font = button_font, bg = "white", width = 20)
    ##### Displaying labels
    tof_mode_value_label <- tklabel(preproc_window, text = tof_mode_value, font = label_font, bg = "white", width = 20)
    transform_data_value_label <- tklabel(preproc_window, text = transform_data_value, font = label_font, bg = "white", width = 20, height = 2)
    smoothing_value_label <- tklabel(preproc_window, text = smoothing_value, font = label_font, bg = "white", width = 20, height = 3)
    baseline_subtraction_value_label <- tklabel(preproc_window, text = baseline_subtraction_value, font = label_font, bg = "white", width = 20, height = 3)
    normalization_value_label <- tklabel(preproc_window, text = normalization_value, font = label_font, bg = "white", width = 20, height = 4)
    spectral_alignment_value_label <- tklabel(preproc_window, text = spectral_alignment_value, font = label_font, bg = "white", width = 20, height = 3)
    #### Geometry manager
    tkgrid(mass_range_label, row = 1, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(mass_range_entry, row = 1, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(tof_mode_label, row = 2, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(tof_mode_entry, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(tof_mode_value_label, row = 2, column = 3, padx = c(10, 10), pady = c(10, 10))
    tkgrid(transform_data_button, row = 3, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(transform_data_value_label, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(smoothing_button, row = 4, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(smoothing_value_label, row = 4, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(baseline_subtraction_button, row = 5, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(baseline_subtraction_iterations_entry, row = 5, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(baseline_subtraction_value_label, row = 5, column = 3, padx = c(10, 10), pady = c(10, 10))
    tkgrid(normalization_button, row = 6, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(normalization_mass_range_entry, row = 6, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(normalization_value_label, row = 6, column = 3, padx = c(10, 10), pady = c(10, 10))
    tkgrid(preprocess_spectra_in_packages_of_label, row = 7, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(preprocess_spectra_in_packages_of_entry, row = 7, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(spectral_alignment_button, row = 8, column = 1, padx = c(10, 10), pady = c(10, 10))
    tkgrid(spectral_alignment_value_label, row = 8, column = 2, padx = c(10, 10), pady = c(10, 10))
    tkgrid(commit_preprocessing_button, row = 9, column = 1, padx = c(10, 10), pady = c(10, 10))
}

##### File type export MATRIX
file_type_export_matrix_choice <- function() {
    # Catch the value from the menu
    file_type_export_matrix <- select.list(c("csv","xlsx","xls"), title = "Choose", multiple = FALSE, preselect = "csv")
    # Default
    if (file_type_export_matrix == "") {
        file_type_export_matrix <- "csv"
    }
    if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
        install_and_load_required_packages("XLConnect")
    }
    # Escape the function
    .GlobalEnv$file_type_export_matrix <- file_type_export_matrix
    # Set the value of the displaying label
    file_type_export_matrix_value_label <- tklabel(window, text = file_type_export_matrix, font = label_font, bg = "white", width = 20)
    tkgrid(file_type_export_matrix_value_label, row = 6, column = 2)
}

##### File type export IMAGES
file_type_export_images_choice <- function() {
    # Catch the value from the menu
    file_type_export_images <- select.list(c("png","jpg","tiff"), title = "Choose", multiple = FALSE, preselect = "png")
    # Default
    if (file_type_export_images == "") {
        file_type_export_images <- "png"
    }
    # Escape the function
    .GlobalEnv$file_type_export_images <- file_type_export_images
    # Set the value of the displaying label
    file_type_export_images_value_label <- tklabel(window, text = file_type_export_images, font = label_font, bg = "white", width = 20)
    tkgrid(file_type_export_images_value_label, row = 6, column = 4)
}

##### Classification mode choice
classification_mode_choice <- function() {
    # Catch the value from the menu
    classification_mode <- select.list(c("pixel", "profile"), title = "Choose", multiple = TRUE, preselect = "pixel")
    # Default
    if (length(classification_mode) == 1 && classification_mode == "") {
        classification_mode <- "pixel"
        classification_mode_value <- "pixel"
    }
    # Choose the pixel grouping
    if ("pixel" %in% classification_mode) {
        ## Choose the pixel grouping
        pixel_grouping <- select.list(c("single", "hca", "moving window average", "graph"), title = "Choose pixel grouping", preselect = "single")
        if (pixel_grouping == "") {
            pixel_grouping <- "single"
        }
        # Escape the function
        .GlobalEnv$pixel_grouping <- pixel_grouping
        if (pixel_grouping == "moving window average") {
            moving_window_size <- as.integer(select.list(c(2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 200, 300, 500, 1000), title = "Choose moving window size", preselect = 5))
            if (moving_window_size == "") {
                moving_window_size <- 5
            }
            # Escape the function
            .GlobalEnv$moving_window_size <- moving_window_size
        } else if (pixel_grouping == "hca") {
            number_of_hca_nodes <- as.integer(select.list(c(2, 3, 4, 5, 10, 15, 20, 25, 50, 100, 200, 300, 500, 1000), title = "Choose number of HC nodes", preselect = 4))
            if (number_of_hca_nodes == "") {
                number_of_hca_nodes <- 4
            }
            # Escape the function
            .GlobalEnv$number_of_hca_nodes <- number_of_hca_nodes
        }
    }
    # Fix the displaying value
    if ("profile" %in% classification_mode && "pixel" %in% classification_mode) {
        classification_mode_value <- "pixel + profile"
    } else if (classification_mode == "profile") {
        classification_mode_value <- "   profile   "
    } else if (classification_mode == "pixel") {
        classification_mode_value <- "  pixel  "
    }
    # Escape the function
    .GlobalEnv$classification_mode <- classification_mode
    .GlobalEnv$classification_mode_value <- classification_mode_value
    # Set the value of the displaying label
    classification_mode_value_label <- tklabel(window, text = classification_mode_value, font = label_font, bg = "white", width = 20)
    tkgrid(classification_mode_value_label, row = 5, column = 2)
}

##### Samples
select_samples_function <- function() {
    # Set the working directory
    setwd(getwd())
    ########## Prompt if a folder has to be selected or a single file
    # Catch the value from the popping out menu
    spectra_input_type <- select.list(c("file","folder"), title = "Folder or file?", multiple = FALSE, preselect = "file")
    # Default
    if (spectra_input_type == "") {
        spectra_input_type <- "file"
    }
    # Folder
    if (spectra_input_type == "folder") {
        filepath_import_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra (stored as imzML files) to be imported: each imzML file should correspond to a patient's spectral dataset, and the imzML files should be put in folders corresponding to the different classes of samples", icon = "info")
        filepath_import <- tclvalue(tkchooseDirectory())
        if (!nchar(filepath_import)) {
            tkmessageBox(message = "No folder selected")
            filepath_import <- NULL
        } else {
            tkmessageBox(message = paste("The sample spectra will be read from:\n\n", filepath_import))
        }
    } else if (spectra_input_type == "file") {
        # File
        filepath_import_select <- tkmessageBox(title = "Sample file", message = "Select the imzML file for the spectra to be imported: the imzML file should correspond to a patient's spectral dataset", icon = "info")
        filepath_import <- tclvalue(tkgetOpenFile(filetypes = "{{imzML files} {.imzML .imzml}}"))
        if (!nchar(filepath_import)) {
            tkmessageBox(message = "No file selected")
            filepath_import <- NULL
        } else {
            tkmessageBox(message = paste("The sample spectra will be read from:\n\n", filepath_import))
        }
    }
    # Set the value for displaying purposes
    filepath_import_value <- filepath_import
    # Exit the function and put the variable into the R workspace
    .GlobalEnv$filepath_import <- filepath_import
}

##### RData file containing the models
select_RData_file_function <- function() {
    # Set the working directory
    setwd(getwd())
    filepath_R_import_select <- tkmessageBox(title = "RData file", message = "Select the RData file containing the model list", icon = "info")
    filepath_R <- tclvalue(tkgetOpenFile(filetypes = "{{RData files} {.RData}}"))
    if (!nchar(filepath_R)) {
        tkmessageBox(message = "No RData file selected")
        filepath_R <- NULL
        RData_file_integrity <- FALSE
        RData_file_integrity_value <- "NO RDATA FILE\nSELECTED"
    } else {
        tkmessageBox(message = paste("The model list will be read from:\n\n", filepath_R))
        # Define the progressbar
        RData_progress_bar <- tkProgressBar(title = "RData integrity verification", label = "", min = 0, max = 1, initial = 0, width = 300)
        # Retrieve the content of the RData workspace
        setTkProgressBar(RData_progress_bar, value = 0, title = NULL, label = "0 %")
        RData_file_integrity <- FALSE
        RData_file_integrity_value <- "INTEGRITY TEST\nFAILED"
        try(RData_variables <- R_workspace_data_retriever(filepath_R)$variable_list)
        setTkProgressBar(RData_progress_bar, value = 0.50, title = NULL, label = "50 %")
        # Check for integrity of the RData file
        try({
            if ("model_list" %in% RData_variables) {
                RData_file_integrity <- TRUE
                RData_file_integrity_value <- "INTEGRITY TEST\nSUCCEDED"
            }
        })
        setTkProgressBar(RData_progress_bar, value = 1, title = NULL, label = "100 %")
        close(RData_progress_bar)
    }
    RData_file_integrity_value_label <- tklabel(window, text = RData_file_integrity_value, font = label_font, bg = "white", width = 20)
    tkgrid(RData_file_integrity_value_label, row = 5, column = 4)
    # Exit the function and put the variable into the R workspace
    .GlobalEnv$filepath_R <- filepath_R
    .GlobalEnv$RData_file_integrity <- RData_file_integrity
    .GlobalEnv$RData_file_integrity_value <- RData_file_integrity_value
}

##### Output
browse_output_function <- function() {
    output_folder <- tclvalue(tkchooseDirectory())
    if (!nchar(output_folder)) {
        # Get the output folder from the default working directory
        output_folder <- getwd()
    }
    tkmessageBox(message = paste("Every file will be saved in:\n\n", output_folder))
    setwd(output_folder)
    # Exit the function and put the variable into the R workspace
    .GlobalEnv$output_folder <- output_folder
}

##### Exit
end_session_function <- function () {
    q(save = "no")
}

##### Peak picking algorithm
peak_picking_algorithm_choice <- function() {
    # Catch the value from the menu
    peak_picking_algorithm <- select.list(c("MAD", "SuperSmoother"), title = "Choose", multiple = FALSE, preselect = "SuperSmoother")
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
    peak_picking_algorithm_value_label <- tklabel(window, text = peak_picking_algorithm_value, font = label_font, bg = "white", width = 20)
    tkgrid(peak_picking_algorithm_value_label, row = 2, column = 2)
    # Escape the function
    .GlobalEnv$peak_picking_algorithm <- peak_picking_algorithm
    .GlobalEnv$peak_picking_algorithm_value <- peak_picking_algorithm_value
}

##### Peaks deisotoping
peaks_deisotoping_choice <- function() {
    # Catch the value from the menu
    peaks_deisotoping <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "NO")
    # Default
    if (peaks_deisotoping == "YES") {
        peaks_deisotoping <- TRUE
    }
    if (peaks_deisotoping == "NO" || peaks_deisotoping == "") {
        peaks_deisotoping <- FALSE
    }
    # Set the value of the displaying label
    if (peaks_deisotoping == TRUE) {
        peaks_deisotoping_value <- "YES"
    } else {
        peaks_deisotoping_value <- "   NO   "
    }
    peaks_deisotoping_value_label <- tklabel(window, text = peaks_deisotoping_value, font = label_font, bg = "white", width = 20)
    tkgrid(peaks_deisotoping_value_label, row = 2, column = 4)
    # Escape the function
    .GlobalEnv$peaks_deisotoping <- peaks_deisotoping
    .GlobalEnv$peaks_deisotoping_value <- peaks_deisotoping_value
}

##### Decision method ensemble
decision_method_ensemble_choice <- function() {
    # Catch the value from the menu
    decision_method_ensemble <- select.list(c("majority"), title = "Choose", multiple = FALSE, preselect = "majority")
    # Default
    if (decision_method_ensemble == "") {
        decision_method_ensemble <- "majority"
    }
    decision_method_ensemble_value_label <- tklabel(window, text = decision_method_ensemble, font = label_font, bg = "white", width = 20)
    tkgrid(decision_method_ensemble_value_label, row = 3, column = 2)
    # Escape the function
    .GlobalEnv$decision_method_ensemble <- decision_method_ensemble
}

##### Vote weights ensemble
vote_weights_ensemble_choice <- function() {
    # Catch the value from the menu
    vote_weights_ensemble <- select.list(c("equal"), title = "Choose", multiple = FALSE, preselect = "equal")
    # Default
    if (vote_weights_ensemble == "") {
        vote_weights_ensemble <- "equal"
    }
    vote_weights_ensemble_value_label <- tklabel(window, text = vote_weights_ensemble, font = label_font, bg = "white", width = 20)
    tkgrid(vote_weights_ensemble_value_label, row = 3, column = 4)
    # Escape the function
    .GlobalEnv$vote_weights_ensemble <- vote_weights_ensemble
}

##### Multicore processing
allow_parallelization_choice <- function() {
    ##### Messagebox
    tkmessageBox(title = "Parallel processing is resource hungry", message = "Parallel processing is resource hungry.\nBy activating it, the computation becomes faster, but the program will eat a lot of RAM, possibly causing your computer to freeze. If you want to play safe, do not enable it.", icon = "warning")
    # Catch the value from the menu
    allow_parallelization <- select.list(c("YES","NO"), title = "Choose", multiple = FALSE, preselect = "NO")
    # Default
    if (allow_parallelization == "YES") {
        allow_parallelization <- TRUE
    }
    if (allow_parallelization == "NO" || allow_parallelization == "") {
        allow_parallelization <- FALSE
    }
    # Set the value of the displaying label
    if (allow_parallelization == TRUE) {
        allow_parallelization_value <- "YES"
    } else {
        allow_parallelization_value <- "  NO  "
    }
    allow_parallelization_value_label <- tklabel(window, text = allow_parallelization_value, font = label_font, bg = "white", width = 20)
    tkgrid(allow_parallelization_value_label, row = 4, column = 2)
    # Escape the function
    .GlobalEnv$allow_parallelization <- allow_parallelization
    .GlobalEnv$allow_parallelization_value <- allow_parallelization_value
}

##### Run the Peaklist Export function
run_patient_classification_function <- function() {
    ######## Run only if all the elements needed are there
    if (!is.null(filepath_import) && RData_file_integrity == TRUE) {
        # Progress bar
        program_progress_bar <- tkProgressBar(title = "Calculating...", label = "", min = 0, max = 1, initial = 0, width = 300)
        setTkProgressBar(program_progress_bar, value = 0, title = NULL, label = "0 %")
        # Choose the legends to plot
        plot_legends <- select.list(c("sample name", "legend", "plot name"), multiple = TRUE, preselect = "legend", title = "Legend to plot")
        setTkProgressBar(program_progress_bar, value = 0.25, title = NULL, label = "25 %")
        ########## Run the classification function
        classification_of_patients <- spectral_classification(spectra_path = filepath_import, filepath_R = filepath_R, model_list_object = "model_list", classification_mode = classification_mode, peak_picking_algorithm = peak_picking_algorithm, deisotope_peaklist = peaks_deisotoping, preprocessing_parameters = list(mass_range = mass_range, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_iterations = baseline_subtraction_iterations, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = spectral_alignment_algorithm, spectral_alignment_reference = spectral_alignment_reference, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), tof_mode = tof_mode, allow_parallelization = allow_parallelization, decision_method_ensemble = decision_method_ensemble, vote_weights_ensemble = vote_weights_ensemble, pixel_grouping = pixel_grouping, moving_window_size = moving_window_size, number_of_hca_nodes = number_of_hca_nodes, number_of_spectra_partitions_graph = 1, partitioning_method_graph = "space", correlation_method_for_adjacency_matrix = "pearson", correlation_threshold_for_adjacency_matrix = 0.95, pvalue_threshold_for_adjacency_matrix = 0.05, max_GA_generations = 50, iterations_with_no_change_GA = 5, seed = 12345, plot_figures = TRUE, plot_graphs = TRUE, plot_legends = plot_legends)
        setTkProgressBar(program_progress_bar, value = 0.75, title = NULL, label = "75 %")
        # Escape the function
        .GlobalEnv$classification_of_patients <- classification_of_patients
        setTkProgressBar(program_progress_bar, value = 1.00, title = NULL, label = "100 %")
        close(program_progress_bar)
        ### Messagebox
        tkmessageBox(title = "Done!", message = "The classification has been performed!", icon = "info")
    } else if (is.null(filepath_import) || RData_file_integrity == FALSE) {
        classification_of_patients <- NULL
        # Escape the function
        .GlobalEnv$classification_of_patients <- classification_of_patients
        ### Messagebox
        tkmessageBox(title = "Something is wrong", message = "Either no proper RData files or spectra files are provided!", icon = "warning")
    }
}

##### Dump the output files
ms_pixel_typer_data_dumper_function <- function() {
    if (!is.null(classification_of_patients)) {
        ########## Dump the files
        ##### Create a subfolder (CLASSIFICATION X)
        ### Go to the working directory
        setwd(output_folder)
        ##### Automatically create a subfolder with all the results
        ## Check if such subfolder exists
        list_of_directories <- list.dirs(output_folder, full.names = FALSE, recursive = FALSE)
        ## Check the presence of a CLASSIFICATION folder
        CLASSIFICATION_folder_presence <- FALSE
        if (length(list_of_directories) > 0) {
            for (dr in 1:length(list_of_directories)) {
                if (length(grep("CLASSIFICATION", list_of_directories[dr], fixed = TRUE)) > 0) {
                    CLASSIFICATION_folder_presence <- TRUE
                }
            }
        }
        ## If it present...
        if (isTRUE(CLASSIFICATION_folder_presence)) {
            ## Extract the number after the STATISTICS
            # Number for the newly created folder
            CLASSIFICATION_new_folder_number <- 0
            # List of already present numbers
            CLASSIFICATION_present_folder_numbers <- integer()
            # For each folder present...
            for (dr in 1:length(list_of_directories)) {
                # If it is named CLASSIFICATION
                if (length(grep("CLASSIFICATION", list_of_directories[dr], fixed = TRUE)) > 0) {
                    # Split the name
                    CLASSIFICATION_present_folder_split <- unlist(strsplit(list_of_directories[dr], "CLASSIFICATION"))
                    # Add the number to the list of CLASSIFICATION numbers
                    try(CLASSIFICATION_present_folder_numbers <- append(CLASSIFICATION_present_folder_numbers, as.integer(CLASSIFICATION_present_folder_split[2])))
                }
            }
            # Sort the STATISTICS folder numbers
            try(CLASSIFICATION_present_folder_numbers <- sort(CLASSIFICATION_present_folder_numbers))
            # The new folder number will be the greater + 1
            try(CLASSIFICATION_new_folder_number <- CLASSIFICATION_present_folder_numbers[length(CLASSIFICATION_present_folder_numbers)] + 1)
            # Generate the new subfolder
            subfolder <- file.path(output_folder, paste("CLASSIFICATION", CLASSIFICATION_new_folder_number))
            # Create the subfolder
            dir.create(subfolder)
        } else {
            # If it not present...
            # Create the folder where to dump the files and go to it...
            subfolder <- file.path(output_folder, paste("CLASSIFICATION", "1"))
            # Create the subfolder
            dir.create(subfolder)
        }
        # Go to the new working directory
        setwd(subfolder)
        ##### Dump the files
        # Determine the number of patients, to create subfolders (try with the msi list, if it is still zero because only the profile has been classified, use the profile list)
        number_of_patients <- length(classification_of_patients$final_result_matrix_msi_list)
        if (number_of_patients == 0) {
            number_of_patients <- length(classification_of_patients$final_result_matrix_profile_list)
        }
        ##### For each patient...
        for (p in 1:number_of_patients) {
            # Retrieve the patient name
            if (length(classification_of_patients$final_result_matrix_msi_list) > 0) {
                patient_name <- names(classification_of_patients$final_result_matrix_msi_list)[[p]]
            } else {
                patient_name <- NULL
            }
            if (is.null(patient_name)) {
                if (length(classification_of_patients$final_result_matrix_profile_list) > 0) {
                    patient_name <- names(classification_of_patients$final_result_matrix_profile_list)[[p]]
                } else {
                    patient_name <- NULL
                }
            }
            if (is.null(patient_name)) {
                patient_name <- p
            }
            # Create the subfolder and go to it
            subfolder_patient <- file.path(subfolder, patient_name)
            dir.create(subfolder_patient)
            setwd(subfolder_patient)
            ### Dump the files
            ## MSI classification matrix
            if (file_type_export_matrix == "csv") {
                try(write.csv(classification_of_patients$final_result_matrix_msi_list[[p]], file = paste("MSI classification matrix", ".", file_type_export_matrix, sep = "")), silent = TRUE)
            } else if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
                try({
                    wb = loadWorkbook(filename = paste("MSI classification matrix", ".", file_type_export_matrix, sep = ""), create = TRUE)
                    createSheet(wb, name = "MSI classification")
                    writeWorksheet(wb, data = classification_of_patients$final_result_matrix_msi_list[[p]], sheet = "MSI classification")
                    saveWorkbook(wb)
                }, silent = TRUE)
            }
            ## MS profile matrix
            if (file_type_export_matrix == "csv") {
                try(write.csv(classification_of_patients$final_result_matrix_profile_list[[p]], file = paste("MS profile classification matrix", ".", file_type_export_matrix, sep = "")), silent = TRUE)
            } else if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
                try({
                    wb = loadWorkbook(filename = paste("MS profile classification matrix", ".", file_type_export_matrix, sep = ""), create = TRUE)
                    createSheet(wb, name = "MS profile classification")
                    writeWorksheet(wb, data = classification_of_patients$final_result_matrix_profile_list[[p]], sheet = "MS profile classification")
                    saveWorkbook(wb)
                }, silent = TRUE)
            }
            ## MSI pixel classification
            # PNG
            if (file_type_export_images == "png") {
                try({
                    for (i in 1:length(classification_of_patients$classification_ms_images_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$classification_ms_images_list[[p]])[i]
                        # Save the plot
                        png(filename = paste("Pixel-by-pixel classification ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                        replayPlot(classification_of_patients$classification_ms_images_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            } else if (file_type_export_images == "tiff") {
                # TIFF
                try({
                    for (i in 1:length(classification_of_patients$classification_ms_images_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$classification_ms_images_list[[p]])[i]
                        # Save the plot
                        tiff(filename = paste("Pixel-by-pixel classification ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                        replayPlot(classification_of_patients$classification_ms_images_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            } else if (file_type_export_images == "jpg" || file_type_export_images == "jpeg") {
                # JPEG
                try({
                    for (i in 1:length(classification_of_patients$classification_ms_images_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$classification_ms_images_list[[p]])[i]
                        # Save the plot
                        jpeg(filename = paste("Pixel-by-pixel classification ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150, quality = 100)
                        replayPlot(classification_of_patients$classification_ms_images_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            }
            ## Ensemble classification MSI matrix
            if (file_type_export_matrix == "csv") {
                try(write.csv(classification_of_patients$classification_ensemble_matrix_msi_all[[p]], file = paste("Ensemble MSI classification matrix", ".", file_type_export_matrix, sep = "")), silent = TRUE)
            } else if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
                try({
                    wb = loadWorkbook(filename = paste("Ensemble MSI classification matrix", ".", file_type_export_matrix, sep = ""), create = TRUE)
                    createSheet(wb, name = "Ensemble MSI classification")
                    writeWorksheet(wb, data = classification_of_patients$classification_ensemble_matrix_msi_all[[p]], sheet = "Ensemble MSI classification")
                    saveWorkbook(wb)
                }, silent = TRUE)
            }
            ## Ensemble classification profile matrix
            if (file_type_export_matrix == "csv") {
                try(write.csv(classification_of_patients$classification_ensemble_matrix_profile_all[[p]], file = paste("Ensemble MS profile classification matrix", ".", file_type_export_matrix, sep = "")), silent = TRUE)
            } else if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
                try({
                    wb = loadWorkbook(filename = paste("Ensemble MS profile classification matrix", ".", file_type_export_matrix, sep = ""), create = TRUE)
                    createSheet(wb, name = "Ensemble MS classification")
                    writeWorksheet(wb, data = classification_of_patients$classification_ensemble_matrix_profile_all[[p]], sheet = "Ensemble MS classification")
                    saveWorkbook(wb)
                }, silent = TRUE)
            }
            ## Ensemble MSI pixel classification
            # PNG
            if (file_type_export_images == "png") {
                try({
                    png(filename = paste("Ensemble Pixel-by-pixel classification", ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                    replayPlot(classification_of_patients$classification_ensemble_ms_image_list[[p]])
                    dev.off()
                }, silent = TRUE)
            } else if (file_type_export_images == "tiff") {
                # TIFF
                try({
                    tiff(filename = paste("Ensemble Pixel-by-pixel classification", ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                    replayPlot(classification_of_patients$classification_ensemble_ms_image_list[[p]])
                    dev.off()
                }, silent = TRUE)
            } else if (file_type_export_images == "jpg" || file_type_export_images == "jpeg") {
                # JPEG
                try({
                    jpeg(filename = paste("Ensemble Pixel-by-pixel classification", ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150, quality = 100)
                    replayPlot(classification_of_patients$classification_ensemble_ms_image_list[[p]])
                    dev.off()
                }, silent = TRUE)
            }
            ## Average spectrum with bars
            # PNG
            if (file_type_export_images == "png") {
                try({
                    for (i in 1:length(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])[i]
                        # Save the plot
                        png(filename = paste("Average spectrum with bars ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                        replayPlot(classification_of_patients$average_spectrum_with_bars_profile_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            } else if (file_type_export_images == "tiff") {
                # TIFF
                try({
                    for (i in 1:length(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])[i]
                        # Save the plot
                        tiff(filename = paste("Average spectrum with bars ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                        replayPlot(classification_of_patients$average_spectrum_with_bars_profile_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            } else if (file_type_export_images == "jpg" || file_type_export_images == "jpeg") {
                # JPEG
                try({
                    for (i in 1:length(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])) {
                        # Retrieve the model name
                        model_name <- names(classification_of_patients$average_spectrum_with_bars_profile_list[[p]])[i]
                        # Save the plot
                        jpeg(filename = paste("Average spectrum with bars ", model_name, ".", file_type_export_images, sep = ""), width = 1920, height = 1080, pointsize = 20, units = "px", res = 150)
                        replayPlot(classification_of_patients$average_spectrum_with_bars_profile_list[[p]][[i]])
                        dev.off()
                    }
                }, silent = TRUE)
            }
        }
        ## Ensemble classification profile matrix ALL
        setwd(subfolder)
        if (file_type_export_matrix == "csv") {
            try(write.csv(classification_of_patients$classification_ensemble_matrix_profile_all, file = paste("Ensemble MS profile classification matrix ALL", ".", file_type_export_matrix, sep = "")), silent = TRUE)
        } else if (file_type_export_matrix == "xls" || file_type_export_matrix == "xlsx") {
            try({
                wb = loadWorkbook(filename = paste("Ensemble MS profile classification matrix ALL", ".", file_type_export_matrix, sep = ""), create = TRUE)
                createSheet(wb, name = "Ensemble MS classification")
                writeWorksheet(wb, data = classification_of_patients$classification_ensemble_matrix_profile_all, sheet = "Ensemble MS classification", header = TRUE, rownames = rownames(classification_of_patients$classification_ensemble_matrix_profile_all))
                saveWorkbook(wb)
            }, silent = TRUE)
        }
        ### Messagebox
        tkmessageBox(title = "Files dumped!", message = "The classification files have been dumped!", icon = "info")
    } else {
        ### Messagebox
        tkmessageBox(title = "No classification found", message = "No classification files have been found!\nRun the classification before dumping the files!", icon = "warning")
    }
}














##################################################################### WINDOW GUI

### Check for updates
check_for_updates_function()



######################## GUI

# Get system info (Platform - Release - Version (- Linux Distro))
system_os = Sys.info()[1]
os_release = Sys.info()[2]
os_version = Sys.info()[3]

### Get the screen resolution
# Windows
if (system_os == "Windows") {
    # Windows 7
    if (length(grep("7", os_release, fixed = TRUE)) > 0) {
        # Get system info
        screen_height <- system("wmic desktopmonitor get screenheight", intern = TRUE)
        screen_width <- system("wmic desktopmonitor get screenwidth", intern = TRUE)
        # Retrieve the values
        screen_height <- as.numeric(screen_height[-c(1, length(screen_height))])
        screen_width <- as.numeric(screen_width[-c(1, length(screen_width))])
    } else if (length(grep("10", os_release, fixed = TRUE)) > 0) {
        # Windows 10
        # Get system info
        screen_info <- system("wmic path Win32_VideoController get VideoModeDescription", intern = TRUE)[2]
        # Get the resolution
        screen_resolution <- unlist(strsplit(screen_info, "x"))
        # Retrieve the values
        screen_height <- as.numeric(screen_resolution[2])
        screen_width <- as.numeric(screen_resolution[1])
    }
} else if (system_os == "Linux") {
    # Get system info
    screen_info <- system("xdpyinfo -display :0", intern = TRUE)
    # Get the resolution
    screen_resolution <- screen_info[which(screen_info == "screen #0:") + 1]
    screen_resolution <- unlist(strsplit(screen_resolution, "dimensions: ")[1])
    screen_resolution <- unlist(strsplit(screen_resolution, "pixels"))[2]
    # Retrieve the wto dimensions...
    screen_width <- as.numeric(unlist(strsplit(screen_resolution, "x"))[1])
    screen_height <- as.numeric(unlist(strsplit(screen_resolution, "x"))[2])
}


### FONTS
# Default sizes (determined on a 1680x1050 screen) (in order to make them adjust to the size screen, the screen resolution should be retrieved)
title_font_size <- 24
other_font_size <- 11
# Windows
if (system_os == "Windows") {
    # Windows 7
    if (length(grep("7", os_release, fixed = TRUE)) > 0) {
        # Determine the font size according to the resolution
        total_number_of_pixels <- screen_width * screen_height
        # Determine the scaling factor (according to a complex formula)
        scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
        scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
        title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
        other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    } else if (length(grep("10", os_release, fixed = TRUE)) > 0) {
        # Windows 10
        # Determine the font size according to the resolution
        total_number_of_pixels <- screen_width * screen_height
        # Determine the scaling factor (according to a complex formula)
        scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
        scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
        title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
        other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    }
    # Define the fonts
    garamond_title_bold = tkfont.create(family = "Garamond", size = title_font_size, weight = "bold")
    garamond_other_normal = tkfont.create(family = "Garamond", size = other_font_size, weight = "normal")
    arial_title_bold = tkfont.create(family = "Arial", size = title_font_size, weight = "bold")
    arial_other_normal = tkfont.create(family = "Arial", size = other_font_size, weight = "normal")
    trebuchet_title_bold = tkfont.create(family = "Trebuchet MS", size = title_font_size, weight = "bold")
    trebuchet_other_normal = tkfont.create(family = "Trebuchet MS", size = other_font_size, weight = "normal")
    trebuchet_other_bold = tkfont.create(family = "Trebuchet MS", size = other_font_size, weight = "bold")
    # Use them in the GUI
    title_font = trebuchet_title_bold
    label_font = trebuchet_other_normal
    entry_font = trebuchet_other_normal
    button_font = trebuchet_other_bold
} else if (system_os == "Linux") {
    # Linux
    # Determine the font size according to the resolution
    total_number_of_pixels <- screen_width * screen_height
    # Determine the scaling factor (according to a complex formula)
    scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
    scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
    title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
    other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    # Ubuntu
    if (length(grep("Ubuntu", os_version, ignore.case = TRUE)) > 0) {
        # Define the fonts
        ubuntu_title_bold = tkfont.create(family = "Ubuntu", size = title_font_size, weight = "bold")
        ubuntu_other_normal = tkfont.create(family = "Ubuntu", size = other_font_size, weight = "normal")
        ubuntu_other_bold = tkfont.create(family = "Ubuntu", size = other_font_size, weight = "bold")
        liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
        liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
        liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
        bitstream_charter_title_bold = tkfont.create(family = "Bitstream Charter", size = title_font_size, weight = "bold")
        bitstream_charter_other_normal = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "normal")
        bitstream_charter_other_bold = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "bold")
        # Use them in the GUI
        title_font = bitstream_charter_title_bold
        label_font = bitstream_charter_other_normal
        entry_font = bitstream_charter_other_normal
        button_font = bitstream_charter_other_bold
    } else if (length(grep("Fedora", os_version, ignore.case = TRUE)) > 0) {
        # Fedora
        # Define the fonts
        cantarell_title_bold = tkfont.create(family = "Cantarell", size = title_font_size, weight = "bold")
        cantarell_other_normal = tkfont.create(family = "Cantarell", size = other_font_size, weight = "normal")
        cantarell_other_bold = tkfont.create(family = "Cantarell", size = other_font_size, weight = "bold")
        liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
        liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
        liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
        # Use them in the GUI
        title_font = cantarell_title_bold
        label_font = cantarell_other_normal
        entry_font = cantarell_other_normal
        button_font = cantarell_other_bold
    } else {
        # Other linux distros
        # Define the fonts
        liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
        liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
        liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
        # Use them in the GUI
        title_font = liberation_title_bold
        label_font = liberation_other_normal
        entry_font = liberation_other_normal
        button_font = liberation_other_bold
    }
} else if (system_os == "Darwin") {
    # macOS
    # Define the fonts
    helvetica_title_bold = tkfont.create(family = "Helvetica", size = title_font_size, weight = "bold")
    helvetica_other_normal = tkfont.create(family = "Helvetica", size = other_font_size, weight = "normal")
    helvetica_other_bold = tkfont.create(family = "Helvetica", size = other_font_size, weight = "bold")
    # Use them in the GUI
    title_font = helvetica_title_bold
    label_font = helvetica_other_normal
    entry_font = helvetica_other_normal
    button_font = helvetica_other_bold
}



# The "area" where we will put our input lines
window <- tktoplevel(bg = "white")
tkpack.propagate(window, FALSE)
tktitle(window) <- "MS PIXEL TYPER"
#### Browse
# Title label
title_label <- tklabel(window, text = "MS PIXEL TYPER", font = title_font, bg = "white")
# Library
select_samples_button <- tkbutton(window, text = "BROWSE\nSPECTRA...", command = select_samples_function, font = button_font, bg = "white", width = 20)
# Output
browse_output_button <- tkbutton(window, text = "BROWSE\nOUTPUT FOLDER...", command = browse_output_function, font = button_font, bg = "white", width = 20)
#### Entries
# Peak picking method
peak_picking_algorithm_entry <- tkbutton(window, text = "PEAK PICKING\nALGORITHM", command = peak_picking_algorithm_choice, font = button_font, bg = "white", width = 20)
# Peaks deisotoping
peaks_deisotoping_entry <- tkbutton(window, text = "PEAK\nDEISOTOPING", command = peaks_deisotoping_choice, font = button_font, bg = "white", width = 20)
# Decision method ensemble
decision_method_ensemble_entry <- tkbutton(window, text = "DECISION METHOD\nENSEMBLE", command = decision_method_ensemble_choice, font = button_font, bg = "white", width = 20)
# Vote weights ensemble
vote_weights_ensemble_entry <- tkbutton(window, text = "VOTE WEIGHTS\nENSEMBLE", command = vote_weights_ensemble_choice, font = button_font, bg = "white", width = 20)
# Classification mode
classification_mode_entry <- tkbutton(window, text = "CLASSIFICATION\nMODE", command = classification_mode_choice, font = button_font, bg = "white", width = 20)
# RData input
select_RData_file_entry <- tkbutton(window, text = "SELECT RData\nWORKSPACE", command = select_RData_file_function, font = button_font, bg = "white", width = 20)
# File type export matrix
file_type_export_matrix_entry <- tkbutton(window, text = "FILE TYPE\nEXPORT\nMATRIX", command = file_type_export_matrix_choice, font = button_font, bg = "white", width = 20)
# File type export images
file_type_export_images_entry <- tkbutton(window, text = "FILE TYPE\nEXPORT\nIMAGES", command = file_type_export_images_choice, font = button_font, bg = "white", width = 20)
# Multicore
allow_parallelization_button <- tkbutton(window, text = "ALLOW\nPARALLEL\nCOMPUTING", command = allow_parallelization_choice, font = button_font, bg = "white", width = 20)
# Spectra preprocessing button
spectra_preprocessing_button <- tkbutton(window, text = "SPECTRA\nPREPROCESSING\nPARAMETERS...", command = preprocessing_window_function, font = button_font, bg = "white", width = 20)
# End session
end_session_button <- tkbutton(window, text = "QUIT", command = end_session_function, font = button_font, bg = "white", width = 20)
# Run the MS Pixel Typer
run_patient_classification_function_button <- tkbutton(window, text = "RUN THE\nMS PIXEL TYPER", command = run_patient_classification_function, font = button_font, bg = "white", width = 20)
# Dump the files
ms_pixel_typer_data_dumper_button <- tkbutton(window, text = "DUMP\nFILES", command = ms_pixel_typer_data_dumper_function, font = button_font, bg = "white", width = 20)
# Updates
download_updates_button <- tkbutton(window, text = "DOWNLOAD\nUPDATE", command = download_updates_function, font = button_font, bg = "white", width = 20)

#### Displaying labels
file_type_export_matrix_value_label <- tklabel(window, text = file_type_export_matrix, font = label_font, bg = "white", width = 20)
file_type_export_images_value_label <- tklabel(window, text = file_type_export_images, font = label_font, bg = "white", width = 20)
peak_picking_algorithm_value_label <- tklabel(window, text = peak_picking_algorithm_value, font = label_font, bg = "white", width = 20)
peaks_deisotoping_value_label <- tklabel(window, text = peaks_deisotoping_value, font = label_font, bg = "white", width = 20)
allow_parallelization_value_label <- tklabel(window, text = allow_parallelization_value, font = label_font, bg = "white", width = 20)
vote_weights_ensemble_value_label <- tklabel(window, text = vote_weights_ensemble_value, font = label_font, bg = "white", width = 20)
classification_mode_value_label <- tklabel(window, text = classification_mode_value, font = label_font, bg = "white", width = 20)
RData_file_integrity_value_label <- tklabel(window, text = RData_file_integrity_value, font = label_font, bg = "white", width = 20)
decision_method_ensemble_value_label <- tklabel(window, text = decision_method_ensemble_value, font = label_font, bg = "white", width = 20)
check_for_updates_value_label <- tklabel(window, text = check_for_updates_value, font = label_font, bg = "white", width = 20)

#### Geometry manager
# Scrollbar
#window_scrollbar <- tkscrollbar(window, command = function(...)tkyview(window,...))
# tkgrid
tkgrid(title_label, row = 1, column = 1, columnspan = 2, padx = c(20, 20), pady = c(20, 20))
tkgrid(select_samples_button, row = 7, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(browse_output_button, row = 7, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(peaks_deisotoping_entry, row = 2, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(peaks_deisotoping_value_label, row = 2, column = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(peak_picking_algorithm_entry, row = 2, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(peak_picking_algorithm_value_label, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(file_type_export_matrix_entry, row = 6, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(file_type_export_matrix_value_label, row = 6, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(file_type_export_images_entry, row = 6, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(file_type_export_images_value_label, row = 6, column = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(decision_method_ensemble_entry, row = 3, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(decision_method_ensemble_value_label, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(vote_weights_ensemble_entry, row = 3, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(vote_weights_ensemble_value_label, row = 3, column = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(classification_mode_entry, row = 5, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(classification_mode_value_label, row = 5, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(select_RData_file_entry, row = 5, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(RData_file_integrity_value_label, row = 5, column = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(allow_parallelization_button, row = 4, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(allow_parallelization_value_label, row = 4, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(spectra_preprocessing_button, row = 4, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(run_patient_classification_function_button, row = 7, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(ms_pixel_typer_data_dumper_button, row = 7, column = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(end_session_button, row = 8, column = 1, columnspan = 4, padx = c(10, 10), pady = c(10, 10))
tkgrid(download_updates_button, row = 1, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(check_for_updates_value_label, row = 1, column = 4, padx = c(10, 10), pady = c(10, 10))



