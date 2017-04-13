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
            #
