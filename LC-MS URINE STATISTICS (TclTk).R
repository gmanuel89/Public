#################### LC-MS URINE STATISTICS (TCL-TK GUI) ####################





### Program version (Specified by the program writer!!!!)
R_script_version <- "2017.02.24.4"
### GitHub URL where the R file is
github_R_url <- "https://raw.githubusercontent.com/gmanuel89/Public-R-UNIMIB/master/LC-MS%20URINE%20STATISTICS%20(TclTk).R"
### Name of the file when downloaded
script_file_name <- "LC-MS URINE STATISTICS.R"
# Change log
change_log <- "1. Correlation plots are now generated and saved\n2. All the outputs are put in a subfolder named 'STATISTICS X'\n3. All the outputs in the subfolder are more organized\n4. The boxplots are generated also without the outliers\n5. The labels in the plots are correctly displayed\n6. New GUI (adaptive fonts also on Windows 10)"





#################### INSTALL AND LOAD THE REQUIRED PACKAGES
##### Install the required packages if not already installed
# FUNCTION
install_and_load_required_packages <- function(required_packages, repository = "http://cran.mirror.garr.it/mirrors/CRAN/") {
    # Update all the packages
    try(update.packages(repos = repository, ask = FALSE), silent = TRUE)
    # Retrieve the installed packages
    installed_packages <- installed.packages()[,1]
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
            install.packages(missing_packages, repos = repository)
        }
    } else {
        if (length(missing_packages) > 0) {
            install.packages(missing_packages)
        }
    }
    # Load the packages
    for (i in 1:length(required_packages)) {
        library(required_packages[i], character.only = TRUE)
    }
}

# RUN THE FUNCTION
install_and_load_required_packages(c("rattle", "phia", "MASS", "ggplot2", "lawstat", "coin", "multcomp", "agricolae", "tcltk", "Hmisc")) # "Rcmdr", "RcmdrPlugin.coin"
# The package lawstat is used for levene test for non parametric anova








###################################### Initialize the variables (default values)
input_file <- ""
output_folder <- getwd()
output_format <- "Comma Separated Values (.csv)"
file_format <- "csv"
image_format <- ".png"
pvalue_expression <- 0.05
pvalue_tests <- 0.05
TestPer_Base <- 0.17    #fair for class
TestPer_Adv <- 0.19        #fair for class
minimum_number_of_patients <- 3
### Data Management
data_record <- TRUE
### Correlation Analysis
correlation_analysis <- FALSE
# Outlier estimation and removal
remove_outliers_correlation_analysis <- FALSE
### Effect Analysis (2-Level Effects)
# One class vs all the others (if more classes are present); one class vs the other class (if only two classes are present) 
two_level_effect_analysis  <- TRUE
# The sampling = TRUE performs the sampling of around 20% of patients/ctrls to be used in R... The other 80% is exported for RapidMiner analysis (with the filtered signals)
sampling <- FALSE
# Outlier estimation in the 2-level effect analysis
remove_outliers_two_level_effect_analysis <- FALSE
# Cumulative class in two level effect
cumulative_class_in_two_level_effect_analysis <- FALSE
### Effect Analysis (Multi-Level Effects) Stima pT (pT4 e pT5 sono uniti) e Grade
multi_level_effect_analysis <- TRUE
remove_outliers_multi_level_effect_analysis <- FALSE
# Age binning
age_binning <- FALSE
age_bins <- 6






################## Values of the variables (for displaying and dumping purposes)
output_file_type_export_value <- "Comma Separated Values (.csv)"
image_file_type_export_value <- "PNG (.png)"
data_record_value <- "YES"
correlation_analysis_value <- "NO"
remove_outliers_correlation_analysis_value <- "NO"
two_level_effect_analysis_value <- "YES"
remove_outliers_two_level_effect_analysis_value <- "NO"
multi_level_effect_analysis_value <- "YES"
remove_outliers_multi_level_effect_analysis_value <- "NO"
cumulative_class_in_two_level_effect_analysis_value <- "NO"
check_for_updates_value <- R_script_version







##################################################### DEFINE WHAT THE BUTTONS DO

##### Check for updates (from my GitHub page) (it just updates the label telling the user if there are updates) (it updates the check for updates value that is called by the label)
check_for_updates_function <- function() {
    ### Initialize the version number
    online_version_number <- NULL
    ### Initialize the variable that says if there are updates
    update_available <- FALSE
    ### Initialize the change log
    online_change_log <- "Bug fixes"
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
        online_version_YYYYMMDDVV <- unlist(strsplit(online_version_number, ".", fixed=TRUE))
        ### Compare with the local version
        local_version_YYYYMMDDVV = unlist(strsplit(R_script_version, ".", fixed = TRUE))
        ### Check the versions
        for (v in 1:length(local_version_YYYYMMDDVV)) {
            if (as.numeric(local_version_YYYYMMDDVV[v]) < as.numeric(online_version_YYYYMMDDVV[v])) {
                update_available <- TRUE
                break
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
    ### Something went wrong: library not installed, retrieving failed, errors in parsing the version number
    if (is.null(online_version_number)) {
        # Update the label
        check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked: connection problems", sep = "")
    }
    # Escape the function
    .GlobalEnv$update_available <- update_available
    .GlobalEnv$online_change_log <- online_change_log
    .GlobalEnv$check_for_updates_value <- check_for_updates_value
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
            download.file(url = github_R_url, destfile = script_file_name, method = "auto")
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

##### Output file type (export)
output_file_type_export_choice <- function() {
    # Catch the value from the menu
    output_format <- select.list(c("Comma Separated Values (.csv)", "Microsoft Excel (.xls)", "Microsoft Excel (.xlsx)"), title = "Choose output file format")
    # Fix the file format
    if (output_format == "Comma Separated Values (.csv)" || output_format == "") {
        file_format <- "csv"
    } else if (output_format == "Microsoft Excel (.xlsx)") {
        file_format <- "xlsx"
        install_and_load_required_packages("XLConnect")
    } else if (output_format == "Microsoft Excel (.xls)") {
        file_format <- "xls"
        install_and_load_required_packages("XLConnect")
    }
    # Set the value of the displaying label
    output_file_type_export_value_label <- tklabel(window, text = output_format, font = label_font)
    tkgrid(output_file_type_export_value_label, row = 2, column = 2)
    # Escape the function
    .GlobalEnv$output_format <- output_format
    .GlobalEnv$file_format <- file_format
}

##### Image file type (export)
image_file_type_export_choice <- function() {
    # Catch the value from the menu
    image_output_format <- select.list(c("JPG (.jpg)", "PNG (.png)", "TIFF (.tiff)"), title = "Choose image format")
    # Fix the file format
    if (image_output_format == "JPG (.jpg)") {
        image_format <- ".jpg"
    } else if (image_output_format == "PNG (.png)" || image_output_format == "") {
        image_format <- ".png"
    } else if (image_output_format == "TIFF (.tiff)") {
        image_format <- ".tiff"
    }
    # Set the value of the displaying label
    image_file_type_export_value_label <- tklabel(window, text = image_output_format, font = label_font)
    tkgrid(image_file_type_export_value_label, row = 3, column = 2)
    # Escape the function
    .GlobalEnv$image_output_format <- image_output_format
    .GlobalEnv$image_format <- image_format
}

##### File import
file_import_function <- function() {
    filepath_import_select <- tkmessageBox(title = "Input file", message = "Select the file containing all the mass spectrometric information for the statistics", icon = "info")
    input_file <- tclvalue(tkgetOpenFile(filetypes = "{{Microsoft Excel files} {.xls .xlsx}} {{Comma Separated Value files} {.csv}}"))
    if (!nchar(input_file)) {
        tkmessageBox(message = "No file selected")
    } else {
        tkmessageBox(message = paste("The following file will be read:", input_file))
    }
    if (input_file != "") {
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
            input_data <- read.csv(input_file, header=TRUE, sep=",")
        }
        ## Rownames
        try({rownames(input_data) <- input_data$No}, silent = TRUE)
        ## Data type
        try({input_data$Age <- as.numeric(input_data$Age)}, silent = TRUE)
        try({input_data$pT <- as.numeric(input_data$pT)}, silent = TRUE)
        try({input_data$pT_2009 <- as.factor(input_data$pT_2009)}, silent = TRUE)
        try({input_data$Dim <- as.numeric(input_data$Dim)}, silent = TRUE)
        try({input_data$Grade <- as.numeric(input_data$Grade)}, silent = TRUE)
        try({input_data$Class <- as.factor(input_data$Class)}, silent = TRUE)
        ## Age binning
        try({
            if (isTRUE(age_binning)) {
                age_bins <- age_bins
                Age_BINNED <- binning(input_data$Age, age_bins, method = "quantile", ordered = TRUE)
                levels_age_binned <- levels(Age_BINNED) <- c(1:age_bins)
                # Replace the Age with the Age_BINNED
                input_data$Age <- cbind(Age_BINNED)
            }
            }, silent = TRUE)
        ##### Separate the mass spectrometric data from the demographic data
        # All features
        feature_vector <- colnames(input_data)
        # Non signals (columns)
        tkmessageBox(title = "Demographical data", message = "Select the demographic data (to be separated from the mass spectrometric data)", icon = "info")
        non_signals <- select.list(feature_vector, title = "Choose the demographical features", multiple = TRUE)
        non_signals_data <- input_data[,(names(input_data) %in% non_signals)]
        signals_data <- input_data[,!(names(input_data) %in% non_signals)]
        number_of_signals <- ncol(signals_data)
        # Matrix of signal names
        signal_name <- data.matrix(colnames(signals_data))
        colnames(signal_name) <- "Signal Name"
        ##### Determine the features of interest
        # Discriminant column
        tkmessageBox(title = "Discriminant feature", message = "Select the discriminant feature", icon = "info")
        discriminant_feature <- select.list(c(non_signals, "NONE"), title = "Choose the discriminant feature", preselect = ifelse("Class" %in% feature_vector, "Class", NULL))
        # Features for correlation analysis
        if (correlation_analysis == TRUE) {
            tkmessageBox(title = "Correlation analysis data", message = "Select the demographic data for correlation analysis", icon = "info")
            non_signals_for_correlation_analysis <- select.list(non_signals[non_signals != discriminant_feature], title = "Choose the features for correlation analysis", multiple = TRUE)
        }
        # Features for two-level effect analysis
        if (two_level_effect_analysis == TRUE) {
            tkmessageBox(title = "Two-level effect analysis data", message = "Select the demographic data for two-level effect analysis", icon = "info")
            two_level_effect_analysis_non_features <- select.list(non_signals[non_signals != discriminant_feature], title = "Choose the features for two-level effect analysis", multiple = TRUE)
        }
        if (multi_level_effect_analysis == TRUE) {
            # Features for multi-level effect analysis
            tkmessageBox(title = "Multi-level effect analysis data", message = "Select the demographic data for multi-level effect analysis", icon = "info")
            multi_level_effect_analysis_non_features <- select.list(non_signals[non_signals != discriminant_feature], title = "Choose the features for multi-level effect analysis", multiple = TRUE)
        }
        ## Class list
        if (discriminant_feature == "NONE") {
            class_list <- discriminant_feature
        } else {
            class_list <- levels(as.factor(input_data[, discriminant_feature]))
        }
        # Escape the function
        .GlobalEnv$input_file <- input_file
        .GlobalEnv$input_data <- input_data
        .GlobalEnv$feature_vector <- feature_vector
        .GlobalEnv$non_signals <- non_signals
        .GlobalEnv$non_signals_data <- non_signals_data
        .GlobalEnv$number_of_signals <- number_of_signals
        .GlobalEnv$signals_data <- signals_data
        .GlobalEnv$signal_name <- signal_name
        .GlobalEnv$discriminant_feature <- discriminant_feature
        .GlobalEnv$non_signals_for_correlation_analysis <- non_signals_for_correlation_analysis
        .GlobalEnv$two_level_effect_analysis_non_features <- two_level_effect_analysis_non_features
        .GlobalEnv$multi_level_effect_analysis_non_features <- multi_level_effect_analysis_non_features
        .GlobalEnv$class_list <- class_list
        tkmessageBox(title = "File imported", message = "The data has been successfully imported from the file!", icon = "info")
    } else {
        # Escape the function
        .GlobalEnv$input_file <- input_file
        tkmessageBox(title = "No input file selected", message = "No input file has been selected!!!\nPlease, select a file to be imported", icon = "warning")
    }
}

##### Output
browse_output_function <- function() {
    output_folder <- tclvalue(tkchooseDirectory())
    if (!nchar(output_folder)) {
        # Get the output folder from the default working directory
        output_folder <- getwd()
    }
    # Go to the working directory
    setwd(output_folder)
    tkmessageBox(message = paste("Every file will be saved in", output_folder))
    tkmessageBox(message = "A sub-directory named 'STATISTICS X' will be created for each run!")
    # Escape the function
    .GlobalEnv$output_folder <- output_folder
}

##### Exit
end_session_function <- function() {
    q(save = "no")
}

##### Data Record
data_record_choice <- function() {
    # Catch the value from the menu
    data_record <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (data_record == "YES" || data_record == "") {
        data_record <- TRUE
    }
    if (data_record == "NO") {
        data_record <- FALSE
    }
    # Set the value of the displaying label
    if (data_record == TRUE) {
        data_record_value <- "YES"
    } else {
        data_record_value <- "NO"
    }
    data_record_value_label <- tklabel(window, text = data_record_value, font = label_font)
    tkgrid(data_record_value_label, row = 4, column = 2)
    # Escape the function
    .GlobalEnv$data_record <- data_record
    .GlobalEnv$data_record_value <- data_record_value
}

##### Correlation analysis
correlation_analysis_choice <- function() {
    # Catch the value from the menu
    correlation_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (correlation_analysis == "YES") {
        correlation_analysis <- TRUE
    }
    if (correlation_analysis == "NO" || correlation_analysis == "") {
        correlation_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (correlation_analysis == TRUE) {
        correlation_analysis_value <- "YES"
    } else {
        correlation_analysis_value <- "NO"
    }
    correlation_analysis_value_label <- tklabel(window, text = correlation_analysis_value, font = label_font)
    tkgrid(correlation_analysis_value_label, row = 5, column = 2)
    # Escape the function
    .GlobalEnv$correlation_analysis <- correlation_analysis
    .GlobalEnv$correlation_analysis_value <- correlation_analysis_value
}

##### Remove outliers correlation analysis
remove_outliers_correlation_analysis_choice <- function() {
    # Catch the value from the menu
    remove_outliers_correlation_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (remove_outliers_correlation_analysis == "YES") {
        remove_outliers_correlation_analysis <- TRUE
    }
    if (remove_outliers_correlation_analysis == "NO" || remove_outliers_correlation_analysis == "") {
        remove_outliers_correlation_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (remove_outliers_correlation_analysis == TRUE) {
        remove_outliers_correlation_analysis_value <- "YES"
    } else {
        remove_outliers_correlation_analysis_value <- "NO"
    }
    remove_outliers_correlation_analysis_value_label <- tklabel(window, text = remove_outliers_correlation_analysis_value, font = label_font)
    tkgrid(remove_outliers_correlation_analysis_value_label, row = 6, column = 2)
    # Escape the function
    .GlobalEnv$remove_outliers_correlation_analysis <- remove_outliers_correlation_analysis
    .GlobalEnv$remove_outliers_correlation_analysis_value <- remove_outliers_correlation_analysis_value
}

##### Two-level effect analysis
two_level_effect_analysis_choice <- function() {
    # Catch the value from the menu
    two_level_effect_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (two_level_effect_analysis == "YES" || two_level_effect_analysis == "") {
        two_level_effect_analysis <- TRUE
    }
    if (two_level_effect_analysis == "NO") {
        two_level_effect_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (two_level_effect_analysis == TRUE) {
        two_level_effect_analysis_value <- "YES"
    } else {
        two_level_effect_analysis_value <- "NO"
    }
    two_level_effect_analysis_value_label <- tklabel(window, text = two_level_effect_analysis_value, font = label_font)
    tkgrid(two_level_effect_analysis_value_label, row = 2, column = 4)
    # Escape the function
    .GlobalEnv$two_level_effect_analysis <- two_level_effect_analysis
    .GlobalEnv$two_level_effect_analysis_value <- two_level_effect_analysis_value
}

##### Cumulative class in Two-level effect analysis
cumulative_class_in_two_level_effect_analysis_choice <- function() {
    # Catch the value from the menu
    cumulative_class_in_two_level_effect_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (cumulative_class_in_two_level_effect_analysis == "YES" || cumulative_class_in_two_level_effect_analysis == "") {
        cumulative_class_in_two_level_effect_analysis <- TRUE
    }
    if (cumulative_class_in_two_level_effect_analysis == "NO") {
        cumulative_class_in_two_level_effect_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (cumulative_class_in_two_level_effect_analysis == TRUE) {
        cumulative_class_in_two_level_effect_analysis_value <- "YES"
    } else {
        cumulative_class_in_two_level_effect_analysis_value <- "NO"
    }
    cumulative_class_in_two_level_effect_analysis_value_label <- tklabel(window, text = cumulative_class_in_two_level_effect_analysis_value, font = label_font)
    tkgrid(cumulative_class_in_two_level_effect_analysis_value_label, row = 7, column = 4)
    # Escape the function
    .GlobalEnv$cumulative_class_in_two_level_effect_analysis <- cumulative_class_in_two_level_effect_analysis
    .GlobalEnv$cumulative_class_in_two_level_effect_analysis_value <- cumulative_class_in_two_level_effect_analysis_value
}

##### Multi-level effect analysis
multi_level_effect_analysis_choice <- function() {
    # Catch the value from the menu
    multi_level_effect_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (multi_level_effect_analysis == "YES" || multi_level_effect_analysis == "") {
        multi_level_effect_analysis <- TRUE
    }
    if (multi_level_effect_analysis == "NO") {
        multi_level_effect_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (multi_level_effect_analysis == TRUE) {
        multi_level_effect_analysis_value <- "YES"
    } else {
        multi_level_effect_analysis_value <- "NO"
    }
    multi_level_effect_analysis_value_label <- tklabel(window, text = multi_level_effect_analysis_value, font = label_font)
    tkgrid(multi_level_effect_analysis_value_label, row = 4, column = 4)
    # Escape the function
    .GlobalEnv$multi_level_effect_analysis <- multi_level_effect_analysis
    .GlobalEnv$multi_level_effect_analysis_value <- multi_level_effect_analysis_value
}

##### Remove outliers two-level effect analysis
remove_outliers_two_level_effect_analysis_choice <- function() {
    # Catch the value from the menu
    remove_outliers_two_level_effect_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (remove_outliers_two_level_effect_analysis == "YES") {
        remove_outliers_two_level_effect_analysis <- TRUE
    }
    if (remove_outliers_two_level_effect_analysis == "NO" || remove_outliers_two_level_effect_analysis == "") {
        remove_outliers_two_level_effect_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (remove_outliers_two_level_effect_analysis == TRUE) {
        remove_outliers_two_level_effect_analysis_value <- "YES"
    } else {
        remove_outliers_two_level_effect_analysis_value <- "NO"
    }
    remove_outliers_two_level_effect_analysis_value_label <- tklabel(window, text = remove_outliers_two_level_effect_analysis_value, font = label_font)
    tkgrid(remove_outliers_two_level_effect_analysis_value_label, row = 3, column = 4)
    # Escape the function
    .GlobalEnv$remove_outliers_two_level_effect_analysis <- remove_outliers_two_level_effect_analysis
    .GlobalEnv$remove_outliers_two_level_effect_analysis_value <- remove_outliers_two_level_effect_analysis_value
}

##### Remove outliers multi-level effect analysis
remove_outliers_multi_level_effect_analysis_choice <- function() {
    # Catch the value from the menu
    remove_outliers_multi_level_effect_analysis <- select.list(c("YES","NO"), title="Choose")
    # Default
    if (remove_outliers_multi_level_effect_analysis == "YES") {
        remove_outliers_multi_level_effect_analysis <- TRUE
    }
    if (remove_outliers_multi_level_effect_analysis == "NO" || remove_outliers_multi_level_effect_analysis == "") {
        remove_outliers_multi_level_effect_analysis <- FALSE
    }
    # Set the value of the displaying label
    if (remove_outliers_multi_level_effect_analysis == TRUE) {
        remove_outliers_multi_level_effect_analysis_value <- "YES"
    } else {
        remove_outliers_multi_level_effect_analysis_value <- "NO"
    }
    remove_outliers_multi_level_effect_analysis_value_label <- tklabel(window, text = remove_outliers_multi_level_effect_analysis_value, font = label_font)
    tkgrid(remove_outliers_multi_level_effect_analysis_value_label, row = 5, column = 4)
    # Escape the function
    .GlobalEnv$remove_outliers_multi_level_effect_analysis <- remove_outliers_multi_level_effect_analysis
    .GlobalEnv$remove_outliers_multi_level_effect_analysis_value <- remove_outliers_multi_level_effect_analysis_value
}

##### Run the statistics
run_statistics_function <- function() {
    # Go to the working directory
    setwd(output_folder)
    if (input_file != "") {
        ##### Automatically create a subfolder with all the results
        ## Check if such subfolder exists
        list_of_directories <- list.dirs(output_folder, full.names = FALSE, recursive = FALSE)
        ## Check the presence of a STATISTICS folder
        STATISTICS_folder_presence <- FALSE
        if (length(list_of_directories) > 0) {
            for (dr in 1:length(list_of_directories)) {
                if (length(grep("STATISTICS", list_of_directories[dr], fixed = TRUE)) > 0) {
                    STATISTICS_folder_presence <- TRUE
                }
            }
        }
        ## If it present...
        if (isTRUE(STATISTICS_folder_presence)) {
            ## Extract the number after the STATISTICS
            # Number for the newly created folder
            STATISTICS_new_folder_number <- 0
            # List of already present numbers
            STATISTICS_present_folder_numbers <- integer()
            # For each folder present...
            for (dr in 1:length(list_of_directories)) {
                # If it is named STATISTICS
                if (length(grep("STATISTICS", list_of_directories[dr], fixed = TRUE)) > 0) {
                    # Split the name
                    STATISTICS_present_folder_split <- unlist(strsplit(list_of_directories[dr], "STATISTICS"))
                    # Add the number to the list of STATISTICS numbers
                    try(STATISTICS_present_folder_numbers <- append(STATISTICS_present_folder_numbers, as.integer(STATISTICS_present_folder_split[2])))
                }
            }
            # Sort the STATISTICS folder numbers
            try(STATISTICS_present_folder_numbers <- sort(STATISTICS_present_folder_numbers))
            # The new folder number will be the greater + 1
            try(STATISTICS_new_folder_number <- STATISTICS_present_folder_numbers[length(STATISTICS_present_folder_numbers)] + 1)
            # Generate the new subfolder
            subfolder <- paste("STATISTICS", STATISTICS_new_folder_number)
            # Create the subfolder
            dir.create(file.path(output_folder, subfolder))
            # Estimate the new output folder
            output_folder <- file.path(output_folder, subfolder)
        } else {
            # If it not present...
            # Create the folder where to dump the files and go to it...
            subfolder <- paste("STATISTICS", "1")
            # Create the subfolder
            dir.create(file.path(output_folder, subfolder))
            # Estimate the new output folder
            output_folder <- file.path(output_folder, subfolder)
        }
        # Go to the new working directory
        setwd(output_folder)
        
        
        ################################################################################
        
        
        # Retrieve the values from the entries
        pvalue_expression <- tclvalue(pvalue_expression)
        pvalue_expression <- as.numeric(pvalue_expression)
        pvalue_expression_value <- as.character(pvalue_expression)
        pvalue_tests <- tclvalue(pvalue_tests)
        pvalue_tests <- as.numeric(pvalue_tests)
        pvalue_tests_value <- as.character(pvalue_tests)
        #TestPer_Base <- tclvalue(TestPer_Base)
        #TestPer_Base <- as.numeric(TestPer_Base)
        #TestPer_Base_value <- as.character(TestPer_Base)
        #TestPer_Adv <- tclvalue(TestPer_Adv)
        #TestPer_Adv <- as.numeric(TestPer_Adv)
        #TestPer_Adv_value <- as.character(TestPer_Adv)
        minimum_number_of_patients <- tclvalue(minimum_number_of_patients)
        minimum_number_of_patients <- as.integer(minimum_number_of_patients)
        minimum_number_of_patients_value <- as.character(minimum_number_of_patients)
        
        
        
        
        ########## Exception handling
        #if (isTRUE(two_level_effect_analysis)) {
        #    ifelse(isTRUE(sampling), print("The Two-Level Effect analysis is performed with sampling for RM"), print("WARNING!!! The Two-Level Effect analysis is performed without sampling!"))
        #}
        
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
        ################################ P-value extractor functions
        
        assumption_p <- function(objt) {
            ifelse(is.object(objt), return(objt$p.value), return(NA))
        }
        assumption_anova <- function(objt) {
            ifelse(is.object(objt), return(objt$Pr[1]), return(NA))
        }
        assumption_permutation <- function(objt) {
            ifelse(is.object(objt), return(pvalue(objt)[1]), return(NA))
        }
        
        
        #################  Post hoc tests
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
            ### Excel
            if (file_format == "xls" || file_format == "xlsx" || file_format == "csv") {
                # Fix the extension
                file_name <- paste(file_name, ".xlsx", sep = "")
                install_and_load_required_packages("XLConnect")
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
                # Fix the extension
                file_name <- paste(file_name, ".", file_format, sep = "")
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
    
        
        
        
        
        
        
        
        
        
        ################ Data from input_data #######################
        if (isTRUE(data_record)) {
            print("########## Data Management ##########")
            # Generate a list of data frames and corresponding filenames to dump
            list_of_dataframes <- list()
            list_of_filenames <- character()
            # If there is only one class...
            if (length(class_list == 1)) {
                # Generate a dataframe with the one class, otherwise dump one dataframe per each class
                list_of_dataframes[[1]] <- input_data[,!(names(input_data) %in% non_signals)]
                list_of_filenames <- "Whole data"
            } else if (length(class_list) == 2) {
            # If there are only two classes...
                # Generate a dataframe with those two classes, otherwise dump one dataframe per each class
                list_of_dataframes[[1]] <- input_data[,!(names(input_data) %in% non_signals)]
                list_of_filenames <- paste("Class", class_list[1], "vs", class_list[2])
            } else if (length(class_list) > 2) {
                # If there are more than two classes...
                for (cl in 1:length(class_list)) {
                    # Put every dataframe in the list
                    list_of_dataframes[[cl]] <- input_data[(input_data[,discriminant_feature] == class_list[cl]), !(names(input_data) %in% non_signals)]
                    list_of_filenames <- append(list_of_filenames, paste(discriminant_feature , class_list[cl]))
                }
            }
            # Create the folder containing all the data frames (and go to it temporarily to dump the files)
            data_rec_subfolder <- file.path(output_folder, "Data")
            dir.create(data_rec_subfolder)
            setwd(data_rec_subfolder)
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
            non_signals_for_correlation_analysis <- sort(non_signals_for_correlation_analysis)
            features_for_correlation_analysis <- append(non_signals_for_correlation_analysis, names(signals_data))
            print("########## CORRELATION ANALYSIS ##########")
            ##### Create the folder for the correlation data
            correlation_subfolder <- file.path(output_folder, "Correlation")
            dir.create(correlation_subfolder)
            setwd(correlation_subfolder)
            ##### Create the folder for the correlation scatter plots
            correlation_plot_subfolder <- file.path(correlation_subfolder, "Plots")
            dir.create(correlation_plot_subfolder)
            ##### Create the folder for the correlation scatter plots (no outliers)
            correlation_plot_subfolder_no_outliers <- file.path(correlation_subfolder, "Plots (No outliers)")
            dir.create(correlation_plot_subfolder_no_outliers)
            ##### Create the folder for the correlation data tables
            correlation_tables_subfolder <- file.path(correlation_subfolder, "Tables")
            dir.create(correlation_tables_subfolder)
            ########## Task #1: -- RCC correlation analysis 
            #Rcc Age VS Signal intensity
            #Rcc Dim VS Signal intensity
            #Rcc pT (four levels) VS Signal intensity
            #Rcc Grade (two levels) VS Signal intensity
            # Temporary data filter for task #1   
            ### For each class...
            for (cl in 1:length(class_list)) {
                # Message
                if (discriminant_feature == "NONE") {
                    print("########## CORRELATION ANALYSIS: Whole data ##########")
                } else {
                    print(paste("########## CORRELATION ANALYSIS, class:", discriminant_feature, class_list[cl], "##########"))
                }
                # Filter the dataframe with only the data for that selected class
                if (discriminant_feature == "NONE") {
                    data_frame_correlation <- input_data[, features_for_correlation_analysis]
                } else {
                    data_frame_correlation <- input_data[input_data[, discriminant_feature] == class_list[cl], features_for_correlation_analysis]
                }
                # Store the original (after it will be modified by the outlier exclusion)
                data_frame_correlation_original <- data_frame_correlation
                ### The correlation list is a list of correlation lists: one list for each correlation, each one of them is a list of the correlation with the signal intensity and the variable (correlation_list[[non_signal1]][[signal1]], correlation_list[[non_signal1]][[signal2]], correlation_list[[non_signal1]][[signal3]], correlation_list[[non_signal2]][[signal1]], correlation_list[[non_signal2]][[signal2]], correlation_list[[non_signal2]][[signal3]], ...)
                correlation_list <- list()
                # List of outliers (one dataframe for each mass)
                outlier_list <- list()
                # List of the correlation matrices
                correlation_matrix_list <- list()
                ### Perform the operations only if a minimum number of patients is present...
                if (nrow(data_frame_correlation) > 1) {
                    ### For each non-signal...
                    for (f in 1:length(non_signals_for_correlation_analysis)) {
                        ns <- non_signals_for_correlation_analysis[f]
                        # Generate the sublist...
                        correlation_sublist <- list()
                        # Isolate the column corresponding to the non-signal...
                        non_signal_column_original <- data_frame_correlation[, ns]
                        ##### Scroll the masses...
                        for (ms in 1:ncol(signals_data)) {
                            m <- colnames(signals_data)[ms]
                            # Isolate the column corresponding to the mass (as a vector)
                            mass_x <- data_frame_correlation[, m]
                            # Remove NAs from the non signal and from the signal
                            mass_x <- mass_x[!is.na(non_signal_column_original)]
                            non_signal_column <- non_signal_column_original[!is.na(non_signal_column_original)]
                            ### Generate a graph before the removal of outliers
                            plot_name <- paste("Correlation", m, "vs", ns)
                            file_name <- sprintf("%s%s", plot_name, image_format)
                            scatter_plot <- qplot(mass_x, non_signal_column, geom = "auto", main = plot_name, ylab = ns, xlab = m)
                            setwd(correlation_plot_subfolder)
                            ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                            ### Outlier detection
                            if (isTRUE(remove_outliers_correlation_analysis)) {
                                # Detect the outliers (fence)
                                outliers <- boxplot(mass_x, plot = FALSE)$out
                                # If there are some...
                                if (length(outliers) > 0) {
                                    # Build the dataframe to be dumped
                                    outliers_dataframe <- data_frame_correlation[mass_x %in% outliers, c(non_signals_for_correlation_analysis, m)]
                                    # Signal name
                                    sign_name <- rep(m, dim(outliers_dataframe)[1])
                                    # Patients
                                    patients <- rownames(outliers_dataframe)
                                    # Append the two columns...
                                    outliers_dataframe <- cbind(outliers_dataframe, patients, sign_name)
                                    # Fix the column names
                                    colnames(outliers_dataframe) <- c(non_signals_for_correlation_analysis,"Intensity","Patient","Mass")
                                    # Store this in the final list of outlier dataframes
                                    outlier_list[[ms]] <- outliers_dataframe
                                    # Replace the intensity in the original dataframe with NA (so that they are excluded)
                                    #data_frame_correlation[mass_x %in% outliers, m] <- NA
                                    mass_x[mass_x %in% outliers] <- NA
                                }
                                # Extract the values without the outliers
                                mass_x_no_outliers <- mass_x[!is.na(mass_x)]
                                non_signal_column_no_outliers <- non_signal_column[!is.na(mass_x)]
                                ### Generate a graph after the removal of outliers
                                plot_name <- paste("Correlation", m, "vs", ns, "(No outliers)")
                                file_name <- sprintf("%s%s", plot_name, image_format)
                                scatter_plot <- qplot(mass_x, non_signal_column, geom = "auto", main = plot_name, ylab = ns, xlab = m)
                                setwd(correlation_plot_subfolder_no_outliers)
                                ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                            }
                            ### Compute the results
                            correlation_result_vector <- c(cor.test(mass_x, as.numeric(non_signal_column), method = "spearman")$estimate, cor.test(mass_x, as.numeric(non_signal_column), method = "spearman")$p.value)
                            ### Correlation (signal m with the non-signal ns)
                            # Append to the final correlation list
                            correlation_list[[ns]][[m]] <- correlation_result_vector
                        }
                        # Restore the original
                        #data_frame_correlation <- data_frame_correlation_original
                    }
                    ### Correlation matrices
                    # For each element of the correlation list...
                    for (l in 1:length(correlation_list)) {
                        # Generate a matrix from the correlation list (one matrix for each non-signal)
                        corr_matrix <- matrix(0, nrow = ncol(signals_data), ncol = 2)
                        corr_matrix <- do.call(rbind, correlation_list[[l]])
                        colnames(corr_matrix) <- c(paste("Correlation", class_list[cl], "vs", non_signals_for_correlation_analysis[l]), "p_value")
                        correlation_matrix_list[[l]] <- corr_matrix
                    }
                    ### Build the final correlation matrix
                    final_correlation_matrix <- as.matrix(cbind(signal_name))
                    for (l in 1:length(correlation_matrix_list)) {
                        final_correlation_matrix <- cbind(final_correlation_matrix, correlation_matrix_list[[l]])
                    }
                    ### Outlier matrices
                    if (isTRUE(remove_outliers_correlation_analysis)) {
                        final_outlier_matrix <- matrix(ncol = (length(non_signals_for_correlation_analysis)) + 3)
                        final_outlier_matrix <- do.call(rbind, outlier_list)
                        colnames(final_outlier_matrix) <- c(non_signals_for_correlation_analysis , "Intensity", "Patient", "Mass")
                    }
                }
                #### Save the outputs
                setwd(correlation_tables_subfolder)
                if (discriminant_feature == "NONE") {
                    write_file(file_name = "Correlation matrix", data = final_correlation_matrix, file_format = file_format)    
                    if (isTRUE(remove_outliers_correlation_analysis)) {
                        write_file(file_name = "Outlier matrix", data = final_outlier_matrix, file_format = file_format)
                    }
                } else {
                    write_file(file_name = paste(discriminant_feature, class_list[cl], "- Correlation matrix"), data = final_correlation_matrix, file_format = file_format)    
                    if (isTRUE(remove_outliers_correlation_analysis)) {
                        write_file(file_name = paste(discriminant_feature, class_list[cl], "- Outlier matrix"), data = final_outlier_matrix, file_format = file_format)
                    }
                }
            }
        }
        
        
        
        #################################################################################
        #################### 2-LEVEL EFFECT ANALYSIS OVER SIGNAL INTENSITY
        if (isTRUE(two_level_effect_analysis)) {
            # Establish the variables for the two-level effect analysis (Append the discriminant column to the previously selected variables)
            if (discriminant_feature != "NONE") {
                two_level_effect_analysis_non_features <- sort(c(discriminant_feature, two_level_effect_analysis_non_features))
            }
            print("########## 2-LEVEL EFFECT ANALYSIS ##########")
            ### Build up the combinations of possibilities
            # Empty vector
            combination_vector <- character()
            # For each non-signal variable to be considered...
            for (t in 1:length(sort(two_level_effect_analysis_non_features))) {
                # Isolate the dataframe column
                data_column <- input_data[, two_level_effect_analysis_non_features[t]]
                # Check if there are more than two level (because if you have 2 levels there is only one combination)
                number_of_levels <- length(levels(as.factor(data_column)))
                if (number_of_levels > 2) {
                    # Extract the (unique) values
                    for (ch in sort(unique(as.character(data_column)))) {
                        # Append the combination ID to the vector...
                        combination <- ch
                        names(combination) <- two_level_effect_analysis_non_features[t]
                        combination_vector <- append(combination_vector, combination)
                    }
                } else {
                    combination <- sort(unique(as.character(data_column)))[1]
                    names(combination) <- two_level_effect_analysis_non_features[t]
                    combination_vector <- append(combination_vector, combination)
                }
            }
            # Generate the combination names (each element of the vector vs the others)
            combination_names <- character()
            for (l in 1:length(combination_vector)) {
                if (isTRUE(cumulative_class_in_two_level_effect_analysis)) {
                    combination_name <- paste(names(combination_vector[l]), combination_vector[l], "vs OTHERS more than", combination_vector[l])
                } else {
                    combination_name <- paste(names(combination_vector[l]), combination_vector[l], "vs OTHERS")
                }
                combination_names <- append(combination_names, combination_name)
            }
            ##### For each combination... (everytime fall on the 0-1 case, by replacing values)
            for (comb in 1:length(combination_names)) {
                # Create the folder in which every file goes
                # Create the folder where to dump the files and go to it...
                combination_subfolder <- file.path(output_folder, combination_names[comb])
                dir.create(combination_subfolder)
                setwd(combination_subfolder)
                # Message
                print(paste("Processing ->", combination_names[comb]))
                ########## Each time put one class with the value 0 and the others with 1
                # Non signal variable
                non_signal_variable <- names(combination_vector[comb])
                # Non signal variable level selected (it will be put = 0)
                l <- combination_vector[comb]
                # Isolate the dataframe with the signals + the selected non-signal variable
                temp_data_frame <- input_data[, c(non_signal_variable, names(signals_data))]
                # Store the original data frame
                temp_data_frame_original <- temp_data_frame
                # Remove NAs from the non signal and from the signal
                temp_data_frame <- temp_data_frame[!is.na(temp_data_frame[, non_signal_variable]), ]
                # Convert the selected non-feature variable to character
                temp_data_frame[, non_signal_variable] <- as.character(temp_data_frame[, non_signal_variable])
                if (isTRUE(cumulative_class_in_two_level_effect_analysis)) {
                    # Replace one level (the current level l) with 0 and the other ones (more than the selected one) with 1 (cumulative option), discarding the others (less than the selected one) (by setting them as NAs and then removing them)
                    #temp_data_frame[temp_data_frame[, non_signal_variable] == l, non_signal_variable] <- "0"
                    try(temp_data_frame[as.numeric(temp_data_frame[, non_signal_variable]) < as.numeric(l), non_signal_variable] <- NA)
                    # Remove NAs from the non signal and from the signal
                    temp_data_frame <- temp_data_frame[!is.na(temp_data_frame[, non_signal_variable]), ]
                    try(temp_data_frame[as.numeric(temp_data_frame[, non_signal_variable]) > as.numeric(l), non_signal_variable] <- paste(">", l))
                    # Generate a corresponding temp data frame but without outliers
                    temp_data_frame_no_outliers <- temp_data_frame
                } else {
                    # Replace one level (the current level l) with 0 and the other ones with 1 (non-cumulative option)
                    #temp_data_frame[temp_data_frame[, non_signal_variable] == l, non_signal_variable] <- "0"
                    temp_data_frame[temp_data_frame[, non_signal_variable] != l, non_signal_variable] <- "OTHERS"
                    # Generate a corresponding temp data frame but without outliers
                    temp_data_frame_no_outliers <- temp_data_frame
                }
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
                    response_variable <- as.factor(temp_data_frame[, non_signal_variable])
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
                            sign_name <- rep(m, dim(outliers_dataframe)[1])
                            # Patients
                            patients <- rownames(outliers_dataframe)
                            # Append the two columns...
                            outliers_dataframe <- cbind(outliers_dataframe, patients, sign_name)
                            # Fix the column names
                            colnames(outliers_dataframe) <- c(non_signal_variable,"Intensity","Patient","Mass")
                            # Store this in the final list of outlier dataframes
                            outlier_list[[ms]] <- outliers_dataframe
                            # Replace the intensity in the original dataframe with NA (so that they are excluded)
                            temp_data_frame_no_outliers[mass_x %in% outliers, m] <- NA
                            mass_x[mass_x %in% outliers] <- NA
                        }
                    }
                    ## Group the dependent variable (the mass column) according to the response variable (for Shapiro test)
                    mass_x_split <- group_dependent_variable(mass_x, as.factor(response_variable))
                    ### Clean the NA values from the mass column and the response variable column (after the grouping, because the NA removal can cause the loss of a response variable, if it contains only the outlier)
                    response_variable <- response_variable[!is.na(mass_x)]
                    mass_x <- mass_x[!is.na(mass_x)]
                    for (x in 1:length(mass_x_split)) {
                        mass_x_split[[x]] <- mass_x_split[[x]][!is.na(mass_x_split[[x]])]
                    }
                    ## If the response variable has only one level, it0s pointless to fodo ther statistics
                    if (length(mass_x_split) > 1) {
                        # Define the groups (response variable = 0, response variable = 1)
                        group_0 <- mass_x_split[[1]]
                        group_1 <- mass_x_split[[2]]
                        # Define the names and the number of patients
                        name_group_0 <- sprintf("%s%s", m, l)
                        name_group_1 <- sprintf("%s%s", m, "_OTHERS")
                        number_of_patients_list[[m]] <- c(length(group_0),length(group_1))
                    } else if (length(mass_x_split) == 1) {
                        # Define the groups (response variable = 0, response variable = 1)
                        group_0 <- mass_x_split[[1]]
                        group_1 <- numeric()
                        # Define the names and the number of patients
                        name_group_0 <- sprintf("%s%s", m, l)
                        name_group_1 <- sprintf("%s%s", m, "_OTHERS")
                        number_of_patients_list[[m]] <- c(length(group_0),length(group_1))
                    }
                    ### If the number of patients in the two lists is more than the minimum number of patients allowed and all the elements of the two groups are not the same...
                    if (length(group_0) >= minimum_number_of_patients && length(group_1) >= minimum_number_of_patients && !all(group_0 == group_0[1]) && !all(group_1 == group_1[1])) {
                        # Checking for normal distributed data in the two groups (class 0 and class 1)
                        if (length(group_0) >= 3 && length(group_0) <= 5000) {
                            shapiro_list_0[[name_group_0]] <- shapiro.test(as.numeric(group_0))
                        } else {
                            shapiro_list_0[[name_group_0]] <- NA
                        }
                        if (length(group_1) >= 3 && length(group_1) <= 5000) {
                            shapiro_list_1[[name_group_1]] <- shapiro.test(as.numeric(group_1))
                        } else {
                            shapiro_list_1[[name_group_1]] <- NA
                        }
                        # Checking variances (Normal data)
                        variance_test_list[[m]] <- var.test(as.numeric(group_0), as.numeric(group_1))
                        # Checking variances (Non-Normal data)
                        Leven_test_list[[m]] <- levene.test(y = as.numeric(mass_x), group = response_variable, location = "mean", kruskal.test = T)
                        # Normal Data, Equal Variances --> Equal Variance t-test
                        EqTTest_list[[m]] <- t.test(as.numeric(group_0), as.numeric(group_1))
                        # Normal Data, Unequal Variances --> Unequal Variance t-test:
                        UneqTTest_list[[m]] <- t.test(as.numeric(group_0), as.numeric(group_1), var.equal = FALSE)
                        # Non-normal Data, Equal Variances --> Mann-Whitney U-test (Wilcoxon)
                        MWUTest_list[[m]] <- wilcox.test(as.numeric(group_0), as.numeric(group_1), paired = FALSE)
                        # Non-normal Data, Unequal Variances --> Kolmogorov-Smirnov test
                        KSTest_list[[m]] <- ks.test(as.numeric(group_0), as.numeric(group_1))
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
                # Message (if the elements in the lists are NA it means that there weren't enough patients per class)
                #if (is.na(vector_of_shapiro_list_0[1])) {
                #    tkmessageBox(title = "Insufficient number of patients in one class", message = "There is an insufficient number of patients in at least one class: the statistics cannot be performed!", icon = "warning")
                #}
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
                non_normal_homoschedastic_data <-    as.data.frame(subset(non_normal_data, vector_of_Leven_test_list > pvalue_tests))
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
                if (length(selected_signals_for_inference) > 0) {
                    print("Differentially expressed signals for inference")
                    print(selected_signals_for_inference)
                } else {
                    print(paste("There are no differentially expressed signals between", names(l), l, "and OTHERS or there is an insuffucient number of patients in at least one class or there is only one class!"))
                }
                # Generate a matrix with the signals for inference and the method used to identify them...
                inference_signals_method_matrix <- rbind(matrix_diff_norm_homo, matrix_diff_norm_hetero, matrix_diff_non_norm_homo, matrix_diff_non_norm_hetero)
                # Select the signals to be used for inference and the non-signal variable
                selected_signals_for_inference_intensity_df <- subset(temp_data_frame, select = c(non_signal_variable, selected_signals_for_inference))
                if (remove_outliers_two_level_effect_analysis == TRUE) {
                    selected_signals_for_inference_intensity_df_no_outliers <- subset(temp_data_frame_no_outliers, select = c(non_signal_variable, selected_signals_for_inference))
                }
                ### data for testing: for diagnostic analysis 
                if (isTRUE(sampling)) TEST <- DIAGTFD[c("No",BaseEffName,non_signal_variable,SelectedSignsForEffect)]
                ############################### Dump the files
                # For each signals of inference...
                ### OUTLIERS
                # Create the folder where to dump the files and go to it...
                plots_two_level_effect_analysis_subfolder <- file.path(combination_subfolder, "Plots")
                dir.create(plots_two_level_effect_analysis_subfolder)
                setwd(plots_two_level_effect_analysis_subfolder)
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
                    jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5), ylab = "Signal intensity", xlab = non_signal_variable)
                    ggsave(jitter_plot, file = file_name, width = 4, height = 4)
                    
                    ##### Box plot
                    plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
                    file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
                    box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot", ylab = "Signal intensity", xlab = non_signal_variable)
                    ggsave(box_plot, file = file_name, width = 4, height = 4)
                    
                    ##### Scatter plot
                    plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
                    # Sort the dataframe rows according to the values of the non-signal variable
                    ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
                    ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
                    colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
                    file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
                    graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
                    scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors, ylab = "Signal intensity", xlab = "Signals (grouped)")
                    ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                }
                ### NO OUTLIERS
                if (remove_outliers_two_level_effect_analysis == TRUE) {
                    # Create the folder where to dump the files and go to it...
                    plots_two_level_effect_analysis_no_outliers_subfolder <- file.path(combination_subfolder, "Plots (No outliers)")
                    dir.create(plots_two_level_effect_analysis_no_outliers_subfolder)
                    setwd(plots_two_level_effect_analysis_no_outliers_subfolder)
                    for (s in selected_signals_for_inference) {
                        # Extract the intensity
                        signal_intensity <- selected_signals_for_inference_intensity_df_no_outliers[[s]]
                        # Extract the non-signal variable as an ordered factor
                        non_signal_as_ordered_factor <- ordered(temp_data_frame_no_outliers[, non_signal_variable])
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
                        jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5), ylab = "Signal intensity", xlab = non_signal_variable)
                        ggsave(jitter_plot, file = file_name, width = 4, height = 4)
                        
                        ##### Box plot
                        plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
                        file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
                        box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot", ylab = "Signal intensity", xlab = non_signal_variable)
                        ggsave(box_plot, file = file_name, width = 4, height = 4)
                        
                        ##### Scatter plot
                        plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
                        # Sort the dataframe rows according to the values of the non-signal variable
                        ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
                        ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
                        colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
                        file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
                        graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
                        scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors, ylab = "Signal intensity", xlab = "Signals (grouped)")
                        ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                    }
                }
                # Create the folder where to dump the files and go to it...
                tables_two_level_effect_analysis_subfolder <- file.path(combination_subfolder, "Tables")
                dir.create(tables_two_level_effect_analysis_subfolder)
                setwd(tables_two_level_effect_analysis_subfolder)
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
                # Restore the original data frame
                #temp_data_frame <- temp_data_frame_original
            }
        }
        
        
        
        
        
        
        
        
        
        #################################################################################
        #################### MULTI-LEVEL EFFECT ANALYSIS OVER SIGNAL INTENSITY
        if (isTRUE(multi_level_effect_analysis)) {
            # Establish the variables for the multi-level effect analysis (Append the discriminant column to the previously selected variables)
            if (discriminant_feature != "NONE") {
                multi_level_effect_analysis_non_features <- sort(c(discriminant_feature, multi_level_effect_analysis_non_features))
            }
            print("########## MULTI-LEVEL EFFECT ANALYSIS ##########")
            ### Build up the combinations of possibilities (the possibilities are simply the name of the non-signal variables, since in each test, the test automatically reads all the levels of the variable)
            combination_vector <- multi_level_effect_analysis_non_features
            ##### For each combination... (everytime read all the levels)
            for (comb in 1:length(combination_vector)) {
                # Create the folder in which every file goes
                # Create the folder where to dump the files and go to it...
                combination_subfolder <- file.path(output_folder, paste(combination_vector[comb], "(multi-level)"))
                dir.create(combination_subfolder)
                setwd(combination_subfolder)
                # Message
                print(paste("Simple effects (multi-level) -> ", combination_vector[comb]))
                # Non signal variable
                non_signal_variable <- combination_vector[comb]
                # Isolate the dataframe with the signals + the selected non-signal variable
                temp_data_frame <- input_data[, c(non_signal_variable, names(signals_data))]
                # Store the original data frame
                temp_data_frame_original <- temp_data_frame
                # Remove NAs from the non signal and from the signal
                temp_data_frame <- temp_data_frame[!is.na(temp_data_frame[, non_signal_variable]), ]
                # Convert the selected non-feature variable to character
                temp_data_frame[, non_signal_variable] <- as.character(temp_data_frame[, non_signal_variable])
                # Generate a corresponding temp data frame but without outliers
                temp_data_frame_no_outliers <- temp_data_frame
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
                    mass_x <- as.numeric(temp_data_frame[, m])
                    response_variable <- temp_data_frame[, non_signal_variable]
                    if (isTRUE(remove_outliers_multi_level_effect_analysis)) {
                        ## Outlier detection
                        # Detect the outliers (fence)
                        outliers <- boxplot(as.numeric(mass_x) ~ response_variable, plot = FALSE)$out
                        # If there are some...
                        if (length(outliers) > 0) {
                            # Build the dataframe to be dumped
                            outliers_dataframe <- temp_data_frame[mass_x %in% outliers, c(non_signal_variable, m)]
                            # Signal name
                            sign_name <- rep(m, dim(outliers_dataframe)[1])
                            # Patients
                            patients <- rownames(outliers_dataframe)
                            # Append the two columns...
                            outliers_dataframe <- cbind(outliers_dataframe, patients, sign_name)
                            # Fix the column names
                            colnames(outliers_dataframe) <- c(non_signal_variable,"Intensity","Patient","Mass")
                            # Store this in the final list of outlier dataframes
                            outlier_list_multi[[ms]] <- outliers_dataframe
                            # Replace the intensity in the original dataframe with NA (so that they are excluded)
                            temp_data_frame_no_outliers[mass_x %in% outliers, m] <- NA
                            mass_x[mass_x %in% outliers] <- NA
                        }
                    }
                    ## Group the dependent variable (the mass column) according to the response variable (for Shapiro test)
                    mass_x_split <- group_dependent_variable(as.numeric(mass_x), as.factor(response_variable))
                    ### Clean the NA values (outliers) from the mass column and the response variable column (after the grouping, because the NA removal can cause the loss of a response variable, if it contains only the outlier)
                    response_variable <- as.factor(response_variable[!is.na(mass_x)])
                    mass_x <- as.factor(mass_x[!is.na(mass_x)])
                    for (x in 1:length(mass_x_split)) {
                        mass_x_split[[x]] <- mass_x_split[[x]][!is.na(mass_x_split[[x]])]
                    }
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
                            if (length(mass_x_split[[x]]) >= 3 && length(mass_x_split[[x]]) <= 5000) {
                                shapiro_list_multi[[x]][[length(shapiro_list_multi[[x]]) + 1]] <- shapiro.test(as.numeric(mass_x_split[[x]]))
                            } else if (length(mass_x_split[[x]]) < 3 || length(mass_x_split[[x]]) > 5000) {
                                shapiro_list_multi[[x]][[length(shapiro_list_multi[[x]]) + 1]] <- NA
                            }
                        }
                        # Do everything only if the factor variable has more than one level
                        if (length(levels(as.factor(response_variable))) > 1) {
                            # Checking variances
                            bartlett_list_multi[[m]] <- bartlett.test(as.numeric(mass_x) ~ response_variable)
                            ### Anova
                            anova_list_multi[[m]] <- anova(lm(formula = as.numeric(mass_x) ~ response_variable))
                            ### Leven
                            leven_list_multi[[m]] <- levene.test(as.numeric(mass_x), response_variable, location = "mean", kruskal.test = T)
                            ### Kruskal
                            kruskal_list_multi[[m]] <- kruskal.test(as.numeric(mass_x) ~ response_variable)
                            ### Welch
                            welch_list_multi[[m]] <- oneway.test(as.numeric(mass_x) ~ response_variable, var.equal = F)
                            ### k-sets Permutation Test for non parametric heteroschedastic data
                            permutation_test_list_multi[[m]] <- oneway_test(as.numeric(mass_x) ~ response_variable, alternative = 'two.sided')
                        } else {
                            bartlett_list_multi[[m]] <- NA
                            anova_list_multi[[m]] <- NA
                            leven_list_multi[[m]] <- NA    
                            kruskal_list_multi[[m]] <- NA
                            welch_list_multi[[m]] <- NA
                            permutation_test_list_multi[[m]] <- NA
                        }
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
                # Message (if the elements in the lists are NA it means that there weren't enough patients per class)
                #if (is.na(vector_of_bartlett_pvalue[1])) {
                #    tkmessageBox(title = "Insufficient number of patients in one class or all patients in one class", message = "There is an insufficient number of patients in at least one class or there is only one class: the statistics cannot be performed!", icon = "warning")
                #}
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
                non_normal_homoschedastic_data <-    as.data.frame(subset(non_normal_data, vector_of_leven_pvalue > pvalue_tests))
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
                if (nrow(diff_normal_homoschedastic_data)) {
                    # For each signal...
                    for (m in rownames(diff_normal_homoschedastic_data)) {
                        ## Extract the column of the signal and the response variable...
                        mass_x <- temp_data_frame[[m]]
                        response_variable <- temp_data_frame[[non_signal_variable]]
                        ### Clean the NA values from the mass column and the response variable column
                        response_variable <- response_variable[!is.na(mass_x)]
                        mass_x <- mass_x[!is.na(mass_x)]
                        # Do not perform if the response variable has only one level...
                        if (length(levels(as.factor(response_variable))) > 1) {
                            ### Anova post hoc analysis: Tukey HSD test
                            anova_test <- aov(lm(formula = mass_x ~ response_variable))
                            # Post Hoc Tukey list creation
                            tukey_post <- TukeyHSD(anova_test)
                            rn <- rownames(tukey_post$response_variable)
                            anova_table <- cbind(rn, tukey_post$response_variable)
                            post_hoc_anova_list[[m]] <- anova_table
                        } else {
                            post_hoc_anova_list[[m]] <- NA
                        }
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
                        # Do not perform if the response variable has only one level...
                        if (length(levels(as.factor(response_variable))) > 1) {
                            ### Kruskal post hoc analysis: Nemenyi-Damico-Wolfe-Dunn test - coin required!
                            signal_data_frame <- data.frame(mass_x, response_variable)
                            # Post hoc test
                            NDWD <- oneway_test(mass_x ~ response_variable, data = signal_data_frame, ytrafo = function(data) trafo(data, numeric_trafo = rank), xtrafo = function(data) trafo(data, factor_trafo = function(mass_x) model.matrix(~mass_x - 1) %*% t(contrMat(table(mass_x), "Tukey"))), teststat = "max", distribution = approximate(B = 90000))
                            post_hoc_pvalue <- pvalue(NDWD, method = "single-step")
                            rn <- rownames(post_hoc_pvalue)
                            post_hoc_table <-  cbind(rn, post_hoc_pvalue)
                            post_hoc_kruskal_list[[m]] <- post_hoc_table
                        } else {
                            post_hoc_kruskal_list[[m]] <- NA
                        }
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
                        # Do not perform if the response variable has only one level...
                        if (length(levels(as.factor(response_variable))) > 1) {
                            ### Welch post hoc analysis: Duncan-Waller post-hoc test - agricolae required!
                            degrees_of_freedom <- df.residual(lm(mass_x ~ response_variable))
                            MSerror <- deviance(lm(mass_x ~ response_variable))/degrees_of_freedom
                            Fc <- anova(lm(mass_x ~ response_variable))[1,4]
                            comparison <- waller.test(mass_x, response_variable, degrees_of_freedom, MSerror, Fc, group = FALSE)
                            rn <- rownames(as.data.frame(comparison[4]))
                            ndf <- cbind(as.data.frame(comparison[4]), rn)
                            post_hoc_welch_list[[m]] <- ndf
                        } else {
                            post_hoc_welch_list[[m]] <- NA
                        }
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
                        # Do not perform if the response variable has only one level...
                        if (length(levels(as.factor(response_variable))) > 1) {
                            ### Permutation post hoc analysis: still follow Nemenyi-Damico-Wolfe-Dunn test - coin required!
                            signal_data_frame <- data.frame(mass_x, response_variable)
                            
                            NDWD <- oneway_test(mass_x ~ response_variable, data = signal_data_frame, ytrafo = function(data) trafo(data, numeric_trafo = rank), xtrafo = function(data) trafo(data, factor_trafo = function(mass_x) model.matrix(~mass_x - 1) %*% t(contrMat(table(mass_x), "Tukey"))), teststat = "max", distribution = approximate(B = 90000))
                            
                            post_hoc_pvalue  <- pvalue(NDWD, method = "single-step")
                            rn <- rownames(post_hoc_pvalue)
                            post_hoc_table <-  cbind(rn, post_hoc_pvalue)
                            post_hoc_permutation_list[[m]] <- post_hoc_table
                        } else {
                            post_hoc_permutation_list[[m]] <- NA
                        }
                    }
                }
                ##### Selected signals for inference
                selected_signals_for_inference <- c(rownames(diff_normal_homoschedastic_data), rownames(diff_normal_heteroschedastic_data), rownames(diff_non_normal_homoschedastic_data), rownames(diff_non_normal_heteroschedastic_data))
                # Print the message...
                if (length(selected_signals_for_inference) > 0) {
                    print("Differentially expressed signals for inference")
                    print(selected_signals_for_inference)
                } else {
                    print(paste("There are no differentially expressed signals among the levels of:", non_signal_variable, "or there is an insuffucient number of patients in at least one class or there is only one class!"))
                }
                # Generate a matrix with the signals for inference and the method used to identify them...
                inference_signals_method_matrix <- rbind(matrix_diff_norm_homo, matrix_diff_norm_hetero, matrix_diff_non_norm_homo, matrix_diff_non_norm_hetero)
                # Select the signals to be used for inference and the non-signal variable
                selected_signals_for_inference_intensity_df <- subset(temp_data_frame, select = c(non_signal_variable, selected_signals_for_inference))
                # Select the signals to be used for inference and the non-signal variable (NO OUTLIERS)
                if (remove_outliers_multi_level_effect_analysis == TRUE) {
                    selected_signals_for_inference_intensity_df_no_outliers <- subset(temp_data_frame_no_outliers, select = c(non_signal_variable, selected_signals_for_inference))
                }
                ########### Dump the files
                # OUTLIERS
                # Create the folder where to dump the files and go to it...
                plots_multi_level_effect_analysis_subfolder <- file.path(combination_subfolder, "Plots")
                dir.create(plots_multi_level_effect_analysis_subfolder)
                setwd(plots_multi_level_effect_analysis_subfolder)
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
                    jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5), ylab = "Signal intensity", xlab = non_signal_variable)
                    ggsave(jitter_plot, file = file_name, width = 4, height = 4)
                    ##### Box plot
                    plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
                    file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
                    box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot", ylab = "Signal intensity", xlab = non_signal_variable)
                    ggsave(box_plot, file = file_name, width = 4, height = 4)
                    ##### Scatter plot
                    plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
                    # Sort the dataframe rows according to the values of the non-signal variable
                    ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
                    ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
                    colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
                    file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
                    graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
                    scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors, ylab = "Signal intensity", xlab = "Signals (grouped)")
                    ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                }
                ### NO OUTLIERS
                if (remove_outliers_multi_level_effect_analysis == TRUE) {
                    # Create the folder where to dump the files and go to it...
                    plots_multi_level_effect_analysis_subfolder <- file.path(combination_subfolder, "Plots (No outliers)")
                    dir.create(plots_multi_level_effect_analysis_subfolder)
                    setwd(plots_multi_level_effect_analysis_subfolder)
                    # For each signals of inference...
                    for (s in selected_signals_for_inference) {
                        # Extract the intensity
                        signal_intensity <- selected_signals_for_inference_intensity_df_no_outliers[[s]]
                        # Extract the non-signal variable as an ordered factor
                        non_signal_as_ordered_factor <- ordered(temp_data_frame_no_outliers[, non_signal_variable])
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
                        jitter_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, geom = "jitter", main = s, alpha = I(1 / 5), ylab = "Signal intensity", xlab = non_signal_variable)
                        ggsave(jitter_plot, file = file_name, width = 4, height = 4)
                        ##### Box plot
                        plot_name <- sprintf("%s%s", non_signal_variable,"_boxplot") 
                        file_name <- sprintf("%s%s%s%s", s, " ", plot_name, image_format)
                        box_plot <- qplot(non_signal_as_ordered_factor, signal_intensity, data = signal_dataframe, main = s, geom = "boxplot", ylab = "Signal intensity", xlab = non_signal_variable)
                        ggsave(box_plot, file = file_name, width = 4, height = 4)
                        ##### Scatter plot
                        plot_name <- sprintf("%s%s", non_signal_variable,"_factor-Ordered_Spectrum__VS__Intensity") 
                        # Sort the dataframe rows according to the values of the non-signal variable
                        ordered_signal_dataframe <- signal_dataframe[order(non_signal_as_ordered_factor),]
                        ordered_signal_dataframe <- cbind(ordered_signal_dataframe, c(1:nrow(signal_dataframe[order(non_signal_as_ordered_factor),])))
                        colnames(ordered_signal_dataframe) <- c("old_id","intensity", non_signal_variable,"id")
                        file_name <- sprintf("%s%s%s%s", " ", s, plot_name, image_format)
                        graph_colors <- as.factor(ordered_signal_dataframe[[non_signal_variable]])
                        scatter_plot <- qplot(id, intensity, data = ordered_signal_dataframe, geom = "line", main = s, color = graph_colors, ylab = "Signal intensity", xlab = "Signals (grouped)")
                        ggsave(scatter_plot, file = file_name , width = 4, height = 4)
                    }
                }
                ### Dump the files
                # Create the folder where to dump the files and go to it...
                tables_multi_level_effect_analysis <- file.path(combination_subfolder, "Tables")
                dir.create(tables_multi_level_effect_analysis)
                setwd(tables_multi_level_effect_analysis)
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
                if (nrow(diff_non_normal_heteroschedastic_data)) {
                    write_posthoc_file(file_name = "PostHoc Permutation", data = post_hoc_permutation_list, file_format = file_format)
                }
                # Go back to the output folder
                setwd(output_folder)
            }
        }
        
        
        ##### FINISH
        finish_message <- tkmessageBox(title = "Operation completed", message = "All the statistical operations have been performed and the files have been dumped", icon = "info")
    } else {
        tkmessageBox(title = "No input file selected", message = "No input file has been selected!!!\nPlease, select a file to be imported", icon = "warning")
    }
}

























###############################################################################




##################################################################### WINDOW GUI

### Check for updates
check_for_updates_function()

########## List of variables, whose values are taken from the entries in the GUI
pvalue_expression <- tclVar("")
pvalue_tests <- tclVar("")
#TestPer_Base <- tclVar("")
#TestPer_adv <- tclVar("")
minimum_number_of_patients <- tclVar("")



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
    }
    # Windows 10
    if (length(grep("10", os_release, fixed = TRUE)) > 0) {
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
    }
    # Windows 10
    if (length(grep("10", os_release, fixed = TRUE)) > 0) {
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
        # Use them in the GUI
        title_font = ubuntu_title_bold
        label_font = ubuntu_other_normal
        entry_font = ubuntu_other_normal
        button_font = ubuntu_other_bold
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
window <- tktoplevel()
tktitle(window) <- "LC-MS PEAK STATISTICS"
#### Browse
select_input_button <- tkbutton(window, text="Import file...", command = file_import_function, font = button_font)
browse_output_button <- tkbutton(window, text="Browse output folder", command = browse_output_function, font = button_font)
#### Entries
output_file_type_export_entry <- tkbutton(window, text="Output file type", command = output_file_type_export_choice, font = button_font)
image_file_type_export_entry <- tkbutton(window, text="Image file type", command = image_file_type_export_choice, font = button_font)
data_record_entry <- tkbutton(window, text="Data REC", command = data_record_choice, font = button_font)
correlation_analysis_entry <- tkbutton(window, text="Correlation analysis", command = correlation_analysis_choice, font = button_font)
remove_outliers_correlation_analysis_entry <- tkbutton(window, text="Remove outliers\nCorrelation analysis", command = remove_outliers_correlation_analysis_choice, font = button_font)
two_level_effect_analysis_entry <- tkbutton(window, text="Two-level effect\nanalysis", command = two_level_effect_analysis_choice, font = button_font)
remove_outliers_two_level_effect_analysis_entry <- tkbutton(window, text="Remove outliers\nTwo-level effect\nanalysis", command = remove_outliers_two_level_effect_analysis_choice, font = button_font)
multi_level_effect_analysis_entry <- tkbutton(window, text="Multi-level effect\nanalysis", command = multi_level_effect_analysis_choice, font = button_font)
remove_outliers_multi_level_effect_analysis_entry <- tkbutton(window, text="Remove outliers\nMulti-level effect\nanalysis", command = remove_outliers_multi_level_effect_analysis_choice, font = button_font)
minimum_number_of_patients_label <- tklabel(window, text="Minimum number of patients", font = label_font)
minimum_number_of_patients_entry <- tkentry(window, width = 10, textvariable = minimum_number_of_patients, font = entry_font)
tkinsert(minimum_number_of_patients_entry, "end", "3")
pvalue_expression_label <- tklabel(window, text="p-value for signal\nexpression difference", font = label_font)
pvalue_expression_entry <- tkentry(window, width = 10, textvariable = pvalue_expression, font = entry_font)
tkinsert(pvalue_expression_entry, "end", "0.05")
pvalue_tests_label <- tklabel(window, text="p-value for significance\nin statistical tests", font = label_font)
pvalue_tests_entry <- tkentry(window, width = 10, textvariable = pvalue_tests, font = entry_font)
tkinsert(pvalue_tests_entry, "end", "0.05")
cumulative_class_in_two_level_effect_analysis_entry <- tkbutton(window, text="Cumulative class in the\ntwo-level effect analysis", command = cumulative_class_in_two_level_effect_analysis_choice, font = button_font)
#TestPer_Base_label <- tklabel(window, text="TestPer_Base", font = label_font)
#TestPer_Base_entry <- tkentry(window, width = 10, textvariable = TestPer_Base, font = entry_font)
#tkinsert(TestPer_Base_entry, "end", "0.17")
#TestPer_Adv_label <- tklabel(window, text="TestPer_Adv", font = label_font)
#TestPer_Adv_entry <- tkentry(window, width = 10, textvariable = TestPer_Adv, font = entry_font)
#tkinsert(TestPer_Adv_entry, "end", "0.19")
# Buttons
download_updates_button <- tkbutton(window, text="DOWNLOAD UPDATE", command = download_updates_function, font = button_font)
run_statistics_function_button <- tkbutton(window, text="RUN STATISTICS", command = run_statistics_function, font = button_font)
end_session_button <- tkbutton(window, text="QUIT", command = end_session_function, font = button_font)
#### Displaying labels
check_for_updates_value_label <- tklabel(window, text = check_for_updates_value, font = label_font)
output_file_type_export_value_label <- tklabel(window, text = output_file_type_export_value, font = label_font)
image_file_type_export_value_label <- tklabel(window, text = image_file_type_export_value, font = label_font)
data_record_value_label <- tklabel(window, text = data_record_value, font = label_font)
correlation_analysis_value_label <- tklabel(window, text = correlation_analysis_value, font = label_font)
remove_outliers_correlation_analysis_value_label <- tklabel(window, text = remove_outliers_correlation_analysis_value, font = label_font)
two_level_effect_analysis_value_label <- tklabel(window, text = two_level_effect_analysis_value, font = label_font)
remove_outliers_two_level_effect_analysis_value_label <- tklabel(window, text = remove_outliers_two_level_effect_analysis_value, font = label_font)
multi_level_effect_analysis_value_label <- tklabel(window, text = multi_level_effect_analysis_value, font = label_font)
remove_outliers_multi_level_effect_analysis_value_label <- tklabel(window, text = remove_outliers_multi_level_effect_analysis_value, font = label_font)
cumulative_class_in_two_level_effect_analysis_value_label <- tklabel(window, text = cumulative_class_in_two_level_effect_analysis_value, font = label_font)
#### Geometry manager
tkgrid(download_updates_button, row = 1, column = 2)
tkgrid(check_for_updates_value_label, row = 1, column = 3)
tkgrid(output_file_type_export_entry, row = 2, column = 1)
tkgrid(image_file_type_export_entry, row = 3, column = 1)
tkgrid(data_record_entry, row = 4, column = 1)
tkgrid(correlation_analysis_entry, row = 5, column = 1)
tkgrid(remove_outliers_correlation_analysis_entry, row = 6, column = 1)
tkgrid(two_level_effect_analysis_entry, row = 2, column = 3)
tkgrid(remove_outliers_two_level_effect_analysis_entry, row = 3, column = 3)
tkgrid(multi_level_effect_analysis_entry, row = 4, column = 3)
tkgrid(remove_outliers_multi_level_effect_analysis_entry, row = 5, column = 3)
tkgrid(output_file_type_export_value_label, row = 2, column = 2)
tkgrid(image_file_type_export_value_label, row = 3, column = 2)
tkgrid(data_record_value_label, row = 4, column = 2)
tkgrid(correlation_analysis_value_label, row = 5, column = 2)
tkgrid(remove_outliers_correlation_analysis_value_label, row = 6, column = 2)
tkgrid(two_level_effect_analysis_value_label, row = 2, column = 4)
tkgrid(remove_outliers_two_level_effect_analysis_value_label, row = 3, column = 4)
tkgrid(multi_level_effect_analysis_value_label, row = 4, column = 4)
tkgrid(remove_outliers_multi_level_effect_analysis_value_label, row = 5, column = 4)
tkgrid(minimum_number_of_patients_label, row = 6, column = 3)
tkgrid(minimum_number_of_patients_entry, row = 6, column = 4)
tkgrid(pvalue_expression_label, row = 7, column = 1)
tkgrid(pvalue_expression_entry, row = 7, column = 2)
tkgrid(pvalue_tests_label, row = 8, column = 1)
tkgrid(pvalue_tests_entry, row = 8, column = 2)
tkgrid(cumulative_class_in_two_level_effect_analysis_entry, row = 7, column = 3)
tkgrid(cumulative_class_in_two_level_effect_analysis_value_label, row = 7, column = 4)
#tkgrid(TestPer_Base_label, row = 7, column = 3)
#tkgrid(TestPer_Base_entry, row = 7, column = 4)
#tkgrid(TestPer_Adv_label, row = 8, column = 3)
#tkgrid(TestPer_Adv_entry, row = 8, column = 4)
tkgrid(browse_output_button, row = 9, column = 1)
tkgrid(select_input_button, row = 9, column = 2)
tkgrid(run_statistics_function_button, row = 9, column = 3)
tkgrid(end_session_button, row = 9, column = 4)
