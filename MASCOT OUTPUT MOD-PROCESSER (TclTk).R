#################### MASCOT OUTPUT MOD-PROCESSER ####################







### Program version (Specified by the program writer!!!!)
R_script_version <- "2017.04.27.0"
### GitHub URL where the R file is
github_R_url <- "https://raw.githubusercontent.com/gmanuel89/Public-R-UNIMIB/master/MASCOT%20OUTPUT%20MOD-PROCESSER.R"
### Name of the file when downloaded
script_file_name <- "MASCOT OUTPUT MOD-PROCESSER"
# Change log
change_log <- "1. New software!\n2. Bugfix"






########## FUNCTIONS
# Check internet connection
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

# Install and load required packages
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

install_and_load_required_packages(c("tcltk", "plyr"), update_packages = TRUE)





###################################### Initialize the variables (default values)
input_file <- ""
output_folder <- getwd()
output_format <- "Comma Separated Values (.csv)"
file_format <- "csv"
image_format <- ".png"





################## Values of the variables (for displaying and dumping purposes)
output_file_type_export_value <- "Comma Separated Values (.csv)"
image_file_type_export_value <- "PNG (.png)"
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
                    check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdate available:\n", online_version_number, sep = "")
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
            tkmessageBox(title = "Updated file downloaded!", message = paste("The updated script, named:\n\n", paste(script_file_name, " (", online_version_number, ").R", sep = ""), "\n\nhas been downloaded to:\n\n", download_folder, "\n\nClose everything, delete this file and run the script from the new file!", sep = ""), icon = "info")
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
    output_file_type_export_value_label <- tklabel(window, text = output_format, font = label_font, bg = "white", width = 30)
    tkgrid(output_file_type_export_value_label, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
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
    image_file_type_export_value_label <- tklabel(window, text = image_output_format, font = label_font, bg = "white", width = 20)
    tkgrid(image_file_type_export_value_label, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    .GlobalEnv$image_output_format <- image_output_format
    .GlobalEnv$image_format <- image_format
}

##### File import
file_import_function <- function() {
    filepath_import_select <- tkmessageBox(title = "Input file", message = "Select the Mascot output file", icon = "info")
    input_file <- tclvalue(tkgetOpenFile(filetypes = "{{Comma Separated Value files} {.csv}}"))
    if (!nchar(input_file)) {
        tkmessageBox(message = "No file selected")
    } else {
        tkmessageBox(message = paste("The following file will be read:\n\n", input_file))
    }
    if (input_file != "") {
        #################### IMPORT THE DATA FROM THE FILE
        ##### XLS or XLSX
        if (length(grep(".xls", input_file, fixed = TRUE)) > 0 || length(grep(".xlsx", input_file, fixed = TRUE)) > 0) {
            input_data <- readWorksheetFromFile(input_file, sheet = 1)
            input_data <- as.data.frame(input_data)
        } else if (length(grep(".csv", input_file, fixed = TRUE)) > 0) {
            ##### CSV
            input_data <- read.csv(input_file, header = TRUE, sep = ",")
        }
        ### Retrieve the input file name
        input_filename <- NULL
        try({
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                input_filename <- unlist(strsplit(input_file, "/"))
                input_filename <- input_filename[length(input_filename)]
                input_filename <- unlist(strsplit(input_filename, ".", fixed = TRUE))[1]
            } else if (Sys.info()[1] == "Windows") {
                input_filename <- unlist(strsplit(input_file, "\\\\"))
                input_filename <- input_filename[length(input_filename)]
                input_filename <- unlist(strsplit(input_filename, ".", fixed = TRUE))[1]
            }
        }, silent = TRUE)
        # Escape the function
        .GlobalEnv$input_file <- input_file
        .GlobalEnv$input_filename <- input_filename
        .GlobalEnv$input_data <- input_data
        tkmessageBox(title = "File imported", message = "The file has been successfully imported!", icon = "info")
    } else {
        # Escape the function
        .GlobalEnv$input_file <- input_file
        .GlobalEnv$input_filename <- NULL
        .GlobalEnv$input_data <- NULL
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
    tkmessageBox(message = "A sub-directory named 'MASCOT X' will be created for each run!")
    # Escape the function
    .GlobalEnv$output_folder <- output_folder
}

##### Exit
end_session_function <- function() {
    q(save = "no")
}

##### Run the mod-processer
run_mascot_output_modprocesser_function <- function() {
    ### Progress bar
    program_progress_bar <- tkProgressBar(title = "Calculating...", label = "", min = 0, max = 1, initial = 0, width = 300)
    setTkProgressBar(program_progress_bar, value = 0, title = NULL, label = "0 %")
    # Go to the working directory
    setwd(output_folder)
    if (input_file != "") {
        ##### Automatically create a subfolder with all the results
        ## Check if such subfolder exists
        list_of_directories <- list.dirs(output_folder, full.names = FALSE, recursive = FALSE)
        ## Check the presence of a MASCOT folder
        MASCOT_folder_presence <- FALSE
        if (length(list_of_directories) > 0) {
            for (dr in 1:length(list_of_directories)) {
                if (length(grep("MASCOT", list_of_directories[dr], fixed = TRUE)) > 0) {
                    MASCOT_folder_presence <- TRUE
                }
            }
        }
        ## If it present...
        if (isTRUE(MASCOT_folder_presence)) {
            ## Extract the number after the STATISTICS
            # Number for the newly created folder
            MASCOT_new_folder_number <- 0
            # List of already present numbers
            MASCOT_present_folder_numbers <- integer()
            # For each folder present...
            for (dr in 1:length(list_of_directories)) {
                # If it is named MASCOT
                if (length(grep("MASCOT", list_of_directories[dr], fixed = TRUE)) > 0) {
                    # Split the name
                    MASCOT_present_folder_split <- unlist(strsplit(list_of_directories[dr], "MASCOT"))
                    # Add the number to the list of MASCOT numbers (if it is NA it means that what was following the MASCOT was not a number)
                    try({
                        if (!is.na(as.integer(MASCOT_present_folder_split[2]))) {
                            MASCOT_present_folder_numbers <- append(MASCOT_present_folder_numbers, as.integer(MASCOT_present_folder_split[2]))
                        } else {
                            MASCOT_present_folder_numbers <- append(MASCOT_present_folder_numbers, as.integer(0))
                        }
                    }, silent = TRUE)
                }
            }
            # Sort the STATISTICS folder numbers
            try(MASCOT_present_folder_numbers <- sort(MASCOT_present_folder_numbers))
            # The new folder number will be the greater + 1
            try(MASCOT_new_folder_number <- MASCOT_present_folder_numbers[length(MASCOT_present_folder_numbers)] + 1)
            # Generate the new subfolder
            subfolder <- paste("MASCOT", MASCOT_new_folder_number)
            # Create the subfolder
            dir.create(file.path(output_folder, subfolder))
            # Estimate the new output folder
            output_folder <- file.path(output_folder, subfolder)
        } else {
            # If it not present...
            # Create the folder where to dump the files and go to it...
            subfolder <- paste("MASCOT", "1")
            # Create the subfolder
            dir.create(file.path(output_folder, subfolder))
            # Estimate the new output folder
            output_folder <- file.path(output_folder, subfolder)
        }
        # Go to the new working directory
        setwd(output_folder)
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.05, title = NULL, label = "5 %")
        
        
        ########## DATA SPLIT: MODIFIED vs NON-MODIFIED
        # Convert the pep_var_mod column to character
        input_data$pep_var_mod <- as.character(input_data$pep_var_mod)
        # The modified peptides have something in this column
        modified_peptides_df <- input_data[input_data$pep_var_mod != "", ]
        # The unmodified peptides do not have anything as modification
        non_modified_peptides_df <- input_data[input_data$pep_var_mod == "", ]
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.10, title = NULL, label = "10 %")
        
        
        ########## IDENTITY
        # Insert a column named "identity", which contains the information about the identity or the homology of the peptides, according to their scores...
        # If pep_score is more than the pep_ident, the peptide is flagged as "identity". If this condition is not satisfied check the pep_homol: if pep_homol is present (not NA or > 0) and the pep_score is more than the pep_homol, flag the peptide as "homology". If pep_homol is not present (NA or 0) or the pep_score is not more than the pep_homol, "discard" the peptide.
        
        non_modified_peptides_df$identity <- non_modified_peptides_df$pep_var_mod
        modified_peptides_df$identity <- modified_peptides_df$pep_var_mod
        
        for (l in 1:length(non_modified_peptides_df$identity)) {
            if (non_modified_peptides_df$pep_score[l] >= non_modified_peptides_df$pep_ident[l]) {
                non_modified_peptides_df$identity[l] <- "identity"
            } else {
                if (!is.na(non_modified_peptides_df$pep_homol[l]) && (non_modified_peptides_df$pep_score[l] > non_modified_peptides_df$pep_homol[l])) {
                    non_modified_peptides_df$identity[l] <- "homology"
                } else {
                    non_modified_peptides_df$identity[l] <- "discard"
                }
            }
        }
        
        for (l in 1:length(modified_peptides_df$identity)) {
            if (modified_peptides_df$pep_score[l] >= modified_peptides_df$pep_ident[l]) {
                modified_peptides_df$identity[l] <- "identity"
            } else {
                if (!is.na(modified_peptides_df$pep_homol[l]) && (modified_peptides_df$pep_score[l] > modified_peptides_df$pep_homol[l])) {
                    modified_peptides_df$identity[l] <- "homology"
                } else {
                    modified_peptides_df$identity[l] <- "discard"
                }
            }
        }
        
        # Discard the peptides to be discarded
        non_modified_peptides_df <- non_modified_peptides_df[non_modified_peptides_df$identity != "discard", ]
        modified_peptides_df <- modified_peptides_df[modified_peptides_df$identity != "discard", ]
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.20, title = NULL, label = "20 %")
        
        
        ########## REMOVE DUPLICATES (NON-MOD)
        non_modified_peptides_df$prot_acc <- as.character(non_modified_peptides_df$prot_acc)
        non_modified_peptides_df$pep_seq <- as.character(non_modified_peptides_df$pep_seq)
        
        # Order the dataframe according to prot_acc, pep_seq, pep_score
        non_modified_peptides_df <- non_modified_peptides_df[order(non_modified_peptides_df$prot_acc, non_modified_peptides_df$pep_seq, -non_modified_peptides_df$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc and pep_seq and only the highest score for prot_acc/pep_seq duplicates)
        non_modified_peptides_df <- ddply(.data = non_modified_peptides_df, .variables = c("prot_acc", "pep_seq"), .fun = head, n = 1)
        
        # Display the 'identity' forst
        non_modified_peptides_df <- non_modified_peptides_df[order(non_modified_peptides_df$pep_score, non_modified_peptides_df$identity, non_modified_peptides_df$prot_acc, decreasing = TRUE), ]
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.40, title = NULL, label = "40 %")
        
        
        ########## REMOVE DUPLICATES (MOD)
        modified_peptides_df$prot_acc <- as.character(modified_peptides_df$prot_acc)
        modified_peptides_df$pep_seq <- as.character(modified_peptides_df$pep_seq)
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df <- modified_peptides_df[order(modified_peptides_df$prot_acc, modified_peptides_df$pep_seq, modified_peptides_df$pep_var_mod, modified_peptides_df$pep_var_mod_pos, -modified_peptides_df$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc, pep_seq, pep_var_mod and pep_var_mod_pos and only the highest score for prot_acc/pep_seq/pep_var_mod/pep_var_mod_pos duplicates)
        modified_peptides_df <- ddply(.data = modified_peptides_df, .variables = c("prot_acc", "pep_seq", "pep_var_mod", "pep_var_mod_pos"), .fun = head, n = 1)
        
        # Display the 'identity' forst
        modified_peptides_df <- modified_peptides_df[order(modified_peptides_df$pep_score, modified_peptides_df$identity, modified_peptides_df$prot_acc, decreasing = TRUE), ]
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.60, title = NULL, label = "60 %")
        
        
        ########## REMOVE SEQUENCE DUPLICATES (MOD)
        modified_peptides_df_sequences <- modified_peptides_df
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df_sequences <- modified_peptides_df_sequences[order(modified_peptides_df_sequences$prot_acc, modified_peptides_df_sequences$pep_seq, -modified_peptides_df_sequences$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc and pep_seq, and only the highest score for prot_acc/pep_seq duplicates)
        modified_peptides_df_sequences <- ddply(.data = modified_peptides_df_sequences, .variables = c("prot_acc", "pep_seq"), .fun = head, n = 1)
        
        # Display the 'identity' forst
        modified_peptides_df_sequences <- modified_peptides_df_sequences[order(modified_peptides_df_sequences$pep_score, modified_peptides_df_sequences$identity, modified_peptides_df_sequences$prot_acc, decreasing = TRUE), ]
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.80, title = NULL, label = "80 %")
        
        
        ########## SAVE FILES
        if (is.null(input_filename)) {
            input_filename <- "Input Data"
        }
        if (file_format == "csv") {
            write.csv(input_data, file = paste(input_filename, ".", file_format, sep = ""), row.names = FALSE)
            write.csv(non_modified_peptides_df, file = paste("Non-modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df, file = paste("Modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df_sequences, file = paste("Modified peptides (unique sequences).", file_format, sep = ""), row.names = FALSE)
        } else if (file_format == "xlsx" || file_format == "xls") {
            wb = loadWorkbook(filename = paste(input_filename, ".", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Input data")
            writeWorksheet(wb, data = input_data, sheet = "Input data", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Non-modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Non-modified peptides")
            writeWorksheet(wb, data = non_modified_peptides_df, sheet = "Non-modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified peptides")
            writeWorksheet(wb, data = modified_peptides_df, sheet = "Modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides (unique sequences).", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified sequences")
            writeWorksheet(wb, data = modified_peptides_df, sheet = "Modified sequences", header = TRUE)
            saveWorkbook(wb)
        }
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 1.00, title = NULL, label = "100 %")
        close(program_progress_bar)
        
        
        tkmessageBox(title = "Success!", message = "The process has successfully finished!", icon = "info")
    } else {
        ########## NO INPUT FILE
        tkmessageBox(title = "No input file selected!", message = "No input file has been selected!", icon = "warning")
    }
}











###############################################################################




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
tktitle(window) <- "MASCOT OUTPUT MOD-PROCESSER"
# Title label
title_label <- tklabel(window, text = "MASCOT\nOUTPUT\nMOD-PROCESSER", font = title_font, bg = "white")
#### Browse
select_input_button <- tkbutton(window, text="IMPORT FILE...", command = file_import_function, font = button_font, bg = "white", width = 20)
browse_output_button <- tkbutton(window, text="BROWSE\nOUTPUT FOLDER...", command = browse_output_function, font = button_font, bg = "white", width = 20)
#### Entries
output_file_type_export_entry <- tkbutton(window, text="Output\nfile type", command = output_file_type_export_choice, font = button_font, bg = "white", width = 20)
# Buttons
download_updates_button <- tkbutton(window, text="DOWNLOAD\nUPDATE...", command = download_updates_function, font = button_font, bg = "white", width = 20)
run_mascot_output_modprocesser_function_button <- tkbutton(window, text = "RUN MASCOT OUTPUT\nMOD-PROCESSER", command = run_mascot_output_modprocesser_function, font = button_font, bg = "white", width = 20)
end_session_button <- tkbutton(window, text="QUIT", command = end_session_function, font = button_font, bg = "white", width = 20)
#### Displaying labels
check_for_updates_value_label <- tklabel(window, text = check_for_updates_value, font = label_font, bg = "white", width = 20)
output_file_type_export_value_label <- tklabel(window, text = output_file_type_export_value, font = label_font, bg = "white", width = 30)

#### Geometry manager
tkgrid(title_label, row = 1, column = 1, padx = c(20, 20), pady = c(20, 20))
tkgrid(download_updates_button, row = 1, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(check_for_updates_value_label, row = 1, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(output_file_type_export_entry, row = 2, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(output_file_type_export_value_label, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(browse_output_button, row = 3, column = 1, padx = c(10, 10), pady = c(10, 10))
tkgrid(select_input_button, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
tkgrid(run_mascot_output_modprocesser_function_button, row = 3, column = 3, padx = c(10, 10), pady = c(10, 10))
tkgrid(end_session_button, row = 4, column = 2, padx = c(10, 10), pady = c(10, 10))

