#################### MASCOT OUTPUT MOD-PROCESSER ####################







### Program version (Specified by the program writer!!!!)
R_script_version <- "2017.03.29.0"
### GitHub URL where the R file is
github_R_url <- ""
### Name of the file when downloaded
script_file_name <- "MASCOT OUTPUT MOD-PROCESSER"
# Change log
change_log <- "1. New software!"






########## FUNCTIONS
# Check internet connection
check_internet_connection <- function(method = "getURL", website_to_ping = "www.google.it") {
    ##### PING
    if (method == "ping") {
        if (Sys.info()[1] == "Linux") {
            # -c: number of packets sent/received (attempts) ; -W timeout in seconds
            there_is_internet <- !as.logical(system(command = paste("ping -c 1 -W 2", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        } else if (Sys.info()[1] == "Windows") {
            # -n: number of packets sent/received (attempts) ; -w timeout in milliseconds
            there_is_internet <- !as.logical(system(command = paste("ping -n 1 -w 2000", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        } else {
            there_is_internet <- !as.logical(system(command = paste("ping", website_to_ping), intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE))
        }
    } else if (method == "getURL") {
        ##### GET URL
        # Install the RCurl package if not installed
        if ("RCurl" %in% installed.packages()[,1]) {
            library(RCurl)
        } else {
            install.packages("RCurl", repos = "http://cran.mirror.garr.it/mirrors/CRAN/", quiet = TRUE, verbose = FALSE)
            library(RCurl)
        }
        there_is_internet <- FALSE
        try({
            there_is_internet <- is.character(getURL(u = website_to_ping, followLocation = TRUE, .opts = list(timeout = 1, maxredirs = 2, verbose = FALSE)))
        }, silent = TRUE)
    }
    return(there_is_internet)
}

# Install and load required packages
install_and_load_required_packages <- function(required_packages, repository = "http://cran.mirror.garr.it/mirrors/CRAN/") {
    ### Check internet connection
    there_is_internet <- check_internet_connection(method = "getURL", website_to_ping = "www.google.it")
    ########## Update all the packages (if there is internet connection)
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
        print("All the packages are up-to-date")
    }
    ##### Load the packages (if there are all the packages)
    if ((length(missing_packages) > 0 && there_is_internet == TRUE) || length(missing_packages) == 0) {
        for (i in 1:length(required_packages)) {
            library(required_packages[i], character.only = TRUE)
        }
    } else {
        print("Packages cannot be loaded... Expect issues...")
    }
}

install_and_load_required_packages(c("tcltk", "plyr"))





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
            download.file(url = github_R_url, destfile = paste(script_file_name, " (", online_version_number, ")", sep = ""), method = "auto")
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
    filepath_import_select <- tkmessageBox(title = "Input file", message = "Select the Mascot output file", icon = "info")
    input_file <- tclvalue(tkgetOpenFile(filetypes = "{{Microsoft Excel files} {.xls .xlsx}} {{Comma Separated Value files} {.csv}}"))
    if (!nchar(input_file)) {
        tkmessageBox(message = "No file selected")
    } else {
        tkmessageBox(message = paste("The following file will be read:", input_file))
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
        # Escape the function
        .GlobalEnv$input_file <- input_file
        .GlobalEnv$input_data <- input_data
        tkmessageBox(title = "File imported", message = "The file has been successfully imported!", icon = "info")
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
    tkmessageBox(message = "A sub-directory named 'MASCOT X' will be created for each run!")
    # Escape the function
    .GlobalEnv$output_folder <- output_folder
}

##### Exit
end_session_function <- function() {
    q(save = "no")
}

##### Run the statistics
run_mascot_output_processer_function <- function() {
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
                    # Add the number to the list of MASCOT numbers
                    try(MASCOT_present_folder_numbers <- append(MASCOT_present_folder_numbers, as.integer(MASCOT_present_folder_split[2])))
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
        non_modified_peptides_df <- non_modified_peptides_df[order(non_modified_peptides_df$identity, decreasing = TRUE), ]
        
        
        ########## REMOVE DUPLICATES (MOD)
        modified_peptides_df$prot_acc <- as.character(modified_peptides_df$prot_acc)
        modified_peptides_df$pep_seq <- as.character(modified_peptides_df$pep_seq)
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df <- modified_peptides_df[order(modified_peptides_df$prot_acc, modified_peptides_df$pep_seq, modified_peptides_df$pep_var_mod, modified_peptides_df$pep_var_mod_pos, -modified_peptides_df$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc, pep_seq, pep_var_mod and pep_var_mod_pos and only the highest score for prot_acc/pep_seq/pep_var_mod/pep_var_mod_pos duplicates)
        modified_peptides_df <- ddply(.data = modified_peptides_df, .variables = c("prot_acc", "pep_seq", "pep_var_mod", "pep_var_mod_pos"), .fun = head, n = 1)
        
        # Display the 'identity' forst
        modified_peptides_df <- modified_peptides_df[order(modified_peptides_df$identity, decreasing = TRUE), ]
        
        
        ########## REMOVE SEQUENCE DUPLICATES (MOD)
        modified_peptides_df_sequences <- modified_peptides_df
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df_sequences <- modified_peptides_df_sequences[order(modified_peptides_df_sequences$prot_acc, modified_peptides_df_sequences$pep_seq, -modified_peptides_df_sequences$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc and pep_seq, and only the highest score for prot_acc/pep_seq duplicates)
        modified_peptides_df_sequences <- ddply(.data = modified_peptides_df_sequences, .variables = c("prot_acc", "pep_seq"), .fun = head, n = 1)
        
        # Display the 'identity' forst
        modified_peptides_df_sequences <- modified_peptides_df_sequences[order(modified_peptides_df_sequences$identity, decreasing = TRUE), ]
        
    }
}
