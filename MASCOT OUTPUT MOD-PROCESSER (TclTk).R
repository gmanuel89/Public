#################### MASCOT OUTPUT MOD-PROCESSER ####################


# Clear the console
#cat("\014")
# Empty the workspace
rm(list = ls())





### Program version (Specified by the program writer!!!!)
R_script_version <- "2017.05.19.0"
### GitHub URL where the R file is
github_R_url <- "https://raw.githubusercontent.com/gmanuel89/Mascot-Output-MOD-Processer/master/MASCOT%20OUTPUT%20MOD-PROCESSER.R"
### Name of the file when downloaded
script_file_name <- "MASCOT OUTPUT MOD-PROCESSER"
# Change log
change_log <- "1. Bugfix"






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

install_and_load_required_packages("tcltk", update_packages = TRUE)





###################################### Initialize the variables (default values)
input_folder <- ""
output_folder <- getwd()
output_format <- "Comma Separated Values (.csv)"
file_format <- "csv"
image_format <- ".png"
class_name <- "CLASS"





################## Values of the variables (for displaying and dumping purposes)
output_file_type_export_value <- "Comma Separated Values (.csv)"
image_file_type_export_value <- "PNG (.png)"






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
                check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked:\nconnection problems", sep = "")
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
    # Folder for the input files
    tkmessageBox(title = "Input files", message = "Select the folder containing the Mascot output files, one relative to the whole proteome/peptidome and the other one relative to the proteome/peptidome which underwent modification.\n\nThe first file's name should start with 'PROT' while the second file's name should start with 'MOD', and both should be in the CSV format", icon = "info")
    input_folder <- tclvalue(tkchooseDirectory())
    setwd(input_folder)
    # Read the CSV files
    input_file_paths <- list.files(path = input_folder, full.names = FALSE, pattern = ".csv", recursive = TRUE)
    # Create a list in which there is the PROT file and the MOD file
    input_file_list <- list()
    input_file_list[["PROT"]] <- input_file_paths[which(startsWith(input_file_paths, "PROT"))]
    input_file_list[["MOD"]] <- input_file_paths[which(startsWith(input_file_paths, "MOD"))]
    # Import data only if a file path is specified
    if (input_folder != "" && length(input_file_list) == 2) {
        #################### IMPORT THE DATA FROM THE FILE
        ### Progress bar
        import_progress_bar <- tkProgressBar(title = "Importing file...", label = "", min = 0, max = 1, initial = 0, width = 300)
        setTkProgressBar(import_progress_bar, value = 0, title = NULL, label = "0 %")
        ### Remove all the lines that are not the desired lines (the useless header)
        setTkProgressBar(import_progress_bar, value = 0.10, title = "Discarding header...", label = "10 %")
        input_file_lines <- list()
        input_file_lines[["PROT"]] <- readLines(input_file_list[["PROT"]])
        input_file_lines[["MOD"]] <- readLines(input_file_list[["MOD"]])
        # Start to read from the matrix header: "prot_hit"
        input_file_header_line_number <- list(PROT = 0, MOD = 0)
        for (l in 1:length(input_file_lines[["PROT"]])) {
            if (startsWith(input_file_lines[["PROT"]][l], "prot_hit")) {
                input_file_header_line_number[["PROT"]] <- l
                break
            }
        }
        for (l in 1:length(input_file_lines[["MOD"]])) {
            if (startsWith(input_file_lines[["MOD"]][l], "prot_hit")) {
                input_file_header_line_number[["MOD"]] <- l
                break
            }
        }
        # Keep only the selected lines
        final_input_file_lines <- list()
        if (input_file_header_line_number[["PROT"]] > 0) {
            final_input_file_lines[["PROT"]] <- input_file_lines[["PROT"]][input_file_header_line_number[["PROT"]] : length(input_file_lines[["PROT"]])]
        } else {
            final_input_file_lines[["PROT"]] <- character()
        }
        if (input_file_header_line_number[["MOD"]] > 0) {
            final_input_file_lines[["MOD"]] <- input_file_lines[["MOD"]][input_file_header_line_number[["MOD"]] : length(input_file_lines[["MOD"]])]
        } else {
            final_input_file_lines[["MOD"]] <- character()
        }
        setTkProgressBar(import_progress_bar, value = 0.40, title = "Reading CSV table...", label = "40 %")
        # Read the CSV file from the lines
        input_data <- list()
        if (length(final_input_file_lines[["PROT"]]) > 0) {
            input_data[["PROT"]] <- read.csv(text = final_input_file_lines[["PROT"]], header = TRUE)
            if (is.null(rownames(input_data[["PROT"]]))) {
                rownames(input_data[["PROT"]]) <- seq(from = 1, to = nrow(input_data[["PROT"]]), by = 1)
            }
        } else {
            input_data[["PROT"]] <- NULL
        }
        if (length(final_input_file_lines[["MOD"]]) > 0) {
            input_data[["MOD"]] <- read.csv(text = final_input_file_lines[["MOD"]], header = TRUE)
            if (is.null(rownames(input_data[["MOD"]]))) {
                rownames(input_data[["MOD"]]) <- seq(from = 1, to = nrow(input_data[["MOD"]]), by = 1)
            }
        } else {
            input_data[["MOD"]] <- NULL
        }
        setTkProgressBar(import_progress_bar, value = 0.95, title = "Setting file name...", label = "95 %")
        ### Retrieve the input file name
        input_filename <- list(PROT = NULL, MOD = NULL)
        try({
            if (Sys.info()[1] == "Linux" || Sys.info()[1] == "Darwin") {
                input_filename[["PROT"]] <- unlist(strsplit(input_file_list[["PROT"]], "/"))
                input_filename[["PROT"]] <- input_filename[["PROT"]][length(input_filename[["PROT"]])]
                input_filename[["PROT"]] <- unlist(strsplit(input_filename[["PROT"]], ".", fixed = TRUE))[1]
                input_filename[["MOD"]] <- unlist(strsplit(input_file_list[["MOD"]], "/"))
                input_filename[["MOD"]] <- input_filename[["MOD"]][length(input_filename[["MOD"]])]
                input_filename[["MOD"]] <- unlist(strsplit(input_filename[["MOD"]], ".", fixed = TRUE))[1]
            } else if (Sys.info()[1] == "Windows") {
                input_filename[["PROT"]] <- unlist(strsplit(input_file_list[["PROT"]], "\\\\"))
                input_filename[["PROT"]] <- input_filename[["PROT"]][length(input_filename[["PROT"]])]
                input_filename[["PROT"]] <- unlist(strsplit(input_filename[["PROT"]], ".", fixed = TRUE))[1]
                input_filename[["MOD"]] <- unlist(strsplit(input_file_list[["MOD"]], "\\\\"))
                input_filename[["MOD"]] <- input_filename[["MOD"]][length(input_filename[["MOD"]])]
                input_filename[["MOD"]] <- unlist(strsplit(input_filename[["MOD"]], ".", fixed = TRUE))[1]
            }
        }, silent = TRUE)
        # Input folder as class name
        class_name <- unlist(strsplit(input_folder, .Platform$file.sep))
        class_name <- class_name[length(class_name)]
        # Escape the function
        .GlobalEnv$input_folder <- input_folder
        .GlobalEnv$class_name <- class_name
        .GlobalEnv$input_file_list <- input_file_list
        .GlobalEnv$input_filename <- input_filename
        .GlobalEnv$input_data <- input_data
        setTkProgressBar(import_progress_bar, value = 1.00, title = "Done!", label = "100 %")
        close(import_progress_bar)
        # Go back to the output folder
        setwd(output_folder)
        # Message
        if (!is.null(input_data[["PROT"]]) && !is.null(input_data[["MOD"]])) {
            tkmessageBox(title = "File imported", message = "The file has been successfully imported!", icon = "info")
        } else {
            tkmessageBox(title = "File not imported", message = "The file could not be imported! Check the file structure and re-import it!", icon = "warning")
        }
    } else {
        # Escape the function
        .GlobalEnv$input_file_list <- input_file_list
        .GlobalEnv$input_filename[["PROT"]] <- NULL
        .GlobalEnv$input_filename[["MOD"]] <- NULL
        .GlobalEnv$input_data[["PROT"]] <- NULL
        .GlobalEnv$input_data[["MOD"]] <- NULL
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
    program_progress_bar <- tkProgressBar(title = "Reading data...", label = "", min = 0, max = 1, initial = 0, width = 300)
    setTkProgressBar(program_progress_bar, value = 0.01, title = "Reading data...", label = "1 %")
    # Go to the working directory
    setwd(output_folder)
    # Check if all the columns are present
    columns_needed <- c("pep_var_mod", "pep_ident", "pep_homol", "pep_score", "prot_acc", "pep_seq", "pep_var_mod_pos")
    all_the_columns_are_present <- list()
    all_the_columns_are_present[["PROT"]] <- all(columns_needed %in% colnames(input_data[["PROT"]]))
    all_the_columns_are_present[["MOD"]] <- all(columns_needed %in% colnames(input_data[["MOD"]]))
    all_the_columns_are_present <- all(unlist(all_the_columns_are_present))
    missing_columns_value <- NULL
    missing_columns <- list()
    missing_columns[["PROT"]] <- columns_needed[columns_needed %in% colnames(input_data[["PROT"]]) == FALSE]
    missing_columns[["MOD"]] <- columns_needed[columns_needed %in% colnames(input_data[["MOD"]]) == FALSE]
    missing_columns <- unlist(missing_columns)
    if (length(missing_columns) > 0) {
        for (m in 1:length(missing_columns)) {
            if (is.null(missing_columns_value)) {
                missing_columns_value <- missing_columns[m]
            } else {
                missing_columns_value <- paste0(missing_columns_value, "\n", missing_columns[m])
            }
        }
    }
    
    
    
    
    
    # Run only if there is data and the data is correct
    if (input_folder != "" && length(input_file_list) == 2 && !is.null(unlist(input_data)) && all_the_columns_are_present == TRUE) {
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
            MASCOT_subfolder <- paste("MASCOT", MASCOT_new_folder_number)
            # Estimate the new output folder
            output_subfolder <- file.path(output_folder, MASCOT_subfolder)
            # Create the subfolder
            dir.create(output_subfolder)
        } else {
            # If it not present...
            # Create the folder where to dump the files and go to it...
            MASCOT_subfolder <- paste("MASCOT", "1")
            # Estimate the new output folder
            output_subfolder <- file.path(output_folder, MASCOT_subfolder)
            # Create the subfolder
            dir.create(output_subfolder)
        }
        # Go to the new working directory
        setwd(output_subfolder)
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.05, title = "Splitting data...", label = "5 %")
        
        
        
        
        
        ########## DATA SPLIT: MODIFIED vs NON-MODIFIED
        # Convert the pep_var_mod column to character
        input_data[["PROT"]]$pep_var_mod <- as.character(input_data[["PROT"]]$pep_var_mod)
        input_data[["MOD"]]$pep_var_mod <- as.character(input_data[["MOD"]]$pep_var_mod)
        # The modified peptides have something in this column
        modified_peptides_df <- list()
        modified_peptides_df[["PROT"]] <- input_data[["PROT"]][input_data[["PROT"]]$pep_var_mod != "", ]
        modified_peptides_df[["MOD"]] <- input_data[["MOD"]][input_data[["MOD"]]$pep_var_mod != "", ]
        # The unmodified peptides do not have anything as modification
        non_modified_peptides_df <- list()
        non_modified_peptides_df[["PROT"]] <- input_data[["PROT"]][input_data[["PROT"]]$pep_var_mod == "", ]
        non_modified_peptides_df[["MOD"]] <- input_data[["MOD"]][input_data[["MOD"]]$pep_var_mod == "", ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.10, title = "Determining homology and identity...", label = "10 %")
        
        
        
        
        
        ########## IDENTITY
        # Insert a column named "identity", which contains the information about the identity or the homology of the peptides, according to their scores...
        # If pep_score is more than the pep_ident, the peptide is flagged as "identity". If this condition is not satisfied check the pep_homol: if pep_homol is present (not NA or > 0) and the pep_score is more than the pep_homol, flag the peptide as "homology". If pep_homol is not present (NA or 0) or the pep_score is not more than the pep_homol, "discard" the peptide.
        
        non_modified_peptides_df[["PROT"]]$identity <- non_modified_peptides_df[["PROT"]]$pep_var_mod
        non_modified_peptides_df[["MOD"]]$identity <- non_modified_peptides_df[["MOD"]]$pep_var_mod
        modified_peptides_df[["PROT"]]$identity <- modified_peptides_df[["PROT"]]$pep_var_mod
        modified_peptides_df[["MOD"]]$identity <- modified_peptides_df[["MOD"]]$pep_var_mod
        
        for (l in 1:length(non_modified_peptides_df[["PROT"]]$identity)) {
            if (non_modified_peptides_df[["PROT"]]$pep_score[l] >= non_modified_peptides_df[["PROT"]]$pep_ident[l]) {
                non_modified_peptides_df[["PROT"]]$identity[l] <- "identity"
            } else {
                if (!is.na(non_modified_peptides_df[["PROT"]]$pep_homol[l]) && (non_modified_peptides_df[["PROT"]]$pep_score[l] > non_modified_peptides_df[["PROT"]]$pep_homol[l])) {
                    non_modified_peptides_df[["PROT"]]$identity[l] <- "homology"
                } else {
                    non_modified_peptides_df[["PROT"]]$identity[l] <- "discard"
                }
            }
        }
        for (l in 1:length(non_modified_peptides_df[["MOD"]]$identity)) {
            if (non_modified_peptides_df[["MOD"]]$pep_score[l] >= non_modified_peptides_df[["MOD"]]$pep_ident[l]) {
                non_modified_peptides_df[["MOD"]]$identity[l] <- "identity"
            } else {
                if (!is.na(non_modified_peptides_df[["MOD"]]$pep_homol[l]) && (non_modified_peptides_df[["MOD"]]$pep_score[l] > non_modified_peptides_df[["MOD"]]$pep_homol[l])) {
                    non_modified_peptides_df[["MOD"]]$identity[l] <- "homology"
                } else {
                    non_modified_peptides_df[["MOD"]]$identity[l] <- "discard"
                }
            }
        }
        
        for (l in 1:length(modified_peptides_df[["PROT"]]$identity)) {
            if (modified_peptides_df[["PROT"]]$pep_score[l] >= modified_peptides_df[["PROT"]]$pep_ident[l]) {
                modified_peptides_df[["PROT"]]$identity[l] <- "identity"
            } else {
                if (!is.na(modified_peptides_df[["PROT"]]$pep_homol[l]) && (modified_peptides_df[["PROT"]]$pep_score[l] > modified_peptides_df[["PROT"]]$pep_homol[l])) {
                    modified_peptides_df[["PROT"]]$identity[l] <- "homology"
                } else {
                    modified_peptides_df[["PROT"]]$identity[l] <- "discard"
                }
            }
        }
        for (l in 1:length(modified_peptides_df[["MOD"]]$identity)) {
            if (modified_peptides_df[["MOD"]]$pep_score[l] >= modified_peptides_df[["MOD"]]$pep_ident[l]) {
                modified_peptides_df[["MOD"]]$identity[l] <- "identity"
            } else {
                if (!is.na(modified_peptides_df[["MOD"]]$pep_homol[l]) && (modified_peptides_df[["MOD"]]$pep_score[l] > modified_peptides_df[["MOD"]]$pep_homol[l])) {
                    modified_peptides_df[["MOD"]]$identity[l] <- "homology"
                } else {
                    modified_peptides_df[["MOD"]]$identity[l] <- "discard"
                }
            }
        }
        
        # Discard the peptides to be discarded
        non_modified_peptides_df[["PROT"]] <- non_modified_peptides_df[["PROT"]][non_modified_peptides_df[["PROT"]]$identity != "discard", ]
        non_modified_peptides_df[["MOD"]] <- non_modified_peptides_df[["MOD"]][non_modified_peptides_df[["MOD"]]$identity != "discard", ]
        modified_peptides_df[["PROT"]] <- modified_peptides_df[["PROT"]][modified_peptides_df[["PROT"]]$identity != "discard", ]
        modified_peptides_df[["MOD"]] <- modified_peptides_df[["MOD"]][modified_peptides_df[["MOD"]]$identity != "discard", ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.20, title = "Discarding non-modified duplicates...", label = "20 %")
        
        
        
        
        
        ########## REMOVE DUPLICATES (NON-MOD)
        non_modified_peptides_df[["PROT"]]$prot_acc <- as.character(non_modified_peptides_df[["PROT"]]$prot_acc)
        non_modified_peptides_df[["MOD"]]$prot_acc <- as.character(non_modified_peptides_df[["MOD"]]$prot_acc)
        non_modified_peptides_df[["PROT"]]$pep_seq <- as.character(non_modified_peptides_df[["PROT"]]$pep_seq)
        non_modified_peptides_df[["MOD"]]$pep_seq <- as.character(non_modified_peptides_df[["MOD"]]$pep_seq)
        
        # Order the dataframe according to prot_acc, pep_seq, pep_score
        non_modified_peptides_df[["PROT"]] <- non_modified_peptides_df[["PROT"]][order(non_modified_peptides_df[["PROT"]]$prot_acc, non_modified_peptides_df[["PROT"]]$pep_seq, -non_modified_peptides_df[["PROT"]]$pep_score), ]
        non_modified_peptides_df[["MOD"]] <- non_modified_peptides_df[["MOD"]][order(non_modified_peptides_df[["MOD"]]$prot_acc, non_modified_peptides_df[["MOD"]]$pep_seq, -non_modified_peptides_df[["MOD"]]$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc and pep_seq and only the highest score for prot_acc/pep_seq duplicates)
        row_ID_to_keep <- list(PROT = vector(), MOD = vector())
        prot_acc_unique <- list()
        non_modified_peptides_df_prot_acc <- list()
        pep_seq_unique <- list()
        non_modified_peptides_df_prot_acc_pep_seq <- list()
        final_df <- list()
        
        # Determine the unique values of the prot_acc column
        prot_acc_unique[["PROT"]] <- unique(non_modified_peptides_df[["PROT"]]$prot_acc)
        prot_acc_unique[["MOD"]] <- unique(non_modified_peptides_df[["MOD"]]$prot_acc)
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["PROT"]])) {
            # Rows with the prot_acc_value
            non_modified_peptides_df_prot_acc[["PROT"]] <- non_modified_peptides_df[["PROT"]][non_modified_peptides_df[["PROT"]]$prot_acc == prot_acc_unique[["PROT"]][pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["PROT"]] <- unique(non_modified_peptides_df_prot_acc[["PROT"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["PROT"]])) {
                # Rows with the pep_seq_value
                non_modified_peptides_df_prot_acc_pep_seq[["PROT"]] <- non_modified_peptides_df_prot_acc[["PROT"]][non_modified_peptides_df_prot_acc[["PROT"]]$pep_seq == pep_seq_unique[["PROT"]][ps],]
                # Keep only the first (with the highest pep_score) (store the row ID)
                row_ID_to_keep[["PROT"]] <- append(row_ID_to_keep[["PROT"]], rownames(non_modified_peptides_df_prot_acc_pep_seq[["PROT"]][1,]))
            }
        }
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["MOD"]])) {
            # Rows with the prot_acc_value
            non_modified_peptides_df_prot_acc[["MOD"]] <- non_modified_peptides_df[["MOD"]][non_modified_peptides_df[["MOD"]]$prot_acc == prot_acc_unique[["MOD"]][pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["MOD"]] <- unique(non_modified_peptides_df_prot_acc[["MOD"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["MOD"]])) {
                # Rows with the pep_seq_value
                non_modified_peptides_df_prot_acc_pep_seq[["MOD"]] <- non_modified_peptides_df_prot_acc[["MOD"]][non_modified_peptides_df_prot_acc[["MOD"]]$pep_seq == pep_seq_unique[["MOD"]][ps],]
                # Keep only the first (with the highest pep_score) (store the row ID)
                row_ID_to_keep[["MOD"]] <- append(row_ID_to_keep[["MOD"]], rownames(non_modified_peptides_df_prot_acc_pep_seq[["MOD"]][1,]))
            }
        }
        # Retrieve the rows to keep
        final_df[["PROT"]] <- non_modified_peptides_df[["PROT"]][rownames(non_modified_peptides_df[["PROT"]]) %in% row_ID_to_keep[["PROT"]], ]
        non_modified_peptides_df[["PROT"]] <- final_df[["PROT"]]
        final_df[["MOD"]] <- non_modified_peptides_df[["MOD"]][rownames(non_modified_peptides_df[["MOD"]]) %in% row_ID_to_keep[["MOD"]], ]
        non_modified_peptides_df[["MOD"]] <- final_df[["MOD"]]
        
        # Display the 'identity' first
        non_modified_peptides_df[["PROT"]] <- non_modified_peptides_df[["PROT"]][order(non_modified_peptides_df[["PROT"]]$identity, non_modified_peptides_df[["PROT"]]$pep_score, non_modified_peptides_df[["PROT"]]$prot_acc, decreasing = TRUE), ]
        non_modified_peptides_df[["MOD"]] <- non_modified_peptides_df[["MOD"]][order(non_modified_peptides_df[["MOD"]]$identity, non_modified_peptides_df[["MOD"]]$pep_score, non_modified_peptides_df[["MOD"]]$prot_acc, decreasing = TRUE), ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.30, title = "Discarding modified duplicates...", label = "30 %")
        
        
        
        
        
        ########## REMOVE DUPLICATES (MOD)
        modified_peptides_df[["PROT"]]$prot_acc <- as.character(modified_peptides_df[["PROT"]]$prot_acc)
        modified_peptides_df[["MOD"]]$prot_acc <- as.character(modified_peptides_df[["MOD"]]$prot_acc)
        modified_peptides_df[["PROT"]]$pep_seq <- as.character(modified_peptides_df[["PROT"]]$pep_seq)
        modified_peptides_df[["MOD"]]$pep_seq <- as.character(modified_peptides_df[["MOD"]]$pep_seq)
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df[["PROT"]] <- modified_peptides_df[["PROT"]][order(modified_peptides_df[["PROT"]]$prot_acc, modified_peptides_df[["PROT"]]$pep_seq, modified_peptides_df[["PROT"]]$pep_var_mod, modified_peptides_df[["PROT"]]$pep_var_mod_pos, -modified_peptides_df[["PROT"]]$pep_score), ]
        modified_peptides_df[["MOD"]] <- modified_peptides_df[["MOD"]][order(modified_peptides_df[["MOD"]]$prot_acc, modified_peptides_df[["MOD"]]$pep_seq, modified_peptides_df[["MOD"]]$pep_var_mod, modified_peptides_df[["MOD"]]$pep_var_mod_pos, -modified_peptides_df[["MOD"]]$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc, pep_seq, pep_var_mod and pep_var_mod_pos and only the highest score for prot_acc/pep_seq/pep_var_mod/pep_var_mod_pos duplicates)
        row_ID_to_keep <- list(PROT = vector(), MOD = vector())
        prot_acc_unique <- list()
        modified_peptides_df_prot_acc <- list()
        pep_seq_unique <- list()
        modified_peptides_df_prot_acc_pep_seq <- list()
        pep_var_mod_unique <- list()
        modified_peptides_df_prot_acc_pep_seq_pep_var_mod <- list()
        pep_var_mod_pos_unique <- list()
        modified_peptides_df_prot_acc_pep_seq_pep_var_mod_pep_var_mod_pos <- list()
        final_df <- list()
        # Determine the unique values of the prot_acc column
        prot_acc_unique[["PROT"]] <- unique(modified_peptides_df[["PROT"]]$prot_acc)
        prot_acc_unique[["MOD"]] <- unique(modified_peptides_df[["MOD"]]$prot_acc)
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["PROT"]])) {
            # Rows with the prot_acc_value
            modified_peptides_df_prot_acc[["PROT"]] <- modified_peptides_df[["PROT"]][modified_peptides_df[["PROT"]]$prot_acc == prot_acc_unique[["PROT"]][pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["PROT"]] <- unique(modified_peptides_df_prot_acc[["PROT"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["PROT"]])) {
                # Rows with the pep_seq_value
                modified_peptides_df_prot_acc_pep_seq[["PROT"]] <- modified_peptides_df_prot_acc[["PROT"]][modified_peptides_df_prot_acc[["PROT"]]$pep_seq == pep_seq_unique[["PROT"]][ps],]
                # Determine the unique values of the pep_var_mod column
                pep_var_mod_unique[["PROT"]] <- unique(modified_peptides_df_prot_acc_pep_seq[["PROT"]]$pep_var_mod)
                # For each pep_var_mod...
                for (pvm in 1:length(pep_var_mod_unique[["PROT"]])) {
                    # Rows with the pep_var_mod_value
                    modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]] <- modified_peptides_df_prot_acc_pep_seq[["PROT"]][modified_peptides_df_prot_acc_pep_seq[["PROT"]]$pep_var_mod == pep_var_mod_unique[["PROT"]][pvm],]
                    # Determine the unique values of the pep_var_mod_pos column
                    modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]]$pep_var_mod_pos <- as.character(modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]]$pep_var_mod_pos)
                    pep_var_mod_pos_unique[["PROT"]] <- unique(modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]]$pep_var_mod_pos)
                    # For each pep_var_mod_pos...
                    for (pvmp in 1:length(pep_var_mod_pos_unique[["PROT"]])) {
                        # Rows with the pep_var_mod_pos_value
                        modified_peptides_df_prot_acc_pep_seq_pep_var_mod_pep_var_mod_pos[["PROT"]] <- modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]][modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["PROT"]]$pep_var_mod_pos == pep_var_mod_pos_unique[["PROT"]][pvmp],]
                        # Keep only the first (with the highest pep_score) (store the row ID)
                        row_ID_to_keep[["PROT"]] <- append(row_ID_to_keep[["PROT"]], rownames(modified_peptides_df_prot_acc_pep_seq_pep_var_mod_pep_var_mod_pos[["PROT"]][1,]))
                    }
                }
            }
        }
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["MOD"]])) {
            # Rows with the prot_acc_value
            modified_peptides_df_prot_acc[["MOD"]] <- modified_peptides_df[["MOD"]][modified_peptides_df[["MOD"]]$prot_acc == prot_acc_unique[["MOD"]][pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["MOD"]] <- unique(modified_peptides_df_prot_acc[["MOD"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["MOD"]])) {
                # Rows with the pep_seq_value
                modified_peptides_df_prot_acc_pep_seq[["MOD"]] <- modified_peptides_df_prot_acc[["MOD"]][modified_peptides_df_prot_acc[["MOD"]]$pep_seq == pep_seq_unique[["MOD"]][ps],]
                # Determine the unique values of the pep_var_mod column
                pep_var_mod_unique[["MOD"]] <- unique(modified_peptides_df_prot_acc_pep_seq[["MOD"]]$pep_var_mod)
                # For each pep_var_mod...
                for (pvm in 1:length(pep_var_mod_unique[["MOD"]])) {
                    # Rows with the pep_var_mod_value
                    modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]] <- modified_peptides_df_prot_acc_pep_seq[["MOD"]][modified_peptides_df_prot_acc_pep_seq[["MOD"]]$pep_var_mod == pep_var_mod_unique[["MOD"]][pvm],]
                    # Determine the unique values of the pep_var_mod_pos column
                    modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]]$pep_var_mod_pos <- as.character(modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]]$pep_var_mod_pos)
                    pep_var_mod_pos_unique[["MOD"]] <- unique(modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]]$pep_var_mod_pos)
                    # For each pep_var_mod_pos...
                    for (pvmp in 1:length(pep_var_mod_pos_unique[["MOD"]])) {
                        # Rows with the pep_var_mod_pos_value
                        modified_peptides_df_prot_acc_pep_seq_pep_var_mod_pep_var_mod_pos[["MOD"]] <- modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]][modified_peptides_df_prot_acc_pep_seq_pep_var_mod[["MOD"]]$pep_var_mod_pos == pep_var_mod_pos_unique[["MOD"]][pvmp],]
                        # Keep only the first (with the highest pep_score) (store the row ID)
                        row_ID_to_keep[["MOD"]] <- append(row_ID_to_keep[["MOD"]], rownames(modified_peptides_df_prot_acc_pep_seq_pep_var_mod_pep_var_mod_pos[["MOD"]][1,]))
                    }
                }
            }
        }
        # Retrieve the rows to keep
        final_df[["PROT"]] <- modified_peptides_df[["PROT"]][rownames(modified_peptides_df[["PROT"]]) %in% row_ID_to_keep[["PROT"]], ]
        modified_peptides_df[["PROT"]] <- final_df[["PROT"]]
        final_df[["MOD"]] <- modified_peptides_df[["MOD"]][rownames(modified_peptides_df[["MOD"]]) %in% row_ID_to_keep[["MOD"]], ]
        modified_peptides_df[["MOD"]] <- final_df[["MOD"]]
        
        # Display the 'identity' first
        modified_peptides_df[["PROT"]] <- modified_peptides_df[["PROT"]][order(modified_peptides_df[["PROT"]]$identity, modified_peptides_df[["PROT"]]$pep_score, modified_peptides_df[["PROT"]]$prot_acc, decreasing = TRUE), ]
        modified_peptides_df[["MOD"]] <- modified_peptides_df[["MOD"]][order(modified_peptides_df[["MOD"]]$identity, modified_peptides_df[["MOD"]]$pep_score, modified_peptides_df[["MOD"]]$prot_acc, decreasing = TRUE), ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.60, title = "Discarding modified sequence duplicates...", label = "60 %")
        
        
        
        
        
        ########## REMOVE SEQUENCE DUPLICATES (MOD)
        modified_peptides_df_sequences <- list()
        modified_peptides_df_sequences[["PROT"]] <- modified_peptides_df[["PROT"]]
        modified_peptides_df_sequences[["MOD"]] <- modified_peptides_df[["MOD"]]
        
        # Order the dataframe according to prot_acc, pep_seq, pep_var_mod, pep_var_mod_pos, pep_score
        modified_peptides_df_sequences[["PROT"]] <- modified_peptides_df_sequences[["PROT"]][order(modified_peptides_df_sequences[["PROT"]]$prot_acc, modified_peptides_df_sequences[["PROT"]]$pep_seq, -modified_peptides_df_sequences[["PROT"]]$pep_score), ]
        modified_peptides_df_sequences[["MOD"]] <- modified_peptides_df_sequences[["MOD"]][order(modified_peptides_df_sequences[["MOD"]]$prot_acc, modified_peptides_df_sequences[["MOD"]]$pep_seq, -modified_peptides_df_sequences[["MOD"]]$pep_score), ]
        
        # Extract the unique values (keep unique prot_acc and pep_seq, and only the highest score for prot_acc/pep_seq duplicates)
        row_ID_to_keep <- list(PROT = vector(), MOD = vector())
        prot_acc_unique <- list()
        modified_peptides_df_sequences_prot_acc <- list()
        pep_seq_unique <- list()
        modified_peptides_df_sequences_prot_acc_pep_seq <- list()
        final_df <- list()
        # Determine the unique values of the prot_acc column
        prot_acc_unique[["PROT"]] <- unique(modified_peptides_df_sequences[["PROT"]]$prot_acc)
        prot_acc_unique[["MOD"]] <- unique(modified_peptides_df_sequences[["MOD"]]$prot_acc)
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["PROT"]])) {
            # Rows with the prot_acc_value
            modified_peptides_df_sequences_prot_acc[["PROT"]] <- modified_peptides_df_sequences[["PROT"]][modified_peptides_df_sequences[["PROT"]]$prot_acc == prot_acc_unique[pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["PROT"]] <- unique(modified_peptides_df_sequences_prot_acc[["PROT"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["PROT"]])) {
                # Rows with the pep_seq_value
                modified_peptides_df_sequences_prot_acc_pep_seq[["PROT"]] <- modified_peptides_df_sequences_prot_acc[["PROT"]][modified_peptides_df_sequences_prot_acc[["PROT"]]$pep_seq == pep_seq_unique[["PROT"]][ps],]
                # Keep only the first (with the highest pep_score) (store the row ID)
                row_ID_to_keep[["PROT"]] <- append(row_ID_to_keep[["PROT"]], rownames(modified_peptides_df_sequences_prot_acc_pep_seq[["PROT"]][1,]))
            }
        }
        # For each prot_acc...
        for (pa in 1:length(prot_acc_unique[["MOD"]])) {
            # Rows with the prot_acc_value
            modified_peptides_df_sequences_prot_acc[["MOD"]] <- modified_peptides_df_sequences[["MOD"]][modified_peptides_df_sequences[["MOD"]]$prot_acc == prot_acc_unique[pa],]
            # Determine the unique values of the pep_seq column
            pep_seq_unique[["MOD"]] <- unique(modified_peptides_df_sequences_prot_acc[["MOD"]]$pep_seq)
            # For each pep_seq...
            for (ps in 1:length(pep_seq_unique[["MOD"]])) {
                # Rows with the pep_seq_value
                modified_peptides_df_sequences_prot_acc_pep_seq[["MOD"]] <- modified_peptides_df_sequences_prot_acc[["MOD"]][modified_peptides_df_sequences_prot_acc[["MOD"]]$pep_seq == pep_seq_unique[["MOD"]][ps],]
                # Keep only the first (with the highest pep_score) (store the row ID)
                row_ID_to_keep[["MOD"]] <- append(row_ID_to_keep[["MOD"]], rownames(modified_peptides_df_sequences_prot_acc_pep_seq[["MOD"]][1,]))
            }
        }
        # Retrieve the rows to keep
        final_df[["PROT"]] <- modified_peptides_df_sequences[["PROT"]][rownames(modified_peptides_df_sequences[["PROT"]]) %in% row_ID_to_keep[["PROT"]], ]
        modified_peptides_df_sequences[["PROT"]] <- final_df[["PROT"]]
        final_df[["MOD"]] <- modified_peptides_df_sequences[["MOD"]][rownames(modified_peptides_df_sequences[["MOD"]]) %in% row_ID_to_keep[["MOD"]], ]
        modified_peptides_df_sequences[["MOD"]] <- final_df[["MOD"]]
        
        # Display the 'identity' first
        modified_peptides_df_sequences[["PROT"]] <- modified_peptides_df_sequences[["PROT"]][order(modified_peptides_df_sequences[["PROT"]]$identity, modified_peptides_df_sequences[["PROT"]]$pep_score, modified_peptides_df_sequences[["PROT"]]$prot_acc, decreasing = TRUE), ]
        modified_peptides_df_sequences[["MOD"]] <- modified_peptides_df_sequences[["MOD"]][order(modified_peptides_df_sequences[["MOD"]]$identity, modified_peptides_df_sequences[["MOD"]]$pep_score, modified_peptides_df_sequences[["MOD"]]$prot_acc, decreasing = TRUE), ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.70, title = "Determining the true modified peptides...", label = "70 %")
        
        
        
        
        
        ########## MERGE MODIFIED PEPTIDE SEQUENCES FROM THE MODIFIED PROTEOME WITH THE MODIFIED PEPTIDE SEQUENCES FROM THE NON-MODIFIED PROTEOME
        # Check the peptide sequence (pep_seq) and then the modification position (pep_var_mod_pos)
        modified_peptides_df[["PROT"]]$pep_var_mod_pos <- as.character(modified_peptides_df[["PROT"]]$pep_var_mod_pos)
        modified_peptides_df[["MOD"]]$pep_var_mod_pos <- as.character(modified_peptides_df[["MOD"]]$pep_var_mod_pos)
        row_ID_to_discard <- vector()
        # Extract the unique sequences from the modified proteome
        pep_seq_unique <- unique(modified_peptides_df[["MOD"]]$pep_seq)
        # For each sequence...
        for (ps in 1:length(pep_seq_unique)) {
            # Rows with the pep_seq value
            modified_peptides_df_pep_seq <- modified_peptides_df[["MOD"]][modified_peptides_df[["MOD"]]$pep_seq == pep_seq_unique[ps], ]
            # Extract the unique modification
            pep_var_mod_pos_unique <- unique(modified_peptides_df_pep_seq$pep_var_mod_pos)
            # For each modification
            for (pvmp in 1:length(pep_var_mod_pos_unique)) {
                # Row with the pep_var_mod_pos value
                modified_peptides_df_pep_seq_pep_var_mod_pos <- modified_peptides_df_pep_seq[modified_peptides_df_pep_seq$pep_var_mod_pos == pep_var_mod_pos_unique[pvmp], ]
                # Search for this in the other file (first for the sequence, and then for the modification)
                same_pep_seq <- modified_peptides_df[["PROT"]][modified_peptides_df[["PROT"]]$pep_seq == pep_seq_unique[ps], ]
                same_pep_var_mod_pos <- same_pep_seq[same_pep_seq$pep_var_mod_pos == pep_var_mod_pos_unique[pvmp], ]
                # Add this to the rows to discard (if there is a match with the other database)
                if (length(same_pep_var_mod_pos) > 0) {
                    row_ID_to_discard <- append(row_ID_to_discard, rownames(modified_peptides_df_pep_seq_pep_var_mod_pos))
                }
            }
        }
        # Retrieve the rows to keep
        true_modified_peptides_df <- modified_peptides_df[["MOD"]][!(rownames(modified_peptides_df[["MOD"]]) %in% row_ID_to_discard), ]
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.85, title = "Merging files...", label = "85 %")
        
        non_modified_peptides_df[["PROT"]]$Type <- rep("NON-MODIFIED PROTEOME", times = nrow(non_modified_peptides_df[["PROT"]]))
        if (nrow(true_modified_peptides_df) > 0) {
            true_modified_peptides_df$Type <- rep("MODIFIED PROTEOME", times = nrow(true_modified_peptides_df))
            non_modified_true_modified_peptides_df <- rbind(true_modified_peptides_df, non_modified_peptides_df[["PROT"]])
        } else {
            modified_peptides_df[["MOD"]]$Type <- rep("MODIFIED PROTEOME", times = nrow(modified_peptides_df[["MOD"]]))
            non_modified_true_modified_peptides_df <- rbind(modified_peptides_df[["MOD"]], non_modified_peptides_df[["PROT"]])
        }
        
        # Add the class
        non_modified_true_modified_peptides_df$Class <- rep(class_name, nrow(non_modified_true_modified_peptides_df))
        
        
        
        
        
        
        
        
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 0.95, title = "Saving files...", label = "95 %")
        
        
        ########## SAVE FILES
        if (is.null(input_filename[["PROT"]])) {
            input_filename[["PROT"]] <- "PROT - Input Data"
        }
        if (is.null(input_filename[["MOD"]])) {
            input_filename[["MOD"]] <- "MOD - Input Data"
        }
        if (file_format == "csv") {
            write.csv(true_modified_peptides_df, file = paste("True modified peptides", ".", file_format, sep = ""), row.names = FALSE)
            write.csv(non_modified_true_modified_peptides_df, file = paste("Non-modified and True modified peptides", ".", file_format, sep = ""), row.names = FALSE)
            PROT_subfolder <- file.path(output_subfolder, "PROTEOME")
            dir.create(PROT_subfolder)
            setwd(PROT_subfolder)
            write.csv(input_data[["PROT"]], file = paste(input_filename[["PROT"]], ".", file_format, sep = ""), row.names = FALSE)
            write.csv(non_modified_peptides_df[["PROT"]], file = paste("Non-modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df[["PROT"]], file = paste("Modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df_sequences[["PROT"]], file = paste("Modified peptides (unique sequences).", file_format, sep = ""), row.names = FALSE)
            MOD_subfolder <- file.path(output_subfolder, "MODIFIED PROTEOME")
            dir.create(MOD_subfolder)
            setwd(MOD_subfolder)
            write.csv(input_data[["MOD"]], file = paste(input_filename[["MOD"]], ".", file_format, sep = ""), row.names = FALSE)
            write.csv(non_modified_peptides_df[["MOD"]], file = paste("Non-modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df[["MOD"]], file = paste("Modified peptides.", file_format, sep = ""), row.names = FALSE)
            write.csv(modified_peptides_df_sequences[["MOD"]], file = paste("Modified peptides (unique sequences).", file_format, sep = ""), row.names = FALSE)
        } else if (file_format == "xlsx" || file_format == "xls") {
            wb = loadWorkbook(filename = paste("True modified peptides", ".", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "True MOD peptides")
            writeWorksheet(wb, data = true_modified_peptides_df, sheet = "True MOD peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Non-modified and True modified peptides", ".", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Non-MOD + True MOD")
            writeWorksheet(wb, data = non_modified_true_modified_peptides_df, sheet = "Non-MOD + True MOD", header = TRUE)
            saveWorkbook(wb)
            
            PROT_subfolder <- file.path(output_subfolder, "PROTEOME")
            dir.create(PROT_subfolder)
            setwd(PROT_subfolder)
            
            wb = loadWorkbook(filename = paste(input_filename[["PROT"]], ".", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Input data")
            writeWorksheet(wb, data = input_data[["PROT"]], sheet = "Input data", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Non-modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Non-modified peptides")
            writeWorksheet(wb, data = non_modified_peptides_df[["PROT"]], sheet = "Non-modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified peptides")
            writeWorksheet(wb, data = modified_peptides_df[["PROT"]], sheet = "Modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides (unique sequences).", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified sequences")
            writeWorksheet(wb, data = modified_peptides_df[["PROT"]], sheet = "Modified sequences", header = TRUE)
            saveWorkbook(wb)
            
            MOD_subfolder <- file.path(output_subfolder, "MODIFIED PROTEOME")
            dir.create(MOD_subfolder)
            setwd(MOD_subfolder)
            wb = loadWorkbook(filename = paste(input_filename[["MOD"]], ".", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Input data")
            writeWorksheet(wb, data = input_data[["MOD"]], sheet = "Input data", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Non-modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Non-modified peptides")
            writeWorksheet(wb, data = non_modified_peptides_df[["MOD"]], sheet = "Non-modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides.", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified peptides")
            writeWorksheet(wb, data = modified_peptides_df[["MOD"]], sheet = "Modified peptides", header = TRUE)
            saveWorkbook(wb)
            
            wb = loadWorkbook(filename = paste("Modified peptides (unique sequences).", file_format, sep = ""), create = TRUE)
            createSheet(wb, name = "Modified sequences")
            writeWorksheet(wb, data = modified_peptides_df[["MOD"]], sheet = "Modified sequences", header = TRUE)
            saveWorkbook(wb)
        }
        
        # Progress bar
        setTkProgressBar(program_progress_bar, value = 1.00, title = "Done!", label = "100 %")
        close(program_progress_bar)
        
        
        # Go to the working directory
        setwd(output_folder)
        
        
        tkmessageBox(title = "Success!", message = "The process has successfully finished!", icon = "info")
    } else {
        if (input_folder == "" || is.null(unlist(input_data))) {
            ########## NO INPUT FILE
            tkmessageBox(title = "No input file selected!", message = "No input file has been selected!", icon = "warning")
        }
        if (all_the_columns_are_present == FALSE) {
            ########## MISSING DATA
            tkmessageBox(title = "Missing data!", message = paste0("Some data is missing from the selected input file!\n\n\n", missing_columns_value), icon = "warning")
        }
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
try({
    # Windows
    if (system_os == "Windows") {
        # Get system info
        screen_info <- system("wmic path Win32_VideoController get VideoModeDescription", intern = TRUE)[2]
        # Get the resolution
        screen_resolution <- unlist(strsplit(screen_info, "x"))
        # Retrieve the values
        screen_height <- as.numeric(screen_resolution[2])
        screen_width <- as.numeric(screen_resolution[1])
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
}, silent = TRUE)



### FONTS
# Default sizes (determined on a 1680x1050 screen) (in order to make them adjust to the size screen, the screen resolution should be retrieved)
title_font_size <- 24
other_font_size <- 11

# Adjust fonts size according to the pixel number
try({
    # Windows
    if (system_os == "Windows") {
        # Determine the font size according to the resolution
        total_number_of_pixels <- screen_width * screen_height
        # Determine the scaling factor (according to a complex formula)
        scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
        scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
        title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
        other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    } else if (system_os == "Linux") {
        # Linux
        # Determine the font size according to the resolution
        total_number_of_pixels <- screen_width * screen_height
        # Determine the scaling factor (according to a complex formula)
        scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
        scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
        title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
        other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    } else if (system_os == "Darwin") {
        # macOS
        print("Using default font sizes...")
    }
}, silent = TRUE)

# Define the fonts
# Windows
if (system_os == "Windows") {
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
    #Linux
    # Ubuntu
    if (length(grep("Ubuntu", os_version, ignore.case = TRUE)) > 0) {
        # Define the fonts
        ubuntu_title_bold = tkfont.create(family = "Ubuntu", size = (title_font_size + 2), weight = "bold")
        ubuntu_other_normal = tkfont.create(family = "Ubuntu", size = (other_font_size + 1), weight = "normal")
        ubuntu_other_bold = tkfont.create(family = "Ubuntu", size = (other_font_size + 1), weight = "bold")
        liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
        liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
        liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
        bitstream_charter_title_bold = tkfont.create(family = "Bitstream Charter", size = title_font_size, weight = "bold")
        bitstream_charter_other_normal = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "normal")
        bitstream_charter_other_bold = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "bold")
        # Use them in the GUI
        title_font = ubuntu_title_bold
        label_font = ubuntu_other_normal
        entry_font = ubuntu_other_normal
        button_font = ubuntu_other_bold
    } else if (length(grep("Fedora", os_version, ignore.case = TRUE)) > 0) {
        # Fedora
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
select_input_button <- tkbutton(window, text="IMPORT FILES...", command = file_import_function, font = button_font, bg = "white", width = 20)
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

