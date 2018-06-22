# Clear the console
cat("\014")
# Empty the workspace
rm(list = ls())


###################### R WORKSPACE DATA DUMPER - 2017.06.13


###################################################### CHECK INTERNET CONNECTION
# This function checks if there is internet connection, by pinging a website. It returns TRUE or FALSE.
# Two methods are available: 'ping' tries to ping the website, while 'getURL' connects to it directly. The default is 'getURL', since it is more reliable than ping.
check_internet_connection <<- function(method = "getURL", website_to_ping = "www.google.it") {
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

##################################################### INSTALL REQUIRED PACKAGES
# This function installs and loads the selected packages
install_and_load_required_packages <<- function(required_packages, repository = "http://cran.mirror.garr.it/mirrors/CRAN/", update_packages = FALSE, print_messages = FALSE) {
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
            if (print_messages == TRUE) {
                cat("\nPackages updated\n")
            }
        } else {
            if (print_messages == TRUE) {
                cat("\nPackages cannot be updated due to internet connection problems\n")
            }
        }
    }
    ##### Retrieve the installed packages
    installed_packages <- installed.packages()[,1]
    ##### Determine the missing packages
    missing_packages <- required_packages[!(required_packages %in% installed_packages)]
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
            if (print_messages == TRUE) {
                cat("\nAll the required packages have been installed\n")
            }
            all_needed_packages_are_installed <<- TRUE
        } else {
            ### If there is NO internet...
            if (print_messages == TRUE) {
                cat("\nSome packages cannot be installed due to internet connection problems\n")
            }
            all_needed_packages_are_installed <<- FALSE
        }
    } else {
        if (print_messages == TRUE) {
            cat("\nAll the required packages are installed\n")
        }
        all_needed_packages_are_installed <<- TRUE
    }
    ##### Load the packages (if there are all the packages)
    if ((length(missing_packages) > 0 && there_is_internet == TRUE) || length(missing_packages) == 0) {
        for (i in 1:length(required_packages)) {
            library(required_packages[i], character.only = TRUE, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
        }
        all_needed_packages_are_installed <<- TRUE
    } else {
        if (print_messages == TRUE) {
            cat("\nPackages cannot be installed/loaded... Expect issues...\n")
        }
        all_needed_packages_are_installed <<- FALSE
    }
}



install_and_load_required_packages("tcltk", update_packages = TRUE, print_messages = TRUE)

###############################################################################





#################### INFO
tkmessageBox(title = "Info", message = "This program explores the content of an R workspace (.RData) file, extracting the list of variables along with their variable type.", icon = "info")





#################### Path where the R workspace file is located (GUI)(FILE RData)
tkmessageBox(title = "Input R workspace", message = "Select the R workspace file to be explored", icon = "info")
filepath_R <- tclvalue(tkgetOpenFile(filetypes = "{{R workspace files} {.RData}}"))
if (!nchar(filepath_R)) {
    tkmessageBox(message = "No R workspace file was selected!")
}    else {
    tkmessageBox(message = paste("The R workspace file selected is", filepath_R))
}





#################### Path where to save the output CSV file
tkmessageBox(title = "Output CSV file", message = "Select where to save the output CSV file", icon = "info")
filepath_export <- tclvalue(tkgetSaveFile(filetypes = "{{Comma Separated Values (CSV) files} {.csv}}"))
# If not specified...
if (filepath_export == "") {
    filepath_export <- paste(getwd(), "/RData content.csv", sep = "")
} else {
    # Fix the extension
    if (length(grep(".csv", filepath_export, fixed = TRUE)) > 0) {
        filepath_export <- filepath_export
    } else {
        filepath_export <- paste(filepath_export, ".csv", sep = "")
    }
}
tkmessageBox(message = paste("The output file will be", filepath_export))





#################### Run only if a file is specified
if (filepath_R != "") {
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
        ## Write the output matrix to the file
        write.csv(output_matrix, file = filepath_export, col.names = TRUE, row.names = FALSE)
        # Done!
        tkmessageBox(title = "R workspace content dumped!", message = paste("The content of the selected R workspace has been dumped to the file\n\n", filepath_export, sep = ""), icon="info")
    } else {
        ### If there are NO variables...
        tkmessageBox(title = "Empty workspace", message = "The selected R workspace is empty!", icon = "warning")
    }
} else {
    ### No R workspaces!!
    tkmessageBox(title = "No R workspaces selected", message = "No R workspaces have been selected!", icon = "warning")
}





