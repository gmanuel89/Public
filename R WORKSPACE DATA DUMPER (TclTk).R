###################### R WORKSPACE DATA DUMPER - 2017.12.21.1




# READ AND DUMP THE CONTENT ON A RDATA WORKSPACE



# Update the packages
update.packages(repos = NULL, ask=FALSE)





##################################################### INSTALL REQUIRED PACKAGES
# This function installs and loads the selected packages
install_and_load_required_packages <- function(required_packages, repository=NULL) {
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



install_and_load_required_packages("tcltk")

###############################################################################





#################### INFO
tkmessageBox(title = "Info", message = "This program explores the content of an R workspace (.RData) file, extracting the list of variables along with their variable type.", icon = "info")





#################### Path where the R workspace file is located (GUI)(FILE RData)
tkmessageBox(title = "Input R workspace", message = "Select the R workspace file to be explored", icon = "info")
filepath_R <- tclvalue(tkgetOpenFile(filetypes = "{{R workspace files} {.RData}}"))
if (!nchar(filepath_R)) {
	tkmessageBox(message = "No R workspace file was selected!")
}	else {
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





