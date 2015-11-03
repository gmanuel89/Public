######################## LC LIBRARY FILE CONVERSION

###INSTALL THE REQUIRED PACKAGES
requiredPackages <- c ("tcltk", "beepr")
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


####### Load the required packages
library ("tcltk")
library ("beepr")







### BEFORE STARTING
tkmessageBox (title = "Before starting", message = "Remember to remove the header first!", icon = "warning")


#### PERCOLATOR DIALOG WINDOW
tt <- tktoplevel()  # Create a new toplevel window
tktitle(tt) <- "Percolator?"
# Declare the percolator variable
percolator <- tclVar(0)
# Define what to do on click
YES.but <- tkbutton (tt, text = "  YES  ",
	command = function() {tclvalue(percolator) <- 1 ; tkdestroy(tt)})
NO.but <- tkbutton (tt, text = "  NO  ",
	command = function() {tclvalue(percolator) <- 0 ; tkdestroy(tt)})
Quit.but <- tkbutton (tt, text = "Cancel", command = function () {q(save="no")})
# Place the buttons
tkgrid (YES.but, NO.but, Quit.but)
tkfocus(tt)
# Do not proceed until a value is set for the variable
tkwait.variable(percolator)
# Create the variable percolatorVal that stores the value of the tcl variable
percolatorVal <- as.numeric (tclvalue(percolator))   # Get and coerce content of a Tcl variable
# Output messages
if (percolatorVal == 1) tkmessageBox(message = "Percolator")
if (percolatorVal == 0) tkmessageBox(message = "No Percolator")






##### DEFINE THE FOLDER IN WHICH THERE IS THE CSV
# Path where the library file is located (GUI)
tkmessageBox (title = "Path", message = "Select the LC library file to be imported", icon = "info")
filepath <- tclvalue(tkgetOpenFile())
if (!nchar(filepath)) {
	tkmessageBox(message = "No file was selected!")
}	else {
	tkmessageBox(message = paste("The file selected is", filepath))
}

## Define the path where to save output files
### Where to save the csv
tkmessageBox (title = "Where to save the output file", message = "Select where to save the processed file", icon = "info")
filepathExport <- tclvalue(tkgetSaveFile())
if (!nchar(filepathExport)) {
    tkmessageBox(message = "No file selected")
}










################################ AUTOMATED SCRIPT
# Run the entire script only if the file paths are selected
if (filepath != "" && filepathExport != "") {


########## IMPORT THE CSV
LClibrary <- read.csv (filepath)


# 1. Remove the empty columns
LClibrary <- LClibrary[, !apply(is.na(LClibrary), 2, all)]

# 2. Sorting according to "pep_calc_mr" (smaller values first), "prot_acc" (smaller values first) and "pep_score" (higher values first)
LClibrary <- LClibrary [order (LClibrary$pep_calc_mr, LClibrary$prot_acc, -LClibrary$pep_score),]

# 3. Create the "score_ord" object
LClibrary$score_ord <- LClibrary$pep_calc_mr

for (i in 2:length(LClibrary$score_ord)) {
	# If pep_calc_mr = pep_calc_mr of the line above
	if (LClibrary$pep_calc_mr[i] == LClibrary$pep_calc_mr[i-1]) {
	# Make that value equal to zero, otherwise do nothing
	LClibrary$score_ord[i] <- 0
	}
}

# 4-5. Create the "acc_ord" object (identical to pep_calc_mr)
LClibrary$acc_ord <- LClibrary$pep_calc_mr

for (i in 2:length(LClibrary$acc_ord)) {
	# If prot_acc = prot_acc of the line above
	if (LClibrary$prot_acc[i] == LClibrary$prot_acc[i-1]) {
	# Make that value equal to zero, otherwise do nothing
	LClibrary$acc_ord[i] <- 0
	}
}

# 6-7. Create the "def" object
LClibrary$def <- LClibrary$pep_calc_mr

for (i in 1:length(LClibrary$def)) {
	if ((LClibrary$acc_ord[i] + LClibrary$score_ord[i]) == 0) {
	LClibrary$def[i] <- 0
	}
}

# 8-9. Erase the rows that have a zero value in the "def" column and srt according to "def"
LClibrary <- LClibrary [order (LClibrary$def),]
LClibrary <- LClibrary [LClibrary$def !=0,]

# 10. Sort according to pep_score and remove the rows with scores less than 13
LClibrary <- LClibrary [order (LClibrary$pep_score),]
LClibrary <- LClibrary [LClibrary$pep_score >= 13,]








########################################################## WITHOUT PERCOLATOR

if (percolatorVal == 0) {
	# 1. Create the "identity" object
	LClibrary$identity <- LClibrary$pep_calc_mr

	for (i in 1:length(LClibrary$ident)) {
		if (LClibrary$pep_score[i] < LClibrary$pep_ident[i]) {
		LClibrary$identity[i] <- 0
		}
	}

	# 2. Sort upon identity and highlight the cells with values
	LClibrary <- LClibrary [order (LClibrary$identity),]

	# 3. On the same column, from the first entry...
	for (i in 1:length(LClibrary$identity)) {
		if (LClibrary$identity[i] == 0) {
			if (LClibrary$pep_score[i] < LClibrary$pep_hom[i]) {
			LClibrary$identity[i] <- 0
			}	else {LClibrary$identity[i] <- LClibrary$pep_calc_mr[i]}
		}
	}

	# 4. Sort identity by colour

	# 5. Create the object "homology"
	LClibrary$homology <- LClibrary$identity

	for (i in 1:length(LClibrary$homology)) {
		if (LClibrary$pep_hom[i] == 0) {
		LClibrary$homology[i] <- 0
		}
	}

	# 6. Erase the rows with 0 in "identity" and with 0 in "homology"
	LClibrary <- LClibrary[LClibrary$identity !=0,]
	LClibrary <- LClibrary[LClibrary$homology !=0,]
}










############ EXPORT THE CSV
# Add the extension if it is not present in the filename
if (length(grep(".csv", filepathExport, fixed=TRUE)) == 1) {
	filename <- filepathExport
}	else {filename <- paste (filepathExport, ".csv", sep="")}
write.csv (LClibrary, file=filename, row.names=FALSE)


beep (3)


### CONVERSION FINISHED
tkmessageBox (title = "Conversion done!", message = "Conversion done! Considering donating 1000$ to the creator! :P", icon = "info")



}
