############## MGF FILE CONVERTER
import sys
from PyQt4.QtGui import *

# Create a PyQt4 Application Object
a = QApplication(sys.argv)

# Before we start
w = QWidget()
selectionMessage = QMessageBox.information (w, "MGF converter", "This application will replace each Cmpd value with the correspondent MSMS value in the selected MGF file")

# Message for selection (open)
w = QWidget()
selectionMessage = QMessageBox.information (w, "MGF selection", "Select the MGF file to be converted")

# CSV file selection
w = QWidget()
inputfile = QFileDialog.getOpenFileName (w, 'Open MGF', '.mgf', 'MGF files (*.mgf)')
print(inputfile)

# Message for selection (save)
w = QWidget()
selectionMessage = QMessageBox.information (w, "MGF selection", "Select where to save the converted MGF file")

# Where to save the CSV file
w = QWidget()
outputfile = QFileDialog.getSaveFileName (w, 'Save MGF', '', 'MGF files (*.mgf)')


### Add the extension to the file automatically
if ".mgf" not in outputfile:
    outputfile = str(outputfile) + ".mgf"





# Conversion script (regular expressions)
with open (inputfile) as mgf:
    # Create the final list of lines
    outputLines = []
    # Scroll the lines
    for line in mgf:
        # If the line starts with ##MSMS
        if line.startswith ("###MSMS"):
            # Write the line to the output file
            outputLines.append (line)
            #### Store the value of the MSMS
            # Break the line at the "###MSMS:" and take the second half
            lineSplitted1 = line.split ("###MSMS:")
            lineSplitted2 = lineSplitted1[1]
            # Break the line at the "/" and take the first half (convert it into a number)
            lineSplitted3 = lineSplitted2.split ("/")
            msmsValue = int (lineSplitted3[0])
        # If there is the TITLE value (the MSMS value has to be passed to TITLE)
        elif line.startswith("TITLE"):
            # Split the line at the commas
            titleLineSplitted = line.split (",")
            # Three pieces are generated: the last two will be joined back together with a comma, the first will have the TITLE Cmpd number replaced by the MSMS value
            titleLineSplitted[0] = "TITLE= Cmpd " + str(msmsValue)
		    # Join back the pieces
            finalTitleLine = titleLineSplitted[0] + "," + titleLineSplitted[1] + "," + titleLineSplitted[2]
		    # Write the final line to the file
            outputLines.append (finalTitleLine)
        # If the line does not start with ##MSMS or TITLE
        else:
            # Write the line to the output file
            outputLines.append (line)


with open (outputfile, 'w') as f:
    for line in outputLines:
        f.writelines (line)


# Conversion done!
w = QWidget()
selectionMessage = QMessageBox.information (w, "MGF conversion done", "The conversion has been successfully performed!")
