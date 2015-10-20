######################## CSV TRANSPOSAL

## CSV LIKE THIS:
# XXXXX
# XXXX
# XXXXXX
# XXXX
#
# RESULT
# XXXX
# XXXX
# XXXX
# XXXX
#  X X
#  X

#### Import the libraries
import sys
from PyQt4.QtGui import *
import csv

# Create a PyQt4 Application Object
a = QApplication(sys.argv)

# Message for selection (open)
w = QWidget()
selectionMessage = QMessageBox.information (w, "CSV selection", "Select the CSV file to be transposed")

# CSV file selection
w = QWidget()
inputfile = QFileDialog.getOpenFileName (w, 'Open CSV', '', 'CSV files (*.csv)')
print(inputfile)

# Message for selection (save)
w = QWidget()
selectionMessage = QMessageBox.information (w, "CSV selection", "Select where to save the transposed CSV file")

# Where to save the CSV file
w = QWidget()
outputfile = QFileDialog.getSaveFileName (w, 'Save CSV', '', 'CSV files (*.csv)')


### Add the extension to the file automatically
if ".csv" not in outputfile:
    outputfile = str(outputfile) + ".csv"

print(outputfile)


############################################# Transposal script
# Open the input file (in a temporary variable f)
with open(inputfile, 'rb') as f:
    # Read the actual csv values
    csvFileRead = csv.reader(f)
    # Create an empty array for storing the values (of the rows) to be stored in columns
    csvColumns = []
    # Scroll the rows of the CSV
    for row in csvFileRead:
        # Add the values in the row to the vector that will become the column
        csvColumns.append(row)

# Open the output file (in a temporary variable f)
with open(outputfile, 'wb') as f:
    # Prepare to write the csv values with the function
    csvFileWrite = csv.writer(f)
    # For each column
    for i in range(len(max(csvColumns, key=len))):
        # Write the values in columns in the output file
        csvFileWrite.writerow([(c[i] if i<len(c) else '') for c in csvColumns])


# Message for completion
w = QWidget()
selectionMessage = QMessageBox.information (w, "CSV transposed!!!", "The CSV file has been succesfully transposed")
