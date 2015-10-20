################## CSV TRANSPOSAL

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
from tkinter import *
import tkinter.filedialog
import tkinter.messagebox
import csv

# Message for selection (open)
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="CSV selection", message="Select the CSV file to be transposed")

# CSV file selection
tkinter.Tk().withdraw()
inputfile = tkinter.filedialog.askopenfilename(filetypes=[('CSV files','.csv'),('TXT files','.txt')])
print(inputfile)

# Message for selection (save)
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="CSV selection", message="Select where to save the transposed CSV file")

# Where to save the CSV file
tkinter.Tk().withdraw()
outputfile = tkinter.filedialog.asksaveasfilename (defaultextension='.csv', filetypes=[('CSV files','.csv'),('TXT files','.txt')])
print(outputfile)


############################################# Transposal script
# Open the input file (in a temporary variable f)
with open(inputfile) as f:
    # Read the actual csv values
    csvFileRead = csv.reader(f)
    # Create an empty array for storing the values (of the rows) to be stored in columns
    csvColumns = []
    # Scroll the rows of the CSV
    for row in csvFileRead:
        # Add the values in the row to the vector that will become the column
        csvColumns.append(row)

# Open the output file (in a temporary variable f)
with open(outputfile, 'w') as f:
    # Prepare to write the csv values with the function
    csvFileWrite = csv.writer(f)
    # max(csvColumns, key=len) iteratively scrolls the rows and returns the longest row, because it measures the length
    # len (max (csvColumns, key=len)) returns the longest row
    for i in range(len(max(csvColumns, key=len))):
        # Write the values in columns in the output file
        csvFileWrite.writerow([(c[i] if i<len(c) else '') for c in csvColumns])


# Message for completion
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="CSV transposed!!!", message="The CSV file has been succesfully transposed")
