############## MGF FILE CONVERTER
from tkinter import *
import tkinter.filedialog
import tkinter.messagebox


# Before we start
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="MGF converter", message="This application will replace each Cmpd value with the correspondent MSMS value in the selected MGF file")

# Message for selection (open)
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="MGF selection", message="Select the MGF file to be converted")

# MGF file selection
tkinter.Tk().withdraw()
inputfile = tkinter.filedialog.askopenfilename(filetypes=[('MGF files','.mgf')])

# Message for selection (save)
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="MGF selection", message="Select where to save the converted MGF file")

# Where to save the CSV file
tkinter.Tk().withdraw()
outputfile = tkinter.filedialog.asksaveasfilename (defaultextension='.mgf', filetypes=[('MGF files','.mgf')])


### Add the extension to the file automatically
# Convert the outputfile into a string
outputfile = str(outputfile)
if not outputfile.endswith(".mgf"):
    outputfile = outputfile + ".mgf"




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
tkinter.Tk().withdraw()
tkinter.messagebox.showinfo(title="MGF conversion done", message="The conversion has been successfully performed!")
