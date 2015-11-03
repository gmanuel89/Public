######################## LC LIBRARY FILE CONVERSION

####### Load the required packages
from Tkinter import *
import tkFileDialog
import tkMessageBox
import pandas
from pandas import DataFrame





### BEFORE STARTING
Tk().withdraw()
tkMessageBox.showinfo (title = "Before starting", message = "Remember to remove the header first!")





##### DEFINE THE FOLDER IN WHICH THERE IS THE CSV
# Path where the library file is located (GUI)
Tk().withdraw()
tkMessageBox.showinfo (title = "Path", message = "Select the LC library file to be imported")
Tk().withdraw()
inputfile = tkFileDialog.askopenfilename(filetypes=[('CSV files','.csv')])

## Define the path where to save output files
### Where to save the csv
Tk().withdraw()
tkMessageBox.showinfo (title = "Where to save the output file", message = "Select where to save the processed file", icon = "info")
Tk().withdraw()
outputfile = tkFileDialog.asksaveasfilename (defaultextension='csv', filetypes=[('CSV files','.csv')])

### Add the extension to the file automatically
# Convert the outputfile into a string
outputfile = str(outputfile)
if not outputfile.endswith(".csv"):
    outputfile = outputfile + ".csv"





#### PERCOLATOR DIALOG WINDOW
# Define the function associated with the command for the buttons
def set_percolator_value(value, window):
	# The variable percolator must exit from the function environment
	global percolator
	# Set the percolator variable value
	percolator = value
	# Close the dialogue window and exit the loop to go on with the code
	window.destroy()
	window.quit()

# Create a new toplevel window
top = Tk()
top.title("Percolator")
top.geometry("300x300")
# Declare the percolator variable
label = Label(top)
label.config(text = "Was the Percolator algorithm used?")
label.pack(expand=YES)
# Define what to do on click
yes_button = Button(top, text = "  YES  ", command = lambda *args: set_percolator_value("YES", top)).pack(expand=YES)
no_button = Button (top, text = "  NO  ", command = lambda *args: set_percolator_value("NO", top)).pack(expand=YES)
# Do not procede until the percolator value has been set
top.mainloop()





# Output messages
if percolator == "YES":
	Tk().withdraw()
	tkMessageBox.showinfo(message = "Percolator")
elif percolator == "NO":
	Tk().withdraw()
	tkMessageBox.showinfo(message = "No Percolator")





################################ AUTOMATED SCRIPT
# Run the entire script only if the file paths are selected
#if (inputfile != "" and outputfile != ""):
########## IMPORT THE CSV
LC_library = pandas.read_csv (inputfile)
# 1. Remove the empty columns
LC_library = LC_library.dropna(axis=1, how='all')
# 2. Sorting according to "pep_calc_mr" (smaller values first), "prot_acc" (smaller values first) and "pep_score" (higher values first)
LC_library = LC_library.sort_values (by=['pep_calc_mr', 'prot_acc', 'pep_score'], axis=0, ascending=[True, True, False])
# 3. Create the "score_ord" object
LC_library['score_ord'] = LC_library['pep_calc_mr']
LC_library = LC_library.reset_index(drop=True)
for i in range(1,len(LC_library['score_ord'])):
	# If pep_calc_mr = pep_calc_mr of the line above
	if LC_library['pep_calc_mr'][i] == LC_library['pep_calc_mr'][i-1]:
		# Make that value equal to zero, otherwise do nothing
		LC_library['score_ord'][i] = 0
	else:
		pass
# 4-5. Create the "acc_ord" object (identical to pep_calc_mr)
LC_library['acc_ord'] = LC_library['pep_calc_mr']
for i in range(len(LC_library['acc_ord'])):
	if i==0:
		pass
	else:
		# If prot_acc = prot_acc of the line above
		if LC_library['prot_acc'][i] == LC_library['prot_acc'][i-1]:
			# Make that value equal to zero, otherwise do nothing
			LC_library['acc_ord'][i] = 0
		else:
			pass
# 6-7. Create the "def" object
LC_library['def'] = LC_library['pep_calc_mr']
for i in range(len(LC_library['def'])):
	if LC_library['acc_ord'][i] + LC_library['score_ord'][i] == 0:
		LC_library['def'][i] = 0
# 8-9. Erase the rows that have a zero value in the "def" column and srt according to "def"
LC_library = LC_library [LC_library['def'] !=0]
LC_library = LC_library.sort_values (by=['def'], axis=0, ascending=False)
LC_library = LC_library.reset_index(drop=True)
# 10. Sort according to pep_score and remove the rows with scores less than 13
LC_library = LC_library.sort_values (by=['pep_score'], axis=0, ascending=False)
LC_library = LC_library [LC_library['pep_score'] >= 13]
LC_library = LC_library.reset_index(drop=True)
########################################################## WITHOUT PERCOLATOR
if percolator == 'NO':
	# 1. Create the "identity" object
	LC_library['identity'] = LC_library['pep_calc_mr']
	for i in range(len(LC_library['ident'])):
		if LC_library['pep_score'][i] < LC_library['pep_ident'][i]:
			LC_library['identity'][i] = 0
	# 2. Sort upon identity and highlight the cells with values
	LC_library = LC_library.sort_values (by=['identity'], axis=0)
	LC_library = LC_library.reset_index(drop=True)
	# 3. On the same column, from the first entry...
	for i in range(len(LC_library['identity'])):
		if LC_library['identity'][i] == 0:
			if LC_library['pep_score'][i] < LC_library['pep_hom'][i]:
				LC_library['identity'][i] = 0
			else:
				LC_library['identity'][i] = LC_library['pep_calc_mr'][i]
	# 4. Sort identity by colour
	# 5. Create the object "homology"
	LC_library['homology'] = LC_library['identity']
	for i in range(len(LC_library['homology'])):
		if LC_library['pep_hom'][i] == 0:
			LC_library['homology'][i] = 0
	# 6. Erase the rows with 0 in "identity" and with 0 in "homology"
	LC_library = LC_library[LC_library['identity'] !=0]
	LC_library = LC_library[LC_library['homology'] !=0]
	LC_library = LC_library.reset_index(drop=True)





############ EXPORT THE CSV
LC_library.to_csv (path_or_buf=outputfile)





### CONVERSION FINISHED
Tk().withdraw()
tkMessageBox.showinfo (title = "Conversion done!", message = "Conversion done! Considering donating 1000$ to the creator! :P")
