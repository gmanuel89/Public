#### Import the libraries
import sys





########### Input values
coordinates_of_spray_x_axis = [-30, 30]
coordinates_of_spray_y_axis = [-80, 80]
height_of_the_needle = 60
distance_between_lines = 1.81
amount_of_solvent = 500
speed_of_movement = 100
matrix_density = 1
number_of_spray_cycles = 10
number_of_valve_rinsing_cycles = 5
number_of_initial_wash_cycles = 5
horizontal_spraying = True





############## Constant values
coordinates_of_spray_z_axis = -80 + abs(height_of_the_needle) # Calculated
amount_of_solvent_per_spray = 11.9760479042 #?????? calculate
max_speed_of_movement = 200
max_height_of_the_needle = 80
max_matrix_density = 5
coordinates_of_washing_x_axis = 0
coordinates_of_washing_y_axis = -110
coordinates_of_washing_z_axis = -50
coordinates_of_z_magnet = 16.5
coordinates_of_z_axis_during_movement = -35
initial_x_coordinates = 0
initial_y_coordinates = 0
initial_z_coordinates = 0
max_x_coordinates = [-60, 60]
max_y_coordinates = [-80, 80]
max_z_coordinates = [-75, 0]





################################################# Avoid that the machine breaks
# Speed
if speed_of_movement > max_speed_of_movement:
    speed_of_movement = max_speed_of_movement

if speed_of_movement <= 0:
    speed_of_movement = 1

# Needle height (Z coordinate)
if coordinates_of_spray_z_axis < max_z_coordinates[0] or coordinates_of_spray_z_axis > max_z_coordinates[1]:
    coordinates_of_spray_z_axis = -35

# Matrix density
if matrix_density > max_matrix_density:
    matrix_density = max_matrix_density

if matrix_density <= 0:
    matrix_density = 1

# X and Y coordinates
if coordinates_of_spray_x_axis[0] < max_x_coordinates[0]:
    coordinates_of_spray_x_axis[0] = max_x_coordinates[0]

if coordinates_of_spray_x_axis[1] > max_x_coordinates[1]:
    coordinates_of_spray_x_axis[1] = max_x_coordinates[1]

if coordinates_of_spray_y_axis[0] < max_y_coordinates[0]:
    coordinates_of_spray_y_axis[0] = max_y_coordinates[0]

if coordinates_of_spray_y_axis[1] > max_y_coordinates[1]:
    coordinates_of_spray_y_axis[1] = max_y_coordinates[1]





######################################################################## Qt4 GUI (requires pyQt4)
from PyQt4.QtGui import *
########### Where to save the output file
# Create a PyQt4 Application Object
a = QApplication(sys.argv)
# Message for selection (open)
w = QWidget()
selectionMessage = QMessageBox.information (w, "File selection", "Select where to save the gcode file")
# Where to save the CSV file
w = QWidget()
outputfile = QFileDialog.getSaveFileName (w, 'Save gcode file', '', 'gcode files (*.gcode)')
### Add the extension to the file automatically
if ".gcode" not in outputfile:
    outputfile = str(outputfile) + ".gcode"





######################################################################## GTK GUI
from Tkinter import *
import tkFileDialog
import tkMessageBox
########### Where to save the output file
Tk().withdraw()
tkMessageBox.showinfo(title="File selection", message="Select where to save the gcode file")
# Where to save the CSV file
Tk().withdraw()
outputfile = tkFileDialog.asksaveasfilename (defaultextension='.gcode', filetypes=[('gcode files','.gcode')])
### Add the extension to the file automatically
if ".gcode" not in outputfile:
    outputfile = str(outputfile) + ".gcode"





################################################################################
################################################################### CODE BLOCKS
go_to_wash_position = [";;;;; go to wash position\n"]
go_to_wash_position_subblock = ["G1 Z%s F%s\n" %(float(coordinates_of_z_axis_during_movement), max_speed_of_movement), "G1 X%s Y%s Z%s F%s\n" %(float(coordinates_of_washing_x_axis), float(coordinates_of_washing_y_axis), float(coordinates_of_z_axis_during_movement), max_speed_of_movement), "G1 Z%s\n" %float(coordinates_of_washing_z_axis)]
# Generate the definitive block
for s in go_to_wash_position_subblock:
    go_to_wash_position.append (s)

# Leave a white line between blocks
go_to_wash_position.append ("\n")

########################################################## Initialisation block
initialisation_block = [";;;;;;;;;; initialisation\n"]
initialisation_subblock = ["G28XYZ\n", "G28P\n", "G90\n"]
# Generate the definitive block
for s in initialisation_subblock:
    initialisation_block.append (s)

# Leave a white line between blocks
initialisation_block.append ("\n")

############################################################## First wash block
first_wash_block = [";;;;;;;;;; first wash\n"]
# Go to wash position first
for s in go_to_wash_position:
    first_wash_block.append (s)

# Define the washing sub-block
first_wash_subblock = ["G1 V3.5 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P4.2 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 V3 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 V0 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P0 F%s\n" %(max_speed_of_movement/2), "G4 S1\n", "\n"]
# Generate the definitive block based on the number of repetitions
for i in range(number_of_initial_wash_cycles):
    for s in first_wash_subblock:
        first_wash_block.append (s)

# After washing
after_washing_block = ["M106\n", "G1 V3 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P2 F%s\n" %max_speed_of_movement, "G1 V1 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 F0.1\n", "G1 P0\n"]
# Generate the definitive block
for s in after_washing_block:
    first_wash_block.append (s)

# Leave a white line between blocks
first_wash_block.append ("\n")

######################################################## Inter-spray wash block
interspray_wash_block = [";;;;;;;;;; interspray wash\n"]
# Go to wash position first
for s in go_to_wash_position:
    interspray_wash_block.append (s)

# Washing sub-block
interspray_wash_subblock = [";;;;; wash\n", "M82\n", "G1 V0 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P0 F%s\n" %max_speed_of_movement, "M106 S0\n", "G4 S0.0\n"]
# Generate the definitive block
for s in interspray_wash_subblock:
    interspray_wash_block.append (s)

# Leave a white line between blocks
interspray_wash_block.append ("\n")

########################################################### Spray rinsing block
spray_rinsing_block = [";;;;; spray rinsing\n"]
spray_rinsing_subblock = ["G1 V1.5 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P4 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 V2 F%s\n" %max_speed_of_movement, "G4 S2\n", "G1 V1 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 F0.1\n", "G1 P0\n", "M106 S0\n"]
# Generate the definitive block
for s in spray_rinsing_subblock:
    spray_rinsing_block.append (s)

# Leave a white line between blocks
spray_rinsing_block.append ("\n")

####################################################### Capillary rinsing block
capillary_rinsing_block = [";;;;; capillary rinsing\n"]
capillary_rinsing_subblock = ["G1 V1.5 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P3 F%s\n" %max_speed_of_movement, "G1 V2 F%s\n" %max_speed_of_movement, "G4 S2\n", "G1 V1 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 F0.1\n", "G1 P0\n", "G1 P1\n"]
# Generate the definitive block
for s in capillary_rinsing_subblock:
    capillary_rinsing_block.append (s)

# Leave a white line between blocks
capillary_rinsing_block.append ("\n")

################################################################## Drying block
drying_block = [";;;;; drying\n"]
drying_subblock = ["M106\n", "G1 P0\n", "G4 S8\n", "M106 S0\n"]
# Generate the definitive block
for s in drying_subblock:
    drying_block.append (s)

# Leave a white line between blocks
drying_block.append ("\n")

########################################################### Parking spray block
parking_spray_block = [";;;;; parking spray\n"]
parking_spray_subblock = ["G1 Z%s F%s\n" %(float(coordinates_of_z_axis_during_movement), max_speed_of_movement), "G1 X%s Y%s Z%s F%s\n" %(float(initial_x_coordinates), float(initial_y_coordinates), float(coordinates_of_z_axis_during_movement), max_speed_of_movement), "G1 Z%s F%s\n" %(float(coordinates_of_z_magnet), max_speed_of_movement)]
# Generate the definitive block
for s in parking_spray_subblock:
    parking_spray_block.append (s)

# Leave a white line between blocks
parking_spray_block.append ("\n")

############################################################## Motors off block
motors_off_block = [";;;;;;;;;; motors off\n", "M18\n"]
# Leave a white line between blocks
#motors_off_block.append ("\n")

######################################################## Syringe emptying block
syringe_emptying_block = [";;;;; syringe emptying\n"]
syringe_emptying_block.append ("G1 P0 F%s\n" %(max_speed_of_movement/2))
# Valve rinsing block (80 microlitres each step)
valve_rinsing_block = [";;;;; valve rinsing\n"]
valve_rinsing_subblock = ["G1 V1.5 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P4.8 F%s\n" %max_speed_of_movement, "G1 V2 F%s\n" %max_speed_of_movement, "G4 S2\n", "G1 V0 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P0 F%s\n" %(max_speed_of_movement/2)]
for s in valve_rinsing_subblock:
    valve_rinsing_block.append (s)

# Generate the definitive block based on the number of repetitions
for i in range(number_of_valve_rinsing_cycles):
    for s in valve_rinsing_block:
        syringe_emptying_block.append (s)

# Leave a white line between blocks
syringe_emptying_block.append ("\n")

############################################################ Cleaning procedure
########### Go to wash + emptying syringe + rinsing valve + rinsing spray + rinsing capillary + drying + parking spray
cleaning_procedure_block = [";;;;;;;;;; cleaning procedure\n"]
cleaning_procedure_subblock = ["M82\n", "G90\n"]
# Generate the definitive block
for s in cleaning_procedure_subblock:
    cleaning_procedure_block.append(s)

# Go to wash
for s in go_to_wash_position:
    cleaning_procedure_block.append (s)

# More code
cleaning_procedure_subblock = ["M106\n", "G1 V0 F%s\n" %max_speed_of_movement, "G4 S1\n"]
for s in cleaning_procedure_subblock:
    cleaning_procedure_block.append(s)

for s in syringe_emptying_block:
    cleaning_procedure_block.append (s)

for s in spray_rinsing_block:
    cleaning_procedure_block.append (s)

for s in capillary_rinsing_block:
    cleaning_procedure_block.append (s)

for s in drying_block:
    cleaning_procedure_block.append (s)

for s in parking_spray_block:
    cleaning_procedure_block.append (s)

################################################################### Spray block
# First part of the block
spray_block = [";;;;;;;;;; spraying phase\n"]
spray_subblock1 = ["M106\n", "M82\n", "G1 V3 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 P%s F%s\n" %(amount_of_solvent_per_spray+2, max_speed_of_movement), "G1 V1 F%s\n" %max_speed_of_movement, "G4 S1\n", "G1 F0.1\n", "G1 P%s\n" %amount_of_solvent_per_spray, "\n", ";;; where to start spraying\n", "G1 Z%s F%s\n" %(float(coordinates_of_z_axis_during_movement), max_speed_of_movement), "G1 X%s Y%s F%s\n" %(float(coordinates_of_spray_x_axis[0]), float(coordinates_of_spray_y_axis[0]), max_speed_of_movement), "G1 Z%s F%s\n" %(float(coordinates_of_spray_z_axis), max_speed_of_movement), "\n", "M83\n", "\n"]
for s in spray_subblock1:
    spray_block.append(s)

################### The X-axis stays fixed
if horizontal_spraying == True:
    ################# Generate the vector of x and y positions
    ###### Y positions
    y_positions = []
    # Calculate the number of Y positions (according to the y length and the distance between lines) (since the path is an S, each y coordinate is reported twice, so the number of positions must be doubled)
    number_of_y_positions = int((abs(coordinates_of_spray_y_axis[0]) + abs(coordinates_of_spray_y_axis[1])) / distance_between_lines) * 2
    # Calculate the other positions (since the path is an S, each y coordinate is reported twice)
    for i in range(number_of_y_positions-1):
        # The first position is where to start spraying (since the path is an S, each y coordinate is reported twice)
        if len(y_positions) == 0:
            y_positions.append (float(coordinates_of_spray_y_axis[0]))
            y_positions.append (float(coordinates_of_spray_y_axis[0]))
        else:
            y_positions.append (float(y_positions[i-1] + distance_between_lines))
    ###### X positions
    x_positions = []
    # Calculate the other positions (since the path is an S, each x coordinate is reported twice)
    for i in range(number_of_y_positions/2):
        # The first position and the last position are on the opposite side of the starting position (where to go spraying, along the x axis) (since the path is an S, each y coordinate is reported twice)
        if len(x_positions) == 0:
            x_positions.append(float(coordinates_of_spray_x_axis[1]))
        elif i%2 != 0:
            x_positions.append(float(coordinates_of_spray_x_axis[0]))
            x_positions.append(float(coordinates_of_spray_x_axis[0]))
        elif i%2 == 0:
            x_positions.append(float(coordinates_of_spray_x_axis[1]))
            x_positions.append(float(coordinates_of_spray_x_axis[1]))
        if len(x_positions) == (len(y_positions)-1):
            x_positions.append(float(coordinates_of_spray_x_axis[1]))
    ###### Values of P
    p_values = []
    for i in range(number_of_y_positions/2):
        p_values.append (float(-amount_of_solvent_per_spray*0.00125))
        p_values.append (float(-amount_of_solvent_per_spray*0.03))
    ###### Generate the final spray ssubblock
    spray_subblock = [";;; spray\n"]
    for i in range(number_of_y_positions):
        # Create the line and append it to the list of the block lines
        line = "G1 X%s Y%s P%s F%s\n" %(x_positions[i], y_positions[i], p_values[i], speed_of_movement)
        spray_subblock.append (line)
################### The Y-axis stays fixed
else:
    pass

spray_subblock.append ("G1 Y%s Z%s F%s\n" %(float(coordinates_of_spray_y_axis[0]), float(coordinates_of_z_axis_during_movement), max_speed_of_movement))
# Attach the spray sub-block to the global spray block
for s in spray_subblock:
    spray_block.append (s)

# Leave a white line between blocks
spray_block.append ("\n")

################################################ Write lines in the output file
with open(outputfile, "w") as f:
    # Initialisation block
    for line in initialisation_block:
        f.writelines(line)
    # First wash block
    for line in first_wash_block:
        f.writelines(line)
    # Spray
    for i in range(number_of_spray_cycles):
        for line in spray_block:
            f.writelines (line)
        for line in interspray_wash_block:
            f.writelines (line)
    # Syringe emptying block
    for line in cleaning_procedure_block:
        f.writelines(line)
    # Motors off
    for line in motors_off_block:
        f.writelines(line)





####################################################################### Qt4 GUI
# Gcode file generated!
w = QWidget()
selectionMessage = QMessageBox.information (w, "gcode file generated", "The gcode file has been successfully generated!")

####################################################################### GTK GUI
# Gcode file generated!
Tk().withdraw()
tkMessageBox.showinfo(title="gcode file generated", message="The gcode file has been successfully generated!")
