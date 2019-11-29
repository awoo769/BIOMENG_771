
import numpy as np
from numpy import matlib
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from scipy import signal

import csv

from write_motion_file import write_motion_file

# This code will open up an existing force plate csv data file and preprocess the forces ready for inverse dynamics in OpenSim

# read in some data from a csv file

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(initialdir = "r",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))

grf_data = []

with open(file_path) as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')

	for row in readCSV:
		grf_data.append(row)

# Close the csv file
csvfile.close()

# First line is the header
grf_data = grf_data[1:]

# Convert to array of floats
grf_data = np.array(grf_data)

for i in range(np.shape(grf_data)[0]):
	for j in range(np.shape(grf_data)[1]):
		if grf_data[i,j] == "": # If it is empty, make 0
			grf_data[i,j] = "0"

grf_data = grf_data.astype(np.float)

# Get data from the csv file data

time = grf_data[:,0] # time

# Force plate 1
force_plate1_forces = grf_data[:,1:4]
force_plate1_COP = grf_data[:,4:7] / 1000
force_plate1_Tz = grf_data[:,7]

# Force plate 2
force_plate2_forces = grf_data[:,8:11]
force_plate2_COP = grf_data[:,11:14] / 1000
force_plate2_Tz = grf_data[:,14]

# Force plate 3
force_plate3_forces = grf_data[:,15:18]
force_plate3_COP = grf_data[:,18:21] / 1000
force_plate3_Tz = grf_data[:,21]

# Cut-off frequency for low pass Butterworth filter
cut_off_frequency = 15 # This should be the same as what you filtered the kinematic data with for dynamic consistancy
analog_rate = 1200

# Filter the forces - do this just for force plate 2 as an example
Wn = cut_off_frequency/(analog_rate/2)

# Describe filter characteristics using 4th order Butterworth filter
b, a = signal.butter(4, Wn)

force_plate2_forces_filt = signal.filtfilt(b, a, force_plate2_forces,axis=0)
force_plate2_Tz_filt = signal.filtfilt(b, a, force_plate2_Tz,axis=0)

# Plot the filtered and unfiltered forces as a visual check
plt.plot(force_plate2_forces,'r')
plt.plot(force_plate2_forces_filt,'b')

plt.show()

# Filter the CoP - note that the CoP should really be re-calculated after filtering the forces
# This is a hack, since we are not going to recalculate the CoP
cut_off_frequency = 50
analog_rate = 1200

Wn = cut_off_frequency/(analog_rate/2)

b, a = signal.butter(2, Wn)
force_plate2_COP_filt = signal.filtfilt(b, a, force_plate2_COP,axis=0)

plt.plot(force_plate2_COP,'r')
plt.plot(force_plate2_COP_filt,'b')

plt.show()

# Find when the vertical forces drop below threshold, and then make all of the forces
# and CoP values 0.0 at these points. This is because the filtering creates residual effects
# at the end points
force_threshold = 20 # Set this to 20 N
force_zero = np.where(force_plate2_forces[:,2] < force_threshold)

force_plate2_forces_filt[force_zero,:] = 0.0
force_plate2_COP_filt[force_zero,:] = 0.0
force_plate2_Tz_filt[force_zero] = 0.0

# Make sure that torque data are in a 3xn array. Note that the x and y torques are zero.
# Only Tz is non-zero
force_plate2_torques = np.zeros((len(force_plate2_Tz_filt), 3))
force_plate2_torques[:,2] = force_plate2_Tz_filt

force_data = np.concatenate((force_plate2_forces_filt, force_plate2_COP_filt, force_plate2_torques),axis=1)
nt, nc = np.shape(force_data)

# First, rotate 90deg about X so that Y is the vertical axis
rot = np.array([(1, 0, 0), (0, 0, 1), (0, -1, 0)]) # Rotation matrix

if nc % 3 != 0:
	print("Error: force columns must have 3 components each")

rotated_force_data = np.zeros(np.shape(force_data))

# Rotate the marker data
for i in range(int(nc/3)):
    rotated_force_data[:,(3*(i+1)-3):3*(i+1)] = np.matmul(rot, force_data[:, (3*(i+1)-3):3*(i+1)].T).T

# Concatenate time with rotated forces and torques into one array
grf_complete = np.concatenate((time[:,np.newaxis], rotated_force_data),axis=1)

# Now output data to a motion file
write_motion_file(grf_complete, file_path)
