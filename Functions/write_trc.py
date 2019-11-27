"""
Create and write a TRC file

Note that you can also export a TRC file using Mocka
File > Export > Motion Analysis Corp. > TRC file 

"""

__author__ = "Nathan Brantly, Alex Woodall"
__version__ = "2.0"
__license__ = "ABI"

import numpy as np
from numpy import matlib
import os
import tkinter as tk
from tkinter import filedialog

def trc_write(marker_labels, marker_data, full_file_name):
	'''
	'''

	# Set default header information values
	frame_rate = 200.0 # 200 Hz, common Vicon mocap camera rate
	units = 'mm' # millimeters, common Vicon mocap marker data units

	[], file_extention = os.path.splitext(full_file_name) # Get the extension of the filename

	if file_extention != '.trc':
		full_file_name = input("Please enter a full file name with the extension '.trc'") # Request a different filename from the user
	
	nframes = np.size(marker_data,axis=1) # Number of frames
	ncols = np.size(marker_data,axis=2) # Number of columns, including frame# and time (nmarkers*3 + 2)
	nmarkers = len(marker_labels) # Number of markers




