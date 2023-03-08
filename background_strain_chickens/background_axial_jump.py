# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 12:15:45 2021

@author: IVitienes
"""


import numpy as np
import os, tkinter.filedialog
import matplotlib.pyplot as plt, tkinter as tk
from scipy import signal
import pandas as pd
import matplotlib

# Written by Isabela Vitienes 2021
# Analysis of background strain data, specifically jumping

############################################################################
# This code analyzes a segment of strain data that contains a jump event.    #
# Code takes as input an .xlsx file with columns for each gauge channel.     #
# Areas of the code that should be tailored to each use are indicated by *** #
# The code is designed for data from 5 channels.                             #
# Outputs are saved in the same location as the input file selected.         #
############################################################################


# LICENSE
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# For more details: https://www.gnu.org/licenses/gpl-3.0.html.



# *** Initial directory for searching for input file
Initialdir = "S:\Projects_Isabela\Current_projects\Chickens\Data"
Title = "Select background strain input file"


# User selects input data
root = tk.Tk()
root.focus_force()
ftypes=[("excel", "*.xlsx")]
file_path = tk.filedialog.askopenfilenames(initialdir=Initialdir, parent=root,title=Title,filetypes=ftypes)
file_names = [os.path.splitext(os.path.basename(file))[0] for file in file_path]
files = [file_path,file_names]
root.withdraw()

filepath1 = file_path[0]
filename1 = file_names[0]


# Load data
df1 = pd.read_excel(filepath1)
df = df1.to_numpy()

# Rosette: R1, R2 (axial element), R3. A1 and A2 axial gauges. 
# *** Adjust indices according to columns in data.
strain_index = df[:,0]
strain_time = df[:,1]
R1 = df[:,4]
R2 = df[:,6]
R3 = df[:,8]
A1 = df[:,5]
A2 = df[:,7] 


# Filter strains
# *** Cutoff frequency selected to ensure that it is above the Nyquist frequency for the activities you have sampled. 
fs = 1000  # Sampling frequency
fc = 40  # Cut-off frequency of the filter
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(4, w, 'low')
R1_filt = signal.filtfilt(b, a, R1)
b, a = signal.butter(4, w, 'low')
R2_filt = signal.filtfilt(b, a, R2)
b, a = signal.butter(4, w, 'low')
R3_filt = signal.filtfilt(b, a, R3)
b, a = signal.butter(4, w, 'low')
A1_filt = signal.filtfilt(b, a, A1)
b, a = signal.butter(4, w, 'low')
A2_filt = signal.filtfilt(b, a, A2)

strains_filt = np.vstack([R1_filt, R2_filt, R3_filt, A1_filt, A2_filt])
strains_filt = np.matrix.transpose(strains_filt)

# Save filtered strains
savepathfilt = filepath1[:-5] + "_filt.xlsx"
dfilt = pd.DataFrame(strains_filt)
dfilt.to_excel(savepathfilt)


# Plots data and prompts user to select points to inform analysis.

# User selects 4 points.
# All that matters is the x-coordinates of the points, not y.
# First two will flank the jump landing. Make sure that peak strains at all gauges of interest are flanked. 
# Next two will select a swing region, either in the steps following the jump or the moment before landing when chicken is in the air. 

L = len(R1_filt)
indices = list(range(0,L))
matplotlib.use('TkAgg')
shift = 500
# In the event that one gauge failed and the output is much higher than outputs from other gauges, it could increase the y-scale such that swing and peak areas are not distinguishable.
# In that case, exclude the failed gauge from the plot below and re-run the routine for this file. 
plt.plot(indices, R1_filt, indices, R2_filt, indices, R3_filt, indices, A1_filt, indices, A2_filt)
plt.gca().legend(('R1', 'R2', 'R3', 'A1', 'A2'))
p = plt.ginput(-1, 300) # Indefinite number of inputs, terminate by 'enter' key
plt.clf()

points = [list(item) for item in p]
pts = np.array(points)
px = pts[:,0] #x-coordinates, indices

# Export ginputs
savepathginput = filepath1[:-5] + "_ginputs.xlsx"
dfginput = pd.DataFrame(pts)
dfginput.to_excel(savepathginput)



# Compute derived parameters
# *** Designed for 5 gauges.

p_pk_l = int(np.round(px[0]))
p_pk_r = int(np.round(px[1]))

pk_strain1_max = np.max(R1_filt[p_pk_l:p_pk_r])
pk_strain2_max = np.max(R2_filt[p_pk_l:p_pk_r])
pk_strain3_max = np.max(R3_filt[p_pk_l:p_pk_r])
pk_strain4_max = np.max(A1_filt[p_pk_l:p_pk_r])
pk_strain5_max = np.max(A2_filt[p_pk_l:p_pk_r])

pk_strain1_min = np.min(R1_filt[p_pk_l:p_pk_r])
pk_strain2_min = np.min(R2_filt[p_pk_l:p_pk_r])
pk_strain3_min = np.min(R3_filt[p_pk_l:p_pk_r])
pk_strain4_min = np.min(A1_filt[p_pk_l:p_pk_r])
pk_strain5_min = np.min(A2_filt[p_pk_l:p_pk_r])

sw_1 = np.mean(R1_filt[int(np.round(px[2])):int(np.round(px[3]))])
sw_2 = np.mean(R2_filt[int(np.round(px[2])):int(np.round(px[3]))])
sw_3 = np.mean(R3_filt[int(np.round(px[2])):int(np.round(px[3]))])
sw_4 = np.mean(A1_filt[int(np.round(px[2])):int(np.round(px[3]))])
sw_5 = np.mean(A2_filt[int(np.round(px[2])):int(np.round(px[3]))])

d_max_1 = abs(pk_strain1_max - sw_1)
d_max_2 = abs(pk_strain2_max - sw_2)
d_max_3 = abs(pk_strain3_max - sw_3)
d_max_4 = abs(pk_strain4_max - sw_4)
d_max_5 = abs(pk_strain5_max - sw_5)

d_min_1 = abs(pk_strain1_min - sw_1)
d_min_2 = abs(pk_strain2_min - sw_2)
d_min_3 = abs(pk_strain3_min - sw_3)
d_min_4 = abs(pk_strain4_min - sw_4)
d_min_5 = abs(pk_strain5_min - sw_5)


# Determine peak strains

# This method is invariant to strain sign (compressive or tensile)
# It will calculate max. and min. within the flanked region and select peak strain as the strain whose magnitude is furthest from the swing strain.

if d_max_1 > d_min_1:
    pk_strain1 = pk_strain1_max
elif d_max_1 < d_min_1:
        pk_strain1 = pk_strain1_min
elif d_max_1 == d_min_1:
    if pk_strain1_max == pk_strain1_min:
        pk_strain1 = pk_strain1_max
    else: 
        pk_strain1 = np.nan

if d_max_2 > d_min_2:
    pk_strain2 = pk_strain2_max
elif d_max_2 < d_min_2:
    pk_strain2 = pk_strain2_min
elif d_max_2 == d_min_2:
    if pk_strain2_max == pk_strain2_min:
        pk_strain2 = pk_strain2_max
    else: 
        pk_strain2 = np.nan
        
if d_max_4 > d_min_4:
    pk_strain4 = pk_strain4_max
elif d_max_4 < d_min_4:
    pk_strain4 = pk_strain4_min
elif d_max_4 == d_min_4:
    if pk_strain4_max == pk_strain4_min:
        pk_strain4 = pk_strain4_max
    else: 
        pk_strain4 = np.nan

if d_max_3 > d_min_3:
    pk_strain3 = pk_strain3_max
elif d_max_3 < d_min_3:
    pk_strain3 = pk_strain3_min
elif d_max_3 == d_min_3:
    if pk_strain3_max == pk_strain3_min:
        pk_strain3 = pk_strain3_max
    else: 
        pk_strain3 = np.nan

if d_max_5 > d_min_5:
    pk_strain5 = pk_strain5_max
elif d_max_5 < d_min_5:
    pk_strain5 = pk_strain5_min
elif d_max_5 == d_min_5:
    if pk_strain5_max == pk_strain5_min:
        pk_strain5 = pk_strain5_max
    else: 
        pk_strain5 = np.nan


# Zero the peak strains with the swing strains
del1 = pk_strain1 - sw_1
del2 = pk_strain2 - sw_2
del3 = pk_strain3 - sw_3
del4 = pk_strain4 - sw_4
del5 = pk_strain5 - sw_5


# Save results
filepath = file_path[0]
results = [pk_strain1, pk_strain2, pk_strain3, pk_strain4, pk_strain5, sw_1, sw_2, sw_3, sw_4, sw_5, del1, del2, del3, del4, del5]
results1 = np.matrix.transpose(np.array(results))
savepath = filepath[:-5] + "_results.xlsx"
dfresults = pd.DataFrame(results1)
dfresults.to_excel(savepath)







