# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:42:12 2021

@author: IVitienes
"""

import numpy as np
import os, tkinter.filedialog
import matplotlib.pyplot as plt, tkinter as tk
from scipy import signal
import pandas as pd
import matplotlib

# Written by Isabela Vitienes 2021

##########################################################################################################
# This code analyzes a segment of strain data that contains a static event such as perching or sitting.  #
# Code takes as input an .xlsx file with columns for each gauge channel.                                 #
# User can select multiple files, and it will analyze each.                                              #
# The .xlsx file should only contain strain during the activity                                          #
# Areas of the code that should be tailored to each use are indicated by ***                             #
# The code is designed for data from 5 channels.                                                         #
# Outputs are saved in the same location as the input file.                                              #
##########################################################################################################


# LICENSE
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# For more details: https://www.gnu.org/licenses/gpl-3.0.html.



# Initial directory    ***
Initialdir = "S:\Projects_Isabela\Current_projects\Chickens\MyStudy\Data\Flock_three\strain_gauging\exported\background"


# User selects ROI file
Title = "Select input file"
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

strain_index = df[:,0]
strain_time = df[:,1]

# Rosette: R1, R2 (axial element), R3. A1 and A2 axial gauges. 
# *** Adjust indices according to columns in data.
R1 = df[:,4]
R2 = df[:,6]
R3 = df[:,8]
A1 = df[:,5]
A2 = df[:,7] 

no_gauges = 5 # *** Adjust according to number of gauges analyzed


# Filter strains
# *** Cutoff frequency selected to ensure that it is above the Nyquist frequency for the activities you have sampled. 

fs = 1000  # Sampling frequency, set when collecting data, in Hz
fc = 40  # Cut-off frequency of the filter, in Hz
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



# Plot data and prompt user to select points to guide analysis.

# All that matters is the x-coordinates of the points, not y.
# First two will flank the swing region (left, then right). Select one swing phase where the strains are steady. 
# Then select two points per step to flank peak strain region. These flanking points must encompass peak strains on all gauges.
# Finally, user selects two points to flank entire segment of consecutive steps. First point at the begining of the onset of step 1. Second at end of swing phase following final step.
# Total points = 4 + 2*number of steps
# Example can be found in Fig. 1 

# In the event that one gauge failed and the output is much higher than outputs from other gauges, it could increase the y-scale such that swing and peak areas are not distinguishable.
# In that case, exclude the failed gauge from the plot below and re-run the routine for this file. 
indices = [x - strain_index[0] for x in strain_index]
matplotlib.use('TkAgg')
shift = 500
plt.plot(indices, R1_filt, indices, R2_filt, indices, R3_filt, indices, A1_filt, indices, A2_filt)
plt.gca().legend(('R1', 'R2', 'R3', 'A1', 'A2'))
p = plt.ginput(-1, 300) # Indefinite number of inputs, terminate by 'enter' key
plt.clf()

points = [list(item) for item in p]
pts = np.array(points)

# Export ginputs
savepathginput = filepath1[:-5] + "_ginputs.xlsx"
dfginput = pd.DataFrame(pts)
dfginput.to_excel(savepathginput)



# Compute derived parameters

results = []

no_steps = int(((len(p)) - 4)/2)

px = pts[:,0] #x-coordinates, indices
px_swing = px[:2]
px_steps = px[2:-2]
px_stepfreq = px[-2:]


# Swing strain per gauge
for i in range(no_gauges):
    
    strain = strains_filt[:,i]
    
    a = int(px_swing[0])
    b = int(px_swing[1])
    
    results.append(np.mean((strain[a:b])))
    

# peak strain per step*gauge
for i in range(no_steps):
    for j in range(no_gauges):
       
        strain = strains_filt[:,j]
        
        a = int(px_steps[2*i])
        b = int(px_steps[2*i + 1])
        
        pk_max = np.max(strain[a:b])
        pk_min = np.min(strain[a:b])
        
        # Determine peak strain. 
        # Generalizable for tensile (-) and compressive (+) strains. Calculate min and max, and choose whichever is furthest from swing strain as the peak.
        d_max = abs(pk_max - results[j])
        d_min = abs(pk_min - results[j])
        if d_max > d_min:
            pk = pk_max
        elif d_max < d_min:
            pk=pk_min
        elif d_max == d_min:
            if pk_max == pk_min:
                pk = pk_max
            else:
                pk = np.nan
        results.append(pk)
        


# Group together swing strains, peak strains, and delta strain
# *** Set up for analyzing 5 gauges. Would need to alter indices for different amount of gauges.
results1 = np.empty([no_steps, 16])
for i in range(no_steps):
    results1[i,0:5] = results[0:5]
    results1[i, 5:10] = results[(i+1)*5:(i+1)*5 + 5]
    results1[i, 10] = results[(i+1)*5] - results[0]
    results1[i, 11] = results[(i+1)*5+1] - results[1]
    results1[i, 12] = results[(i+1)*5 + 2] - results[2]
    results1[i, 13] = results[(i+1)*5 + 3] - results[3]
    results1[i, 14] = results[(i+1)*5 + 4] - results[4]
    results1[i, 15] = results[-1]

results1 = np.matrix.transpose(np.array(results1))
    
# Each column is data from one step.
# Rows are: swing 1, ... swing 5, peak 1, ..., peak 5, delta 1, ..., delta 5, i.e. 3 x number of gauges
savepath = filepath1[:-5] + "_results_per_step.xlsx"
dfresults = pd.DataFrame(results1)
dfresults.to_excel(savepath)



# Average and standard deviation of delta strains across all steps within this segment
# Include variance across steps in a segment as this can be included in statistical model
rps = np.empty([1, 12])

rps[0,0] = np.mean(results1[:,10])
rps[0,5] = np.stdev(results1[:,10])
rps[0,1] = np.mean(results1[:,11])
rps[0,6] = np.stdev(results1[:,11])
rps[0,2] = np.mean(results1[:,12])
rps[0,7] = np.stdev(results1[:,12])
rps[0,3] = np.mean(results1[:,13])
rps[0,8] = np.stdev(results1[:,13])
rps[0,4] = np.mean(results1[:,14])
rps[0,9] = np.stdev(results1[:,14])


# Step frequency. 
# Peak strains are a function of locomotion speed, and therefore this is important to record.
step_freq = no_steps/((strain_index[int(px_stepfreq[1])] - strain_index[int(px_stepfreq[0])])/1000)

rps[0,10] = step_freq
rps[0,11] = no_steps


savepath = filepath1[:-5] + "_results_per_seg.xlsx"
dfrps = pd.DataFrame(rps)
dfrps.to_excel(savepath)
        
        
    
    
    
    
