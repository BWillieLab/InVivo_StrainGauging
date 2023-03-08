# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 10:39:12 2021

@author: IVitienes
"""



import numpy as np
import os, tkinter.filedialog
import tkinter as tk
from scipy import signal
import pandas as pd


# Written by Isabela Vitienes 2021

##########################################################################################################
# This code analyzes a segment of strain data that contains a static event such as perching or sitting.  #
# Code takes as input an .xlsx file with columns for each gauge channel.                                 #
# User can select multiple files, and it will analyze each.                                              #
# The .xlsx file should only contain strain during the activity                                          #
# Areas of the code that should be tailored to each use are indicated by *                               #
# The code is designed for data from 5 channels.                                                         #
# Outputs are saved in the same location as the input file.                                              #
##########################################################################################################


# LICENSE
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# For more details: https://www.gnu.org/licenses/gpl-3.0.html.




# *** Initial directory
Initialdir = "S:\Projects_Isabela\Current_projects\Chickens\MyStudy\Data\Flock_three\strain_gauging\exported\background"   # *
Title = "Select background strain ROI file"


# User selects ROI file(s)

root = tk.Tk()
root.focus_force()
ftypes=[("excel", "*.xlsx")]
file_path = tk.filedialog.askopenfilenames(initialdir=Initialdir, parent=root,title=Title,filetypes=ftypes)
file_names = [os.path.splitext(os.path.basename(file))[0] for file in file_path]
files = [file_path,file_names]
root.withdraw()


# Analyze each input file

L = len(file_names)
for i in range(L):
    
    filepath1 = file_path[i]
    filename1 = file_names[i]


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

    # Calculate mean strains for each gauge
    mean_strain1 = np.mean(R1_filt)
    sd_strain1 = np.std(R1_filt)
    med_strain1 = np.median(R1_filt)
    max_strain1 = np.max(R1_filt)
    min_strain1 = np.min(R1_filt)
    
    mean_strain2 = np.mean(R2_filt)
    sd_strain2 = np.std(R2_filt)
    med_strain2 = np.median(R2_filt)
    max_strain2 = np.max(R2_filt)
    min_strain2 = np.min(R2_filt)
    
    mean_strain3 = np.mean(R3_filt)
    sd_strain3 = np.std(R3_filt)
    med_strain3 = np.median(R3_filt)
    max_strain3 = np.max(R3_filt)
    min_strain3 = np.min(R3_filt)
    
    mean_strain4 = np.mean(A1_filt)
    sd_strain4 = np.std(A1_filt)
    med_strain4 = np.median(A1_filt)
    max_strain4 = np.max(A1_filt)
    min_strain4 = np.min(A1_filt)
    
    mean_strain5 = np.mean(A2_filt)
    sd_strain5 = np.std(A2_filt)
    med_strain5 = np.median(A2_filt)
    max_strain5 = np.max(A2_filt)
    min_strain5 = np.min(A2_filt)
    
    
    # Compile and save results
    results = [mean_strain1, sd_strain1, med_strain1, min_strain1, max_strain1, mean_strain2, sd_strain2, med_strain2, min_strain2, max_strain2, mean_strain3, sd_strain3, med_strain3, min_strain3, max_strain3, mean_strain4, sd_strain4, med_strain4, min_strain4, max_strain4, mean_strain5, sd_strain5, med_strain5, min_strain5, max_strain5]
    results1 = np.matrix.transpose(np.array(results))
    savepath = filepath1[:-5] + "_results.xlsx"
    dfresults = pd.DataFrame(results1)
    dfresults.to_excel(savepath)
    
