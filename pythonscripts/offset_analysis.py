# Feb 2024
# calculate offsets.

import os, glob, time
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from scipy.optimize import curve_fit
import pathlib

def append_to_file(file_name, data):
    """
    appends measurement data for flipping experiment to existing csv file, or create new csv file if it
    is the first measurement in campaign.
    
    Parameters
    ----------
    file_name : string
        name of file to sppend to.
    data : dict
        dictionary of data.

    Returns
    -------
    None, writes to a file

    """
    # if file exists, append
    file_exists = os.path.isfile(file_name)
    fieldnames = ["# time",
                  "# coil Bx (pT) mean",
                  "# coil Bx (pT) std",
                  "# coil By (pT) mean",
                  "# coil By (pT) std",
                  "# coil Bz (pT) mean",
                  "# coil Bz (pT) std",
                  "# Orientation",
                  "# cell Bx (pT) mean",
                  "# cell Bx (pT) std",
                  "# cell By (pT) mean",
                  "# cell By (pT) std",
                  "# cell Bz (pT) mean",
                  "# cell Bz (pT) std",
                  ]  
    if file_exists:
        with open(file_name, 'a') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(data)
    else:
       df_file = pd.DataFrame(data, index=[0])
       df_file.to_csv(file_name, mode='a', index=False)

def find_start_idx(series, sr, stab_cond):
    '''
    given time series of compensation coils in pT, find idx when it stablized (<stab_cond pT/s change).
    '''      
    arr = np.asarray(series)
    last_idx = -1
    while True:
        new_idx = last_idx - int(sr)
        if abs(arr[new_idx] - arr[last_idx])/(int(sr)/sr) <= stab_cond:
            last_idx = new_idx
            continue
        else:
            return last_idx
    
#%% get mean and std of coil readings + cell readings compiled into separate csv file

sensor = 'AAL9'
axis = 'y'
date = 'jun06'
trial = 1
#%%
write_to_file_path = '../data/Offset/processed/{}/{}/'.format(sensor, date)
write_to_file_name = 'around_{}axis_trial{}.csv'.format(axis, trial)

plt.figure()
for file in glob.glob('..//data//Offset//raw//{}//around_{}//*COIL*.csv'.format(date, axis)):
    coil_df = pd.read_csv(file)
    coil_x = coil_df['x'].iloc[:-10]
    coil_y = coil_df['y'].iloc[:-10]
    coil_z = coil_df['z'].iloc[:-10]
    cell_df = pd.read_csv(file[:52]+'CELL'+file[56:])
    cell_x = cell_df['x']
    cell_y = cell_df['y']
    cell_z = cell_df['z']
    
    if file[-5] == 'p':
        orientation = 'up'
    elif file[-5] == 'n':
        orientation = 'down'
    
    # convert coil index to cell index
    coil_sr = 7.5
    cell_sr = 100
    
    x_start = len(coil_x) + find_start_idx(coil_x, coil_sr, 100)
    x_end = -1
    y_start = len(coil_y) + find_start_idx(coil_y, coil_sr, 100)
    y_end = -1
    z_start = len(coil_z) + find_start_idx(coil_z, coil_sr, 100)
    z_end = -1
    
    x_second_start = int(x_start / 7.5 * cell_sr)
    if x_end != -1: x_second_end = int(x_end / 7.5 * cell_sr)
    if x_end == -1: x_second_end = -1
    y_second_start = int(y_start / 7.5 * cell_sr)
    if y_end != -1: y_second_end =int(y_end / 7.5 * cell_sr)
    if y_end == -1: y_second_end = -1
    z_second_start = int(z_start / 7.5 * cell_sr)
    if z_end != -1: z_second_end =int(z_end / 7.5 * cell_sr)
    if z_end == -1: z_second_end = -1
    
    # setup dataframe madeup of relevant data of current measurement, to be appended to file.
    data = {"# time":file[28:33],
            "# coil Bx (pT) mean":np.mean(coil_x[x_start:x_end]),
            "# coil Bx (pT) std":(np.max(coil_x[x_start:x_end])-np.min(coil_x[x_start:x_end]))/2,
            "# coil By (pT) mean":np.mean(coil_y[y_start:y_end]),
            "# coil By (pT) std":(np.max(coil_y[y_start:y_end])-np.min(coil_y[y_start:y_end]))/2,
            "# coil Bz (pT) mean":np.mean(coil_z[z_start:z_end]),
            "# coil Bz (pT) std":np.std(coil_z[z_start:z_end]),
            "# Orientation":orientation,
            "# cell Bx (pT) mean":np.mean(cell_x[x_second_start:x_second_end]),
            "# cell Bx (pT) std":np.std(cell_x[x_second_start:x_second_end]),
            "# cell By (pT) mean":np.mean(cell_y[y_second_start:y_second_end]),
            "# cell By (pT) std":np.std(cell_y[y_second_start:y_second_end]),
            "# cell Bz (pT) mean":np.mean(cell_z[z_second_start:z_second_end]),
            "# cell Bz (pT) std":np.std(cell_z[z_second_start:z_second_end]),
            }
    
    plt.plot(coil_y)
    plt.vlines(y_start, -10000, 10000)
    # append to file
    append_to_file(write_to_file_path+write_to_file_name, data)
    

#%%
# perform offset calculations 
sensor = 'AAY4'
date = 'may31'
trial = 1

# which offset to get
axis = 'x'

if axis == 'x':
    allresults = pd.read_csv('../data/Offset/processed/{}/{}/around_{}axis_trial{}.csv'.format(sensor, date, 'y', trial))
    x_offsets = []
    for i in range(allresults.shape[0] - 1):
        x_offsets.append(1/2*(allresults['# coil Bx (pT) mean'].iloc[i]+allresults['# coil Bx (pT) mean'].iloc[i+1]))
    x_offsets = np.asarray(x_offsets)
    print('mean x offset = {}, std = {}'.format(np.mean(x_offsets), np.std(x_offsets)))

if axis == 'y':
    allresults = pd.read_csv('../data/Offset/processed/{}/{}/around_{}axis_trial{}.csv'.format(sensor, date, 'x', trial))
    y_offsets = []
    for i in range(allresults.shape[0] - 1):
        y_offsets.append(1/2*(allresults['# coil By (pT) mean'].iloc[i]+allresults['# coil By (pT) mean'].iloc[i+1]))
    y_offsets = np.asarray(y_offsets)
    print('mean y offset = {}, std = {}'.format(np.mean(y_offsets), np.std(y_offsets)))

if axis == 'z':
    allresults1 = pd.read_csv('../data/Offset/processed/{}/{}/around_{}axis_trial{}.csv'.format(sensor, date, 'x', trial))
    allresults2 = pd.read_csv('../data/Offset/processed/{}/{}/around_{}axis_trial{}.csv'.format(sensor, date, 'y', trial))
    z_offsets1, z_offsets2 = [], []
    for i in range(allresults1.shape[0] - 1):
        z_offsets1.append(1/2*(allresults1['# coil Bz (pT) mean'].iloc[i]+allresults1['# coil Bz (pT) mean'].iloc[i+1]))
    for i in range(allresults2.shape[0] - 1):
        z_offsets2.append(1/2*(allresults2['# coil Bz (pT) mean'].iloc[i]+allresults2['# coil Bz (pT) mean'].iloc[i+1]))    
    z_offsets1, z_offsets2 = np.asarray(z_offsets1), np.asarray(z_offsets2)
    print('mean z offset = {}, std = {}'.format(np.mean([np.mean(z_offsets1), np.mean(z_offsets2)]), 
                                                np.sqrt(np.std(z_offsets1)**2 + np.std(z_offsets2)**2)))


