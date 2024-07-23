# Feb 2024
# to import, fine tune and analysis offset of x and z of quspin based on data collected on Feb 20 2024.

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
axis = 'x'
date = 'jun06'
trial = 2
#%%
write_to_file_path = '..//results//{}//{}//'.format(sensor, date)
write_to_file_name = 'around_{}axis_trial{}.csv'.format(axis, trial)


for file in glob.glob('..//data//{}//{}_around_{}axis//*COIL*.csv'.format(date, sensor, axis)):
    coil_df = pd.read_csv(file)
    coil_x = coil_df['x']
    coil_y = coil_df['y']
    coil_z = coil_df['z']
    cell_df = pd.read_csv(file[:48]+'CELL'+file[52:])
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

    # append to file
    append_to_file(write_to_file_path+write_to_file_name, data)
    

# perform offset calculations 

allresults = pd.read_csv('../results/{}/{}/around_{}axis_trial{}.csv'.format(sensor, date, axis, trial))

#%%

if axis == 'z':
    x_offsets, y_offsets = [], []
    x_err, y_err = [], []
    
    for i in range(allresults.shape[0] - 1):
        x_offsets.append(1/2*(allresults['# coil Bx (pT) mean'].iloc[i]+allresults['# coil Bx (pT) mean'].iloc[i+1]))
        y_offsets.append(1/2*(allresults['# coil By (pT) mean'].iloc[i]+allresults['# coil By (pT) mean'].iloc[i+1]))
        x_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil Bx (pT) std'].iloc[i]**2+allresults['# coil Bx (pT) std'].iloc[i+1]**2)]))
        y_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil By (pT) std'].iloc[i]**2+allresults['# coil By (pT) std'].iloc[i+1]**2)]))
    x_offsets, y_offsets = np.asarray(x_offsets), np.asarray(y_offsets)
    x_err, y_err = np.asarray(x_err), np.asarray(y_err)
    print('mean x offset = {}, std = {}'.format(np.mean(x_offsets), np.std(x_offsets)))
    print('mean y offset = {}, std = {}'.format(np.mean(y_offsets), np.std(y_offsets)))
elif axis == 'y':
    x_offsets, z_offsets = [], []
    x_err, z_err = [], []
    
    for i in range(allresults.shape[0] - 1):
        x_offsets.append(1/2*(allresults['# coil Bx (pT) mean'].iloc[i]+allresults['# coil Bx (pT) mean'].iloc[i+1]))
        z_offsets.append(1/2*(allresults['# coil Bz (pT) mean'].iloc[i]+allresults['# coil Bz (pT) mean'].iloc[i+1]))
        x_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil Bx (pT) std'].iloc[i]**2+allresults['# coil Bx (pT) std'].iloc[i+1]**2)]))
        z_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil Bz (pT) std'].iloc[i]**2+allresults['# coil Bz (pT) std'].iloc[i+1]**2)]))
    x_offsets, z_offsets = np.asarray(x_offsets), np.asarray(z_offsets)
    x_err, z_err = np.asarray(x_err), np.asarray(z_err)
    print('mean x offset = {}, std = {}'.format(np.mean(x_offsets), np.std(x_offsets)))
    print('mean z offset = {}, std = {}'.format(np.mean(z_offsets), np.std(z_offsets)))
elif axis == 'x':
    y_offsets, z_offsets = [], []
    y_err, z_err = [], []
    
    for i in range(allresults.shape[0] - 1):
        y_offsets.append(1/2*(allresults['# coil By (pT) mean'].iloc[i]+allresults['# coil By (pT) mean'].iloc[i+1]))
        z_offsets.append(1/2*(allresults['# coil Bz (pT) mean'].iloc[i]+allresults['# coil Bz (pT) mean'].iloc[i+1]))
        y_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil By (pT) std'].iloc[i]**2+allresults['# coil By (pT) std'].iloc[i+1]**2)]))
        z_err.append(np.max([np.sqrt(2)*100,
                             1/2*np.sqrt(allresults['# coil Bz (pT) std'].iloc[i]**2+allresults['# coil Bz (pT) std'].iloc[i+1]**2)]))
    y_offsets, z_offsets = np.asarray(y_offsets), np.asarray(z_offsets)
    y_err, z_err = np.asarray(y_err), np.asarray(z_err)
    print('mean y offset = {}, std = {}'.format(np.mean(y_offsets), np.std(y_offsets)))
    print('mean z offset = {}, std = {}'.format(np.mean(z_offsets), np.std(z_offsets)))
#%%
a = np.asarray(allresults['# coil Bx (pT) mean'].iloc[1::2]).astype(float)
b = np.asarray(allresults['# coil Bx (pT) mean'].iloc[::2]).astype(float)


plt.figure(figsize=(16,9))
for i in range(allresults.shape[0] - 1):
    if i == 0:
        plt.vlines(i, 
                   np.minimum(allresults['# coil Bx (pT) mean'].iloc[i], allresults['# coil Bx (pT) mean'].iloc[i+1]),
                   np.maximum(allresults['# coil Bx (pT) mean'].iloc[i], allresults['# coil Bx (pT) mean'].iloc[i+1]),
                   linestyles='dashed')
        plt.plot(i, 1/2*(allresults['# coil Bx (pT) mean'].iloc[i]+allresults['# coil Bx (pT) mean'].iloc[i+1]), 'r.', label='x offset = (up+down)/2')
        plt.plot(i, allresults['# coil Bx (pT) mean'].iloc[i], 'b.', label='Orientation Up')
        plt.plot(i, allresults['# coil Bx (pT) mean'].iloc[i+1], 'g.', label='Orientation Down')
    else:
        plt.vlines(i, 
                   np.minimum(allresults['# coil Bx (pT) mean'].iloc[i], allresults['# coil Bx (pT) mean'].iloc[i+1]),
                   np.maximum(allresults['# coil Bx (pT) mean'].iloc[i], allresults['# coil Bx (pT) mean'].iloc[i+1]),
                   linestyles='dashed')
        plt.plot(i, 1/2*(allresults['# coil Bx (pT) mean'].iloc[i]+allresults['# coil Bx (pT) mean'].iloc[i+1]), 'r.')
        plt.plot(i, allresults['# coil Bx (pT) mean'].iloc[i], 'b.')
        plt.plot(i, allresults['# coil Bx (pT) mean'].iloc[i+1], 'g.')
plt.plot(np.arange(allresults.shape[0] -1 ), np.full((allresults.shape[0] -1 ), np.mean(x_offsets)), 'k--', label='Averaged Offset')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
plt.ylabel('B (pT)')
plt.xlabel('Set Number')
plt.show()
plt.savefig('xoffsets.pdf')

#%%
plt.figure(figsize=(16,9))
for i in range(allresults.shape[0] // 2):
    if i == 0:
        plt.vlines(i, 
                   np.minimum(allresults['# coil By (pT) mean'].iloc[2*i], allresults['# coil By (pT) mean'].iloc[2*i+1]),
                   np.maximum(allresults['# coil By (pT) mean'].iloc[2*i], allresults['# coil By (pT) mean'].iloc[2*i+1]),
                   linestyles='dashed')
        plt.plot(i, 1/2*(allresults['# coil By (pT) mean'].iloc[2*i]+allresults['# coil By (pT) mean'].iloc[2*i+1]), 'r.', label='y offset = (up+down)/2')
        plt.plot(i, allresults['# coil By (pT) mean'].iloc[2*i], 'b.', label='Orientation Up')
        plt.plot(i, allresults['# coil By (pT) mean'].iloc[2*i+1], 'g.', label='Orientation Down')
    else:
        plt.vlines(i, 
                   np.minimum(allresults['# coil By (pT) mean'].iloc[2*i], allresults['# coil By (pT) mean'].iloc[2*i+1]),
                   np.maximum(allresults['# coil By (pT) mean'].iloc[2*i], allresults['# coil By (pT) mean'].iloc[2*i+1]),
                   linestyles='dashed')
        plt.plot(i, 1/2*(allresults['# coil By (pT) mean'].iloc[2*i]+allresults['# coil By (pT) mean'].iloc[2*i+1]), 'r.')
        plt.plot(i, allresults['# coil By (pT) mean'].iloc[2*i], 'b.')
        plt.plot(i, allresults['# coil By (pT) mean'].iloc[2*i+1], 'g.')
plt.plot(np.arange(allresults.shape[0] // 2), np.full((allresults.shape[0] // 2), np.mean(y_offsets)), 'k--', label='Averaged Offset')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
plt.ylabel('B (pT)')
plt.xlabel('Set Number')
plt.show()
plt.savefig('yoffsets.pdf')


#%% box plot
x_a_z = x_offsets
#%%
x_a_y = x_offsets
#%%
#  data
x_data = [x_a_y,x_a_z]
xlabels = [ 'Rotate around y','Rotate around z']

#%%
y_a_x = y_offsets
#%%
y_a_z = y_offsets 
#%%
y_data = [y_a_x, y_a_z]
ylabels = [ 'Rotate around x','Rotate around z']


#%%
z_a_x = z_offsets
#%%
z_a_y = z_offsets
#%%
z_data = [z_a_x, z_a_y]
zlabels = [ 'Rotate around x','Rotate around y']



#%%
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))

# rectangular box plot
bplot1 = ax1.boxplot(x_data,
                     vert=True,  # vertical box alignment
                     patch_artist=False,  # fill with color
                     labels=xlabels)  # will be used to label x-ticks
ax1.set_title('x Offsets')



# adding horizontal grid lines

ax1.yaxis.grid(True)
ax1.set_xlabel('')
ax1.set_ylabel(r'B (pT)')

# rectangular box plot
bplot2 = ax2.boxplot(y_data,
                     vert=True,  # vertical box alignment
                     patch_artist=False,  # fill with color
                     labels=ylabels)  # will be used to label x-ticks
ax2.set_title('y Offsets')



# adding horizontal grid lines

ax2.yaxis.grid(True)
ax2.set_xlabel('')

# rectangular box plot
bplot3 = ax3.boxplot(z_data,
                     vert=True,  # vertical box alignment
                     patch_artist=False,  # fill with color
                     labels=zlabels)  # will be used to label x-ticks
ax3.set_title('z Offsets')



# adding horizontal grid lines

ax3.yaxis.grid(True)
ax3.set_xlabel('')

fig.suptitle('Offsets Measured in 2 Ways')
plt.tight_layout()
plt.show()


