# mar 2024
# analysis script for QZFM crosstalk

import os, glob, time, sys
import matplotlib.pyplot as plt
import numpy as np
from time import time
import time as time_
import pandas as pd
import pathlib
import scipy
from scipy import signal
from scipy import fft
from scipy import optimize
from scipy.signal import butter, sosfilt, sosfreqz

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y
    
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

plt.rcParams.update({'font.size': 14})
#%% Histogram plot
direction = 'y'
freq = 35
tc_idx = 4000
highcut = 36
lowcut = 34
sr = 1000
reject_threshold = 0.004
for fp in glob.glob('..//data//Crosstalk//gain_change_raw//x_axis//2cm//trial2//'):
    
    single_pp = []
    duo_pp = []
    
    sct = 0
    for file in glob.glob(fp+'*single*.csv'):
        df = pd.read_csv(file)
        
        B = df[direction]
        
        # find index of maximums and minimums
        max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
        min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
        

        
        period = sr / freq
        
        maximas = []
        minimas = []
        for idx in max_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    
                    maximas.append(max(yn(x)))
                except IndexError:
                    continue
        for idx in min_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    
                    
                    minimas.append(min(yn(x)))
                except IndexError:
                    continue
    
        
        
        extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
        for i in range(len(extremas)-1):
            single_pp.append(abs(extremas[i]-extremas[i+1]))
        
        sct += 1
      
    dct = 0
    for file in glob.glob(fp+'*duo*.csv'):
        
        df = pd.read_csv(file)
        
        B = df[direction]
        
        # find index of maximums and minimums
        max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
        min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
        
        
        period = sr / freq
        
        maximas = []
        minimas = []
        for idx in max_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    
                    maximas.append(max(yn(x)))
                except IndexError:
                    continue
        for idx in min_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    
                    minimas.append(min(yn(x)))
                except IndexError:
                    continue
        
        extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
        for i in range(len(extremas)-1):
            duo_pp.append(abs(extremas[i]-extremas[i+1]))  
        dct += 1



fig, axs = plt.subplots(1, 1, figsize=(6, 6))

legends = ['Perturbation Off', 'Perturbation On']     
for i, j in enumerate([single_pp, duo_pp]):
    data = np.divide(j, 1) 
    binwidth = 0.001
    bins = np.arange(min(data), max(data) + binwidth, binwidth)
    alpha = 0.75
    axs.hist(data, bins, alpha=alpha, density=True, label=legends[i])

axs.legend()
axs.set_xlabel('Peak-to-Peak Amplitude (nT)')
axs.set_ylabel('Density')
axs.set_xlim(6.35, 6.6)
axs.grid()
plt.tight_layout()
fig.show()



#%%
'''
Crosstalk results
'''
amp_dict = {'0cm':[],
            '1cm':[],
            '2cm':[],
            '3cm':[],
            '4cm':[],
            '5cm':[],
            '6cm':[],
            '7cm':[]}
 
lowcut = 76
highcut = 78
freq = 77
sr = 1000
tc_idx = 2000
 
def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos
 
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y
 
for fp in glob.glob('..//data//Crosstalk//StefanMayer//*cm//'):
    distance = fp[-4:-1]
    for file in glob.glob(fp+'*.csv'):
        df = pd.read_csv(file)
        B = df['SM']
        B = butter_bandpass_filter(B, lowcut, highcut, sr)[tc_idx:]
        # find index of maximums and minimums
        max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
        min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
        # plt.plot(B)
        # plt.plot(max_idx[0], B[max_idx[0]])
        period = sr / freq
        maximas = []
        minimas = []
        for idx in max_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    # plt.figure()
                    # plt.plot(x, y)
                    # plt.plot(x, yn(x))
                    # plt.show()
                    # time_.sleep(1)
                    #print(max(yn(x)))
                    maximas.append(max(yn(x)))
                except IndexError:
                    continue
        for idx in min_idx[0]:
            if idx > period // 4:
                try:
                    y = B[int(idx-period//4):int(idx+period//4)]
                    x = np.arange(len(y))
                    p = np.polyfit(x, y, 4)
                    yn = np.poly1d(p)
                    # plt.figure()
                    # plt.plot(x, y)
                    # plt.plot(x, yn(x))
                    # plt.show()
                    # time_.sleep(0.2)
                    minimas.append(min(yn(x)))
                except IndexError:
                    continue
        extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
        for i in range(len(extremas)-1):
            amp_dict[distance].append(abs(extremas[i]-extremas[i+1]) - 0.0825)
 
def inverse_cube(x, C):
    return C / x**3 

width = 0.7
height = 1.3
 
fig, ax = plt.subplots(2, 1, figsize=(width*6.4, height*4.8))
 
distances = []
mean_amps = []
# inset
x1, x2, y1, y2 = 2, 7, -0.1, 0.5  # subregion of the original image
axins = ax[0].inset_axes(
    [0.4, 0.4, 0.47, 0.47],
    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
 
 
for key, value in amp_dict.items():
    distance = int(key[0])+0.5
    mean_amp = np.mean(value)
    binwidth = 0.03
    bins = np.arange(min(value), max(value) + binwidth, binwidth)
    alpha = 0.65
    ax[0].errorbar(distance, mean_amp, np.std(value), fmt='*', markersize=12, color='blue')
    axins.errorbar(distance, mean_amp, np.std(value), fmt='*', markersize=12, color='blue')
    distances.append(distance)
    mean_amps.append(mean_amp)
# best fit inverse cube
popt, pcov = optimize.curve_fit(inverse_cube, np.array(distances)[2:], np.array(mean_amps)[2:])
ax[0].plot(np.linspace(0.5, 7.5, 100), inverse_cube(np.linspace(0.5, 7.5, 100), *popt), 'g--', label='Inverse Cubed Field')
ax[0].set_ylim(-0.2, 3)
 
 
axins.plot(np.linspace(0.5, 7.5, 100), inverse_cube(np.linspace(0.5, 7.5, 100), *popt), 'g--')
 
ax[0].indicate_inset_zoom(axins, edgecolor="black", linewidth=2)
 
 
ax[0].legend(loc='upper center')
ax[0].grid()
ax[0].set_xlim(left=0.2, right=7)
axins.grid()
 
#ax[0].set_xlabel('Distance (cm)')
ax[0].set_ylabel('Perturbation Amplitude (nT)')
 
 
directions = ['x', 'y', 'z']
 
holding_arr = np.empty((3, 6))
 
def fitfunc(x, a):
    return a*x**(-2)
 
ct = 0
 
shape_dict={'x':'*',
            'y':'^',
            'z':'v'}
 
for direction in directions:
    df = pd.read_csv('..//data//Crosstalk//Gain_change_processed//{}_direction.csv'.format(direction))
    x = np.asarray(df['distance'])
    y = np.asarray(np.divide(df['shift mean'], df['shift mean'].iloc[0]))
    holding_arr[ct, :] = y
    ax[1].plot(x, y, shape_dict[direction],label=direction, markersize=10)
    ct += 1
    # plt.show()
y = np.mean(holding_arr, axis=0)    
p = np.polyfit(x, y, 4)
yn = np.poly1d(p)
 
popt, pcov = optimize.curve_fit(fitfunc, np.asarray(df['distance']), np.mean(holding_arr, axis=0) )
ax[1].plot(np.linspace(1, 7, 100), fitfunc(np.linspace(1, 6, 100), *popt), '--')
#plt.plot(np.linspace(1, 6, 100), yn(np.linspace(1, 6, 100)), '--')    
ax[1].legend()
ax[1].set_xlim(left=0.2, right=7)
ax[1].grid()
ax[1].set_xlabel('Distance (cm)')
ax[1].set_ylabel(r'Gain Change Factor $\rho$')
 
plt.tight_layout()
fig.show()




