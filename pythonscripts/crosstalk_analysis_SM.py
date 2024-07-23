# apr 2024
# crosstalk analysis with SM

import os, glob, time, sys
import matplotlib.pyplot as plt
from labjack import ljm 
from datetime import datetime
import ctypes
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
    
plt.rcParams.update({'font.size': 14})

#%%
# 1 cm 
df_off = pd.read_csv('..//data//apr08//StefanMayer//QZFMoff/240408T130649.csv')
df_2cm = pd.read_csv('..//data//apr08//StefanMayer//1cm//240408T122553.csv')

lowcut = 76
highcut = 78
sr = 1000

#%%
plt.figure()
a, b = signal.periodogram(df_off['SM'], sr)
c, d = signal.periodogram(df_2cm['SM'], sr)

plt.semilogy(c, np.sqrt(d), label='QZFM On')
plt.semilogy(a, np.sqrt(b), label='QZFM Off')
plt.legend()
plt.grid()
plt.show()

#%%
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


for fp in glob.glob('..//data//apr08//StefanMayer//*cm//'):
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

#%%
def inverse_cube(x, C):
    return C / x**3 



fig, ax = plt.subplots(1, 2, figsize=(1.75*6.4, 4.8))
distances = []
mean_amps = []
# inset
x1, x2, y1, y2 = 2, 7.7, -0.1, 0.5  # subregion of the original image
axins = ax[1].inset_axes(
    [0.4, 0.4, 0.47, 0.47],
    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

for key, value in amp_dict.items():
    distance = int(key[0])+0.5
    mean_amp = np.mean(value)
    binwidth = 0.03
    bins = np.arange(min(value), max(value) + binwidth, binwidth)
    alpha = 0.65
    ax[0].hist(value, bins, alpha=alpha, density=True, label=str(distance)+' cm')
    ax[1].errorbar(distance, mean_amp, np.std(value), fmt='*', markersize=12, color='blue')
    axins.errorbar(distance, mean_amp, np.std(value), fmt='*', markersize=12, color='blue')
    distances.append(distance)
    mean_amps.append(mean_amp)
# best fit inverse cube
popt, pcov = optimize.curve_fit(inverse_cube, np.array(distances)[2:], np.array(mean_amps)[2:])
ax[1].plot(np.linspace(0.5, 7.5, 100), inverse_cube(np.linspace(0.5, 7.5, 100), *popt), 'g--', label='Best Fit Inverse Cubed Field')
ax[1].set_ylim(-0.2, 3)



axins.plot(np.linspace(0.5, 7.5, 100), inverse_cube(np.linspace(0.5, 7.5, 100), *popt), 'g--')

ax[1].indicate_inset_zoom(axins, edgecolor="black", linewidth=5)



ax[0].legend(loc='upper center', ncols=3)
ax[1].legend()
ax[1].grid()
axins.grid()
ax[0].set_xlabel('Measured Amplitude (nT)')
ax[0].set_ylabel('Probability Density')
ax[1].set_xlabel('Distance From Cell (cm)')
ax[1].set_ylabel('Measured Amplitude (nT)')
plt.tight_layout()
fig.show()


#%% check what is the 77 Hz detected when QZFM is completely off
off_amp = []
for file in glob.glob('..//data//apr08//StefanMayer//QZFMoff//*'):
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
        off_amp.append(abs(extremas[i]-extremas[i+1]))








