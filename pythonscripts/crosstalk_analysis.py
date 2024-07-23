# mar 2024
# analysis script for QZFM crosstalk

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
    
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

plt.rcParams.update({'font.size': 14})

#%%
# single = pd.read_csv('..//data//mar25//AAL9//sensitive_axis_change//single1.csv')
# duo = pd.read_csv('..//data//mar25//AAL9//sensitive_axis_change//duo.csv')
# single2 = pd.read_csv('..//data//mar25//AAL9//sensitive_axis_change//single2.csv')
freq = 35
tc_idx = 4000
highcut = 36
lowcut = 34
sr = 1000
reject_threshold = 0.004


fp = '..//data//apr04//AAL9//sensitive_axis_change//z_axis//1cm//trial1//'

def center_around_zero(arr):
    return arr - np.mean(arr)

def lin_func(x, m, b):
    return m * x + b

for direction in ['z']:
    
    directiondata = []
    directionlabel = []
    
    single_pp = []
    duo_pp = []
    single_again_pp = []
    
    sct = 0
    for file in glob.glob(fp+'*single*.csv'):

        df = pd.read_csv(file)
        # head, tail = os.path.split(file)
        
        B = df[direction]
        B = butter_bandpass_filter(B, lowcut, highcut, sr, 3)[tc_idx:]
        
        # find index of maximums and minimums
        max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
        min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
        
        # plt.figure()
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
        
        # get interweaved peak values
        # maximas = B[max_idx[0]]
        # minimas = B[min_idx[0]]
        
        # envelope check
        popt, pcov = optimize.curve_fit(lin_func, np.arange(len(maximas))/len(maximas), maximas)
        if abs(popt[0]) >= reject_threshold:
            print('rejected single file {}'.format(sct))
            sct += 1
            continue
        
        
        extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
        for i in range(len(extremas)-1):
            single_pp.append(abs(extremas[i]-extremas[i+1]))
        
        sct += 1
      
    dct = 0
    for file in glob.glob(fp+'*duo*.csv'):
        
        df = pd.read_csv(file)
        # head, tail = os.path.split(file)
        
        B = df[direction]
        B = butter_bandpass_filter(B, lowcut, highcut, sr, 3)[tc_idx:]
        
        # find index of maximums and minimums
        max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
        min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
        
        # plt.figure()
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
                    # time_.sleep(0.2)
                    
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
        # get interweaved peak values
        # maximas = B[max_idx[0]]
        # minimas = B[min_idx[0]]
        
        popt, pcov = optimize.curve_fit(lin_func, np.arange(len(maximas))/len(maximas), maximas)
        if abs(popt[0]) >= reject_threshold:
            print('rejected duo file {}'.format(dct))
            dct += 1
            continue
        
        extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
        for i in range(len(extremas)-1):
            duo_pp.append(abs(extremas[i]-extremas[i+1]))  
        dct += 1

ori_mean = np.mean(single_pp)
ori_std = np.std(single_pp)

new_mean = np.mean(duo_pp)

num_std_away = abs(new_mean-ori_mean) / ori_std
#print(num_std_away)
#%%
plt.figure(figsize=(8, 4))
legends = ['Only AAL9', 'Both AAL9 and AAY4']
for i, j in enumerate([single_pp, duo_pp]):
    data = np.divide(j, 2) 
    binwidth = 0.0001
    bins = np.arange(min(data), max(data) + binwidth, binwidth)
    alpha = 0.75
    plt.hist(data, bins, alpha=alpha, density=True, label=legends[i])
#plt.title('Sensor Taped Together')
plt.xlabel('Measured B (nT)')
plt.ylabel('PDF')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()


#%% calculate number of stds away from original mean
ori_mean = np.mean(single_pp)
ori_std = np.std(single_pp)

new_mean = np.mean(duo_pp)

num_std_away = abs(new_mean-ori_mean) / ori_std
#print(num_std_away)
print(abs(new_mean-ori_mean))
#%% ks test
print(scipy.stats.kstest(single_pp, duo_pp, alternative='less', method='asymp'))


#%%
'''
'''
freq = 35
tc_idx = 4000
highcut = 36
lowcut = 34
sr = 1000
reject_threshold =  0.005

def center_around_zero(arr):
    return arr - np.mean(arr)

def lin_func(x, m, b):
    return m * x + b

shifts_1 = {'x':[],
          'y':[],
          'z':[]}
shifts_2 = {'x':[],
          'y':[],
          'z':[]}    
shifts_3 = {'x':[],
          'y':[],
          'z':[]}


#%%

for direction in ['y']:
    for fp in glob.glob('..//data//mar27//AAL9//sensitive_axis_change//x_axis//2cm//trial2//'):
       
        shift_dict = shifts_1
        
        single_pp = []
        duo_pp = []
        
        sct = 0
        for file in glob.glob(fp+'*single*.csv'):
            df = pd.read_csv(file)
            # head, tail = os.path.split(file)
            
            B = df[direction]
            #B = butter_bandpass_filter(B, lowcut, highcut, sr, 3)[tc_idx:]
            
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
            
            # get interweaved peak values
            # maximas = B[max_idx[0]]
            # minimas = B[min_idx[0]]
            
            # envelope check
            # popt, pcov = optimize.curve_fit(lin_func, np.arange(len(maximas))/len(maximas), maximas)
            # if abs(popt[0]) >= reject_threshold:
            #     #print('rejected this single file {}'.format(sct))
            #     sct += 1
            #     continue
            
            
            extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
            for i in range(len(extremas)-1):
                single_pp.append(abs(extremas[i]-extremas[i+1]))
            
            sct += 1
          
        dct = 0
        for file in glob.glob(fp+'*duo*.csv'):
            
            df = pd.read_csv(file)
            # head, tail = os.path.split(file)
            
            B = df[direction]
            #B = butter_bandpass_filter(B, lowcut, highcut, sr, 3)[tc_idx:]
            
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
                        # time_.sleep(0.2)
                        
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
            # get interweaved peak values
            # maximas = B[max_idx[0]]
            # minimas = B[min_idx[0]]
            
            # popt, pcov = optimize.curve_fit(lin_func, np.arange(len(maximas))/len(maximas), maximas)
            # if abs(popt[0]) >= reject_threshold:
            #     #print('rejected this duo file {}'.format(dct))
            #     dct += 1
            #     continue
            
            extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
            for i in range(len(extremas)-1):
                duo_pp.append(abs(extremas[i]-extremas[i+1]))  
            dct += 1




#%%


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
axs.grid()
plt.tight_layout()
fig.show()


#%%
ori_mean = np.mean(single_pp)
ori_std = np.std(single_pp)

new_mean = np.mean(duo_pp)

num_std_away = abs(new_mean-ori_mean) / ori_std

#print(abs(new_mean-ori_mean))

shift_dict[direction].append(scipy.stats.ttest_ind(single_pp, duo_pp, equal_var=False)[1])
#print(scipy.stats.kstest(single_pp, duo_pp, alternative='less', method='asymp'))

#print(scipy.stats.ttest_ind(single_pp, duo_pp, equal_var=False, alternative='greater')[1])
#print(scipy.stats.kstest(single_pp, duo_pp, alternative='less', method='asymp')[0])
#print(scipy.stats.kstest(single_pp, duo_pp, alternative='greater', method='asymp')[1])
print(scipy.stats.kstest(single_pp[::20], duo_pp[::20], alternative='two-sided', method='asymp')[1])
# # see if new distribution shifted to left or to right
# if np.median(single_pp) <= np.median(duo_pp):
#shifts_3[direction].append(scipy.stats.kstest(single_pp, duo_pp, alternative='less', method='asymp')[0])
# else:
#     shifts_3[direction].append(scipy.stats.kstest(single_pp, duo_pp, alternative='greater', method='asymp')[1])
        
#%% for x and z
plt.figure()
use_dicts = [shifts_1[direction], shifts_2[direction], shifts_3[direction]]


for direction in ['z']:
    
    
    plt.plot(np.arange(len(use_dicts[0]))+1,
             np.mean(use_dicts,
                     axis=0) , 
             '*--',
             label=direction)
plt.yscale('log')
plt.legend()
plt.grid()
plt.xlabel('Distance Apart (cm)')
plt.ylabel('Amplitude Difference (nT)')
plt.tight_layout()
plt.show()


#%% save results 

direction = 'y'
name = '{}_direction'.format(direction)
starting_distance = 2

save_dict = {'distance':np.arange(len(use_dicts[0]))+starting_distance,
             'shift mean':np.mean(use_dicts, axis=0),
             'shift std':np.std(use_dicts, axis=0),
           }

save_df = pd.DataFrame(save_dict)
save_df.to_csv('..//results//{}.csv'.format(name))


#%%
plt.figure()
plt.plot(range(len(single_pp)), single_pp, '*')
plt.plot(range(+len(single_pp), len(single_pp)+len(duo_pp)), duo_pp, '*')
plt.show()
#%%
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
axs[0].plot(range(len(single_pp)), single_pp, '*', label='Only AAL9 Powered On')
axs[0].plot(range(+len(single_pp), len(single_pp)+len(duo_pp)), duo_pp, '*', label='Both AAL9 AAY4 Powered On')
axs[0].set_ylabel('Measured Peak-to-Peak Amplitude (nT)')

axs[0].grid()
axs[0].legend()
legends = ['Only AAL9', 'Both AAL9 and AAY4']     
for i, j in enumerate([single_pp, duo_pp]):
    data = np.divide(j, 1) 
    binwidth = 0.001
    bins = np.arange(min(data), max(data) + binwidth, binwidth)
    alpha = 0.75
    axs[1].hist(data, bins, alpha=alpha, density=True, label=legends[i])
axs[1].set_xlabel('Measured Peak-to-Peak Amplitude (nT)')
axs[1].set_ylabel('PDF')
axs[1].grid()
plt.tight_layout()
fig.show()


#%%
directions = ['x', 'y', 'z']

plt.figure(figsize=(1.2*6.4, 1.2*4.8))
holding_arr = np.empty((3, 6))

def fitfunc(x, a):
    return a*x**(-2)

ct = 0
for direction in directions:
    df = pd.read_csv('..//results//{}_direction.csv'.format(direction))
    x = np.asarray(df['distance'])
    y = np.asarray(np.divide(df['shift mean'], df['shift mean'].iloc[0]))
    holding_arr[ct, :] = y
    plt.plot(x, y, '*',label=direction, markersize=10)
    ct += 1
    
    # plt.show()
y = np.mean(holding_arr, axis=0)    
p = np.polyfit(x, y, 4)
yn = np.poly1d(p)

popt, pcov = optimize.curve_fit(fitfunc, np.asarray(df['distance']), np.mean(holding_arr, axis=0) )
plt.plot(np.linspace(1, 6, 100), fitfunc(np.linspace(1, 6, 100), *popt), '--')
#plt.plot(np.linspace(1, 6, 100), yn(np.linspace(1, 6, 100)), '--')    
plt.legend()
plt.grid()
plt.xlabel('Baseline Distance (cm)')
plt.ylabel(r'$\rho$')
plt.tight_layout()
plt.show()






