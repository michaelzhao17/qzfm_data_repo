# mar 2024
# measurement script for QZFM crosstalk

from QZFM import QZFM
import os, glob, time, sys
import matplotlib.pyplot as plt
import queue
from threading import Thread
from labjack import ljm 
from datetime import datetime
import ctypes
import numpy as np
from time import time
import time as time_
import pandas as pd
from labjackmeasure import labjack_measure
import pathlib
from scipy import signal
from scipy import fft
from SiglentDevices import DG1032Z
from scipy.signal import butter, sosfilt, sosfreqz

def make_folder(fp, folder_name):
    '''
    Parameters
    ----------
    fp : str
        file path to location where folder is to be created
    axis : str, x|y|z
        axis of rotation.
    sensor : str
        Name of sensor being measured e.g., 'AAY4'.

    Returns
    -------
    folder_name : str
        name of folder
    
    creates folder and then returns folder name
    '''
    # check if folder already exists, if not create it
    pathlib.Path(fp+folder_name).mkdir(parents=True, exist_ok=True) 
    return folder_name

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
    

#%%
# Digital Function Generator
awg = DG1032Z(hostname='USB0::0x1AB1::0x0642::DG1ZA232603182::INSTR')
awg.query('*IDN?')

#%%
AAL9 = QZFM("COM3")  

#%%
AAY4 = QZFM("COM4")
#%%
AAL9.auto_start(zero_calibrate=False)
#%%
AAY4.auto_start(zero_calibrate=False)
#%%
gain = '0.33x'
#%%
AAL9.set_gain(gain)
#%%
AAY4.set_gain(gain)
#%%
# sample at 1kHz
sr = 1000
zero_time = 15
measure_time = 10
save = True
gain_dict = {'0.1x':0.27,
             '0.33x':0.9,
             '1x':2.7,
             '3x':8.1}
# primary sensor
sensor = AAL9
#%%


for trial in range(3):

    fp = '..//data//apr08//'
    folder_name = 'StefanMayer'
    folder_name = make_folder(fp, folder_name)
    
    phases = ['single', 'duo']
    
    n_measurement = 5
    for phase in phases:
        print('currently on {} phase'.format(phase))
        for i in range(n_measurement):
            print('performing measurement no.{} out of {}'.format(i+1, n_measurement))
            if i == 0:
                # zero
                sensor.field_zero(True, True, False)
                for i in range(zero_time):
                    time_.sleep(1)
                sensor.field_zero(False)
                time_.sleep(40)
            # turn on AWG
            awg.set_ch_state(1, state=True)
            time_.sleep(0.1)
            # measure
            out = labjack_measure(measure_time, sr, ["AIN0", "AIN1", "AIN2"], [gain_dict[gain], gain_dict[gain], gain_dict[gain]], 
                                                      [10.0, 10.0, 10.0], False)
            
            # turn off AWG
            awg.set_ch_state(1, state=False)
            
            out_df = pd.DataFrame(out.T)
            out_df.columns = ['Epoch Time', 'x', 'y', 'z']
            out_df.set_index("Epoch Time", inplace=True)
            
            time = datetime.now().strftime('%y%m%dT%H%M%S')
            # save
            if save:
                out_df.to_csv(fp+folder_name+'//'+time+phase+".csv")
        input('Just finished {} phase, turn on/off second sensor'.format(phase))
        

#%%
lowcut = 76
highcut = 78


for axis in range(1):
    plt.figure()
    plt.plot(out[axis+1, :])
    plt.plot(butter_bandpass_filter(out[axis+1, :], lowcut, highcut, sr, 5))
    plt.grid()
    plt.show()
#%%
plt.figure()
for axis in range(1):
    a, b = signal.periodogram(out[axis+1, :], sr)
    plt.semilogy(a, np.sqrt(b), label=axis)
plt.legend()
plt.grid()
plt.show()

#%%
no_q = out

#%%

plt.figure()

a, b = signal.periodogram(out[1, :], sr)
plt.semilogy(a, np.sqrt(b), label='Quspin Off')

c, d = signal.periodogram(no_q[1, :], sr)
plt.semilogy(c, np.sqrt(d), label='Quspin On')
plt.legend()
plt.grid()
plt.show()



#%% measure with SM
fp = '..//data//apr08//StefanMayer//'
folder_name = '0cm'
folder_name = make_folder(fp, folder_name)
save = True
n_measurements = 5
# measure
for i in range(n_measurements):
    out = labjack_measure(measure_time, sr, ["AIN3"], [0.001], 
                                              [1.0], False)
    
    out_df = pd.DataFrame(out.T)
    out_df.columns = ['Epoch Time', 'SM']
    out_df.set_index("Epoch Time", inplace=True)
    
    time = datetime.now().strftime('%y%m%dT%H%M%S')
    # save
    if save:
        out_df.to_csv(fp+folder_name+'//'+time+".csv")
  






