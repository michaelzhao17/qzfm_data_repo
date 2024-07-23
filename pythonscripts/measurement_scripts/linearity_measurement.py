# mar 2024
# to fully automate linearity measurements

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
from tqdm import tqdm
from scipy import signal
from scipy.signal import butter, lfilter
import pyvisa 
from SiglentDevices import DG1032Z
import pathlib
from DMM6500 import DMM6500
from DMM6500_SCPI import Function
from labjackmeasure import labjack_measure
#%% initialize instruments as objects
rm = pyvisa.ResourceManager()
#print(rm.list_resources())

# Digital Multimeter
DMM = rm.open_resource('USB0::0x05E6::0x6500::04550534::INSTR')
DMM.read_termination = '\n'
DMM.write_termination = '\n'
DMM.query('*IDN?')
mm = DMM6500(DMM)
mm.function = Function.AC_VOLTAGE
mm.timeout = 25000


# Digital Function Generator
awg = DG1032Z(hostname='USB0::0x1AB1::0x0642::DG1ZA232603182::INSTR')
awg.query('*IDN?')

#%%
# QZFM
q = QZFM('COM3')

#%%

q.set_master(True)
#%%



q.auto_start(zero_calibrate=False)

#%% auto folder creation function
def make_folder(fp, axis, freq):
    '''
    Parameters
    ----------
    fp : str
        file path to location where folder is to be created
    axis : str, x|y|z
        axis being measured.
    freq : float
        frequency being measured.

    Returns
    -------
    folder_name : str
        name of folder
    
    creates folder and then returns folder name
    '''
    folder_name = '{}axis_{}Hz'.format(axis, freq)
    # check if folder already exists, if not create it
    pathlib.Path(fp+folder_name).mkdir(parents=True, exist_ok=True) 
    return folder_name

#%% initial configuration
# AWG settings
ch = 1
freq = 35
offset = 0
vpp = 1
waveform = 'SIN'
hp = True

# circuit setting
r = 3390 # resistance value of r2 (the one NOT apart of highpass)

# QZFM settings
gain = '0.1x' # possible gains are 0.1x|0.33x|1x|3x
# corresponding conversion V/nT for each gain
gain_dict = {'0.1x':0.27,
             '0.33x':0.9,
             '1x':2.7,
             '3x':8.1}
q.set_gain(gain)
zero_t = 8
t = 5 # number of seconds to record
sr = 10000 # sr
axis = 'x' # axis being measured

# save file settings
save = True # save as csv if True
fp = '..//data//AAL9//may21//'

# dictionary of axis and corresponding labjack channel
ljch = {'x':'AIN1',
        'y':'AIN0',
        'z':'AIN2'}
#%% main script
if __name__ == '__main__':
    # iterate over frequencies
    for freq in [35]:
        if hp:
            r1 = 5062
            r2 = r
            res = 1/ (1/r1 + 1/r2)
            cap = 8.67e-6
            imp = np.sqrt(res**2+1/(2*np.pi*freq*cap)**2)
        if not hp:
            imp = r
        awg.set_ch_state(ch, False)
        awg.set_impedance(ch, imp)
        awg.set_wave(ch, waveform, freq, vpp, offset, phase=0)
        # iterable of Vpp values to output
        vpps = np.linspace(0.05, 17, 200, endpoint=False)
        i = 39
        # make folder and get folder name 
        folder_name = make_folder(fp, axis, freq)
        # zero 
        q.field_zero(True, show=False)
        # sleep 10 seconds
        for i in range(zero_t):
            time_.sleep(1)
        q.field_zero(False)
        
        
        for idx, vpp in enumerate(vpps):
            awg.set_wave(ch, waveform, freq, vpp, offset, phase=0)
            
            
            # turn on AWG
            awg.set_ch_state(ch, state=True)
            time_.sleep(2)
            # measure rms voltage from DMM            
            v_rms_list = []
            instr = 'DMM'
            j = 0
            while j < 10:
                try:
                    ret = mm.measure()
                    v_rms_list.append(ret)
                    j += 1
                except Exception:
                    pass
            # calculate average vrms measured and convert to vpp
            v_rms_meas = np.mean(v_rms_list)
            v_pp_meas = round(v_rms_meas * 2 * np.sqrt(2), 2)
            strV = str(v_pp_meas)
            print('DMM measures {} Vpp'.format(strV))
            if len(strV) == 4:
                pass
            else:
                strV = strV + '0'
          

            # labjack measure
            out = labjack_measure(t, sr, [ljch[axis]], [gain_dict[gain]], 
                                                      [10.0], False)

            
            # turn off AWG
            awg.set_ch_state(ch, state=False)
            # save data
            out_df = pd.DataFrame(out.T)
            out_df.columns = ['Epoch Time', axis]
            out_df.set_index("Epoch Time", inplace=True)
            if save:
                current_time = datetime.now().strftime('%y%m%dT%H%M%S')
                out_df.to_csv(fp+folder_name+'//'+strV+'-'+datetime.now().strftime('%y%m%dT%H%M%S')+'.csv')
            i += 1
    
    


#%% testing
vpp = 5
awg.set_wave(ch, waveform, freq, vpp, offset, phase=0)

# zero 
q.field_zero(True, show=False)
# sleep 10 seconds
for i in range(zero_t):
    time_.sleep(1)
q.field_zero(False)

# turn on AWG
awg.set_ch_state(ch, state=True)
time_.sleep(2)
# measure rms voltage from DMM            
v_rms_list = []
instr = 'DMM'
j = 0
while j < 10:
    try:
        ret = mm.measure()
        v_rms_list.append(ret)
        j += 1
    except Exception:
        pass
# calculate average vrms measured and convert to vpp
v_rms_meas = np.mean(v_rms_list)
v_pp_meas = round(v_rms_meas * 2 * np.sqrt(2), 2)
strV = str(v_pp_meas)
print('DMM measures {} Vpp'.format(strV))
if len(strV) == 4:
    pass
else:
    strV = strV + '0'


# labjack measure
out = labjack_measure(t, sr, ["AIN0", "AIN1", "AIN2"], [gain_dict[gain], gain_dict[gain], gain_dict[gain]], 
                                          [10.0, 10.0, 10.0], False)


# turn off AWG
awg.set_ch_state(ch, state=False)

#%%

plt.figure()
# for i in range(1):
#     plt.plot(out[i+1, :])
plt.plot(out[2])
plt.show()






