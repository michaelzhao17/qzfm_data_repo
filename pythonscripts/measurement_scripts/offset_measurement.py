# Jan 2024
# To measure from Labjack and QZFM module serial roughly simultaneously 

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
import pandas as pd
from labjackmeasure import labjack_measure
import pathlib

q = QZFM("COM4")

# make folder function
def make_folder(fp, sensor, axis):
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
    folder_name = '{}_around_{}axis'.format(sensor, axis)
    # check if folder already exists, if not create it
    pathlib.Path(fp+folder_name).mkdir(parents=True, exist_ok=True) 
    return folder_name
#%%
q.auto_start(zero_calibrate=False)


#%%
q.set_master(True)


#%%
if __name__ == '__main__':
    t = 40
    orientation = 1
    fp = '..//data//jun06//'
    sensor = 'AAY4'
    axis = 'x' # of rotaton
    folder_name = make_folder(fp, sensor, axis)
    save = True
    flip_num = 10
    
    for i in range(flip_num): 
        print('currently on measurement no.{} out of {}'.format(i+1, flip_num))
        q.field_zero(True, True, False)
        # initialize queue object for saving outputs of functions
        output_queue = queue.Queue()
         
        # turn on field zeroing
        # t0 = Thread(target=q.field_zero, args=(True, True, False))
        # t0.start()
        # t0.join()
        # print('zeroing started')
        
        # initialize thread of QZFM module for coil offset reading
        # q = QZFM('COM3')
        t1 = Thread(target=q.read_offsets_custom, args=(int(t*7.5), output_queue))
        print('QZFM thread defined')
        
        # initialize thread of labjack for cell reading
        t2 = Thread(target=labjack_measure, args=(t, 100, ["AIN0", "AIN1", "AIN2"], [0.0027, 0.0027, 0.0027], 
                                                  [1.0, 1.0, 1.0], True, output_queue))
        print("Labjack thread defined")
                                                                   
        # start the threads
        t1.start()
        print('QZFM thread started')
        t2.start()
        print('Labjack thread started')
        
     
        # join the threads
        t1.join()
        print("p1 joined")
        t2.join()
        print("p2 joined")
        
        q.field_zero(False)
        q.field_reset()
        
        cnt = 0
        out = []
        while cnt < 2:
            try:
                out.append(output_queue.get())
                cnt +=1
            except Exception:
                break
        
        for item in out:
            if isinstance(item, pd.DataFrame):
                coil_offsets = item
            elif isinstance(item, np.ndarray):
                cell_readings = item 
    
        cell_readings = pd.DataFrame(cell_readings.T)
        cell_readings.columns = ['Epoch Time', 'x', 'y', 'z']
        cell_readings.set_index("Epoch Time", inplace=True)
        
        time = datetime.now().strftime('%y%m%dT%H%M%S')
        
        if save:
            if orientation == 1:
                coil_offsets.to_csv(fp+folder_name+'//'+time+"COIL"+"up"+".csv")
                cell_readings.to_csv(fp+folder_name+'//'+time+"CELL"+"up"+".csv")
                orientation *= -1
            elif orientation == -1:
                coil_offsets.to_csv(fp+folder_name+'//'+time+"COIL"+"dn"+".csv")
                cell_readings.to_csv(fp+folder_name+'//'+time+"CELL"+"dn"+".csv")
                orientation *= -1
                
        q.field_zero(False)
        input('press enter to continue\n')
    # turn off field zeroing
    q.field_zero(False)
    print('zeroing turned off')
#%%

plt.figure()
for axis in ['x', 'y', 'z']: 
    
    plt.plot(cell_readings.index-cell_readings.index[0], cell_readings[axis], label="Cell Reading {}".format(axis))
    plt.plot((coil_offsets.index-cell_readings.index[0])[:], coil_offsets[axis].iloc[:], label="Coil Reading {}".format(axis))
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("B field [pT]")
    plt.title("Simultaneous coil offset field and cell field in {} direction".format(axis))

plt.grid()
plt.show()
