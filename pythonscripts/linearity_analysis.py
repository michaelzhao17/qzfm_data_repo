# mar 2024
# linearity data analysis

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import glob
from scipy import signal
from scipy import optimize
import ntpath

def quad_fit(x, a, b, c):
    return a*x**2 + b*x + c

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def aLorentfit(x, A, freq, tau):
    return A*(7*x*tau)/(1+(7*x*tau)**2)

# curve fit to determine proportionality constant k between input voltage and input B

# fitting function
def fitfunc(V, a, b):
    return a * V + b

# def fitfunc(V, a):
#     return a * V 

# function to get convert factor from input voltage axis to input magnetic field 
def conv_factor_v2b(df, cutoff_idx):
    popt, pcov = optimize.curve_fit(fitfunc, df['Voltage p2p (V)'][:cutoff_idx],
                           df['Magnetic Field p2p (nT)'][:cutoff_idx], sigma=df['Magnetic Uncertainty (nT)'][:cutoff_idx],
                           absolute_sigma=True)
    return popt


#%% read data and find peak to peak with simple peak finding method
# sensor
sensor = 'AAL9'
# axis parallel to applied B
axis = 'x'
# driving frequency of applied B
freq = 35
# sampling frequency
sr = 10000

# bandpass
bandpass = False
 

results = {'Voltage p2p (V)':[],
           'Magnetic Field p2p (nT)':[],
           'Magnetic Uncertainty (nT)':[]}

    
for file in glob.glob('../data/Linearity/{}/{}axis_{}Hz/*.csv'.format(sensor, axis, freq)):
    df = pd.read_csv(file)
    # Vin
    filename = path_leaf(file)

    
    V = float(filename[:filename.index('-')]) 
    results['Voltage p2p (V)'].append(V)
    
    # time
    t = df['Epoch Time'] - df['Epoch Time'].iloc[0]
    # B along axis 
    B = df[axis].to_numpy()
  
    # find index of maximums and minimums
    max_idx = signal.find_peaks(B, distance=int(0.9*sr/freq))
    min_idx = signal.find_peaks(-B, distance=int(0.9*sr/freq))
    
    # ensure both arrays have same number of elements
    if len(max_idx) == len(min_idx):
        pass
    elif len(max_idx) == len(min_idx) + 1:
        max_idx = max_idx[:-1]
    elif len(max_idx) == len(min_idx) - 1:
        min_idx = min_idx[:-1]
    
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
    extremas = np.asarray([val for pair in zip(maximas, minimas) for val in pair])
  
    # plt.figure()
    # plt.plot(B)
    # plt.plot(max_idx[0], B[max_idx[0]], "x")
    # plt.plot(min_idx[0], B[min_idx[0]], "x")
    # plt.show()
    
    #time_.sleep(4)
    # calculate peak to peak 
    pp = []
    for i in range(len(extremas)-1):
        pp.append(abs(extremas[i]-extremas[i+1]))
    results['Magnetic Field p2p (nT)'].append(np.mean(pp))
    results['Magnetic Uncertainty (nT)'].append(np.std(pp))

#%% save data
results_df = pd.DataFrame(results)
input('Sure you want to save? If not careful, WILL overwrite existing file at {}'.format('../data/{}/processed/{}axis_035Hz.csv'.format(sensor, axis)))
results_df.to_csv('../data/Linearity/{}/processed/{}axis_035Hz.csv'.format(sensor, axis))

#%% v2
width = 0.7
height = 1.5

marker_dict = {'x':'*',
            'y':'^',
            'z':'v'}


fig = plt.figure(figsize=(width*6.4, height*4.8))

outer = gridspec.GridSpec(2, 1, height_ratios = [1, 4]) 
#make nested gridspecs
gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[1], hspace = .05)


# gs = fig.add_gridspec(8,1)
# # gs.update(wspace=1, hspace=1.5) # set the spacing between axes. 
ax = fig.add_subplot(gs2[0, 0])
ax1 = fig.add_subplot(gs2[1, 0])

ax2 = fig.add_subplot(gs1[:, 0])

cutoff_idx = {'x':10, 'y':10, 'z':10}

# inset 
x1, x2, y1, y2 = 2, 6, 2, 6  # subregion of the original image
axins = ax.inset_axes(
    [0.02, 0.57, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2))
axins.set_yticklabels([])
axins.set_xticklabels([])

x1, x2, y1, y2 = 0, 3, -2, 6  # subregion of the original image
axins1 = ax1.inset_axes(
    [0.07, 0.5, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2))

colorct = 0
for file in glob.glob('..//data//Linearity//AAL9//processed//*035Hz.csv'):
    filename = path_leaf(file)
    
    
    freq = filename[6:9]
    axis = filename[0]
    df = pd.read_csv(file)
    pctdev = abs(1 -  df['Magnetic Field p2p (nT)'] / 
                 (np.multiply(df['Voltage p2p (V)'], conv_factor_v2b(df, cutoff_idx[axis])[0])  + conv_factor_v2b(df, cutoff_idx[axis])[1]))
    
    
    
    
    
    
    ax1.plot(np.multiply(df['Voltage p2p (V)'], conv_factor_v2b(df, cutoff_idx[axis])[0]+
                   conv_factor_v2b(df, cutoff_idx[axis])[1])[::2], 
                 100*pctdev[::2], 
                 marker_dict[axis],
                 markersize=2,
                 label='{} Hz {} axis'.format(freq, axis))
    axins1.plot(np.multiply(df['Voltage p2p (V)'], conv_factor_v2b(df, cutoff_idx[axis])[0]+
                    conv_factor_v2b(df, cutoff_idx[axis])[1])[::2], 
                  100*pctdev[::2], 
                  marker_dict[axis],
                  markersize=2,
                  label='{} Hz {} axis'.format(freq, axis))
    
    popt, pcov = optimize.curve_fit(fitfunc, df['Voltage p2p (V)'][:cutoff_idx[axis]],
                           df['Magnetic Field p2p (nT)'][:cutoff_idx[axis]], sigma=df['Magnetic Uncertainty (nT)'][:cutoff_idx[axis]],
                           absolute_sigma=True)
    residuals = (fitfunc(df['Voltage p2p (V)'], *popt) - df['Magnetic Field p2p (nT)'])[:cutoff_idx[axis]]
    ax2.errorbar(x=df['Voltage p2p (V)'][:cutoff_idx[axis]], y=1000*residuals, 
                    yerr=1000*df['Magnetic Uncertainty (nT)'][:cutoff_idx[axis]], fmt=marker_dict[axis], 
                    capsize=3,)
    
    
    df.sort_values('Voltage p2p (V)', inplace=True)
    axis = file[-15]
    lin_idx = {'x':10, 'y':10, 'z':10}
    x = np.multiply(df['Voltage p2p (V)'], conv_factor_v2b(df, lin_idx[axis])[0]+
                   conv_factor_v2b(df, lin_idx[axis])[1])
    cutoff = np.argmin(abs(x-20))
    popt, pcov = optimize.curve_fit(aLorentfit, 
                                    xdata=x[:cutoff],
                           ydata=df['Magnetic Field p2p (nT)'].iloc[:cutoff],
                           sigma=df['Magnetic Uncertainty (nT)'].iloc[:cutoff],
                           absolute_sigma=True)
    ax.errorbar(x[:cutoff:4], 
                 y=df['Magnetic Field p2p (nT)'][:cutoff:4], 
                 yerr=df['Magnetic Uncertainty (nT)'][:cutoff:4],
                 capsize=1,
                 fmt=marker_dict[axis],
                 markersize=2,
                 label='{} axis'.format(axis),
                 color="C{}".format(colorct))
    ax.plot(x[:cutoff],
             aLorentfit(x[:cutoff], *popt), '-', color="C{}".format(colorct), alpha=0.5)
    
    axins.errorbar(x[:cutoff:4], 
                 y=df['Magnetic Field p2p (nT)'][:cutoff:4], 
                 yerr=df['Magnetic Uncertainty (nT)'][:cutoff:4],
                 capsize=1,
                 fmt=marker_dict[axis],
                 markersize=2,
                 label='{} axis'.format(axis),
                 color="C{}".format(colorct))
    axins.plot(x[:cutoff],
             aLorentfit(x[:cutoff], *popt), '-', color="C{}".format(colorct), alpha=0.5)
    
    colorct += 1
    
    
#ax.set_xlabel('Environment Peak-to-Peak Amplitude (nT)')
ax.set_ylabel('Response (nT)')
ax1.set_xlabel('Peak-to-Peak Amplitude (nT)')
ax1.set_ylabel('Percent Deviation')
ax2.set_xlabel('Driving Peak-to-Peak Voltage (V)')
ax2.set_ylabel('Residuals (pT)')
    
ax.text(16, 13.3, 'y', color='C1')
ax.text(16, 8.5, 'z', color='C2')
ax.text(16, 4, 'x', color='C0')


axins1.plot(np.arange(5), np.full((5, 1), 2), 'r--')
ax.indicate_inset_zoom(axins, edgecolor="black", linewidth=1)
ax1.indicate_inset_zoom(axins1, edgecolor="black", linewidth=1)
ax.plot(np.linspace(0, 20, 100), np.linspace(0, 20, 100),
         'k--', label='Ideal Linearity')
axins.plot(np.linspace(0, 20, 100), np.linspace(0, 20, 100),
         'k--', label='Ideal Linearity')

ax2.axhline(0, -1, 11, linewidth=2, color='Black')
#ax.legend()
ax.set_xlim(left=-1, right=20)
ax.set_ylim(top=15)
ax.set_xticklabels([])
ax1.set_xlim(left=-1, right=20)
ax.grid() 
ax1.grid()
ax2.grid()
axins.grid()
axins1.grid()
plt.tight_layout()  
plt.show()



