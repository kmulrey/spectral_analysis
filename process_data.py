import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import sys
import glob
import process_files
import re
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.optimize import curve_fit

def filter(times,traces, fmin, fmax): # 2d-> time,data
    if traces.ndim == 2:
        traces = np.expand_dims(traces, axis=0)

    filetered_traces = np.full_like(traces, 0)

    nTraces=len(traces)
    tstep=times[0][1]-times[0][0]
    ldata=len(traces[1][0])

    rate=1/tstep
    nyq=0.5*rate
    lowcut=fmin*1e6
    highcut=fmax*1e6
    low = lowcut/nyq
    high = highcut/nyq
    order=1
    b, a = signal.butter(order, [low, high], btype='band')

    for i in np.arange(nTraces):
        filetered_traces[i].T[0]= signal.filtfilt(b, a, traces[i].T[0])  # this is data in the time domain
        filetered_traces[i].T[1]= signal.filtfilt(b, a, traces[i].T[1])  # this is data in the time domain
        filetered_traces[i].T[2]= signal.filtfilt(b, a, traces[i].T[2])  # this is data in the time domain

    return np.squeeze(filetered_traces)
