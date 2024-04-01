import pickle
import numpy as np

from linien_client.connection import LinienClient
import time
from matplotlib import pyplot as plt
from time import sleep
from scipy.signal import savgol_filter as sf
def remove_background(w): # remove filtered background
    x = np.arange(0,len(w))
    sf_fit = sf(w, 150, 3)
    poly = np.polyfit(x,w,1)
    poly_y = np.poly1d(poly)(x)
    return w - poly_y

def get_waveform(c: LinienClient): #averaged over 10
    #returns reflection signal while sweeping, length is 2048
    waveform = np.zeros(2048)
    for _ in range(1): # couldn't average because piezo voltage is shifting?
        time.sleep(0.1)
        plot_data = pickle.loads(c.parameters.to_plot.value)
        wf = np.array(plot_data['error_signal_1'])
        wf = remove_background(wf)
        wf = -1 * wf
        waveform = waveform + wf
    return waveform

def set_scan_range(c: LinienClient, sweep_center, sweep_amplitude) -> None:
    c.parameters.sweep_center.value = sweep_center
    c.parameters.sweep_amplitude.value = sweep_amplitude
    c.connection.root.write_registers() # apply changes
