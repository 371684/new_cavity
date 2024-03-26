import pickle
import numpy as np

from linien_client.connection import LinienClient
import time
from matplotlib import pyplot as plt
from time import sleep

def remove_background(w): # remove linear background
    x = np.arange(0,len(w))
    poly = np.polyfit(x,w,1)
    poly_y = np.poly1d(poly)(x)
    return w - poly_y

def get_waveform(c: LinienClient): #averaged over 10
    #returns reflection signal while sweeping, length is 2048
    waveform = np.zeros(2048)
    for _ in range(10):
        time.sleep(1.0/20)
        plot_data = pickle.loads(c.parameters.to_plot.value)
        wf = np.array(plot_data['error_signal_1'])
        wf = remove_background(wf)
        wf = -1 * wf
        waveform = waveform + wf/10
    return waveform

def set_scan_range(c: LinienClient, sweep_center, sweep_amplitude) -> None:
    c.parameters.sweep_center.value = sweep_center
    c.parameters.sweep_amplitude.value = sweep_amplitude
    c.connection.root.write_registers() # apply changes
