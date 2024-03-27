import pickle
import numpy as np
import time
from linien_client.connection import LinienClient

from matplotlib import pyplot as plt
from time import sleep

def remove_background(w): # remove linear background
    x = np.arange(0,len(w))
    poly = np.polyfit(x,w,1)
    poly_y = np.poly1d(poly)(x)
    return w - poly_y

def get_waveform(c: LinienClient):
    #returns reflection signal while sweeping, length is 2048
    wf = []
    for _ in range(5):
        time.sleep(0.05)
        plot_data = pickle.loads(c.parameters.to_plot.value)
        waveform = np.array(plot_data['error_signal_1'])
        waveform = remove_background(waveform)
        waveform = -1 * waveform
        wf.append(waveform)
    waveform = np.average(wf, axis=0)
    return waveform

def set_scan_range(c: LinienClient, sweep_center, sweep_amplitude) -> None:
    c.parameters.sweep_center.value = sweep_center
    c.parameters.sweep_amplitude.value = sweep_amplitude
    c.connection.root.write_registers() # apply changes
