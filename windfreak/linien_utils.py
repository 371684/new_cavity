import pickle
import numpy as np

from linien_client.connection import LinienClient

from matplotlib import pyplot as plt
from time import sleep

def get_waveform(c: LinienClient):
    #returns reflection signal while sweeping, length is 2048
    plot_data = pickle.loads(c.parameters.to_plot.value)
    waveform = np.array(plot_data['error_signal_1'])
    waveform = -1 * waveform + np.average(waveform)
    return waveform

def set_scan_range(c: LinienClient, sweep_center, sweep_amplitude) -> None:
    c.parameters.sweep_center.value = sweep_center
    c.parameters.sweep_amplitude.value = sweep_amplitude
    c.connection.root.write_registers() # apply changes
