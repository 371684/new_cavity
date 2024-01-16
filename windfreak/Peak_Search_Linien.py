import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks
from linien_utils import LinienClient, set_scan_range, get_waveform
import time

FSR = 970e6  # FSR of the cavity is measured to be 970MHz
MAX_SCAN_RANGE = 2 # 2V scan range

def ini_peak_search(c: LinienClient, threshold:float = 10):
    # find the peaks with minimium required *prominence* and distance between two peaks
    # to avoid get local maximum peaks or noise 
    set_scan_range(c, 0, 1)
    time.sleep(0.5)
    ini_rev_data = get_waveform(c)

    data = np.sort(ini_rev_data)
    height = data[-10]
    ini_peaks, _ = find_peaks(ini_rev_data,height=height,threshold=threshold,prominence=100,distance=500)

    # determine the frequency resolution of the data (the frequency corrsponding to distance between two poins)
    # MHz/point at max scanning range
    if len(ini_peaks) == 2:
        FSR_points = ini_peaks[1] - ini_peaks[0]
        reso = FSR / FSR_points
        # print(2*ini_peaks[0]/len(ini_rev_data))
        # print(2*ini_peaks[1]/len(ini_rev_data))
        # print(2*FSR_points/len(ini_rev_data))
        # raise KeyboardInterrupt
    else:
        plt.plot(ini_rev_data)
        raise RuntimeError("not enough/too many peaks detected")
    return ini_peaks, reso

def mod_peak_search(c: LinienClient):
    set_scan_range(c, 0, 1)
    # search for the peaks whrn turn on the EOM modulaiton signal
    mod_scope_data = get_waveform(c)
    if max(mod_scope_data) > 200:
        mod_peaks,_ =  find_peaks(mod_scope_data, threshold=15,prominence=50, distance=20)
    elif max(mod_scope_data) > 100:
        mod_peaks, _ =  find_peaks(mod_scope_data, threshold=5,prominence=10, distance=20)
    else:
        mod_peaks,_ =  find_peaks(mod_scope_data, threshold=1,prominence=10, distance=20)
    # failed to find peak, likely due to redpitaya not yet got the data
    if len(mod_peaks) == 0:
        time.sleep(0.1)
        mod_scope_data, mod_peaks = mod_peak_search(c)
    return mod_scope_data, mod_peaks

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def zero_measure(
    mod_scope_data,
    ini_peaks,
    mod_peaks,
    reso: float
    ):
    # re-measure the first TEM00 peak's position in case it shifted between two measurements
    zero_posi_est1 =  ini_peaks[0]  # roughly determine the position of the 0th order peak from the first TEM00 peak
    zero_posi_est2 =  ini_peaks[1]  # roughly determine the second position of the 0th order peak from the second TEM00 peak
    zero_posi_mea1 = find_nearest(mod_peaks, zero_posi_est1) # find the nearest peak position 
    zero_posi_mea2 = find_nearest(mod_peaks, zero_posi_est2) # find the nearest peak position 
    # print(mod_peaks)
    if np.abs(zero_posi_mea1 - zero_posi_est1) * reso < 20e6:
        return zero_posi_mea1, zero_posi_mea2, mod_scope_data[zero_posi_mea1]
    else:
        #zeroth peak is too small 
        zero_amp = 0
        return zero_posi_est1, zero_posi_est2, zero_amp

def fir_measure(
    c: LinienClient,
    mod_scope_data,
    ini_peaks,
    mod_peaks, 
    mod_fre: float, 
    reso: float,
    ):

    fir_posi_est = int((mod_fre%1e9) / reso + ini_peaks[0])  # roughly determine the position of the 1st order peak feom the second TEM00 peak
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    # re-measure the height of this peak by decreasing the piezo scanning range
    # print(abs(fir_posi_est-fir_posi_mea))
    if abs(fir_posi_est-fir_posi_mea) > 30:
        plt.plot(mod_scope_data)
        plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
        plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
        plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "8", color='orange')
        raise RuntimeError("peak is too far from expected")
        # fir_posi_mea = fir_posi_est
        # plt.plot(mod_scope_data)
        # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
        # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
        # plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "o")
        # raise KeyboardInterrupt
    # plt.plot(mod_scope_data)
    # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
    # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
    # plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "8", color='orange')
    # raise KeyboardInterrupt
    sweep_amplitude = 0.5
    sweep_center = 0
    max_scan_range = MAX_SCAN_RANGE
    while max_scan_range > 0.1:
        sweep_center = sweep_center + (-1/2  + fir_posi_mea/len(mod_scope_data)) * max_scan_range # find the voltage corresponds to 1st order peak
        set_scan_range(c, sweep_center, sweep_amplitude)    
        mod_scope_data_new, mod_peaks_new = mod_peak_search(c)
        fir_posi_est = len(mod_scope_data_new)/2 
        fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) 
        MAX_SCAN_RANGE = sweep_amplitude
        sweep_amplitude = sweep_amplitude / 2
    fir_posi_est = len(mod_scope_data_new) / 2
    fir_posi_mea = find_nearest(mod_peaks_new, fir_posi_est)
    fir_amp = mod_scope_data_new[fir_posi_mea]
    
    # plt.plot(mod_scope_data_new)
    # plt.plot(mod_peaks, mod_scope_data_new[mod_peaks], "x")
    # plt.plot(fir_posi_mea, mod_scope_data_new[fir_posi_mea], "o", color='red')
    # raise KeyboardInterrupt
    return fir_posi_mea, fir_amp

