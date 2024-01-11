import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks
import linien_utils
import time

FSR = 970e6  # FSR of the cavity is measured to be 970MHz
MAX_SCAN_RANGE = 2 # 2V scan range

def ini_peak_search(c: linien_utils.LinienClient, threshold:float = 10):
    # find the peaks with minimium required *prominence* and distance between two peaks
    # to avoid get local maximum peaks or noise 
    linien_utils.set_scan_range(c, 0, 1)
    time.sleep(0.3)
    ini_rev_data = linien_utils.get_waveform(c)

    data = np.sort(ini_rev_data)
    height = data[-20]
    ini_peaks, _ = find_peaks(ini_rev_data,height=height,threshold=threshold,prominence=40,distance=500)

    # determine the frequency resolution of the data (the frequency corrsponding to distance between two poins)
    # MHz/point at max scanning range
    if len(ini_peaks) == 2:
        FSR_points = ini_peaks[1] - ini_peaks[0]
        reso = FSR / FSR_points
    else:
        plt.plot(ini_rev_data)
        raise RuntimeError("not enough/too many peaks detected")
    return ini_peaks, reso

def mod_peak_search(mod_scope_data, threshold:float = 10):
    # search for the peaks whrn turn on the EOM modulaiton signal
    mod_peaks,_ =  find_peaks(mod_scope_data, threshold=threshold, prominence=40, distance=20)
    return mod_peaks

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
    c: linien_utils.LinienClient,
    mod_scope_data,
    ini_peaks,
    mod_peaks, 
    mod_fre: float, 
    reso: float
    ):

    fir_posi_est = int((mod_fre%1e9) / reso + ini_peaks[0])  # roughly determine the position of the 1st order peak feom the second TEM00 peak
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    # re-measure the height of this peak by decreasing the piezo scanning range
    sweep_center = (-1/2  + fir_posi_mea/len(mod_scope_data)) * MAX_SCAN_RANGE # find the voltage corresponds to 1st order peak
    sweep_amplitude = 0.05
    # sweep_amplitude = 100e6/reso/len(mod_scope_data) * MAX_SCAN_RANGE # scan range 1st order peak +/- 20MHz
    # fix bug in linien scanning range
    sweep_center = sweep_center - (-0.20301619*sweep_amplitude + 0.18310253)

    linien_utils.set_scan_range(c, sweep_center, sweep_amplitude)
    time.sleep(0.3)
    mod_scope_data = linien_utils.get_waveform(c)
    mod_peaks = mod_peak_search(mod_scope_data)
    fir_posi_est = len(mod_scope_data)/2 # expected at center
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 

    linien_utils.set_scan_range(c, 0, 1) # restore

    if np.abs(fir_posi_mea-fir_posi_est)*reso > 40e6:
        # failed to find first order peak
        fir_amp = 0    
    else:
        fir_amp = mod_scope_data[fir_posi_mea]
    
    return fir_posi_mea, fir_amp

