import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks
import linien_utils
import time

FSR = 970e6  # FSR of the cavity is measured to be 970MHz
MAX_SCAN_RANGE = 2 # 2V scan range

def ini_peak_search(ini_scope_data: np.array(), prominence:float = 0.05):
    # find the peaks with minimium required *prominence* and distance between two peaks
    # to avoid get local maximum peaks or noise 
    ini_peaks, _ = find_peaks(ini_scope_data,prominence=prominence,distance=500)
    
    return ini_peaks

def reso_det(ini_peaks: np.array()):
    # determine the frequency resolution of the data (the frequency corrsponding to distance between two poins)
    # MHz/point at max scanning range
    peaks_dis = np.diff(ini_peaks)
    FSR_points = peaks_dis[0]
    reso = FSR / FSR_points
    return reso

def mod_peak_search(mod_scope_data: np.array(), prominence:float = 0.05):
    # search for the peaks whrn turn on the EOM modulaiton signal
    mod_peaks,_ =  find_peaks(mod_scope_data, prominence=prominence, distance=20)
    return mod_peaks

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def zero_measure(
    mod_scope_data: np.array(),
    ini_peaks: np.array(),
    mod_peaks: np.array(),
    reso: float
    ):
    # re-measure the first TEM00 peak's position in case it shifted between two measurements
    zero_posi_est1 =  ini_peaks[0]  # roughly determine the position of the 0th order peak from the first TEM00 peak
    zero_posi_est2 =  ini_peaks[1]  # roughly determine the second position of the 0th order peak from the second TEM00 peak
    zero_posi_mea1 = find_nearest(mod_peaks, zero_posi_est1) # find the nearest peak position 
    zero_posi_mea2 = find_nearest(mod_peaks, zero_posi_est2) # find the nearest peak position 
    if np.abs(zero_posi_mea1 - zero_posi_est1) * reso < 20e6:
        return zero_posi_mea1, zero_posi_mea2, mod_scope_data[zero_posi_mea1]
    else:
        #zeroth peak is too small 
        zero_amp = 0
        return zero_posi_est1, zero_posi_est2, zero_amp

def fir_measure(
    c: linien_utils.LinienClient,
    mod_scope_data: np.array(),
    ini_peaks: np.array(),
    mod_peaks: np.array(), 
    mod_fre: float, 
    reso: float
    ):

    fir_posi_est = int((mod_fre%1e9) / reso + ini_peaks[0])  # roughly determine the position of the 1st order peak feom the second TEM00 peak
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    
    # re-measure the height of this peak by decreasing the piezo scanning range
    sweep_center = (-1/2  + fir_posi_mea/len(mod_scope_data)) * MAX_SCAN_RANGE # find the voltage corresponds to 1st order peak
    sweep_amplitude = 20e6/FSR * reso/len(mod_scope_data) * MAX_SCAN_RANGE # scan range 1st order peak +/- 20MHz
    linien_utils.set_scan_range(c, sweep_center, sweep_amplitude)
    time.sleep(0.1)
    mod_scope_data = linien_utils.get_waveform(c)
    mod_peaks = mod_peak_search(mod_scope_data)
    fir_posi_est = len(mod_scope_data)/2 # expected at center
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 

    if np.abs(fir_posi_mea-fir_posi_est)*reso > 40e6:
        # failed to find first order peak
        fir_amp = 0    
    else:
        fir_amp = mod_scope_data[fir_posi_mea]
    
    return fir_posi_mea, fir_amp

