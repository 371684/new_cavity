import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks

FSR = 970e6  # FSR of the cavity is measured to be 970MHz

def ini_peak_search(ini_scope_data: npt.NDArray[np.float_], min_height:float = 0.05):
    # find the peaks with minimium required height and distance between two peaks
    # to avoid get local maximum peaks
    ini_peaks, _ = find_peaks(ini_scope_data,min_height,distance=10)

    # # piezo scan voltage set with 4.5V amplitude, at oscilloscope usually have 4 TEM00 peaks. 
    # while True:
    #     if len(ini_peaks) > 4:
    #         min_height += 0.01
    #         ini_peaks, _ = find_peaks(ini_scope_data,min_height,distance=20)
    #         # to make sure two peaks we will use are the peaks in the middle, total peaks
    #         # number should > 4
    #     elif len(ini_peaks) < 4:
    #         print("didn't detect enough peaks, check piezo voltage")
    #         break
    #     else:
    #         break
    
    return ini_peaks

def reso_det(ini_peaks: npt.NDArray[np.float_]):
    # determine the frequency resolution of the data (the frequency corrsponding to distance between two poins)
    peaks_dis = np.diff(ini_peaks)
    FSR_points = peaks_dis[0]
    reso = FSR / FSR_points
    return reso

def mod_peak_search(mod_scope_data: npt.NDArray[np.float_],min_height:float = 0.05):
    # search for the peaks whrn turn on the EOM modulaiton signal
    mod_peaks,_ =  find_peaks(mod_scope_data,min_height,distance=20)
    return mod_peaks

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def zero_measure(
    mod_scope_data: npt.NDArray[np.float_],
    ini_peaks: npt.NDArray[np.float_],
    mod_peaks: npt.NDArray[np.float_],
    reso: float
    ):
    # re-measure the first TEM00 peak's position incase it shifted between two measurements
    zero_posi_est1 =  ini_peaks[0]  # roughly determine the position of the 0th order peak from the first TEM00 peak
    zero_posi_est2 =  ini_peaks[1]  # roughly determine the second position of the 0th order peak from the second TEM00 peak
    zero_posi_mea1 = find_nearest(mod_peaks, zero_posi_est1) # find the nearest peak position 
    zero_posi_mea2 = find_nearest(mod_peaks, zero_posi_est2) # find the nearest peak position 
    if np.abs(zero_posi_mea1-zero_posi_est1)*reso < 20e6:
        return zero_posi_mea1,zero_posi_mea2, mod_scope_data[zero_posi_mea1]
    else:
        #('zeroth peak is too small or has moved')
        zero_amp = 0
        return zero_posi_est1,zero_posi_est2, zero_amp

def fir_measure(
    mod_scope_data: npt.NDArray[np.float_],
    ini_peaks: npt.NDArray[np.float_],
    mod_peaks: npt.NDArray[np.float_], 
    mod_fre: float, 
    reso: float
    ):

    fir_posi_est = int((mod_fre%1e9) / reso + ini_peaks[0])  # roughly determine the position of the 1st order peak feom the second TEM00 peak
    fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    
    if np.abs(fir_posi_mea-fir_posi_est)*reso > 40e6:
        fir_amp = 0.05
    # failed to find first order peak
    else:
        fir_amp = mod_scope_data[fir_posi_mea]
    
    return fir_posi_mea, fir_amp

