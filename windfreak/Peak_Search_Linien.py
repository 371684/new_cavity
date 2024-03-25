import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks
from linien_utils import LinienClient, set_scan_range, get_waveform
import time

def find_nearest(array, value):
    try:
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    except Exception:
        return 0
FSR = 970e6  # FSR of the cavity is measured to be 970MHz
# MAX_SCAN_RANGE = 2 # 2V scan range
class LinienDevice():
    def __init__(self, host="10.10.222.36", user="root", password="root") -> None:
        super().__init__()
        self.client = LinienClient(
            host=host,
            user=user,
            password=password
        )
        self.client.connect(autostart_server=True, use_parameter_cache=False)
        self.ini_peaks = []
        self.mod_scope_data = []
        self.mod_peaks = []
        self.zero_posi_mea1 = 0
        self.zero_posi_mea2 = 0
        self.mod_zero_amp = 0
        self.sweep_center = 0
        self.fir_amp = 0       
        self.reso = 0


    def ini_peak_search(self, threshold:float = 10):
        # find the peaks with minimium required *prominence* and distance between two peaks
        # to avoid get local maximum peaks or noise 
        set_scan_range(self.client, 0, 0.8)
        # determine the frequency resolution of the data (the frequency corrsponding to distance between two poins)
        # MHz/point at max scanning range
        ini_peaks = []
        while len(ini_peaks) != 2:
            time.sleep(0.1)
            ini_rev_data = get_waveform(self.client)
            data = np.sort(ini_rev_data)
            height = data[-10]
            ini_peaks, _ = find_peaks(ini_rev_data, height=height, threshold=threshold, prominence=100, distance=500)                           
        FSR_points = ini_peaks[1] - ini_peaks[0]
        reso = FSR / FSR_points
        self.ini_peaks = ini_peaks
        self.reso = reso
    
    def mod_peak_search(self, first_search = False):
        # search for the peaks when turn on the EOM modulaiton signal
        start = time.time()
        mod_peaks = [] 
        mod_scope_data = np.array([0])
        while len(mod_peaks) == 0:  # failed to find peak, likely due to redpitaya not yet got the data
            time.sleep(0.1)
            if first_search == False: # search for peaks without modulation, so need to find peaks large enough
                mod_scope_data = get_waveform(self.client)
                mod_peaks, _ =  find_peaks(mod_scope_data, threshold=3,prominence=10, distance=10)
                finish = time.time()
                if finish - start > 5:
                    mod_peaks = []
                    break
            else:
                mod_scope_data = get_waveform(self.client)
                mod_peaks, _ =  find_peaks(mod_scope_data, threshold=3,prominence=10, distance=10)
                finish = time.time()
                if finish - start > 5:
                    mod_peaks = []
                    break
        self.mod_scope_data = mod_scope_data
        self.mod_peaks = mod_peaks

    def zero_measure(self):
        # re-measure the first TEM00 peak's position in case it shifted between two measurements
        zero_posi_est1 =  self.ini_peaks[0]  # roughly determine the position of the 0th order peak from the first TEM00 peak
        zero_posi_est2 =  self.ini_peaks[1]  # roughly determine the second position of the 0th order peak from the second TEM00 peak
        zero_posi_mea1 = find_nearest(self.mod_peaks, zero_posi_est1) # find the nearest peak position 
        zero_posi_mea2 = find_nearest(self.mod_peaks, zero_posi_est2) # find the nearest peak position 
        if np.abs(zero_posi_mea1 - zero_posi_est1) * self.reso < 20e6 and np.abs(zero_posi_mea2 - zero_posi_est2) * self.reso < 20e6:
            self.zero_posi_mea1 = zero_posi_mea1
            self.zero_posi_mea2 = zero_posi_mea2
            self.mod_zero_amp = (self.mod_scope_data[zero_posi_mea1] + self.mod_scope_data[zero_posi_mea2])/2
        else:
            #zeroth peak is too small 
            self.mod_zero_amp = 0
            self.zero_posi_mea1 = zero_posi_est1
            self.zero_posi_mea2 = zero_posi_est2  

    def fir_measure(self, mod_fre: float):
        if self.sweep_center == 0:
            # roughly determine the position of the 1st order peak feom the second TEM00 peak
            # using cavity scanning the peak shifts to lower voltage
            fir_posi_est = int(-(mod_fre%1e9) / self.reso + self.ini_peaks[1])  
            
            fir_posi_mea = find_nearest(self.mod_peaks, fir_posi_est) # find the nearest peak position 
            # re-measure the height of this peak by decreasing the piezo scanning range
            sweep_amplitude = 0.8
            max_scan_range = 0.8*2
            for sweep_amplitude in [0.6, 0.4, 0.32, 0.25, 0.18, 0.12, 0.12]:
                self.sweep_center = self.sweep_center + (-1/2  + fir_posi_mea/len(self.mod_scope_data)) * max_scan_range # find the voltage corresponds to 1st order peak
                max_scan_range = sweep_amplitude*2
                set_scan_range(self.client, self.sweep_center, sweep_amplitude)   
                time.sleep(0.5)
                self.mod_peak_search()
                if len(self.mod_peaks) == 0: # peak too small
                    self.fir_amp = 0
                if sweep_amplitude > 0.18:
                    fir_posi_est = len(self.mod_scope_data)/2 
                    fir_posi_mea = find_nearest(self.mod_peaks, fir_posi_est) 
                else:
                    fir_posi_mea = np.argmax(self.mod_scope_data)
            self.fir_amp = self.mod_scope_data[fir_posi_mea]
        else:
            self.mod_peak_search()
            fir_posi_est = len(self.mod_scope_data)/2 
            fir_posi_mea = find_nearest(self.mod_peaks, fir_posi_est) 
            self.fir_amp = self.mod_scope_data[fir_posi_mea]

    # def mod_peak_search(c: LinienClient, first_search = False):
    #     # search for the peaks whrn turn on the EOM modulaiton signal
    #     start = time.time()
    #     mod_peaks = [] 
    #     mod_scope_data = np.array([0])
    #     while len(mod_peaks) == 0:  # failed to find peak, likely due to redpitaya not yet got the data
    #         time.sleep(0.1)
    #         if first_search == False: # search for peaks without modulation, so need to find peaks large enough
    #             mod_scope_data = get_waveform(c)
    #             if max(mod_scope_data) > 400:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=30,prominence=300, distance=10)
    #             elif max(mod_scope_data) > 150:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=15,prominence=100, distance=10)
    #             else:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=3,prominence=10, distance=10)
    #             finish = time.time()
    #             if finish - start > 5:
    #                 mod_peaks = []
    #                 break
    #         else:
    #             mod_scope_data = get_waveform(c)
    #             if max(mod_scope_data) > 200:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=20,prominence=50, distance=10)
    #             elif max(mod_scope_data) > 100:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=10,prominence=25, distance=10)
    #             else:
    #                 mod_peaks, _ =  find_peaks(mod_scope_data, threshold=3,prominence=10, distance=10)
    #             finish = time.time()
    #             if finish - start > 5:
    #                 mod_peaks = []
    #                 break


    #     return mod_scope_data, mod_peaks

  
    # def fir_measure(
    #     c: LinienClient,
    #     mod_scope_data,
    #     ini_peaks,
    #     mod_peaks, 
    #     mod_fre: float, 
    #     reso: float,
    #     ):
    #     # roughly determine the position of the 1st order peak feom the second TEM00 peak
    #     # using cavity scanning the peak shifts to lower voltage
    #     fir_posi_est = int(-(mod_fre%1e9) / reso + ini_peaks[1])  
        
    #     fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    #     # re-measure the height of this peak by decreasing the piezo scanning range
    #     # print(abs(fir_posi_est-fir_posi_mea))
    #     start = time.time()
        
    #     # while abs(fir_posi_est-fir_posi_mea) > 50:
    #     #     p=10
    #     #     mod_peaks, _ = find_peaks(mod_scope_data, threshold=2,prominence=p, distance=10)
    #     #     fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 
    #     #     end = time.time()
    #     #     if end - start > 5:
    #     #         fir_posi_mea = fir_posi_est # "peak is too far from expected"
    #     #         # plt.plot(mod_scope_data)
    #     #         # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
    #     #         # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red', markersize=1)
    #     #         # plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "8", color='orange')
    #     #         # print(fir_posi_est)
    #     #         # print(mod_peaks)
    #     #         # raise RuntimeError("peak is too far from expected")
    #     #     p=p-1
    #     #     # fir_posi_mea = fir_posi_est
    #     #     # plt.plot(mod_scope_data)
    #     #     # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
    #     #     # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
    #     #     # plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "o")
    #     #     # raise KeyboardInterrupt
    #     # # plt.plot(mod_scope_data)
    #     # # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
    #     # # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
    #     # # plt.plot(fir_posi_est, mod_scope_data[fir_posi_est], "8", color='orange')
    #     # # raise KeyboardInterrupt

    #     sweep_amplitude = 0.8
    #     sweep_center = 0
    #     max_scan_range = 0.8*2
    #     #sweep_center = sweep_center + (-1/2  + fir_posi_mea/len(mod_scope_data)) * max_scan_range # find the voltage corresponds to 1st order peak
    #     for sweep_amplitude in [0.5, 0.25, 0.18, 0.1]:
    #     #for sweep_amplitude in  [ 0.2, 0.1, 0.08]:
    #         sweep_center = sweep_center + (-1/2  + fir_posi_mea/len(mod_scope_data)) * max_scan_range # find the voltage corresponds to 1st order peak
    #         max_scan_range = sweep_amplitude*2
    #         set_scan_range(c, sweep_center, sweep_amplitude)   
    #         time.sleep(0.5)
    #         mod_scope_data, mod_peaks = mod_peak_search(c)
    #         if len(mod_peaks) == 0: # peak too small
    #             return 0, sweep_center
    #         if sweep_amplitude > 0.18:
    #             fir_posi_est = len(mod_scope_data)/2 
    #             fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) 
    #         else:
    #             fir_posi_mea = np.argmax(mod_scope_data)
    #     fir_amp = mod_scope_data[fir_posi_mea]
    #     # fir_amps = []
    #     # sweep_center = sweep_center + (-1/2  + fir_posi_mea/len(mod_scope_data)) * max_scan_range # find the voltage corresponds to 1st order peak
    #     # for offset in [0]:
    #     #     sweep_amplitude = 0.05
    #     #     set_scan_range(c, sweep_center+offset, sweep_amplitude)   
    #     #     time.sleep(1)
    #     #     mod_scope_data, mod_peaks = mod_peak_search(c)
    #     #     if len(mod_peaks) == 0:
    #     #         fir_amp = 0
    #     #     else:
    #     #         fir_amp = np.max(mod_scope_data[mod_peaks])
    #     #     fir_amps.append(fir_amp)
    #     # plt.plot(mod_scope_data)
    #     # plt.plot(mod_peaks, mod_scope_data[mod_peaks], "x")
    #     # plt.plot(fir_posi_mea, mod_scope_data[fir_posi_mea], "o", color='red')
    #     return fir_amp, sweep_center

    # def fir_measure_new(
    #     c: LinienClient,        
    #     mod_scope_data,
    #     ini_peaks,
    #     mod_peaks, 
    #     mod_fre: float, 
    #     reso: float,
    #     ):
    #     # roughly determine the position of the 1st order peak feom the second TEM00 peak
    #     # using cavity scanning the peak shifts to lower voltage
    #     fir_posi_est = int(-(mod_fre%1e9) / reso + ini_peaks[1])  
        
    #     fir_posi_mea = find_nearest(mod_peaks, fir_posi_est) # find the nearest peak position 

    #     fir_amp = mod_scope_data[fir_posi_mea]

    #     return fir_amp, fir_posi_mea