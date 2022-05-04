# SignalBrowser
## Description 
SignalBrowser is a graphic user interface software tool written in Matlab to load, process and analyze psychophysiological data (ECG, respiration and skin conductance). This software can do the following:
-	Filter and process ECG, respiration (RESP), and skin conductance (SCL)
-	Automatically detect and subsequently manually edit ECG R-peaks to compute the RR or interbeat-interval (IBI) signal 
-	Compute ECG derived respiratory signal (EDR)
-	Compute several Heart Rate Variability (HRV) measures of IBI in the time (Mean RR, RMSSD…) and frequency (ultra low, low and high relative and absolute frequency power) domain.  
-	Compute several frequency domain measures of RESP and EDR 
-	Compute Galvanic skin conductance level (SCL) or responses (SCR) 
-	Compute the time frequency power spectrum using wavelet (Morlet, Dog, Paul) or Fourier Transforms  
-	Compute the time course of the low (LF), high (HF) and low to high ratio (LF/HF) of the HRV frequency spectrum and related frequency of the peak respiratory sinus arrhythmia (RSA; frequency of the high frequency HRV power spectrum) and corresponding respiratory signal (RESP, EDR) as a measure of the sympathetic (LF) and parasympathetic (HF) balance of the autonomic nervous system.
-	The HRV or RESP/EDR power spectrum analyses can be performed using various computational methods (Fourier transform, Wavelet, Welch, Welch-Spargle, Autoregressive based power spectrum density calculation) and several conditioning window filtering (Hamming, Hanning, Square, Barlett, Blackman…)
-	Export outcome measures in a file for each of the defined analysis “Sample” data segments

The RR editing view allows you to principally load ECG and other signals and compute/edit ECG R-peaks to extract the RR or interbeat-interval signal (IBI). It can also allow you to define 
1.	Some “Sample” data segments (manually or interactively with a click and drag of the mouse) from which outcome measures are computed 
2.	Some noisy “Signal” data segments (manually or with the mouse) which are excluded from the computation of outcome measures. 

The All Signal view allows you to scroll through and visualize all signals (IBI, RESP, EDR, SCL), as well as outcome measures for each of the “Sample” data segments. It also computes and displays the time-frequency power spectrum of IBI and the time course of the autonomic nervous system (ANS) represented by the low, high and low to high ratio of the HRV frequency spectrum and related time course of the frequency of the peak respiratory sinus arrhythmia (RSA; frequency of the high frequency HRV power spectrum) and corresponding respiratory signal (RESP, EDR)

## Disclaimer
This software has been developed organically over 9 yrs ago for my research. This is not the best example of a clean and organized coding :-). I have recently made some udpates to fix some bugs, add a couple of features and make it available with the hope that it can be useful to someone. Improvements to the software organizations, architecture and annotations is a work in progress. Please send me an email if you encounter any bugs or issues. I will try to address them as quickly as I can, with the understanding that this can only be done in my spare time.

The software is provided "as is", without warranty of any kind, express or implied. Please refer to the license document.

## A Short User Manual

![RREditingView1](/Pictures/RREditingView1.jpg)

![RREditingView2](/Pictures/RREditingView2.jpg)

![AllSignalView](/Pictures/AllSignalsView.jpg)
