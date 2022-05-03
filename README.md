# SignalBrowser
SignalBrowser is a psychophysiology Signal Processing and Analysis Toolbox written in Matlab to load, process and analyze psychophysiological data (ECG, respiration and skin conductance). This software can do the following:
-	Filter and process ECG, respiration (RESP), and skin conductance (SCL)
-	Automatically detect and subsequently manually edit ECG R-peaks to compute the RR or interbeat-interval (IBI) signal 
-	Compute ECG derived respiratory signal (EDR)
-	Compute several Heart Rate Variability measures of IBI in the time (Mean RR, RMSSD…) and frequency (ultra low, low and high relative and absolute frequency power) domains.  
-	Compute several frequency domain measures of RESP and EDR 
-	Compute Galvanic skin conductance level (SCL) or responses (SCR) 
-	Compute the time frequency power spectrum using wavelet (Morlet, Dog, Paul) or Fourier Transforms (second screenshots on next page, bottom plot) 
-	HRV or RESP/EDR power spectrum analyses can be performed using various computational methods (Fourier transform, Wavelet, Welch, Welch-Spargle, Autoregressive based power spectrum density calculation) and several conditioning window filtering (Hamming, Hanning, Square, Barlett, Blackman…)

The RR editing view allows you to principally load ECG and other signals and compute/edit ECG R-peaks to extract the RR or interbeat-interval signal (IBI). It can also allow you to define 
1.	Some “Sample” data segments (manually or interactively with a click and drag of the mouse) from which outcome measures are computed 
2.	Some noisy “Signal” data segments (manually or with the mouse) which are excluded from the computation of outcome measures. 

The All Signal view allows you to scroll through and visualize all signals (IBI, RESP, EDR, SCL), as well as outcome measures for each of the “Sample” data segments. It also computes and displays the time-frequency power spectrum of IBI and the time course of the autonomic nervous system (ANS) represented by the low, high and low to high ratio of the HRV frequency spectrum and related time course of the frequency of the peak respiratory sinus arrhythmia (RSA; frequency of the high frequency HRV power spectrum) and corresponding respiratory signal (RESP, EDR)

![RREditingView1](/Pictures/RREditingView1.jpg)

![RREditingView2](/Pictures/RREditingView2.jpg)

![AllSignalView](/Pictures/AllSignalsView.jpg)
