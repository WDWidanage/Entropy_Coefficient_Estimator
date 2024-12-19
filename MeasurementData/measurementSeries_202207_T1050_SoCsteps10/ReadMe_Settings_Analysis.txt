Settings and Analysis measurement data
======================================================================================

Analysis of measurement data...

- ... used reference profiles for measurements
	-> potentiometric: ./refSig/potentiometric_T10T50.csv
	-> freq controlled: ./refSig/refSig_freq_T10T50.csv

- ... are relaxation used data sufficient?
	-> Relaxation criteria
	- t_min = 30 min (no t_max!)
	- dT/dt = 20 mK/h
	- dU/dt = 10 mV/h

- ... do we need a pause between potentiometric and frequency controlled measurements?
	-> atm no pause but I could wait 1 or 2 hours after potentiometric measurement (until it is at T30 again) before starting potentiometric measurement 
	
- ... are written interval and length of data sets as needed?

- ... temperature distribution over cell. T at
	-> surfaces
	-> tabs
	-> ...
	=> as we discussed last Friday (2022/08/05)
	
- ... some further information about the measurement series:
	- at T25: 100% SoC achieved via data sheet recommendation (CC with current C/2, CV until C/40 reached) => by definition SoC 100 %
	- afterwards T change to T30 -> serves as "reference temperature"
	- discharge steps of 10 % SoC (Ah based derived from C_std determination -> C/2 CC discharge at T25 before bringing cell to SoC 100 %)
	- for SoC change: above SoC 30 % 3 h relaxation time; from SoC 30 % (and below): 4 h