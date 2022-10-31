"""
Calculate transfer function 

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians

"""

#Input frequency spectrum 
Xs = []

#Output frequency spectrum

