"""
Smoothing routines for Python 

Routines: 
savgol_filter (savgol_filter)
sliding box average (sliding_box_avg)
"""


"""
savgol_filter
https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way

e.g. 
x = np.linspace(0,2*np.pi,100)
y = np.sin(x) + np.random.random(100) * 0.2
ysmoo = savgol_filter(y, window=50, poly=3)

"""
def savgol_filter(data,window=50,poly=3):
    from scipy.signal import savgol_filter    
    #savgol_filter(data, window_size, poly_order)
    return savgol_filter(data, window, poly) 



"""
Sliding box average filter
https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way

E.g.:
x = np.linspace(0,2*np.pi,100)
y = np.sin(x) + np.random.random(100) * 0.8
box_pts = 3
ysmoo = sliding_box_avg(y,box_pts=10)
"""

def sliding_box_avg(data, box_pts=3):
    import numpy as np
    box = np.ones(box_pts)/box_pts
    return np.convolve(data, box, mode='same')






