
"""
Makes a "movie" hodogram. 

xvals, yvals - must be the same size. 

pauseT = number of seconds to pause b/t each iteration
npts = number of data points to plot in each step. 



"""


def plot_hodogram_dynamic(xvals, yvals, npts=3, plot_kwargs={}, pauseT=0.1):

    import matplotlib.pyplot as plt

    plt.scatter(xvals,yvals,color='blue',marker='+')
    plt.gca().set_aspect(1)

    if 'xlim' in plot_kwargs.keys():
        plt.xlim(plot_kwargs['xlim'])
        plt.ylim(plot_kwargs['xlim'])
    if 'xlabel' in plot_kwargs.keys():
        plt.xlabel(plot_kwargs['xlabel'])
    if 'ylabel' in plot_kwargs.keys():
        plt.ylabel(plot_kwargs['ylabel'])
    if 'title' in plot_kwargs.keys():
        plt.title(plot_kwargs['title'])


    for f in range(len(xvals)):
        points = plt.scatter(xvals[f:npts+f],yvals[f:npts+f],color='orange')
        plt.pause(pauseT)
        points.remove()

