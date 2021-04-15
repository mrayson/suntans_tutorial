"""
Plot a vertical slice of the data
"""

from datetime import datetime
import numpy as np

from sfoda.suntans.sunslice import SliceEdge

import matplotlib
#matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Set some default parameters
matplotlib.rcParams['text.color']='white'
matplotlib.rcParams['savefig.facecolor']='black'
matplotlib.rcParams['savefig.edgecolor']='black'
matplotlib.rcParams['figure.facecolor']='black'
matplotlib.rcParams['figure.edgecolor']='black'
matplotlib.rcParams['axes.facecolor']='black'
matplotlib.rcParams['axes.edgecolor']='white'
matplotlib.rcParams['axes.labelcolor']='white'
matplotlib.rcParams['xtick.color']='white'
matplotlib.rcParams['ytick.color']='white'
matplotlib.rcParams['font.family']='serif'

import pdb

def main(ncfile):
    #
    ###
    # Inputs

    #ncfile = 'rundata/IWave_nosponge_0000.nc'

    xpt = np.array([0.5e2, 3.19e5])
    #xpt = np.array([1e1, 1.5e5])
    ypt = np.array([1., 1.])
    
    #clevs = np.linspace(0,2.4,12)
    clevs = np.arange(0,27,0.25)

    varname = 'uc'
    clim = [-0.8, 0.8]

    #varname = 'w'
    #clim = [-0.10, 0.10]

    #varname = 'salt'
    #clim = None

    #varname = 'temp'
    #clim = None


    cmap = 'RdBu'

    t0 = 0
    ######

    # Load the slice object
    sun = SliceEdge(ncfile, xpt=xpt, ypt=ypt, MAXITER=20000,\
       edgemethod=0, abortedge=False)
    sun.tstep = [t0]

    if clim is not None:
        sun.clim = clim

    # This removes masked values
    cmap = getattr(matplotlib.cm, cmap)
    cmap.set_bad('k')

    #print sun.j

    # Load the data
    data = sun.loadData(variable = varname, method='max')


    # Build the figure
    fig = plt.figure(figsize = (12,8), num = 'SUNTANS Vertical Slice Tool')

    # Time slider axes
    axtime = plt.subplot2grid((7,3), (6,0), colspan=2, rowspan=1)

    ax = plt.subplot2grid((7,3), (0,0), colspan=3, rowspan=5)
    #sun.contourslice(data,t=0)
    h1, axcb, title = sun.pcolorslice(data, cmap = cmap,\
        bathyoverlay=True)

    # Load the density
    rho = sun.loadData(variable ='rho', method='max')*1000.
    sun.clim = [clevs[0], clevs[-1]] 
    h2 = sun.contourslice(rho, clevs=clevs,\
        filled=False, colors='k',linewidths=0.5)

    # Create the time slider
    valstr = ' of %d'%(sun.Nt-1)
    ts = Slider(axtime, 'Time', 0, sun.Nt-1, valinit=t0, valfmt='%d'+valstr,facecolor='0.5',)

    # On change of slider: load new data and set the plot objects accordingly
    def update_slider(val):
        t = int(np.floor(val))
        if not sun.tstep[0] == t:
            sun.tstep = [t]


            ax.collections=[]

            # Update the pcolor plot
            data = sun.loadData(variable = varname, method='max')
            #h1.set_array(data[:-1,:-1].ravel())
            sun.clim = clim
            h1, axcb, title = sun.pcolorslice(data, cmap = cmap,\
                bathyoverlay=True, titlestr='', colorbar=False )

            # Update the contour plot
            rho = sun.loadData(variable ='rho', method='max')*1000.
            sun.clim = [clevs[0], clevs[-1]] 
            h2 = sun.contourslice(rho, clevs=clevs,titlestr='',\
                filled=False, colors='0.5',linewidths=0.5,)

            title.set_text('%s [%s]\n%s'%(sun.long_name, sun.units,\
                datetime.strftime(sun.time[t], '%Y-%m-%d %H:%M:%S')))
            #fig.canvas.draw_idle()

    ts.on_changed(update_slider)

    plt.tight_layout()

    plt.show()

if __name__=='__main__':
    import sys
    ncfile = sys.argv[1]
    main(ncfile)
