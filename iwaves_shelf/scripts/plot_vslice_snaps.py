"""
Plot a vertical slice of the data
"""

from datetime import datetime
import numpy as np

from sfoda.suntans.sunslice import SliceEdge

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.animation as animation


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

def main(ncfile, outpath):
    #
    ###
    # Inputs

    #ncfile = 'rundata/IWave_nosponge_0000.nc'

    xpt = np.array([0.5e2, 3.19e5])
    #xpt = np.array([1e1, 1.5e5])
    ypt = np.array([1., 1.])
    
    #clevs = np.linspace(0,2.4,12)
    clevs = np.arange(0,27,0.125)

    varname = 'uc'
    clim = [-0.8, 0.8]

    #varname = 'w'
    #clim = [-0.10, 0.10]

    #varname = 'salt'
    #clim = None

    #varname = 'temp'
    #clim = None


    cmap = 'RdBu'

    #outpath = 'FIGURES/iwave_ridge'

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
    fig = plt.figure(figsize = (12,8), num = 'SUNTANS Vertical Slice')

    # Time slider axes

    ax = plt.subplot2grid((7,3), (0,0), colspan=3, rowspan=6)
    #sun.contourslice(data,t=0)
    h1, axcb, title = sun.pcolorslice(data, cmap = cmap,\
        bathyoverlay=True)

    # Load the density
    rho = sun.loadData(variable ='rho', method='max')*1000.
    sun.clim = [clevs[0], clevs[-1]] 
    h2 = sun.contourslice(rho, clevs=clevs,\
        filled=False, colors='k',linewidths=0.5)

    plt.tight_layout()

    def initanim():
        return h1, h2, title

    def update_slider(val):
        t = int(np.floor(val))
        ax.collections=[]
 
        sun.tstep = [t]
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

        return h1, h2, title


    timedays = sun.time.shape[0]

    anim = animation.FuncAnimation(fig, update_slider,\
        init_func=initanim, frames=timedays,  blit=False)
            #init_func=init, frames=[1,2], interval=50, blit=True)

    anim.save("%s.gif"%outpath, fps=12, dpi=90)

    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=12, bitrate=14400, codec='mpeg4')
    #anim.save("%s.avi"%outpath, writer=writer)
    print('Saved to %s.gif'%outpath)
    #plt.savefig(outfile)
    #print('Figure saved to the output file...')

if __name__=='__main__':
    import sys
    ncfile = sys.argv[1]
    outpath = sys.argv[2]
    main(ncfile, outpath)
