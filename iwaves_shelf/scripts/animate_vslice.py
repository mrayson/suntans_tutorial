"""
Plot a vertical slice of the data
"""

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.animation as animation
import types

from sfoda.suntans.sunslice import SliceEdge

import matplotlib
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

def fix_quadcontour(c1, ax1, fig):
    #################################################################
    ## Bug fix for Quad Contour set not having attribute 'set_visible'
    def setvisible(self,vis):
        for c in self.collections:    c.set_visible(vis)
    def setanimated(self,vis):
        for c in self.collections:    c.set_animated(vis)
    def setdraw(self,vis):
        for c in self.collections:    c.draw(vis)


    c1.set_visible = types.MethodType(setvisible,c1)
    c1.set_animated = types.MethodType(setanimated,c1)
    c1.draw = types.MethodType(setdraw,c1)

    c1.axes = ax1
    c1.figure=fig

    return c1


def main(ncfile, outfile):

    #######
    # Inputs

    #ncfile = 'SCENARIOS/IWS_NH_qp2_Ny1_Nz100_Nx3000_dt2_hmix30_mode1/Soliton_0000.nc'

    #outfile = 'MOVIES/ISW_linear_hmix30'

    #xpt = np.array([1e5, 1.5e5])
    #xpt = np.array([1e1, 3.0e5])
    xpt = np.array([0.5e5,1.5e5])
    ypt = np.array([1., 2.])

    zmin = -300.

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

    dpi=150
    #########

    # Load the slice object
    sun = SliceEdge(ncfile, xpt=xpt, ypt=ypt, Npt=1000)
    sun.tstep = [t0]

    if clim is not None:
        sun.clim = clim

    # This removes masked values
    cmap = getattr(matplotlib.cm, cmap)
    cmap.set_bad('k')

    #print sun.j

    # Load the data
    data = sun.loadData(variable = varname, method='max')

    ##########
    # Build the figure
    fig = plt.figure(figsize = (8,5))
    ax = plt.subplot(111)

    #sun.contourslice(data,t=0)
    h1, axcb, title = sun.pcolorslice(data, cmap = cmap,\
        bathyoverlay=True, titlestr='' )

    # Load the density
    rho = sun.loadData(variable ='rho', method='max')*1000.
    sun.clim = [clevs[0], clevs[-1]] 
    h2, = sun.contourslice(rho, clevs=clevs,\
        filled=False, colors='k',linewidths=0.5)

    # Hack to fix the quadcontour object
    h2 = fix_quadcontour(h2, ax, fig)
    plt.tight_layout()

    def init():
        return h1, h2#, title

    # On change of slider: load new data and set the plot objects accordingly
    def update_slider(t):
       sun.tstep = [t]

       ax.collections=[]

       # Update the pcolor plot (Not working)
       data = sun.loadData(variable = varname, method='max')
       sun.clim = clim
       h1, axcb, title = sun.pcolorslice(data, cmap = cmap,\
           bathyoverlay=True, titlestr='', colorbar=False )

       #h1.set_array(data[:-1,:-1].ravel())

       # Update the contour plot
       rho = sun.loadData(variable ='rho', method='max')*1000.
       sun.clim = [clevs[0], clevs[-1]] 
       h2, = sun.contourslice(rho, clevs=clevs,titlestr='',\
           filled=False, colors='0.5',linewidths=0.5,)
       # Hack to fix the quadcontour object
       h2 = fix_quadcontour(h2, ax, fig)

       ax.set_ylim(zmin,0)


       #title.set_text('%s [%s]\n%s'%(sun.long_name, sun.units,\
       #    datetime.strftime(sun.time[t], '%Y-%m-%d %H:%M:%S')))

       return h1, h2#, title


    anim = animation.FuncAnimation(fig, update_slider,\
        init_func=init, frames=range(t0, sun.Nt), interval=10, blit=True)

    Writer = animation.writers['mencoder']
    writer = Writer(fps=3, metadata=dict(artist='Matt Rayson'), bitrate=3600)

    anim.save("%s.mp4"%outfile, writer=writer)

    #anim.save("%s.mp4"%outfile, writer='mencoder', fps=6, bitrate=3600)
    #anim.save("%s.mp4"%outfile, writer='ffmpeg', fps=6, bitrate=3600)
    #anim.save("%s.gif"%outfile,writer='imagemagick',dpi=dpi)

    print( 'Saved to %s.gif'%outfile)


    #plt.show()

if __name__=='__main__':
    import sys
    ncfile = sys.argv[1]
    outfile = sys.argv[2]
    main(ncfile, outfile)
