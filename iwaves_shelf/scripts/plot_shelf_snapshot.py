"""
Generates a series of animations
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


from soda.dataio.suntans.sunpy import Spatial
from soda.dataio.suntans.sunprofile import Profile
from soda.dataio.suntans.sunslice import SliceEdge
from soda.utils import timeseries
from soda.utils.maptools import plotmap,readShpPointLine
from soda.utils.otherplot import axcolorbar

from mycurrents import oceanmooring as om

import pdb

##########

basedir = 'SCENARIOS/data_Apr_Nk75_U110'
basedir = 'data'
ncfile = '%s/IWaveRidge_00*nc'%basedir
scenario = 'Apr_Nk75_U110'

tstarts = [
        #'20000104.0400',\
        #'20000104.0530',\
        '20000104.0700',\
        '20000104.0830',\
        '20000104.1000',\
        '20000104.1130',\
        '20000104.1300',\
        '20000104.1430',\
        '20000104.1600',\
        '20000104.1730',\
        ]

tstarts = [
        '20000105.0700',\
        '20000105.0830',\
        '20000105.1000',\
        '20000105.1130',\
        '20000105.1300',\
        '20000105.1430',\
        '20000105.1600',\
        '20000105.1730',\
        ]



#tstarts = [
#        '20000104.0800',\
#        '20000104.0830',\
#        '20000104.0900',\
#        '20000104.0930',\
#        '20000104.1000',\
#        '20000104.1030',\
#        '20000104.1100',\
#        '20000104.1130',\
#        ]



#labels = ['a','b','c','d','e','f','g','h']
labels = ['b','c','d','e','f','g','h','i']

RHO0 = 1000.


ylim = [-300,0.]
xpt = np.array([5e4, 1.5e5])
ypt = np.array([10, 10])

# SCR400
xyloc = [110000., 10.]

outfile = 'FIGURES/suntans_2d_Apr_snaps'


#clevs = np.arange(8,32,1.0).tolist()
clevs = np.arange(0,30,0.25).tolist()
clevsu = np.arange(-1.0,1.1,0.1).tolist()

clevsbold = [12., 15., 18., 21.]
##########

t = [datetime.strptime(tstart,'%Y%m%d.%H%M') for tstart in tstarts]
dt = timedelta(1./24.)

# Read the shapefile with all of the lines
# Load the spatial object and a slice object
#sun = Spatial(ncfile,klayer=[2],variable='salt')
slice = SliceEdge(ncfile,Npt=500,xpt=xpt,ypt=ypt)

####
## Load the Profile data for a time series
sunfile = '%s/%s_Profile.nc'%(basedir,scenario)
sunTS = Profile(sunfile)

Smod_ts = sunTS(xyloc[0], xyloc[1], None, 'salt')
Smod_ts = om.from_sunprofile(Smod_ts)
#
#Umod_ts = sunTS(xyloc[0], xyloc[1], None, 'uc')
#Umod_ts = om.from_sunprofile(Umod_ts) # Convert to an ocean mooring object

####
# Plot a figure of the timeseries

fig = plt.figure(figsize=(12,14))
ax=plt.subplot2grid((5,2),(0,0), colspan=2)


#Umod_ts.contourf(clevsu, \
#    cbar=False, filled=True, cmap='RdBu_r', extend='both')
Smod_ts.contourf(clevs, \
    cbar=False, filled=False)
Smod_ts.contourf(clevsbold, \
    cbar=False, filled=False, linewidths=1.)


plt.xlim(t[0]-dt, t[-1]+dt)
plt.ylim(-250,0)
#plt.xlabel('Time [hours since %s]'%(tstarts[0][0:8]))
plt.xlabel('Time [mm-dd HH]')
plt.ylabel('Depth [m]')
plt.xticks(rotation=0)

for ll,tt in zip(labels, t):
    plt.plot([tt,tt],[-400,0],'--',c='0.5')
    plt.text(tt, -350, '(%s)'%ll)

#plt.show()
####

row=-1
ii=0
def get_rowcol(ii):
    return ii//2,ii%2

for ll, tstart in zip(labels,tstarts):
    row,col = get_rowcol(ii)
    row+=1
    ii+=1
    # Get the time steps
    tstep = slice.getTstep(tstart,tstart)
    print( tstep)
    slice.tstep=[tstep[0]]

    #Load the data and calculate the mean

    print('Loading the sliced data...')
    temp = slice.loadData(variable='rho')*RHO0
    #temp = slice.loadData(variable='temp')
    uc = slice.loadData(variable='uc')
    #S_s = S_s.mean(axis=0)



    #clevs=map(float,range(0,35))
    #clevs = np.arange(20,32,0.25).tolist()
    
    ax=plt.subplot2grid((5,2),(row,col))
    slice.clim = [-1.2,1.2]
    #h2 = slice.contourslice(uc,xaxis='xslice',clevs=clevsu,\
    #    colorbar=False,titlestr='',outline=True,filled=True,\
    #    cmap='RdBu_r',extend='both')

    slice.clim = [0,34]
    #h1 = slice.contourslice(rho,xaxis='xslice',clevs=clevs,\
    #    colorbar=False,titlestr='',outline=True,filled=False,\
    #    #cmap='gist_ncar',\
    #    colors='k',linewidths=0.4)
    h1 = slice.contourslice(temp,xaxis='xslice',clevs=clevs,\
        colorbar=False,titlestr='',outline=True,filled=False,\
        #cmap='gist_ncar',\
        colors='k',linewidths=0.2)

    #h1 = slice.contourslice(temp,xaxis='xslice',clevs=clevsbold,\
    #    colorbar=False,titlestr='',outline=True,filled=False,\
    #    #cmap='gist_ncar',\
    #    colors='k',linewidths=1.0)

    ax.tick_params(direction='out')

    ax.set_ylim(ylim)
    if col==1:
        ax.set_yticklabels([])
        plt.ylabel('')

    if not row == 4:
        ax.set_xticklabels([])
        plt.xlabel('')

    # Plot the mooring location
    plt.plot([xyloc[0],xyloc[0]],[-400,0], '--', lw=2, c='0.5')

    #if ii == 1:
    #    cb=axcolorbar(h2[0],ax=ax,pos=[0.45,0.12,0.4,0.07],ticks=[-1.2,0.,1.2])

    #h2 = slice.contourslice(S_s,xaxis='distslice',clevs=majclevs,\
    #    colorbar=False,titlestr='',outline=True,filled=False,\
    #    #cmap='gist_ncar',\
    #    colors='k',linewidths=0.8)
    #plt.clabel(h2[0],fmt='%2.0f',inline_spacing=6)

    #if row<2:
    #    ax.set_xticklabels([])
    #    plt.xlabel('')
    #ax.tick_params(direction='out')
 
    #txtlabel = datetime.strftime(sun.time[sun.tstep[0]],'%B %Y')#Month YYYY
    plt.text(0.92,0.05,'(%s)'%ll,transform=ax.transAxes,\
        fontsize=16,fontstyle='italic',zorder=1e7)

#cbarpos = [0.63,0.13,0.15,0.02]
#cbaxes = fig.add_axes(cbarpos)
#cb = fig.colorbar(h2[0],cax = cbaxes,orientation='horizontal')
#cb.ax.set_title('Salinity [psu]')
#cb.set_ticks([0.,9.0,18.0,34.0])
#cb.set_ticklabels([0.,9.,18.,34.])
#plt.subplots_adjust(hspace=0.05,wspace=0.05)
#plt.subplots_adjust(hspace=0.05,wspace=0.05,right=0.96,left=0.1,top=0.96,bottom=0.1)
plt.tight_layout()

plt.savefig('%s.pdf'%outfile, dpi=150)
plt.savefig('%s.png'%outfile, dpi=150)
plt.show()
