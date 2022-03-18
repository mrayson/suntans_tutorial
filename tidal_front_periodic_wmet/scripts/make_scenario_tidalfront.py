# -*- coding: utf-8 -*-
"""
Create a one cell quad grid

Created on Mon Nov 18 10:23:26 2013

@author: mrayson
"""

import os
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import operator
from scipy.interpolate import interp1d

from sfoda.ugrid.ugridgen import cartesian_ugrid_gen
import periodic_ugrid_gen as ugridgen
from sfoda.suntans.sunboundary import modifyBCmarker, Boundary, InitialCond
from sfoda.suntans.metfile import SunMet
from sfoda.suntans.sunpy import Grid

PI = np.pi
GRAV=9.81

###
# Bathymetry generation functions
def lambda_factory(ab):
    return lambda x:x*ab[0]+ab[1]

def broken_line(x, x0, y0):
    cl = []
    fl = []
    for i in range(len(x0)-1):
        ab = np.polyfit(x0[i:i+2], y0[i:i+2], 1)
        # Compute and append a "condition" interval
        cl.append(np.logical_and(x >= x0[i], x <= x0[i+1]))
        # Create a line function for the interval
        fl.append(lambda_factory(ab))
    return(np.piecewise(x, condlist=cl, funclist=fl))

def gauss_kernel(L, dx):
    X = np.arange(-2*L,2*L+dx,dx)
    Y = np.exp(- (X/L)**2 )
    return Y/np.sum(Y)

def gaussian(x,x0,a,b,pow=2.):
    return b*np.exp(-(x-x0)**pow/(2*a**pow))

def make_suntans(suntanspath):
    ####################################################
    # Inputs

    # Boundary forcing parameters
    wave_period=12.42*3600.
    U0 = 0.25        # Barotropic velocity [m/s]

    # Temperature intitial condition parameters
    T0 = 23.0
    dTdz = 0.07

    # Size of domain
    L = 370e3  # Domain length [m]
    W = 50e3 # Domain width [m]
    H = -250 # Maximum depth

    # Number of grid cells along each dimensions
    ny = 250 # 32
    nx = 200 #192
    #gridfolder=''

    nz = 40
    rstretch = 1.0

    # Bathymetry parameters
    x0 = np.array([0, 50e3, 60e3, 180e3, 300e3, 310e3, L])
    y0 = np.array([H,H, H+130, H+190, H+130,H,H])


    # Boundary condition time
    starttime = '20000601.000000'
    endtime = '20000701.000000'
    dt = 3600.

    # Filenames initial and boundary condition file names
    icfile = 'TFront_IC.nc'
    bcfile = 'TFront_BC.nc'
    #####################################################

    ######
    # Generate the bathymetry
    Lsmooth = 6000
    dx = L/nx
    #x = np.arange(0, L+dx, dx)
    x = np.linspace(0, L, nx)
    #nx = x.shape[0]
    h =  broken_line(x, x0, y0)
    K = gauss_kernel(Lsmooth, dx)
    Nk = K.shape[0]//2
    h_smooth = np.zeros_like(h)
    h_smooth[Nk:-Nk+1] = np.convolve(h, K, mode='valid')
    h_smooth[:Nk]=H
    h_smooth[-Nk:]=H
    
    # Bathymetry generation function
    # Depths should be positive
    F_bathy = interp1d(x, -h_smooth, kind=3)

    #######
    # Compute wave parameters to get the domain size etc
    k_z = PI / H        
    omega = 2*PI/wave_period
    W = ny*dx

    #####
    # Create the grid
    if not os.path.isdir(suntanspath):
        print('Creating new directory: %s'%suntanspath)
        os.mkdir(suntanspath)
        copyfile('rundata/suntans.dat','%s/suntans.dat'%suntanspath)

    xlims = [0,L]
    ylims = [0,W]

    # Create the grid
    #grd= cartesian_ugrid_gen(xlims, ylims, dx, suntanspath=suntanspath)
    #xgrd = np.arange(xlims[0],xlims[1]+1.0*dx,dx)
    #ygrd = np.arange(ylims[0],ylims[1]+1.0*dx,dx)
    #X,Y = np.meshgrid(xgrd,ygrd)
    #grd= ugridgen.periodic_ugrid_gen(X, Y, suntanspath=suntanspath,
    #    periodic_y=True, periodic_x=False)
    # Copy the grid files across

    # Load the grid
    grd = Grid(suntanspath)

    #grd.dv = H*np.ones_like(grd.xv)
    grd.dv = F_bathy(grd.xv)
    grd.saveBathy('%s/depth.dat-voro'%suntanspath)

    grd.dz = grd.calcVertSpace(nz, rstretch, H)
    grd.saveVertspace('%s/vertspace.dat'%suntanspath)
    grd.setDepth(grd.dz)

    # Temperature field
    temp = T0 + dTdz*grd.z_r
    print(grd.z_r, temp)

    # Create the boundary conditions

    ##########
    # Modify the boundary markers and create the boundary condition file
    ##########
    # This changes the edge types based on the polygons in the shapefile
    #modifyBCmarker(suntanspath,bcpolygonfile)

    # Modify the left and right edges and convert to type 2
    #hgrd = grd.convert2hybrid()

    grd.mark[grd.mark>0]=1 # reset all edges to type-1

    ## convert edges +/- half a grid cell from the edge
    #dx = hgrd.dg.max()
    xmin = grd.xv.min()-dx/4.0
    xmax = grd.xv.max()+dx/4.0

    grd.calc_edgecoord()

    indleft = operator.and_(grd.mark==1, grd.xe < xmin) # all boundaries
    indright = operator.and_(grd.mark==1, grd.xe > xmax) # all boundaries

    ## Free-surface boundaries
    ##grd.mark[indleft]=3
    ##grd.mark[indright]=3
    #
    ## Velocity bouundaries
    grd.mark[indleft]=2
    grd.mark[indright]=2
    ##grd.edge_id[indleft]=1
    ##grd.edge_id[indright]=2
    #
    edgefile = suntanspath+'/edges.dat'
    grd.saveEdges(edgefile)
    #print 'Updated markers written to: %s'%(edgefile)
    #grd.write2suntans(suntanspath)

    #hgrd.dv=grd.dv
    #hgrd.z_r = grd.z_r
    #grd = hgrd

    #Load the boundary object from the grid
    #   Note that this zeros all of the boundary arrays
    bnd = Boundary(suntanspath,(starttime,endtime,dt))

    bnd.setDepth(grd.dv)

    t = bnd.tsec-bnd.tsec[0]
    
    # Velocity boundary
    for k in range(bnd.Nk):
       for ii, eptr in enumerate(bnd.edgep.tolist()):
           if indleft[eptr]:
               #bnd.boundary_u[:,k,ii] = phi0*np.cos(k_z*grd.z_r[k])*np.sin(omega*t)
               bnd.boundary_u[:,k,ii] = U0*np.sin(omega*t)
               bnd.boundary_T[:,k,ii] = temp[k]# - phi0*np.sin(k_z*grd.z_r[k])*np.sin(omega*t)
           elif indright[eptr]:
               bnd.boundary_u[:,k,ii] = U0*np.sin(omega*t)
               bnd.boundary_T[:,k,ii] = temp[k]
               #bnd.boundary_u[:,k,ii] = u*nx[ii]
               #bnd.boundary_v[:,k,ii] = u*ny[ii]
               #bnd.boundary_h[:,ii] = h
    #
   # Write the boundary file
    bnd.write2NC(suntanspath+'/'+bcfile)

    #########
    # Create the initial conditions file
    #########
    IC = InitialCond(suntanspath,starttime)

    #seiche=Seiche(L,W,H,A)
    #u,IC.h[:] = seiche(IC.xv,0)
    IC.h[:] = 0

    IC.S[:,nz//2:, :] = 1.0 # set S = 1 in the lower water column
    IC.T[:,:,:] = temp[np.newaxis,:,np.newaxis]

    # Write the initial condition file
    IC.writeNC(suntanspath+'/'+icfile,dv=grd.dv)

    #grd.plotmesh()
    #plt.show()
    makewinds = 1
    if makewinds:
        ####
        # Atmospheric input options
        ####

        metfile = 'metfile.nc'

        Uwind = 4.0
        Vwind = 0.0
        RH = 40.0
        Tair = 20.0
        Pair = 1016.0
        cloud = 0.
        ####
        # End of options
        ####


        xpt = grd.xv.mean()
        ypt = grd.yv.mean()
        zpt = 10.0
        met = SunMet(xpt,ypt,zpt,(starttime,endtime,dt))
        
        # Set a cyclic air temperature and all other variables constant (rain and cloud =0)
        #omegaT = 2*PI/(24.0 * 3600.0) # diurnal frequency
        #for ii in range(met.Tair.Npt):
        #    met.Tair.data[:,ii] = Tair_mean + Tair_amp * np.sin(omegaT * met.nctime)
            
        met.Tair.data[:] = Tair
        met.Pair.data[:] = Pair
        met.Uwind.data[:] = Uwind
        met.Vwind.data[:] = Vwind
        met.RH.data[:] = RH 
        met.cloud.data[:] = cloud
        met.rain.data[:] = 0.0

        met.write2NC(suntanspath+'/'+metfile)

if __name__=='__main__':
    import sys
    sunpath = sys.argv[1]

    make_suntans(sunpath)

