
# coding: utf-8

# # Fit vertical mode to SUNTANS harmonic data
# 
# $$
# b(z,t) = A(t)N^2(z)\phi(z)
# $$

# In[1]:



#from soda.dataio.suntans.suntides import suntides
from soda.dataio.suntans.sunpy import Spatial
from soda.dataio.suntans.suntans_ugrid import ugrid
from soda.utils import mynumpy as mynp
from soda.utils import othertime

from iwaves.utils.isw import iwave_modes_uneven

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg as la

from glob import glob
import pdb

def calc_alpha(phi, c, N2, z):
    """
    Nonlinearity parameter
    """
    #phi_z = np.gradient(phi,-np.abs(dz))
    z2d = np.repeat(z[:,np.newaxis],phi.shape[1], axis=1)
    phi_z = mynp.grad_z(phi, z2d)

    #phi_z = mynp.grad_z(phi, z, axis=-1)
    num = 3*c*np.trapz( phi_z**3., x=z2d, axis=0)
    den = 2*np.trapz( phi_z**2., x=z2d, axis=0)

    return num/den


# In[2]:


###
nmodes = 5
mindepth = 100.
###


# In[3]:


def main(sunfile, outfile):
    #Constants
    RHO0=1024.
    g=9.81

    # Load the object
    sun=Spatial(sunfile, klayer=[-99])
    #Nf=sun.frq.shape[0]
    Nf=sun.time.shape[0]

    # Re-calculate the grid variables
    #sun.reCalcGrid()
    sun.calc_def()


    ## calculate the buoyancy frequecny
    sun.tstep=0
    rho = sun.loadData(variable='rho')*1000.0
    #rho = sun.Mean['rho']*1000

    # Get the mask array
    mask = rho >= 999999.
    mask3d = np.tile(mask,(Nf,1,1))

    rho[mask]=0

    print ('Calculating N^2...')
    N2 = -g/RHO0 * sun.gradZ(rho)
    N2[N2<1e-5]=1e-5
    N2[mask]=0

    # Put everything in complex form
    #print 'Converting arrays into complex form...'
    #rho = sun.Amp['rho']*1000.0*np.cos(sun.Phs['rho']) + 1j * sun.Amp['rho']*1000.0*np.sin(sun.Phs['rho'])
    print('Loading all rho...')
    sun.tstep =range(sun.Nt)
    b = sun.loadData(variable='rho')*1000.0
    b -= rho[np.newaxis,...]

    # Mask these arrays
    b[mask3d]=0

    # Calculate the buoyancy perturbation
    b *=  -g/RHO0 

    nf, nz, nc = b.shape

    def mid_extrap(A, nk):
        B = np.zeros((nk+1,), dtype=A.dtype)
        B[1:-1] = 0.5*A[0:-1]+0.5*A[1::]
        B[0] = B[1]
        B[-1] = B[-2]
        return B

    # Create the output arrays
    cnall = np.zeros((nmodes, nc))
    alphaall = np.zeros((nmodes, nc))
    phiall = np.zeros((nz+1, nc, nmodes)) # Put modes last so we can array multiply
    ampall = np.zeros((nf, nc, nmodes))
    N2out = np.zeros((nz+1, nc)) # Put modes last so we can array multiply


    maxamp = 0.
    for ii in range(nc):
        if ii%10 == 0:
            print( "%d of %d (max amp = %3.1f m)..."%(ii,nc, maxamp))
            maxamp = 0.

        # Extract the N2 profile at a point
        nk = sun.Nk[ii]
        
        if nk < 10:
            continue


        myN2 = N2[0:nk,ii]
        z = -sun.z_r[0:nk]
        zw = -sun.z_w[0:nk+1]

        # Interpolate N2 onto the edge points
        N2mid = mid_extrap(myN2, nk)

        # Calculate the mode shapes
        phi, cn = iwave_modes_uneven(N2mid, zw)

        alpha_n = calc_alpha(phi[:,0:nmodes], cn[0:nmodes], N2mid, zw)

        # for each frequency
        for ff in range(nf):
        #for ff in range(100,101):
            # Get buoyancy at cell-edges
            myb = b[ff,0:nk,ii]
            bmid = mid_extrap(myb, nk)

            # Compute the LHS
            L = phi[:,0:nmodes] * N2mid[:,np.newaxis]

            myamp,_,_,_ = la.lstsq(L, bmid)

            ampall[ff,ii,:] = myamp
            
            if myamp.max() > maxamp:
                maxamp = myamp.max()*1

        # Output the variables
        cnall[:,ii] = cn[0:nmodes]
        alphaall[:,ii] = alpha_n
        phiall[0:nk+1,ii,:] = phi[:,0:nmodes]
        N2out[0:nk+1,ii] = N2mid

    #plt.subplot(121)
    #plt.plot(myN2,z)
    #plt.plot(N2mid,zw)
    #
    #plt.subplot(122)
    #plt.plot(phi[:,1],zw)
    #
    #plt.show()

    ####
    # Write the output to netcdf
    print( 'Writing the output to netcdf...')

    sun.writeNC(outfile)

    nc = Dataset(outfile,'a')
    nc.Title = 'SUNTANS harmonic modal amplitudes'
    #nc.Constituent_Names = ' '.join(sun.frqnames)

    # Add another dimension
    #nc.createDimension('Ntide', nf)
    #nc.createDimension('time', nf)
    nc.createDimension('Nmode', nmodes)
    nc.close()

    coords = 'xv yv time'

    #sun.create_nc_var(outfile,'omega', ('Ntide',),\
    #        {'long_name':'frequency','units':'rad s-1'})
    sun.create_nc_var(outfile,'time', ugrid['time']['dimensions'], ugrid['time']['attributes'])

    sun.create_nc_var(outfile,'modes', ('Nmode',),\
            {'long_name':'Vertical mode number','units':''})


    sun.create_nc_var(outfile, 'cn', ('Nmode','Nc'),\
            {'long_name':'Baroclinic phase speed','units':'m s-1',\
            'coordinates':coords})

    sun.create_nc_var(outfile, 'alpha_n', ('Nmode','Nc'),\
            {'long_name':'Nonlinearity Parameter','units':'m-1',\
            'coordinates':coords})


    sun.create_nc_var(outfile, 'phi', ('Nkw','Nc', 'Nmode'),\
            {'long_name':'Vertical eigenfunction',\
              'units':'',\
              'coordinates':coords})

    sun.create_nc_var(outfile, 'N2', ('Nkw','Nc'),\
            {'long_name':'Buoyancy frequency squared',\
              'units':'s-2',\
              'coordinates':coords})


    sun.create_nc_var(outfile, 'amp_b_re', ('time','Nc', 'Nmode'),        {'long_name':'Buoyancy real amplitude',         'units':'m s-2',         'coordinates':coords})

    #sun.create_nc_var(outfile, 'amp_b_im', ('Ntide','Nc', 'Nmode'),        {'long_name':'Buoyancy imaginary amplitude',         'units':'m s-2',         'coordinates':coords})

    # Write the data

    print( 'Writing the variable data to netcdf...')

    nc = Dataset(outfile,'a')
    #nc.variables['omega'][:]=sun.frq
    tsec = othertime.SecondsSince(sun.time)
    nc.variables['time'][:]=tsec
    nc.variables['modes'][:]=range(nmodes)
    nc.variables['cn'][:]=cnall
    nc.variables['alpha_n'][:]=alphaall
    nc.variables['phi'][:]=phiall
    nc.variables['N2'][:]=N2out
    nc.variables['amp_b_re'][:]=np.real(ampall)
    #nc.variables['amp_b_im'][:]=np.imag(ampall)

    nc.close()

    print('Done.')


    # In[ ]:





#ncfile = 'SCENARIOS/data_Apr_Nk75/IWaveRidge_0*.nc'
#outfile = 'SCENARIOS/data_Apr_Nk75/IWaveRidge_ModeAmp.nc'

ncfile = 'SCENARIOS/data_Mar_Nk75/IWaveRidge_0*.nc'
outfile = 'SCENARIOS/data_Mar_Nk75/IWaveRidge_ModeAmp.nc'


print( '\n', 72*'#', '\n')
print(outfile)

main(ncfile, outfile)



