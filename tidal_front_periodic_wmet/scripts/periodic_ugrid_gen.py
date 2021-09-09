
from sfoda.ugrid.hybridgrid import HybridGrid
import numpy as np

def periodic_ugrid_gen(X, Y, periodic_x=False, periodic_y=False, suntanspath=None):
    """
    Creates a curvilinear mesh from grid corners points stored
    in arrays X and Y.
    """

    ny,nx = X.shape

    XY = np.vstack((X.ravel(),Y.ravel())).T

    cells=[]
    neighs = []
    xp = []
    yp = []

    def pntindx(j,i,nrows):
        return j*nrows + i

    def neighidx(j, i, nrows, ncols):
        """
        Return index for the cell neighbours
        """
        N = nrows-1
        M = ncols-1

        neigh = []

        if periodic_x and periodic_y:
            assert 'Error - only set one periodic direction (x OR y)'

        if (i < 0) | (i == M) | (j < 0) | (j == N):
             neigh1 = -1
        else:
             neigh1 = j*M+i

        if periodic_x:
            if i  < 0: # First column --> go to last
                neigh1 = (j+1)*M-1 
            elif i  >= M: # last column --> go to first
                neigh1 = j*M + 0
        if periodic_y:
            if j < 0: # First row --> go to last
                neigh1 = N*M - M + i
            elif j >= N: # last row --> go to first
                print(j,i,N,M)
                neigh1 = i
        #    else:
        #        neigh1 = j*N + i

        #elif periodic_y:
        #    if j < 0: # First row --> go to last
        #        neigh1 = M*N+i
        #    elif j  > N: # Last row --> back to first row
        #        neigh1 = 0*N+i
        #    else:
        #        neigh1 = j*N + i
        #else:
        #    neigh1 = j*N+i

        return neigh1   

    for jj in range(ny):
        for ii in range(nx):
            #if mask[jj,ii]:
            #xp.append(xgrd[ii])
            #yp.append(ygrd[jj])
            xp.append(X[jj,ii])
            yp.append(Y[jj,ii])

    for jj in range(ny-1):
        for ii in range(nx-1):
            cells.append([pntindx(jj,ii,nx), pntindx(jj+1,ii,nx),\
                pntindx(jj+1,ii+1,nx),pntindx(jj,ii+1,nx)])
            neighs.append([neighidx(jj,ii-1,ny,nx), neighidx(jj+1,ii,ny,nx),\
                neighidx(jj,ii+1,ny, nx),neighidx(jj-1,ii,ny,nx)])
            

    Nc = len(cells)
    nfaces = 4*np.ones((Nc,),np.int)

    cells = np.array(cells,dtype=np.int32)
    neighs = np.array(neighs, dtype=np.int32)
    # Convert to a suntans grid
    grd = HybridGrid(np.array(xp),np.array(yp),cells,nfaces=nfaces, neigh=neighs)
    #neighs = None
    #grd = HybridGrid(np.array(xp),np.array(yp),cells,nfaces=nfaces, neigh=neighs)

    if not suntanspath is None:
        print('Writing grid...')
        grd.write2suntans(suntanspath)

    return grd

####
# Test
dx = 500.
ny = 3
nx = 2
xlims = [0,nx*dx]
ylims = [0,ny*dx]
xgrd = np.arange(xlims[0],xlims[1]+1.0*dx,dx)
ygrd = np.arange(ylims[0],ylims[1]+1.0*dx,dx)

# Create a mask polygon
X,Y = np.meshgrid(xgrd,ygrd)

grd = periodic_ugrid_gen(X,Y, periodic_x=True, periodic_y=False)
#grd = periodic_ugrid_gen(X,Y, periodic_x=False, periodic_y=True)

print(grd.mark)
#print(grd.neigh)
