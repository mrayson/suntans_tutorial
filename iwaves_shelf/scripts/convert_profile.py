"""
Convert the profile binary files into netcdf format
"""
from soda.dataio.suntans.sunprofile import save_profile_nc


def main(basedir, outfile):
    basetime = '20000101.0000'
    save_profile_nc(basedir, outfile, basetime, varnames=['temp','salt'])


if __name__=='__main__':
    import sys
    #basedir = 'sunhex100m'
    #outfile = '../PROFILES/ScottReef3D_hex100m_Profile.nc'
    basedir = sys.argv[1]
    outfile = sys.argv[2]
    print( outfile)
    main(basedir, outfile)
#
#p = Profile('%s/%s'%(basedir, outfile))
#
#temp = p(5e5, 6e7, 100, 'temp')
