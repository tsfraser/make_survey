import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import fitsio
from pathlib import Path
import healpy as hp
from cosmoprimo.cosmology import Cosmology
#plt.style.use(['enrique-science', 'bright'])
def read_cmass_weights(filename):
    data = fitsio.read(filename)
    weights = data['WEIGHT_SYSTOT'] * (data['WEIGHT_CP'] + data['WEIGHT_NOZ'] - 1)
    return weights
def read_patchy_weights(filename):
    data = np.genfromtxt(filename)
    weights = data[:, 6] * data[:, 7]
    return weights

def spl_nofz(zarray, fsky, cosmo, zmin, zmax, Nzbins=100, weights=None):
    """Compute the number density of galaxies as a function of redshift"""
    zbins = np.linspace(zmin, zmax, Nzbins+1)
    Nz, zbins = np.histogram(zarray, zbins, weights=weights)
    zmid = zbins[0:-1] + (zmax-zmin)/Nzbins/2.0
    # set z range boundaries to be zmin and zmax and avoid the interpolation error
    zmid[0], zmid[-1] = zbins[0], zbins[-1]
    rmin = cosmo.comoving_radial_distance(zbins[0:-1])
    rmax = cosmo.comoving_radial_distance(zbins[1:])
    vol = fsky * 4./3*np.pi * (rmax**3.0 - rmin**3.0)
    nz_array = Nz/vol
    spl_nz = InterpolatedUnivariateSpline(zmid, nz_array)
    return spl_nz

zmin = 0.46
zmax = 0.6
cosmo = Cosmology(
    Omega_m=0.307115,
    Omega_b=0.048,
    sigma8=0.8288,
    h=0.677,
    engine='class'
)

fsky = 6851 / 41253  # taken from Reid at al. 2015​
# CMASS data

data_dir = '/home/epaillas/projects/rrg-wperciva/MOCKS/SDSS-DATA/BOSS/data'
data_fn = Path(data_dir, 'galaxy_DR12v5_CMASS_North.fits.gz')
data_cmass = fitsio.read(data_fn)
weights_cmass = read_cmass_weights(data_fn)

z_cmass = data_cmass['Z']
nz_cmass = spl_nofz(z_cmass, fsky, cosmo, zmin, zmax,
                    Nzbins=25, weights=weights_cmass)​
# Plot n(z)

fig, ax = plt.subplots(figsize=(4, 4))
z = np.linspace(zmin, zmax, 50)
ax.plot(z, nz_cmass(z)*1e4, label='CMASS')
ax.legend()
ax.set_xlabel('z')
ax.set_ylabel('n(z)')
plt.savefig('cmass_nz.png', dpi=300)