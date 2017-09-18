
#written by Keith Tirimba


'''
The purpose of this script is to acquire magnitudes from your spectrum for the purpose of observing color changes as a function of 
redshift.

This first portion of the script assumes you are using the filter curves available in the speclite package. If you're using your 
own curves, use the code below.'''

from astropy.table import Table
import speclite
import numpy as np
import numpy.ma as ma
from astropy.io import fits

files=glob.glob('*.fits')

z_list=np.arange(0,0.25,0.05) # range of redshifts for your spectrum
for k in files:    
    for rs in z_list:
        
        hdu=fits.open(k)
        header=hdu[1].header
        data=hdu[1].data
        
        
        a_rs=hdu[2].data['Z'][0]
        
        wave=(10**data['loglam']/(1+a_rs))*(1+rs)

        flux=data['flux']*1e-17
        
        
        const=np.median(flux[0:200])/((wave[0])**(-2.0))  #median value of the UV spec
        UV_wave=np.arange(2000,wave[0],0.5) #set your left boundary for the UV. 
        UV_normalized_flux=(const*UV_wave**-2.0)  # model of UV. Read (Wilkins et. al 2012) for justification on used values
        wave_full=np.concatenate((UV_wave,wave),axis=0)  # combining the UV
        flux_full=np.concatenate((UV_normalized_flux,flux),axis=0)
        f_error=ma.array((1/np.sqrt(data['ivar']))*1e-17)
        
       #you can load your filters here if they are available in the speclite package.
        #Check here for available curves. It also shows how you can implement your own filter curves
        '''http://speclite.readthedocs.io/en/latest/filters.html'''
    
        filters = speclite.filters.load_filters('sdss2010-g','sdss2010-r','sdss2010-i')  
        
        #Acquires the magnitudes given the data you provide
        mags = filters.get_ab_magnitudes(flux_full, wave_full,mask_invalid=True) 
        
        ##Example
        g_r= (mags['sdss2010-g']-mags['sdss2010-r'])[0]  #values at a specific redshift
        r_i= (mags['sdss2010-r']-mags['sdss2010-i'])[0]
    
########################################################
#ALTERNATIVE IF YOU'D LIKE TO USE YOUR OWN CURVES
from astropy.table import Table
import numpy as np
from scipy import interpolate
from astropy.io import fits

u = Table.read('u.dat.txt', format='ascii')  # , delimiter=' ', guess=False)
g = Table.read('g.dat.txt', format='ascii')  # , delimiter=' ', guess=False)
r = Table.read('r.dat.txt', format='ascii')  # , delimiter=' ', guess=False)
i = Table.read('i.dat.txt', format='ascii')  # , delimiter=' ', guess=False)
z = Table.read('z.dat.txt', format='ascii')  # , delimiter=' ', guess=False)


hdu = fits.open('*.fits')
header = hdu[1].header
data = hdu[1].data

a_rs = hdu[2].data['Z'][0]
z_list = np.linspace(0, 0.1, 100)


# Acquiring the "limits" of the filter in wavelength array
def fil(filt):
    if filt['col1'][0] == 2980:
        val1 = min(wave, key=lambda x: abs(x - 2980))  # locating w_0 in wave
        val2 = min(wave, key=lambda x: abs(x - 4130))
        b = 1.4e-10  # softening parameter
    if filt['col1'][0] == 3630:
        val1 = min(wave, key=lambda x: abs(x - 3630))  # locating w_0 in wave
        val2 = min(wave, key=lambda x: abs(x - 5830))
        b = 0.9e-10
    if filt['col1'][0] == 5380:
        val1 = min(wave, key=lambda x: abs(x - 5380))  # locating w_0 in wave
        val2 = min(wave, key=lambda x: abs(x - 7230))
        b = 1.2e-10
    if filt['col1'][0] == 6430:
        val1 = min(wave, key=lambda x: abs(x - 6430))  # locating w_0 in wave
        val2 = min(wave, key=lambda x: abs(x - 8630))
        b = 1.8e-10
    if filt['col1'][0] == 7730:
        val1 = min(wave, key=lambda x: abs(x - 7730))  # locating w_0 in wave
        val2 = min(wave, key=lambda x: abs(x - 11230))
        b = 7.4e-10
    return val1, val2, b


for rs in z_list:
    u_g=[]
    g_r=[]
    r_i=[]

    wave = (10 ** data['loglam'] / (1 + a_rs)) * (1 + rs)  ###shift wave as a function of redshift

    flux = data['flux'] * 1e-17

    const = np.median(flux[0:200]) / ((wave[0]) ** (-2.0))
    UV_wave = np.arange(2000, wave[0], 0.5)
    UV_normalized_flux = (const * UV_wave ** -2.0)
    wave = np.concatenate((UV_wave, wave), axis=0)
    flux = np.concatenate((UV_normalized_flux, flux), axis=0)
    f0 = (2.99792e-13 * 3631) / ((wave / 1e4) ** 2)  ##f0 is the flux of an object with conventional magnitude of zero.

    filters = [u, g, r, i]
    mags = []
    for c in filters:
        ind1 = np.where(wave == fil(c)[0])[0][0]
        ind2 = np.where(wave == fil(c)[1])[0][0]
        wave = wave[ind1:ind2]
        flux = flux[ind1:ind2]
        tck = interpolate.splrep(c['col1'], c['col2'], s=0)  # take the SDSS filter data and interpolate
        transmissions = interpolate.splev(wave, tck, der=0)  # spits out the interpolated data

        f0 = f0[ind1:ind2]
        fxt = transmissions * flux
        fxt_0 = transmissions * f0
        d_lambda = []
        for lam in np.arange(len(fxt)):
            try:
                d_lambda.append(wave[lam + 1] - wave[lam])
            except IndexError:
                d_lambda.append(wave[lam] - wave[lam - 1])

        numer = np.sum(fxt * d_lambda)
        denom = np.sum(transmissions * d_lambda)
        f = numer / denom

        f_zero = np.sum(fxt_0 * d_lambda) / np.sum(transmissions * d_lambda)

        m = (-2.5 / np.log(10)) * (np.arcsinh((f / f_zero) / (2 * fil(c)[2])) + np.log(fil(c)[2]))
        # this is the SDSS asinh magnitude

        mags.append(m)# list of all 4 magnitudes (u,g,r,i) ..Index for specific one

        wave = (10 ** data['loglam'] / (1 + a_rs)) * (1 + rs)
        flux = data['flux'] * 1e-17
        wave = np.concatenate((UV_wave, wave), axis=0)
        flux = np.concatenate((UV_normalized_flux, flux), axis=0)
        f0 = (2.99792e-13 * 3631) / ((wave / 1e4) ** 2)

    print mags
    u_g.append(mags[0] - mags[1])
    g_r.append(mags[1] - mags[2])
    r_i.append(mags[2] - mags[3])

