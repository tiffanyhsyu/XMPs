from astropy.io import fits
import matplotlib
matplotlib.use('Qt4Agg')
import numpy as np
import matplotlib.pyplot as plt

import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import numpy.ma as ma


Hb=4862.721
oiii_4959=4960.295
oiii_5007=5008.239
Ha = 6564.614
nii_6584=6585.27
Hg=4341.68
Hd=4102.89


def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))


def LF(f,el):
    # give a file and emission line
    hdu=fits.open(f)
    try:
        wave=10**hdu[1].data['loglam']
    except KeyError:
        wave=hdu[1].data['LAMBDA'][0]
    try:
        flux=hdu[1].data['flux']
    except KeyError:
        flux=hdu[1].data['SPEC'][0]

    error=(1/np.sqrt(ma.array(hdu[1].data['ivar'])))
    if el=='oiii_4363':
        el=4364.44
    if el=='Hb':
        el=4862.721
    if el == 'oiii_4959':
        el= 4960.295
    if el == 'oiii_5007':
        el=5008.239
    if el ==  'Ha':
        el= 6564.614
    if el == 'nii_6584':
        el = 6585.27
    if el == 'Hg':
        el =4341.68
    if el == 'Hd':
        el=4102.89

    z=hdu[2].data['Z']
    l=abs(wave-(z+1)*el)
    pos=np.argmin(l)
    pos2=np.argmin(abs(wave-(z+1)*6300))

    peak_flux=flux[pos-8:pos+8]#-np.median(flux)
    peak_error=error[pos-8:pos+8]
    oiii_wave=wave[pos-8:pos+8]
    error_list=error[pos-8:pos+8]

    mean_oiii = sum(oiii_wave * peak_flux) / sum(peak_flux)
    sigma_oiii = np.sqrt(sum(peak_flux * (oiii_wave - mean_oiii)**2) / sum(peak_flux))


    popt,pcov = curve_fit(gaus,oiii_wave,peak_flux,p0=[np.max(peak_flux),mean_oiii,sigma_oiii])
    error_val=np.sqrt(np.diag(pcov))

    perr = np.sqrt(np.diag(pcov))
    oiii_flux=popt[2]*popt[0]*np.sqrt(2*np.pi)
    cont=np.median(flux[pos2-50:pos2+50])
    err_cont=np.median(error[pos2-50:pos2+50])
    EW=[]
    EW_error=[]
    for p in np.arange(len(peak_flux)):
        EW.append(1-(peak_flux[p]/cont))
        EW_error.append(1-(peak_error[p]/err_cont))
    EW=np.sum(EW)
    EW_error=np.sqrt(perr[1]**2+ perr[2]**2)

    plt.plot(oiii_wave,peak_flux,'k',linewidth=2.0,drawstyle='steps-mid')
    plt.plot(oiii_wave,gaus(oiii_wave,*popt),linewidth=4.0,label='fit')
    plt.plot(oiii_wave,error_list,'r',label='flux error',linewidth=1.0,drawstyle='steps-mid')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (10$^{-17}$erg/s/cm$^2$/$\AA$)')
    plt.xlim(oiii_wave[0]-0.5,oiii_wave[-1]+0.5)
    plt.legend(loc='best')
    plt.show()

    print 'EW:',np.abs(EW), 'EW_err:',np.abs(EW_error), 'Line Flux:',oiii_flux
    return

LF('spec-7027-56448-0036.fits','Hb')