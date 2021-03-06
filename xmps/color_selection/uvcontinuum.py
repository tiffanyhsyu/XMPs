from astropy.io import fits
import numpy as np

###simple addition of UV component to spectrum. Modeled as a f propto lambda**-2.3. This value is acceptable within 0.3 of 2.0 
#for star forming galaxies. You should compare resulting u band spectroscopic magnitudes to photometric magnitudes to get a sense of what value to use.
def add_UV_continuum(file,wave_lim):
  hdu=fits.open(file)
  data=hdu[1].data
  wave=10**data['loglam']
  flux=data['flux']*1e-17

  const=np.median(flux[0:50])/((wave[0])**(-2.3))
  #print const
  UV_wave=np.arange(wave_lim,wave[0],0.5)
  UV_normalized_flux=(const*UV_wave**-2.3)
  wave=np.concatenate((UV_wave,wave),axis=0)
  flux=np.concatenate((UV_normalized_flux,flux),axis=0)
  return wave,flux  

