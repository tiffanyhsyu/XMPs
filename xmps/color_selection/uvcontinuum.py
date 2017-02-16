from astropy.io import fits

hdu=fits.open(file)
data=hdu[1].data
wave=10**data['loglam']
flux=data['flux']*1e-17

const=np.median(flux[0:50])/((wave[0])**(-2.3))
#print const
UV_wave=np.arange(2000,wave[0],0.5)
UV_normalized_flux=(const*UV_wave**-2.3)
wave=np.concatenate((UV_wave,wave),axis=0)
flux=np.concatenate((UV_normalized_flux,flux),axis=0)
