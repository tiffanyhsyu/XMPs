
from astropy.table import Table
u=Table.read('u.dat.txt', format='ascii')#, delimiter=' ', guess=False) 
g=Table.read('g.dat.txt', format='ascii')#, delimiter=' ', guess=False) 
r=Table.read('r.dat.txt', format='ascii')#, delimiter=' ', guess=False) 
i=Table.read('i.dat.txt', format='ascii')#, delimiter=' ', guess=False) 
z=Table.read('z.dat.txt', format='ascii')#, delimiter=' ', guess=False) 

# Almeida
#spec-0596-52370-0581
#spec-0964-52646-0351
#random spec
#spec-0596-52370-0589
hdu=fits.open('i_zw_18.fits') #spec-7027-56448-0036
header=hdu[1].header
data=hdu[1].data
plt.plot(10**data['loglam'],data['flux'])
plt.show()
max_flux_ind=np.where(data['flux']== np.max(data['flux']))
a_rs=(10**data['loglam'][max_flux_ind[0][0]]/6563)-1  ##calculates actual redshift by assuming H alpha as strongest line
print a_rs
#print a_rs
z_list=np.linspace(0,0.25,100)
u_g=[]
g_r=[]
r_i=[]

for rs in z_list:
    
    wave=(10**data['loglam']/(1+a_rs))*(1+rs)  


    #print data['sky']

    flux=data['flux']*1e-17

    const=np.median(flux[0:50])/((wave[0])**(-2.3))
    #print const
    UV_wave=np.arange(2000,wave[0],0.5)
    UV_normalized_flux=(const*UV_wave**-2.3)
    wave=np.concatenate((UV_wave,wave),axis=0)
    flux=np.concatenate((UV_normalized_flux,flux),axis=0)
    f0=(2.99792e-13*3631)/((wave/1e4)**2)




    def fil(filt):
        if filt['col1'][0]==2980:
            val1=min(wave,key=lambda x:abs(x-2980))  # locating w_0 in wave
            val2=min(wave,key=lambda x:abs(x-4130))
            b=1.4e-10
        if filt['col1'][0]==3630:
            val1=min(wave,key=lambda x:abs(x-3630))  # locating w_0 in wave
            val2=min(wave,key=lambda x:abs(x-5830))
            b=0.9e-10
        if filt['col1'][0]==5380:
            val1=min(wave,key=lambda x:abs(x-5380))  # locating w_0 in wave
            val2=min(wave,key=lambda x:abs(x-7230))
            b=1.2e-10
        if filt['col1'][0]==6430:
            val1=min(wave,key=lambda x:abs(x-6430))  # locating w_0 in wave
            val2=min(wave,key=lambda x:abs(x-8630))
            b=1.8e-10
        if filt['col1'][0]==7730:
            val1=min(wave,key=lambda x:abs(x-7730))  # locating w_0 in wave
            val2=min(wave,key=lambda x:abs(x-11230))
            b=7.4e-10
        return val1,val2,b

    filters=[u,g,r,i]
    mags=[]
    for c in filters:
        ind1=np.where(wave==fil(c)[0])[0][0]    
        ind2=np.where(wave==fil(c)[1])[0][0]
        wave=wave[ind1:ind2]
        tck = interpolate.splrep(c['col1'], c['col2'],s=0) #take the SDSS filter data
        transmissions = interpolate.splev(wave,tck, der=0) #spits out the interpolated data

        flux=flux[ind1:ind2]

        #plt.plot(wave,transmissions*1e-15)
        f0=f0[ind1:ind2]
        fxt=transmissions*flux
        fxt_0=transmissions*f0
        d_lambda=[]
        for lam in np.arange(len(fxt)):
            try:
                d_lambda.append(wave[lam+1]-wave[lam])
            except IndexError:
                d_lambda.append(wave[lam]-wave[lam-1])

        numer=np.sum(fxt*d_lambda)
        denom=np.sum(transmissions*d_lambda)
        f=numer/denom

        f_zero=np.sum(fxt_0*d_lambda)/np.sum(transmissions*d_lambda)
        #print (f/f_zero)*1e9

        #print f 
        m = (-2.5/np.log(10)) * (np.arcsinh((f/f_zero)/(2*fil(c)[2])) + np.log(fil(c)[2]))
        #print m
        m2 = 22.5 - 2.5 *np.log10((f/f_zero)*1e9)
        #print m2
        wave=(10**data['loglam']/(1+a_rs))*(1+rs)
        flux=data['flux']*1e-17
        wave=np.concatenate((UV_wave,wave),axis=0)
        flux=np.concatenate((UV_normalized_flux,flux),axis=0)
    #     wave=np.concatenate((UV_wave,wave),axis=0)
    #     flux=np.concatenate((UV_normalized_flux,flux),axis=0)
        f0=(2.99792e-13*3631)/((wave/1e4)**2)
    #         print f
        
        mags.append(m)
    #print mags
    # print mags
    # print hdu[2].data['CMODELMAG'][0]
    # #print hdu[2].data['CMODELMAG'][0][:4]-mags
    # print mags[0]-mags[1]
    # print mags[1]-mags[2]
    # print mags[2]-mags[3]
    # # #plt.scatter(mags[0]-mags[1],mags[2]-mags[3])
    u_g.append(mags[0]-mags[1])
    g_r.append(mags[1]-mags[2])
    r_i.append(mags[2]-mags[3])
cm=plt.get_cmap('jet')
plt.scatter(u_g,np.linspace(0,0.25,len(u_g)),s=3,lw=0,c=z_list,cmap=cm)#label='SDSS data')
plt.colorbar()
#plt.title('Non_XMP')
plt.xlabel('u-g')
plt.ylabel('z')
 

# scale=1e-14
# plt.plot(u['col1'],u['col2']*scale)
# plt.plot(g['col1'],g['col2']*scale)
# plt.plot(r['col1'],r['col2']*scale)
# plt.plot(i['col1'],i['col2']*scale)
# plt.plot(z['col1'],z['col2']*scale)
# plt.plot(wave,flux)
# #plt.plot(UV_wave,UV_normalized_flux*0.9*np.mean(flux[:20]),'b')
# #plt.xlim(4500,7000)
plt.show()
