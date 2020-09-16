from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from numpy import *
from astropy.nddata import Cutout2D
from reproject import reproject_interp, reproject_exact
import sys
import time



def master_cutout(path, filename, np=4400):
    '''
    Function to make master cutout image from the master stack. path: directory where the master stack is, filename: master stack image filename, np is the dimension of master cutout image.
    '''


    filename = filename
    filename_nosuffix = filename.split('.')[0]
    wtname = filename_nosuffix + '.wt.fits'

    master = fits.open(path + filename)
    masterwt = fits.open(path + wtname)

    cp1 = master[0].header['CRPIX1']
    cp2 = master[0].header['CRPIX2']


    #header = master[0].header

    #master[0].header['NAXIS1'] = np
    #master[0].header['NAXIS2'] = np

    master[0].header['CTYPE1'] = 'RA---TAN'
    master[0].header['CTYPE2'] = 'DEC--TAN'

    ra = master[0].header['CRVAL2']
    dec = master[0].header['CRVAL1']
    cd11 = master[0].header['CD1_1']
    cd12 = master[0].header['CD1_2']
    cd21 = master[0].header['CD2_1']
    cd22 = master[0].header['CD2_2']

    master[0].header['CRVAL1'] = ra
    master[0].header['CRVAL2'] = dec
    master[0].header['CD1_1'] = cd21
    master[0].header['CD1_2'] = cd22
    master[0].header['CD2_1'] = cd11
    master[0].header['CD2_2'] = cd12


    masterwt[0].header['CTYPE1'] = 'RA---TAN'
    masterwt[0].header['CTYPE2'] = 'DEC--TAN'

    ra = masterwt[0].header['CRVAL2']
    dec = masterwt[0].header['CRVAL1']
    cd11 = masterwt[0].header['CD1_1']
    cd12 = masterwt[0].header['CD1_2']
    cd21 = masterwt[0].header['CD2_1']
    cd22 = masterwt[0].header['CD2_2']

    masterwt[0].header['CRVAL1'] = ra
    masterwt[0].header['CRVAL2'] = dec
    masterwt[0].header['CD1_1'] = cd21
    masterwt[0].header['CD1_2'] = cd22
    masterwt[0].header['CD2_1'] = cd11
    masterwt[0].header['CD2_2'] = cd12


    #header_co = master[0].header
    #header_co['NAXIS1'] = np
    #header_co['NAXIS2'] = np
    print('Making cutout of master_stack......')
    t0 = time.time()
    dat_co = zeros((np,np))

    
    hdr = master[0].header.copy()

    wcs = WCS(hdr)
    center = SkyCoord(ra, dec, unit='deg')

    co = Cutout2D(master[0].data, center, np, wcs=wcs, mode='partial', fill_value=0)
    co_wt = Cutout2D(masterwt[0].data, center, np, wcs=wcs, mode='partial', fill_value=0)

    hdr['CRPIX1'] = co.wcs.wcs.crpix[0]
    hdr['CRPIX2'] = co.wcs.wcs.crpix[1]

    fits.writeto(path + filename_nosuffix + '_cutout.fits',co.data, hdr)
    fits.writeto(path + filename_nosuffix + '_cutout.wt.fits',co_wt.data, hdr)
    t1 = time.time()

    print(filename_nosuffix+'_cutout.fits and weighted image created. Time = ', t1-t0)

    master.close()
    masterwt.close()


def reproject2mastercutout(path, filename, master_cutout):
    '''
    Function to reproject individual stacks to the master stack cutout image so that they all align and have same physical dimensions. path: directory of the individual stack image, \
    filename: filename of individual stack image, master_cutout: filename of the master cutout image.
    '''


    filename = filename
    filename_nosuffix = filename.split('.')[0]
    wtname = filename_nosuffix + '.wt.fits'
    master_cutout = master_cutout

    dat = fits.open(path + filename)
    datwt = fits.open(path + wtname)
    master = fits.open(path + master_cutout)

    gain0 = dat[0].header['gain']
    satr0 = dat[0].header['saturate']

    ra_dat = dat[0].header['CRVAL2']
    dec_dat = dat[0].header['CRVAL1']
    cd11_dat = dat[0].header['CD1_1']
    cd12_dat = dat[0].header['CD1_2']
    cd21_dat = dat[0].header['CD2_1']
    cd22_dat = dat[0].header['CD2_2']

    ra_datwt = datwt[0].header['CRVAL2']
    dec_datwt = datwt[0].header['CRVAL1']
    cd11_datwt = datwt[0].header['CD1_1']
    cd12_datwt = datwt[0].header['CD1_2']
    cd21_datwt = datwt[0].header['CD2_1']
    cd22_datwt = datwt[0].header['CD2_2']

    #ra_master = master[0].header['CRVAL2']
    #dec_master = master[0].header['CRVAL1']
    #cd11_master = master[0].header['CD1_1']
    #cd12_master = master[0].header['CD1_2']
    #cd21_master = master[0].header['CD2_1']
    #cd22_master = master[0].header['CD2_2']

    dat[0].header['CTYPE1'] ='RA---TAN'
    dat[0].header['CTYPE2'] ='DEC--TAN'
    dat[0].header['CRVAL1'] = ra_dat
    dat[0].header['CRVAL2'] = dec_dat
    dat[0].header['CD1_1'] = cd21_dat
    dat[0].header['CD1_2'] = cd22_dat
    dat[0].header['CD2_1'] = cd11_dat
    dat[0].header['CD2_2'] = cd12_dat

    datwt[0].header['CTYPE1'] ='RA---TAN'
    datwt[0].header['CTYPE2'] ='DEC--TAN'
    datwt[0].header['CRVAL1'] = ra_datwt
    datwt[0].header['CRVAL2'] = dec_datwt
    datwt[0].header['CD1_1'] = cd21_datwt
    datwt[0].header['CD1_2'] = cd22_datwt
    datwt[0].header['CD2_1'] = cd11_datwt
    datwt[0].header['CD2_2'] = cd12_datwt

    #master[0].header['CTYPE1'] ='RA---TAN'
    #master[0].header['CTYPE2'] ='DEC--TAN'
    #master[0].header['CRVAL1'] = ra_master
    #master[0].header['CRVAL2'] = dec_master
    #master[0].header['CD1_1'] = cd21_master
    #master[0].header['CD1_2'] = cd22_master
    #master[0].header['CD2_1'] = cd11_master
    #master[0].header['CD2_2'] = cd12_master


    print('Making reprojected image.....')
    t0 = time.time()
    # array1, footprint1 = reproject_interp(dat[0], master[0].header)
    # array2, footprint2 = reproject_interp(datwt[0], master[0].header)

    # array1, footprint1 = reproject_exact(dat[0], master[0].header, parallel=True)
    # array2, footprint2 = reproject_exact(datwt[0], master[0].header, parallel=True)
    array1, footprint1 = reproject_exact(dat[0], master[0].header, parallel=False)
    array2, footprint2 = reproject_exact(datwt[0], master[0].header, parallel=False)

    # name = file0.split('_')[2][0:3]
    # cluster = file0.split('_')[1]

    mh = master[0].header
    newh = mh
    newh['gain'] = gain0
    newh['saturate'] = satr0


    fits.writeto(path + 'cutout_' + filename, array1, newh)
    fits.writeto(path + 'cutout_' + wtname, array2, newh)
    
    print('Reproject cutout image: cutout_{} and cutout_{} created.'.format(filename, wtname))
    t1 = time.time()
    print('Total time: ', (t1-t0))
    print('\n')

    master.close()
    dat.close()
    datwt.close()




def cutout2mastercutout(path, filename, master_cutout):
    '''
    Function to cutout individual stacks using crval1, crval2 of master stack cutout image as center. path: directory of the individual stack image, \
    filename: filename of individual stack image, master_cutout: filename of the master cutout image.
    '''


    filename = filename
    filename_nosuffix = filename.split('.')[0]
    wtname = filename_nosuffix + '.wt.fits'
    master_cutout = master_cutout

    dat = fits.open(path + filename)
    datwt = fits.open(path + wtname)
    master = fits.open(path + master_cutout)

    ra0 = master[0].header['CRVAL1']
    dec0 = master[0].header['CRVAL2']
    np = master[0].header['NAXIS1']

    center_coord = SkyCoord(ra0, dec0, unit='deg')

    dhdr = dat[0].header.copy()
    wthdr = datwt[0].header.copy()

    dhdr['CRVAL1'] = dat[0].header['CRVAL2']
    dhdr['CRVAL2'] = dat[0].header['CRVAL1']
    dhdr['CTYPE1'] = dat[0].header['CTYPE2']
    dhdr['CTYPE2'] = dat[0].header['CTYPE1']
    dhdr['CD1_1'] = dat[0].header['CD2_1']
    dhdr['CD1_2'] = dat[0].header['CD2_2']
    dhdr['CD2_1'] = dat[0].header['CD1_1']
    dhdr['CD2_2'] = dat[0].header['CD1_2']

    wthdr['CRVAL1'] = datwt[0].header['CRVAL2']
    wthdr['CRVAL2'] = datwt[0].header['CRVAL1']
    wthdr['CTYPE1'] = datwt[0].header['CTYPE2']
    wthdr['CTYPE2'] = datwt[0].header['CTYPE1']
    wthdr['CD1_1'] = datwt[0].header['CD2_1']
    wthdr['CD1_2'] = datwt[0].header['CD2_2']
    wthdr['CD2_1'] = datwt[0].header['CD1_1']
    wthdr['CD2_2'] = datwt[0].header['CD1_2']

    wcsd = WCS(dhdr)
    wcswt = WCS(wthdr)

    print('Making individual cutout image.....')
    t0 = time.time()

    cutout = Cutout2D(dat[0].data, center_coord, np, wcs=wcsd, mode='partial', fill_value=0.0)
    cutoutwt = Cutout2D(datwt[0].data, center_coord, np, wcs=wcswt, mode='partial', fill_value=0.0)

    # Update header
    dhdr['NAXIS1'] = np
    dhdr['NAXIS2'] = np
    dhdr['CRPIX1'] = cutout.wcs.to_header()['CRPIX1']
    dhdr['CRPIX2'] = cutout.wcs.to_header()['CRPIX2']

    wthdr['NAXIS1'] = np
    wthdr['NAXIS2'] = np
    wthdr['CRPIX1'] = cutoutwt.wcs.to_header()['CRPIX1']
    wthdr['CRPIX2'] = cutoutwt.wcs.to_header()['CRPIX2']

    fits.writeto(path + 'cutout_' + filename, cutout.data, dhdr)
    fits.writeto(path + 'cutout_' + wtname, cutoutwt.data, wthdr)
    
    print('Cutout image: cutout_{} and cutout_{} created.'.format(filename, wtname))
    t1 = time.time()
    print('Total time: ', (t1-t0))
    print('\n')

    master.close()
    dat.close()
    datwt.close()
