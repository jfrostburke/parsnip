import time
import numpy as np
from astropy.table import Table
from photutils.psf import DAOGroup, BasicPSFPhotometry, EPSFModel
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.io.fits import getdata, getheader
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.exceptions import AstropyUserWarning
import warnings
import matplotlib.pyplot as plt

from utils import load_pickle, create_or_update_pickle


def psf_photometry(image, show):

    start = time.time()
    
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    warnings.simplefilter('ignore', category=AstropyUserWarning)

    print(f'Working on {image.filename}')

    full_filepath = image.filepath + image.filename
    image_data = getdata(full_filepath)
    hdr = getheader(full_filepath)
    fwhm = hdr['L1FWHM']/hdr['PIXSCALE']

    metadata = load_pickle(image.filename)
    epsf_data = np.array(metadata['epsf'])
    epsf = EPSFModel(epsf_data, fwhm=fwhm, oversampling=2)
    daogroup = DAOGroup(2.0*fwhm)
    bkg = MMMBackground()
    fitter = LevMarLSQFitter()
    fitshape = 25
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=bkg,
                                    psf_model=epsf,
                                    fitter=fitter,
                                    fitshape=fitshape,
                                    aperture_radius=fitshape)

    psfmags = []

    print('\tExtracting other stars . . .')
    counter = 0
    for star in metadata['psf_fitted_stars']:
        counter += 1
        (x, y) = star
        psfmags = _do_phot(x, y, image_data, fitshape, photometry, psfmags)
        print(f'\t\tStars extracted: {counter}/{len(metadata["psf_fitted_stars"])}', end='\r')
    print()

    print('\tExtracting supernova . . .')
    x, y = _get_sn_xy(image)
    psfmags = _do_phot(x, y, image_data, fitshape, photometry, psfmags)

    create_or_update_pickle(image.filename, key='psfmags', val=psfmags)

    end = time.time()
    print(f'Time to perform photometry (s): {end-start:.3f}')

    if show:
        checkmag(image_data, photometry.get_residual_image(), x, y, fitshape)
    print()

    # TODO: pickle residual image around sn? can't save photometry


def _get_sn_xy(image):

    full_filepath = image.filepath + image.filename
    hdr = getheader(full_filepath)
    ra = hdr['CAT-RA']
    dec = hdr['CAT-DEC']
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    (ra0, dec0) = (c.ra.degree, c.dec.degree)
    wcs = WCS(hdr)
    x, y = wcs.wcs_world2pix(ra0, dec0, 1)
    return float(x), float(y)


def _do_phot(x, y, image_data, fitshape, photometry, psfmags):

    data_fit = image_data[
        int(y)-fitshape:int(y)+fitshape, int(x)-fitshape:int(x)+fitshape
    ]
    init_guesses = Table(np.array([fitshape, fitshape]), names=['x_0', 'y_0'])
    results = photometry(image=data_fit, init_guesses=init_guesses)

    psfflux = results['flux_fit'][0]
    dpsfflux = results['flux_unc'][0]

    psfmag = -2.5*np.log10(psfflux)
    dpsfmag = abs((-2.5*dpsfflux)/(psfflux*np.log(10)))

    star_measurement = {
        'x': x,
        'y': y,
        'psfmag': psfmag,
        'dpsfmag': dpsfmag,
    }

    psfmags.append(star_measurement)
    return psfmags


def checkmag(image_data, residual_data, x, y, fitshape):

    sn_data = image_data[
        int(y)-fitshape:int(y)+fitshape, int(x)-fitshape:int(x)+fitshape
    ]

    plt.ion()
    plt.clf()

    plt.subplot(1, 2, 1)
    vmin = np.percentile(sn_data, 5)
    vmax = np.percentile(sn_data, 95)
    plt.imshow(sn_data, vmin=vmin, vmax=vmax,
               cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower')
    plt.title('Data')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

    vmin = np.percentile(residual_data, 5)
    vmax = np.percentile(residual_data, 95)
    plt.subplot(1, 2, 2)
    plt.imshow(residual_data, vmin=vmin, vmax=vmax,
               cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower')
    plt.title('Residuals')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

    plt.show()

    response = input('Does this look good?')
