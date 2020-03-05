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
import warnings
import matplotlib.pyplot as plt

from utils import load_pickle


def psf_photometry(image):

    start = time.time()

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

    ra = hdr['CAT-RA']
    dec = hdr['CAT-DEC']
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    (ra0, dec0) = (c.ra.degree, c.dec.degree)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    wcs = WCS(hdr)
    x, y = wcs.wcs_world2pix(ra0, dec0, 1)
    init_guesses = Table(np.array([x, y]), names=['x_0', 'y_0'])

    print('\tExtracting supernova . . .')
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=bkg,
                                    psf_model=epsf,
                                    fitter=fitter,
                                    fitshape=fitshape,
                                    aperture_radius=fitshape)
    results = photometry(image=image_data, init_guesses=init_guesses)
    end = time.time()
    print(f'Time to extract supernova (s): {end-start:.3f}')

    # TODO: add show option for plotting
    #   good but hella slow rn, cut image to just SN region to optimize?
    #   save results, think about how to save psfmags for stars
    #   run on all objects? unsure about zcat stage

    plt.ion()
    plt.clf()

    plt.subplot(1, 2, 1)
    x = int(results['x_fit'][0])
    y = int(results['y_fit'][0])
    data_display = image_data[y-fitshape:y+fitshape, x-fitshape:x+fitshape]
    vmin = np.percentile(data_display, 5)
    vmax = np.percentile(data_display, 95)
    plt.imshow(data_display, vmin=vmin, vmax=vmax,
               cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower')
    plt.title('Data')

    residual_image = photometry.get_residual_image()
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.subplot(1, 2, 2)
    plt.imshow(residual_image[y-fitshape:y+fitshape, x-fitshape:x+fitshape],
               vmin=vmin, vmax=vmax,
               cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower')
    plt.title('Residuals')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.show()

    response = input('Does this look good?')
    print()
