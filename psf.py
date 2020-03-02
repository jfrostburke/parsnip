from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch, simple_norm
from photutils import find_peaks
from photutils.psf import extract_stars
from astropy.nddata import NDData
from photutils import EPSFBuilder
from astroscrappy import detect_cosmics
import time
from astropy.coordinates import SkyCoord, match_coordinates_3d
import astropy.units as u


def make_psf(images, show):

    plt.ion()

    for image in images:

        t_start = time.time()
        print(f'Working on {image.filename}')
        full_filepath = image.filepath + image.filename
        hdul = fits.open(full_filepath)
        data = hdul[0].data
        banzai = hdul[1].data

        print('\tDetecting cosmic rays . . .')
        (crmask, data) = detect_cosmics(data)
        cr_coords = np.argwhere(crmask)
        # TODO: fix this jank
        #   double for loop faster than this? just want to not pick CRs
        #   separate CR stage? Add tables to image itself?
        #   KDTree from https://stackoverflow.com/questions/39107896/efficiently-finding-the-closest-coordinate-pair-from-a-set-in-python
        # TODO: add filtering like lscloop.py --banzai
        # also add nstars, figure out how to save/pass objects, etc.
        # TODO: check PEP8 throughout new code!
        cr_cat = SkyCoord(x=cr_coords[:,0], y=cr_coords[:,1],
            z=[1e6]*cr_coords.shape[0], unit='kpc', representation_type='cartesian')
        banzai_cat = SkyCoord(x=banzai['x'], y=banzai['y'],
            z=[1e6]*banzai.shape[0], unit='kpc', representation_type='cartesian')

        idx, d2d, d3d = match_coordinates_3d(banzai_cat, cr_cat)
        min_sep = 50*u.kpc # at least 50 px from CR
        sep_constraint = d3d > min_sep
        stars_tbl = Table(banzai[sep_constraint][:12])
        size = 25

        """
        # copied from https://photutils.readthedocs.io/en/stable/epsf.html
        print('\tDetecting stars . . .')
        peaks_tbl = find_peaks(data, threshold=500.)
        hsize = (size - 1) / 2
        x = peaks_tbl['x_peak']  
        y = peaks_tbl['y_peak']  
        # Exclude stars close to the image edge
        mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) &
            (y > hsize) & (y < (data.shape[0] -1 - hsize)))
        stars_tbl = Table()
        stars_tbl['x'] = x[mask]  
        stars_tbl['y'] = y[mask]
        stars_tbl = stars_tbl[:12]
        """
        nddata = NDData(data=data)
        stars = extract_stars(nddata, stars_tbl, size=size)
        print(f'\t\t{len(stars)} stars detected')

        print('\tBuilding psf . . .')
        epsf_builder = EPSFBuilder(oversampling=2, maxiters=3,
            progress_bar=False)
        epsf, fitted_stars = epsf_builder(stars)

        hdul.close()
    
        t_end = time.time()
        print(f'Time to generate psf (s): {int(t_end-t_start)}')

        if show:
            check_psf(data, epsf.data, stars_tbl, size)

    return


def check_psf(data_img, data_psf, stars_tbl, size):

    # TODO: pretty this up

    plt.clf()
    gs = gridspec.GridSpec(3, 2)
    ax_img = plt.subplot(gs[:2, :2])
    ax_psf_2d = plt.subplot(gs[2, 0])
    ax_psf_3d = plt.subplot(gs[2, 1], projection='3d')

    (vmin, vmax) = ZScaleInterval().get_limits(data_img)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    im = ax_img.imshow(data_img, norm=norm, origin='lower', cmap='gray')
        
    ax_img.scatter(stars_tbl['x'], stars_tbl['y'], s=size, edgecolor='red',
        facecolor=(0, 0, 0, 0))

    norm = simple_norm(data_psf, 'linear', percent=99.)
    im_psf = ax_psf_2d.imshow(data_psf, norm=norm, origin='lower', cmap='viridis')

    x = np.arange(0, data_psf.shape[0])
    y = np.arange(0, data_psf.shape[1])

    X, Y = np.meshgrid(x, y)
    ax_psf_3d.plot_wireframe(X, Y, data_psf, rcount=size, ccount=size)

    response = input('good psf [[y]/n] or [b] bad quality image?')

