import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.io.fits import getdata
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch, simple_norm
from photutils import EPSFBuilder
from photutils.psf import extract_stars
from scipy import spatial
from scipy.stats import sigmaclip

from utils import load_pickle, create_or_update_pickle


def make_psf(filepath, filename, show, datamax, nstars, ncores):

    # TODO: add stage checking generally to stop errors

    t_start = time.time()

    if not nstars:
        nstars = 12

    print(f'Working on {filename}')
    full_filepath = filepath + filename
    data = getdata(full_filepath, 0)
    banzai = getdata(full_filepath, 1)
    banzai_coords = _filter_banzai(banzai, datamax)
    print(f'\t{len(banzai_coords)} stars detected')

    size = 25
    metadata = load_pickle(filename)
    if 'cr_coords' in metadata.keys():
        cr_coords = metadata['cr_coords']
        print(f'\tIgnoring {len(cr_coords)} cosmic rays . . .')

        tree = spatial.KDTree(cr_coords)
        distances = tree.query(banzai_coords)[0]
        sep_constraint = distances > size
        psf_stars = Table(banzai_coords[sep_constraint][:nstars],
                          names=['x', 'y'])
        banzai_coords = Table(banzai_coords[sep_constraint],
                              names=['x', 'y'])
    else:
        print('\t-s cosmic not run, not ignoring cosmic rays . . .')
        psf_stars = Table(banzai_coords[:nstars],
                          names=['x', 'y'])
        banzai_coords = Table(banzai_coords,
                              names=['x', 'y'])

    nddata = NDData(data=data)
    stars = extract_stars(nddata, psf_stars, size=size)

    print('\tBuilding psf . . .')
    epsf_builder = EPSFBuilder(oversampling=2, maxiters=3,
                               progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)
    # can only pickle python built-ins?
    create_or_update_pickle(filename=filename, key='epsf',
                            val=epsf.data.tolist())
    create_or_update_pickle(filename=filename, key='psf_fitted_stars',
                            val=np.array(banzai_coords).tolist())
    #create_or_update_pickle(filename=filename, key='psf_fitted_stars', val=fitted_stars)

    t_end = time.time()
    print(f'Time to generate psf (s): {t_end-t_start:.2f}')

    if show and ncores == 1:
        check_psf(data, filename, epsf.data, banzai_coords, psf_stars, size)
    elif show and ncores != 1:
        print("Can't display psfs while using multiple cores")

    print()


def check_psf(data_img, filename, data_psf, stars_all, stars_psf, size):

    # TODO: make checkpsf stage work, probably change this to view_psf and add function

    plt.ion()
    plt.clf()

    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
    ax_img = plt.subplot(gs[0, 0])
    ax_psf = plt.subplot(gs[0, 1])

    (vmin, vmax) = ZScaleInterval().get_limits(data_img)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    ax_img.imshow(data_img, norm=norm, origin='lower', cmap='gray')
    ax_img.scatter(stars_all['x'], stars_all['y'], s=size, edgecolor='blue',
                   facecolor=(0, 0, 0, 0))
    ax_img.scatter(stars_psf['x'], stars_psf['y'], s=size, edgecolor='red',
                   facecolor=(0, 0, 0, 0))
    ax_img.set_title('Image: red stars selected for PSF')

    norm = simple_norm(data_psf, 'linear', percent=99.)
    ax_psf.imshow(data_psf, norm=norm, origin='lower', cmap='viridis')
    ax_psf.set_title('PSF')

    plt.gcf().suptitle(f'PSF for {filename}')
    plt.gcf().set_size_inches(10, 5)

    response = input('good psf [[y]/n] or [b] bad quality image?')
    if 'n' in response.lower():
        print('\tUpdating PSF to bad')
        # TODO: add db update when we get there


def _filter_banzai(banzai, datamax):

    # TODO: add edge padding

    if datamax:
        mask = (banzai['PEAK'] < datamax)
        banzai = banzai[mask]

    sigma = 3.0
    mask = np.array([True]*len(banzai))
    for param in ['ELLIPTICITY', 'BACKGROUND', 'FWHM']:
        filtered_vals, _, _ = sigmaclip(banzai[param], low=sigma, high=sigma)
        mask &= (banzai[param] >= min(filtered_vals))
        mask &= (banzai[param] <= max(filtered_vals))
    banzai_filtered = banzai[mask]

    sorted_indices = banzai_filtered['PEAK'].argsort()
    x = banzai_filtered['x'][sorted_indices[::-1]]
    y = banzai_filtered['y'][sorted_indices[::-1]]
    banzai_coords = np.array([coord for coord in zip(x, y)])

    return banzai_coords
