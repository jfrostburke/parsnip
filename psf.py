import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch, simple_norm
from photutils import EPSFBuilder
from photutils.psf import extract_stars
from astroscrappy import detect_cosmics
from scipy import spatial
from scipy.stats import sigmaclip


def make_psf(image, show, datamax, nstars):

    t_start = time.time()

    if not nstars:
        nstars = 12

    print(f'Working on {image.filename}')
    full_filepath = image.filepath + image.filename
    hdul = fits.open(full_filepath)
    data = hdul[0].data

    banzai = hdul[1].data
    banzai_filtered = _filter_banzai(banzai, datamax)
    banzai_coords = np.array([coord for coord in zip(banzai_filtered['x'],
                                                     banzai_filtered['y'])])

    print('\tIgnoring cosmic rays . . .')
    (crmask, data) = detect_cosmics(data)
    cr_coords = np.argwhere(crmask)
    print(f'\t\t{len(cr_coords)} cosmic rays detected')
    # TODO: Think about structure more
    #   CR discovery is by far the slowest part, separate step and save?
    #   Best way to save PSF/python objects?

    size = 25
    tree = spatial.KDTree(cr_coords)
    distances = tree.query(banzai_coords)[0]
    sep_constraint = distances > size
    print(f'\t\t{len(banzai_coords)} stars detected')
    psf_stars = Table(banzai_coords[sep_constraint][:nstars],
                      names=['x', 'y'])
    banzai_coords = Table(banzai_coords[sep_constraint],
                          names=['x', 'y'])

    nddata = NDData(data=data)
    stars = extract_stars(nddata, psf_stars, size=size)
    # extract_stars is fast on banzai_coords, will need them for psfmag stage
    # refactor this to have a find_stars function to use in psfmag stage?

    print('\tBuilding psf . . .')
    epsf_builder = EPSFBuilder(oversampling=2, maxiters=3,
                               progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)

    hdul.close()

    t_end = time.time()
    print(f'Time to generate psf (s): {int(t_end-t_start)}')

    if show:
        check_psf(data, image.filename, epsf.data, banzai_coords, psf_stars, size)

    return epsf, fitted_stars


def check_psf(data_img, filename, data_psf, stars_all, stars_psf, size):

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

    if datamax:
        mask = (banzai['PEAK'] < datamax)
        banzai = banzai[mask]

    sigma = 3.0
    mask = np.array([True]*len(banzai))
    for param in ['ELLIPTICITY', 'BACKGROUND', 'FWHM']:
        filtered_vals, _, _ = sigmaclip(banzai[param], low=sigma, high=sigma)
        mask &= (banzai[param] >= min(filtered_vals))
        mask &= (banzai[param] <= max(filtered_vals))
    # Sort by S/N?
    banzai_filtered = banzai[mask]

    return banzai_filtered
