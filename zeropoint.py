import warnings
import matplotlib.pyplot as plt
import numpy as np
from sqlalchemy import Table as SQLTable
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io.fits import getheader
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from utils import load_pickle


def get_catalog(db_session, metadata, base, image):

    # TODO: generalize to using astroquery? just returning apass to start

    class Targets(base):
        __table__ = SQLTable('targets', metadata, autoload=True)

    catalog_filename = db_session.query(Targets).filter(
        (Targets.id == image.targetid)
    ).first().apass_cat
    catalog_filepath = '/supernova/github/lcogtsnpipe/trunk/src/lsc/standard/cat/apass/'
    catalog = catalog_filepath + catalog_filename

    return catalog


def make_zeropoint(image, catalog, show):

    if image.filter not in ['B', 'V', 'gp', 'rp', 'ip']:
        print(f'Skipping {image.filename} (cannot calibrate filter)')

    print(f'Working on {image.filename}')

    # TODO: add color term
    #   will need to update this to nightly zps to get batches of images

    # TODO: put this into utils and replace in code
    full_filepath = image.filepath + image.filename
    hdr = getheader(full_filepath)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    wcs = WCS(hdr)

    metadata = load_pickle(image.filename)
    psfmag_stars = np.array(metadata['psfmags'])
    coords_psfmag_stars_xy = [[star['x'], star['y']] for star in psfmag_stars]
    coords_psfmag_stars = wcs.wcs_pix2world(coords_psfmag_stars_xy, 1)
    coords_psfmag_stars = SkyCoord(coords_psfmag_stars, unit=u.degree)

    catalog_stars = ascii.read(catalog,
                               names=['ra', 'dec', 'id', 'B', 'Berr', 'V', 'Verr',
                                      'g', 'gerr', 'r', 'rerr', 'i', 'ierr'
                               ]
    )
    coords_catalog_stars = SkyCoord(catalog_stars['ra'], catalog_stars['dec'], unit=u.degree)

    idx, sep_angle, _ = coords_psfmag_stars.match_to_catalog_sky(coords_catalog_stars)
    max_sep = 2.0*u.arcsec
    max_mag = 25
    filt_name = image.filter[0]
    sep_constraint = sep_angle < max_sep
    mag_constraint = catalog_stars[idx][filt_name] < max_mag
    mask = sep_constraint & mag_constraint
    print(f'\tLimiting to {len(np.argwhere(mask))}/{len(psfmag_stars)} with known mag')
    psfmag_matches = psfmag_stars[mask]
    catalog_matches = catalog_stars[idx[mask]]

    psfmag_matches = Table(np.array([list(star.values()) for star in psfmag_matches]),
                           names=list(psfmag_matches[0].keys())
                     )

    zeropoints = catalog_matches[filt_name] - psfmag_matches['psfmag']
    dzeropoints = np.sqrt(catalog_matches[f'{filt_name}err']**2
                        + psfmag_matches['dpsfmag']**2)
    zeropoint, dzeropoint = _weighted_avg_and_std(zeropoints, 1/(dzeropoints**2))

    # TODO: make real mag stage
    sn = psfmag_stars[-1]
    sn_mag = sn['psfmag'] + zeropoint
    dsn_mag = np.sqrt(sn['dpsfmag']**2 + dzeropoint**2)
    print(f'\t{sn_mag:.3f} ± {dsn_mag:.3f}')

    if show:
        plot_zeropoint(zeropoints, dzeropoints, zeropoint)

    print()


def plot_zeropoint(zeropoints, dzeropoints, zeropoint):
    
    plt.ion()
    plt.clf()

    x = np.arange(len(zeropoints))
    plt.errorbar(x=x, y=zeropoints, yerr=dzeropoints, marker='o', linestyle='none')
    plt.axhline(y=zeropoint, color='r', lw=2)
    plt.xlabel('Star ID')
    plt.ylabel('Catalog - instrumental mag')
    response = input('Enter to continue')


def _weighted_avg_and_std(values, weights):

    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)

    return average, np.sqrt(variance)

