import warnings
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io.fits import getheader
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning


def catalog_crossmatch(filepath, filename, stars, catalog, filt_name, max_sep, max_mag):

    full_filepath = filepath + filename
    hdr = getheader(full_filepath)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    wcs = WCS(hdr)

    coords_stars_xy = [[star['x'], star['y']] for star in stars]
    coords_stars = wcs.wcs_pix2world(coords_stars_xy, 1)
    coords_stars = SkyCoord(coords_stars, unit=u.degree)

    catalog_stars = ascii.read(catalog,
                               names=['ra', 'dec', 'id', 'B', 'Berr', 'V', 'Verr',
                                      'g', 'gerr', 'r', 'rerr', 'i', 'ierr'
                               ]
    )
    coords_catalog_stars = SkyCoord(catalog_stars['ra'], catalog_stars['dec'], unit=u.degree)

    idx, sep_angle, _ = coords_stars.match_to_catalog_sky(coords_catalog_stars)
    max_sep = max_sep*u.arcsec
    sep_constraint = sep_angle < max_sep
    if filt_name in ['B', 'V', 'g', 'r', 'i']:
        mag_constraint = catalog_stars[idx][filt_name] < max_mag
    else:
        print('\tSorry, that filter does not exist in this catalog')
        mag_constraint = True
    mask = sep_constraint & mag_constraint
    print(f'\t\tLimiting to {len(np.argwhere(mask))}/{len(stars)} catalog stars')
    star_matches = stars[mask]
    catalog_matches = catalog_stars[idx[mask]]

    if type(star_matches) == Table:
        star_matches = Table(np.array([list(star) for star in star_matches]),
                               names=list(star_matches.keys())
                         )
    else:
        star_matches = Table(np.array([list(star.values()) for star in star_matches]),
                               names=list(star_matches[0].keys() )
                         )

    return star_matches, catalog_matches

