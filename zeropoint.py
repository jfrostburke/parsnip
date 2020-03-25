import matplotlib.pyplot as plt
import numpy as np

from utils import load_pickle
from crossmatch import catalog_crossmatch


def make_zeropoint(image, catalog, show):

    if image.filter not in ['B', 'V', 'gp', 'rp', 'ip']:
        print(f'Skipping {image.filename} (cannot calibrate filter)')
        return
        

    print(f'Working on {image.filename}')

    # TODO: add color term
    #   will need to update this to nightly zps to get batches of images
    filt_name = image.filter[0]
    metadata = load_pickle(image.filename)
    psfmag_stars = np.array(metadata['psfmags'])

    psfmag_matches, catalog_matches = catalog_crossmatch(image.filepath, image.filename, psfmag_stars, catalog, 
                                                         filt_name, max_sep=2, max_mag=25)

    zeropoints = catalog_matches[filt_name] - psfmag_matches['psfmag']
    dzeropoints = np.sqrt(catalog_matches[f'{filt_name}err']**2
                        + psfmag_matches['dpsfmag']**2)
    zeropoint, dzeropoint = _weighted_avg_and_std(zeropoints, 1/(dzeropoints**2))

    # TODO: make real mag stage
    sn = psfmag_stars[-1]
    sn_mag = sn['psfmag'] + zeropoint
    dsn_mag = np.sqrt(sn['dpsfmag']**2 + dzeropoint**2)
    print(f'\tparsnip:     {sn_mag:.3f} ± {dsn_mag:.3f}')
    print(f'\tlcogtsnpipe: {image.mag:.3f} ± {image.dmag:.3f}')

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

