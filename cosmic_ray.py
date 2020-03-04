import numpy as np
from astropy.io import fits
from astroscrappy import detect_cosmics

from utils import create_or_update_pickle, get_pickle_filename


def get_cosmic_rays(image):

    print(f'Working on {image.filename}')
    full_filepath = image.filepath + image.filename
    with fits.open(full_filepath) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
        gain = hdr['gain']
        readnoise = hdr['rdnoise']
        satlevel = hdr['saturate']

    print('\tDetecting cosmic rays . . .')
    (crmask, _) = detect_cosmics(data, gain=gain, readnoise=readnoise,
                                 satlevel=satlevel)
    cr_coords = np.argwhere(crmask)
    print(f'\t\t{len(cr_coords)} cosmic rays detected')

    print(f'\tSaving cosmic ray coordinates to {get_pickle_filename(image.filename)}')
    print()
    create_or_update_pickle(filename=image.filename, key='cr_coords', val=cr_coords)

    return
