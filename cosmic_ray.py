import numpy as np
from astropy.io.fits import getdata, getheader
from astroscrappy import detect_cosmics

from utils import create_or_update_pickle, get_pickle_filename


def get_cosmic_rays(filepath, filename):

    print(f'Working on {filename}')
    full_filepath = filepath + filename
    data = getdata(full_filepath)
    hdr = getheader(full_filepath)
    gain = hdr['gain']
    readnoise = hdr['rdnoise']
    satlevel = hdr['saturate']

    print('\tDetecting cosmic rays . . .')
    (crmask, _) = detect_cosmics(data, gain=gain, readnoise=readnoise,
                                 satlevel=satlevel)
    cr_coords = np.argwhere(crmask)
    print(f'\t\t{len(cr_coords)} cosmic rays detected')

    print(f'\tSaving cosmic ray coordinates to {get_pickle_filename(filename)}')
    print()
    create_or_update_pickle(filename=filename, key='cr_coords', val=cr_coords)
