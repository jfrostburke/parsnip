import time
import sys
import os
import argparse

try:
    import utils
except (SyntaxError, ImportError) as e:
    parsnip_error = '\n\033[1m\033[91mERROR: \033[0m' \
                    'Need to be in parsnip conda environment:\n' \
                    '\tconda activate parsnip\n' \
                    'If conda environment not installed yet, run:\n' \
                    '\tconda env create -f parsnip.yaml\n'

    print(parsnip_error)
    sys.exit()


if __name__ == "__main__":

    start = time.time()

    parser = argparse.ArgumentParser(
        description='PARSNIP: Pipeline for the Analysis and Reduction of SN '
        'Images and Photometry.')
    parser.add_argument('-n', '--name', type=str, help='name of SN to reduce', required=True)
    parser.add_argument('-e', '--epoch', type=str, help='epoch range to reduce')
    parser.add_argument('-f', '--filt', type=str, help='filters to reduce')
    parser.add_argument('--datamax', type=int, help='max value to pick during psf generation')
    parser.add_argument('--nstars', type=int, help='number of stars to pick for psf generation')
    parser.add_argument('--filetype', type=int, choices=(1, 3, 4),
                        default=1, const=1, nargs='?',
                        help='filetype to reduce: 1 (default), 3 (subtracted), 4 (template)')
    parser.add_argument('-s', '--stage', type=str, choices=('cosmic', 'psf', 'psfmag'),
                        help='stage to run')
    parser.add_argument('--show', action='store_true',
                        help='displays plots for visualization and analysis')
    args = parser.parse_args()

    db_name = os.getenv('pipeline_db_name')
    db_user = os.getenv('pipeline_db_user')
    db_host = os.getenv('pipeline_db_host')
    db_pwd = os.getenv('pipeline_db_pwd')

    db_dic = utils.get_db_session(db_name, db_user, db_pwd, db_host)
    images = utils.get_images(**db_dic, epoch=args.epoch, name=args.name, 
                              filt=args.filt, filetype=args.filetype)

    if args.stage == 'cosmic':
        from cosmic_ray import get_cosmic_rays
        print('Detecting cosmic rays . . .')
        print()
        for image in images:
            get_cosmic_rays(image)
        print()

    if args.stage == 'psf':
        from psf import make_psf  # saves >1 second of import on startup
        print('Generating psfs for images . . .')
        print()
        for image in images:
            make_psf(image, show=args.show, datamax=args.datamax,
                     nstars=args.nstars)
        print()

    if args.stage == 'psfmag':
        from photometry import psf_photometry
        print('Running psf photometry . . .')
        print()
        for image in images:
            psf_photometry(image)
        print()

    db_dic['db_session'].close()
    end = time.time()
    print(f'Elapsed time (s): {end-start:.3f}')
