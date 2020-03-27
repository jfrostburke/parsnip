import time
import sys
import os
import argparse
from joblib import Parallel, delayed

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
    parser.add_argument('-s', '--stage', type=str, 
                        choices=('cosmic', 'psf', 'psfmag', 'zcat', 'all'),
                        help='stage to run')
    parser.add_argument('--field', type=str, 
                        choices=('apass'),
                        help='catalog to reference in choosing psf stars ' \
                             'or calibrating zeropoints')
    parser.add_argument('--show', action='store_true',
                        help='displays plots for visualization and analysis')
    parser.add_argument('--plot_mag', action='store_true',
                        help='displays and compares mags during zcat stage')
    parser.add_argument('--ncores', type=int,
                        default=1, const=1, nargs='?',
                        help='number of cores to run on')
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
        Parallel(n_jobs = args.ncores)(
            delayed(get_cosmic_rays)
                   (filepath=image.filepath, filename=image.filename) 
            for image in images
        )
        print()

    if args.stage == 'psf' or args.stage == 'all':
        from psf import make_psf
        print('Generating psfs for images . . .')
        print()
        catalog_filepath = (utils.get_catalog(**db_dic, image=images[0]) 
                            if args.field == 'apass' else '')
        Parallel(n_jobs = args.ncores)(
            delayed(make_psf)
                   (filepath=image.filepath, filename=image.filename, filt_name=image.filter[0],
                    show=args.show, datamax=args.datamax, nstars=args.nstars,
                    ncores=args.ncores, catalog_filepath=catalog_filepath
                   )
            for image in images
        )
        print()

    if args.stage == 'psfmag' or args.stage == 'all':
        from photometry import psf_photometry
        print('Running psf photometry . . .')
        print()
        sn_ra, sn_dec = utils.get_sn_radec(**db_dic, image=images[0])
        Parallel(n_jobs = args.ncores)(
            delayed(psf_photometry)
                   (filepath=image.filepath, filename=image.filename,
                    show=args.show, sn_ra=sn_ra, sn_dec=sn_dec
                   )
            for image in images
        )
        print()

    if args.stage == 'zcat' or args.stage == 'all':
        from zeropoint import make_zeropoint
        import matplotlib.pyplot as plt
        print('Generating zeropoints . . .')
        print()
        fig, ax = plt.subplots()
        catalog = utils.get_catalog(**db_dic, image=images[0])
        for image in images:
            make_zeropoint(image=image, catalog=catalog, show=args.show, 
                           plot_mag=args.plot_mag, ax=ax)
        if args.plot_mag:
            ax.legend()
            plt.show()
            plt.close()
        print()

    db_dic['db_session'].close()
    end = time.time()
    print(f'Elapsed time (s): {end-start:.3f}')
