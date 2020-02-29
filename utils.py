import sys

try:
    from sqlalchemy import create_engine, MetaData, Table
except ImportError:
    alch_error = '\n\033[1m\033[91mERROR: \033[0m' \
                 'Need to be in parsnip conda environment:\n' \
                 '\tconda activate parsnip\n' \
                 'If conda environment not installed yet, run:\n' \
                 '\tconda env create -f parsnip.yaml\n'

    print(alch_error)
    sys.exit()
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base


def get_db_session(db_name, db_user, db_pwd, db_host):

    db_engine = create_engine(
        f'mysql://{db_user}:{db_pwd}@{db_host}/{db_name}'
        '?charset=utf8&use_unicode=1')

    db_engine.connect()
    base = declarative_base()
    metadata = MetaData(bind=db_engine)

    db_session = sessionmaker(bind=db_engine)
    session = db_session()

    db_dic = {
        'db_session': session,
        'metadata': metadata,
        'base': base
    }

    return db_dic


def get_images(epoch, name, db_session, metadata, base):

    class TargetNames(base):
        __table__ = Table('targetnames', metadata, autoload=True)

    class PhotLCO(base):
        __table__ = Table('photlco', metadata, autoload=True)

    targetid = db_session.query(TargetNames).filter(
        (TargetNames.name == name)
    ).first().targetid

    (epoch_start, epoch_end) = epoch.split('-')
    images = db_session.query(PhotLCO).filter(
        (PhotLCO.targetid == targetid) &
        (PhotLCO.dateobs >= epoch_start) &
        (PhotLCO.dateobs <= epoch_end)
    )

    print_images(images)

    return images


def print_images(images):

    x_bad = '\033[1m\033[91mX\033[0m'
    y_good = '\033[1m\033[92mY\033[0m'

    print('########################')
    print('#image name\t\t\t\tobject\t\tfilter\tWCS\tPSF\tpsfmag\tapmag\tzcat\t\tmag')
    zcat_width = len(max([x.zcat for x in images], key=len))
    length = 0
    for image in images:

        length += 1

        filename = image.filename.split('.fits')[0]

        # TODO: replace this with helper function so copy less code
        if 'psf.fits' in image.psf:
            psf = y_good
        else:
            psf = x_bad

        if int(image.wcs) == 0:
            wcs = y_good
        else:
            wcs = x_bad

        col_width = 6
        if int(image.psfmag) != 9999:
            psfmag = str(image.psfmag)[:col_width]
        else:
            psfmag = x_bad.ljust(col_width)

        if int(image.apmag) != 9999:
            apmag = str(image.apmag)[:col_width]
        else:
            apmag = x_bad.ljust(col_width)

        if int(image.mag) != 9999:
            mag = str(image.mag)[:col_width]
        else:
            mag = x_bad.ljust(col_width)

        if 'cat' in image.zcat:
            zcat = image.zcat.ljust(zcat_width)
        else:
            zcat = x_bad.ljust(zcat_width)

        print(
            f'{filename}\t\t{image.objname}\t{image.filter}\t'
            f'{wcs}\t{psf}\t{psfmag}\t{apmag}\t'
            f'{zcat}\t{mag}'
        )

    print()
    print(f'Number of images: {length}')
    print()

    return
