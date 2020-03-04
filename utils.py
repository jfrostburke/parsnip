import os
import pickle
from sqlalchemy import create_engine, MetaData, Table
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


def get_images(db_session, metadata, base, epoch, name, filt, filetype):

    class TargetNames(base):
        __table__ = Table('targetnames', metadata, autoload=True)

    class PhotLCO(base):
        __table__ = Table('photlco', metadata, autoload=True)

    targetid = db_session.query(TargetNames).filter(
        (TargetNames.name == name)
    ).first().targetid
    criteria = (PhotLCO.targetid == targetid)

    if epoch:
        (epoch_start, epoch_end) = epoch.split('-') if '-' in epoch else (epoch, epoch)
        criteria &= (PhotLCO.dayobs >= epoch_start)
        criteria &= (PhotLCO.dayobs <= epoch_end)

    if filt:
        if filt == 'sloan':
            sloan_crit = (PhotLCO.filter == 'up')
            for sloan_filt in ['gp', 'rp', 'ip', 'zs']:
                sloan_crit |= (PhotLCO.filter == sloan_filt)
            criteria &= sloan_crit
        elif filt == 'landolt':
            landolt_crit = (PhotLCO.filter == 'U')
            for landolt_filt in ['B', 'V', 'R', 'I']:
                landolt_crit |= (PhotLCO.filter == landolt_filt)
            criteria &= landolt_crit
        else:
            criteria &= (PhotLCO.filter.like(f'%{filt}%'))

    criteria &= (PhotLCO.quality == 127)  # good images
    criteria &= (PhotLCO.filetype == filetype)

    images = db_session.query(PhotLCO).filter(criteria).order_by(PhotLCO.filename)

    # Need width of longest zcat filename to format printout nicely
    zcats = [image.zcat for image in images.distinct(PhotLCO.zcat)]

    if not zcats:
        print()
        print('No images selected')
        print()
    else:
        zcat_width = len(max(zcats, key=len))
        _print_images(images, zcat_width)

    return images


def _print_images(images, zcat_width):
    
    number_of_images = 0
    
    print()
    print('filename\t\t\t\tobject\t\tfilter\tWCS\tPSF\tpsfmag\tapmag\tzcat\t\t\tmag')

    for image in images:

        number_of_images += 1

        filename = image.filename.split('.fits')[0]
        psf = _format_output(image.psf, good_value='psf.fits')
        wcs = _format_output(image.wcs, good_value=0.0)
        psfmag = _format_output(image.psfmag, bad_value=9999)
        apmag = _format_output(image.apmag, bad_value=9999)
        mag = _format_output(image.mag, bad_value=9999)
        zcat = _format_output(image.zcat, pad=zcat_width, good_value='cat')

        print(
            f'{filename}\t\t{image.objname}\t{image.filter}\t'
            f'{wcs}\t{psf}\t{psfmag}\t{apmag}\t'
            f'{zcat}\t{mag}'
        )

    print()
    print(f'Number of images: {number_of_images}')
    print()

    return


def _format_output(param, pad=6, bad_value=None, good_value=None):
    
    bad = '\033[1m\033[91mX\033[0m'.ljust(pad)
    good = '\033[1m\033[92mY\033[0m'.ljust(pad)

    if bad_value:
        output = str(param)[:pad] if int(param) != bad_value else bad
    elif type(good_value) == str:
        if good_value in param:
            output = param.ljust(pad) if 'cat' in good_value else good
        else:
            output = bad
    elif type(good_value) == float:
        output = good if int(param) == good_value else bad
    else:
        output = bad

    return output


def create_or_update_pickle(filename, key, val):

    pickle_filename = get_pickle_filename(filename)
    if os.path.isfile(pickle_filename):
        metadata = pickle.load(open(pickle_filename, 'rb'))
        metadata[key] = val
        pickle.dump(metadata, open(pickle_filename, 'wb'))
    else:
        metadata = {
            key: val,
        }
        pickle.dump(metadata, open(pickle_filename, 'wb'))


def load_pickle(filename):

    pickle_filename = get_pickle_filename(filename)
    if os.path.isfile(pickle_filename):
        metadata = pickle.load(open(pickle_filename, 'rb'))
    else:
        metadata = {}
    return metadata


def get_pickle_filename(filename):

    # TODO: change pickle location to /supernova/data when the time comes
    pickle_filename = './pickles/' + os.path.splitext(filename)[0] + '.pickle'
    return pickle_filename
