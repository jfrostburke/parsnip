import sys
import os
import argparse

try:
    import utils
except:
    alch_error = '\n\033[1m\033[91mERROR: \033[0m' \
                 'Need to be in parsnip conda environment:\n' \
                 '\tconda activate parsnip\n' \
                 'If conda environment not installed yet, run:\n' \
                 '\tconda env create -f parsnip.yaml\n'

    print(alch_error)
    sys.exit()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='PARSNIP: Pipeline for the Analysis and Reduction of SN '
        'Images and Photometry.')
    parser.add_argument('-n', '--name', type=str, help='name of SN to reduce', required=True)
    parser.add_argument('-e', '--epoch', type=str, help='epoch range to reduce')
    parser.add_argument('-f', '--filt', type=str, help='filters to reduce')
    parser.add_argument('--filetype', type=int, choices=(1,3,4),
        default=1, const=1, nargs='?',
        help='filetype to reduce: 1 (default), 3 (subtracted), 4 (template)')
    args = parser.parse_args()

    db_name = os.getenv('pipeline_db_name')
    db_user = os.getenv('pipeline_db_user')
    db_host = os.getenv('pipeline_db_host')
    db_pwd = os.getenv('pipeline_db_pwd')

    db_dic = utils.get_db_session(db_name, db_user, db_pwd, db_host)
    images = utils.get_images(**db_dic, epoch=args.epoch, name=args.name, 
        filt=args.filt, filetype=args.filetype)

    db_dic['db_session'].close()

