import os
import argparse

import utils


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='PARSNIP: Pipeline for the Analysis and Reduction of SN '
        'Images and Photometry.')
    parser.add_argument('-n', '--name', type=str, help='name of SN to reduce', required=True)
    parser.add_argument('-e', '--epoch', type=str, help='epoch range to reduce')
    args = parser.parse_args()

    db_name = os.getenv('pipeline_db_name')
    db_user = os.getenv('pipeline_db_user')
    db_host = os.getenv('pipeline_db_host')
    db_pwd = os.getenv('pipeline_db_pwd')

    db_dic = utils.get_db_session(db_name, db_user, db_pwd, db_host)
    images = utils.get_images(epoch=args.epoch, name=args.name, **db_dic)

    db_dic['db_session'].close()
