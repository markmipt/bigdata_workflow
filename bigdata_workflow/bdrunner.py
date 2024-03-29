from __future__ import division
import argparse
import logging
import pkg_resources
from . import main
import os


def check_input_dir(dir_path):
    if not os.path.isdir(dir_path):
        raise argparse.ArgumentTypeError(
            "readable_dir:{0} is not a valid path".format(dir_path))
    if not os.access(dir_path, os.R_OK):
        raise argparse.ArgumentTypeError(
            "readable_dir:{0} is not a readable dir".format(dir_path))


def run():
    parser = argparse.ArgumentParser(
        description='wrapper for running bigdata workflow',
        epilog='''

  Example usage
  -------------
  $ bdrunner input.raw/input_folder

  Available modes:
  1 - convert RAW to mzML
  2 - canonical search
  3 - variant search
  4 - brute-force search
  5 - variant table generation
  6 - group FDR filtering
  7 - Prosit MS/MS prediction
  8 - DeepLC RT prediction
  9 - final table composition
  -------------
    ''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('indir', help='input directory')
    parser.add_argument(
        '-m',
        '--mode',
        help='mode for workflow. Can be 1-4; 1,2,5; 1; 1-;\
             read help for mode description',
        default='1-')
    parser.add_argument(
        '-ao',
        '--overwrite',
        help='overwrite existing files',
        action='store_true')
    # parser.add_argument(
    #     '-fmods',
    #     '--fixed_mods',
    #     help='fixed modifications.\
    #          in mass1@aminoacid1,mass2@aminoacid2 format')
    parser.add_argument(
        '-ff',
        '--featurefinder',
        help='path to biosaur or dinosaur bin file',
        default='/usr/bin/biosaur2')
    parser.add_argument(
        '-idpy',
        '--identipy',
        help='path to identipy',
        default='/usr/bin/identipy')
    parser.add_argument(
        '-scav',
        '--scavager',
        help='path to scavager',
        default='/usr/bin/scavager')
    parser.add_argument(
        '-wdb',
        '--wdatabase',
        help='path to wild fasta file',
        required=True)
    parser.add_argument(
        '-cdb',
        '--cdatabase',
        help='path to combined peptide fasta file',
        required=True)
    parser.add_argument(
        '-cfg',
        '--config',
        help='path to IdentiPy config file',
        required=True)
    parser.add_argument(
        '-gmp',
        '--genemap',
        help='path to strange genes',
        default='')
    parser.add_argument(
        '-prosit',
        help='path to prosit folder',
        default='')
    parser.add_argument(
        '-prosit_model',
        help='path to prosit fragmentation model',
        default='')
    parser.add_argument(
        '-prosit_irt_model',
        help='path to prosit iRT model',
        default='')
    parser.add_argument(
        '-corm',
        '--cormap',
        help='path to prosit correlation map',
        default='')
    parser.add_argument('-deeplc', help='path to deeplc', default='')
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debugging output')
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%s' % (pkg_resources.require("bigdata_workflow")[0], ))
    args = vars(parser.parse_args())
    logging.basicConfig(
        format='%(levelname)9s: %(asctime)s %(message)s',
        datefmt='[%H:%M:%S]',
        level=[logging.INFO, logging.DEBUG][args['debug']])
    logger = logging.getLogger(__name__)
    logger.debug('Starting with args: %s', args)
    check_input_dir(args['indir'])

    main.process_folder(args)


if __name__ == '__main__':
    run()
