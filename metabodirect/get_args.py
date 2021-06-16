# get_args.py

import argparse
import os
import metabodirect


# --------------------------------------------------
def get_args():
    """Get command-line arguments."""

    parser = argparse.ArgumentParser(
        description='Program for running all the MetaboDirect analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('data_file',
                        help='Name of the file with the Direct Injection MS data in .csv format',
                        metavar='DATA',
                        type=str)

    parser.add_argument('metadata_file',
                        help='Name of the file with the sample information (metadata) in .csv format',
                        metavar='METADATA',
                        type=str)

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory',
                        metavar='OUTDIR',
                        type=str,
                        default='MetaboDirect_output')

    parser.add_argument('-m',
                        '--mass_filter',
                        help='Range to filter m/z data (min_mz, max_mz). The pipeline will not filter '
                             'm/z values by default',
                        metavar='INT',
                        type=float,
                        nargs=2)

    parser.add_argument('-g',
                        '--group',
                        help='Grouping variables for coloring and faceting figures (Max 2)',
                        metavar='STR',
                        type=str,
                        nargs='+',
                        required=True)

    parser.add_argument('-f',
                        '--filter_by',
                        help='Filter samples based on metadata. First enter the name of the feature,'
                             'followed by the values associated with the samples you want to keep in the analysis.'
                             '(Example -f Habitat Bog,Palsa)',
                        metavar='STR',
                        type=str,
                        nargs=2)

    parser.add_argument('-b',
                        '--biochem_key',
                        help='File with the biochemical key to use for the transformation network',
                        metavar='STR',
                        type=str,
                        default='Default key')

    parser.add_argument('-v',
                        '--version',
                        action='version',
                        version='%(prog)s {version}'.format(version=metabodirect.__version__))

    parser.add_argument('-k',
                        '--kegg_annotation',
                        help='Set this option to perform annotation of the molecular formulas using'
                             'the KEGG database',
                        action='store_true',
                        default=False)

    norm = parser.add_argument_group('Normalization methods',
                                     'Options to define how data normalization will be carried out')

    norm.add_argument('-n',
                      '--norm_method',
                      help="Available methods to normalize data are: 'mean', 'median', 'zscore', "
                           "'sum', 'max', 'minmax', 'binary', 'none'",
                      metavar='STR',
                      type=str,
                      choices=['mean', 'median', 'zscore', 'sum', 'max', 'minmax', 'binary', 'none'],
                      default='max')

    norm.add_argument('--norm_subset',
                      help='Subset of the data to use for normalization purpouses. '
                           'Available subset methods: ALL, LOS, PPP'
                           'LOS uses peaks in the top L '
                           'order statistics, PPP uses peaks having a minimum percentage of observed values.',
                      metavar='STR',
                      type=str,
                      choices=['ALL', 'LOS', 'PPP'],
                      default='ALL'
                      )

    norm.add_argument('--subset_parameter',
                      help='If using a sample subset for nomalization, this parameter defines the subsample of peaks '
                           'that will be used for normalization.'
                           'If not defined, the default values will be 0.3 for LOS and 0.5 for PPP',
                      metavar='FLOAT',
                      type=float
                      )

    norm.add_argument('--log_transform',
                      help='Set this option to log transform the data. (Program will fail if there are peaks with '
                           'intensities of 0. Consider tranforming this values into 1 if log transformation is desired',
                      action='store_true',
                      default=False
                      )

    transf = parser.add_argument_group('Transformation network options',
                                       'Options to control wheter transformations will be calculated and if '
                                       'networks will be constructed')

    transf.add_argument('-t',
                        '--calculate_transformations',
                        help='Set this option to calculate transformations based on biochemical key',
                        action='store_true',
                        default=False)

    transf.add_argument('-c',
                        '--create_networks',
                        help='Set this option to build transformation networks based on transfomations calculated'
                             'with the biochemical key (this options turns -t automatically)',
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    # Check that command line arguments are specified properly

    if not os.path.isfile(args.data_file):
        parser.error(f'File {args.data_file} does not exist')

    if not os.path.isfile(args.metadata_file):
        parser.error(f'File {args.metadata_file} does not exist')

    if len(args.group) > 2:
        parser.error('Incorrect number of variables for grouping. Please enter exactly two values')

    if args.mass_filter:
        if args.mass_filter[1] < args.mass_filter[0]:
            parser.error('Incorrect order of values for filtering. Please enter min_mz first (e.g. -f 200 900)')

    if args.norm_subset != 'ALL':
        if args.subset_parameter > 1:
            parser.error(f'Incorrect subset parameter. Subset parameter must be between 0 and 1 (proportion).')
        if not args.subset_parameter:
            args.subset_parameter = 0.3 if args.norm_subset == 'LOS' else 0.5

    if args.create_networks:
        args.calculate_transformations = True

    return args
