#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2021-04-10
Purpose: Program to run MetaboDirect scripts
"""

import argparse
import os
import datetime
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW
from functions_metabodirect import *


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Program will run all the MetaboDirect analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('data_file',
                        help='Name of the file with the DI-MS data in .csv format',
                        metavar='DATA',
                        type=str)

    parser.add_argument('sampleinfo_file',
                        help='Name of the fie with the sample information (metadata) in tabular format',
                        metavar='SAMPLEINFO',
                        type=str)

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory',
                        metavar='OUTDIR',
                        type=str,
                        default='MetaboDirect_output')

    parser.add_argument('-f',
                        '--filter',
                        help='Range to filter m/z data (min_mz, max_mz). The pipeline will not filter '
                             'm/z values by default',
                        metavar='mz',
                        type=float,
                        nargs=2,
                        default= None)

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error(f'File {args.data_file} does not exist')

    if not os.path.isfile(args.sampleinfo_file):
        parser.error(f'File {args.sampleinfo_file} does not exist')

    #if not len(args.filter) == 2:
    #    parser.error(f'Incorrect number of values for filtering m/z. Please enter exactly two values')

    if args.filter:
        if args.filter[1] < args.filter[0]:
            parser.error(f'Incorrect order of values for filtering. Please enter min_mz first (e.g. -f 200 900')

    return args


# --------------------------------------------------
def make_directories(outdir):
    project_name = outdir
    preprocess_dir = '1_preprocessing_output'
    diagnostics_dir = '2_diagnostics'
    exploratory_dir = '3_exploratory'
    stats_dir = '4_statistics'
    transf_dir = '5_transformations'

    list_dir = [
        preprocess_dir, diagnostics_dir, exploratory_dir, stats_dir, transf_dir
    ]

    list_dir = [os.path.join(project_name, x) for x in list_dir]

    if not os.path.exists(project_name):
        os.makedirs(project_name)

    for d in list_dir:
        if not os.path.exists(d):
            os.makedirs(d)


# --------------------------------------------------
def data_filtering(df_file, filter_values):
    df = pd.read_csv(df_file)
    print(f'Number of m/z in provided file: {df.shape[0]}')

    filt_df, n_excluded = filter_C13(df)
    print(f'Number of masses excluded because of C13 isotope: {n_excluded}')
    print(f'Number of m/z remaining after isotope-filtering: {filt_df.shape[0]}')

    filt_df = filter_mz(filt_df, filter_values[0], filter_values[1]) if filter_values else filt_df
    print(f'Number of m/z remaining after m/z-filtering: {filt_df.shape[0]}')

    filt_df = filter_error_ppm(filt_df, err_range=0.5)
    print(f'Number of m/z after error filtering (0.5 ppm): {filt_df.shape[0]}')

    return filt_df


# --------------------------------------------------
def main():
    """Main body to run all MetaboDirect scripts"""

    args = get_args()

    print('========================\nWelcome to MetaboDirect\n========================\n')
    print('Analysis starting on {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")))
    print('Results will be saved in directory: {}\n'.format(os.path.abspath(args.outdir)))

    print(f'Data file is {args.data_file}')
    print(f'Metadata file is {args.sampleinfo_file}')

    make_directories(args.outdir)
    filt_df = data_filtering(args.data_file, args.filter)


# --------------------------------------------------
if __name__ == '__main__':
    main()
