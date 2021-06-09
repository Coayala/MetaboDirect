#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2021-04-10
Purpose: Program to run MetaboDirect scripts
"""

import argparse
import os
import datetime
import subprocess
import pandas as pd
from shutil import which
from platform import system


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Program for running all the MetaboDirect analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('data_file',
                        help='Name of the file with the DI-MS data in .csv format',
                        metavar='DATA',
                        type=str)

    parser.add_argument('metadata_file',
                        help='Name of the file with the sample information (metadata) in tabular format',
                        metavar='METADATA',
                        type=str)

    parser.add_argument('group',
                        help='Grouping variables to test for normalization significance',
                        metavar='GROUP',
                        type=str)

    parser.add_argument('-f',
                        '--filter_by',
                        help='Filter samples based on metadata. First enter the name of the feature,'
                             'followed by the values associated with the samples you want to keep in the analysis.'
                             '(Example -f Habitat Bog,Palsa)',
                        metavar='STR',
                        type=str,
                        nargs=2)

    parser.add_argument('--log_transform',
                        help='Set this if you plan to log transform your data before normalization',
                        action='store_true',
                        default=False)

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error(f'File {args.data_file} does not exist')

    if not os.path.isfile(args.metadata_file):
        parser.error(f'File {args.metadata_file} does not exist')

    return args


# --------------------------------------------------
def sample_filtering(df, metadata, filter_by):
    """Filter samples based on selected features and values."""

    # Get the variable a values specified for sample filtering
    filter_col = filter_by[0]
    filter_values = filter_by[1].split(sep=',')

    # Saving a new metadata file containing only the samples remaining after filtering
    filt_metadata = pd.DataFrame()
    for i in filter_values:
        filt_metadata = filt_metadata.append(metadata[metadata[filter_col] == i])
    filt_metadata_file = 'filtered_metadata_norm_test.csv'
    filt_metadata.to_csv(filt_metadata_file, index=False)

    # Saving a new input file containing only the samples remaining after filtering
    from_formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]
    col_df = filt_metadata['SampleID'].to_list()
    from_formularity.extend(col_df)
    filt_df = df[from_formularity]
    filt_df_file = 'filtered_input_norm_test.csv'
    filt_df.to_csv(filt_df_file, index=False)

    return filt_df_file, filt_metadata_file


# --------------------------------------------------
def write_r_script(rscript, data_file, metadata_file, group, log):
    """Function to write R scripts based on the provided Rscripts templates"""
    current_dir = os.getcwd().replace('\\', '/')
    data_file = data_file.replace('\\', '/')
    metadata_file = metadata_file.replace('\\', '/')
    metabo_home = os.path.split(os.path.realpath(__file__))[0].replace('\\', '/')
    r_in = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'R_scripts_templates', rscript)
    r_in = open(r_in)
    r_file = os.path.join(current_dir, rscript.replace('_template', ''))
    r_fh = open(r_file, 'w')
    for line in r_in:
        r_fh.write(line.replace('%currentdir%', current_dir)
                   .replace('%metadata%', metadata_file)
                   .replace('%group%', group)
                   .replace('%data%', data_file)
                   .replace('%Metabo_HOME%', metabo_home)
                   .replace('%log%', 'TRUE' if log else 'FALSE')
                   )
    r_fh.close()

    return r_file


# --------------------------------------------------
def run_r(rscript):
    """Function to execute R scripts"""
    rscript_exe = 'Rscript.exe' if system() == 'Windows' else 'Rscript'

    if which(rscript_exe) is not None:
        r_cmd = [rscript_exe, rscript]
        p = subprocess.Popen(
            r_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = p.communicate()
        if out:
            print(out)

        if err:
            print(err)
    else:
        print('Rscript not accesible')


# --------------------------------------------------
def main():
    """Main body to run test normalization analysis"""

    args = get_args()

    print('========================\nWelcome to MetaboDirect\n========================\n')
    print('\n------------------------\nNormalization methods testing\n------------------------\n')
    print('Analysis starting on {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")))
    print('Results will be saved in directory: {}\n'.format(os.path.abspath(os.getcwd())))

    print(f'Data file is {args.data_file}')
    print(f'Metadata file is {args.metadata_file}')

    if args.filter_by:
        df = pd.read_csv(args.data_file)
        metadata = pd.read_csv(args.metadata_file)
        print(f'Filtering samples based on {args.filter_by[0]} = {args.filter_by[1]}')
        df_file, metadata_file = sample_filtering(df, metadata, args.filter_by)
        args.data_file = df_file
        args.metadata_file = metadata_file

    else:
        print('Option -f not detected, all samples will be used')

    spans_score_script = write_r_script('SPANS_score_template.R', args.data_file, args.metadata_file,
                                        args.group, args.log_transform)

    print(f'Running R script: {spans_score_script}')
    run_r(spans_score_script)
    print(f'Find results and R script in the directory: {os.path.abspath(os.getcwd())}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
