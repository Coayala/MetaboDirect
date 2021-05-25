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
    """Main body to run all MetaboDirect scripts"""

    args = get_args()

    print('========================\nWelcome to MetaboDirect\n========================\n')
    print('\n------------------------\nNormalization methods testing\n------------------------\n')
    print('Analysis starting on {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")))
    print('Results will be saved in directory: {}\n'.format(os.path.abspath(os.getcwd())))

    print(f'Data file is {args.data_file}')
    print(f'Metadata file is {args.metadata_file}')

    spans_score_script = write_r_script('SPANS_score_template.R', args.data_file, args.metadata_file,
                                        args.group, args.log_transform)

    print(f'Running R script: {spans_score_script}')
    run_r(spans_score_script)
    print(f'Find results and R script in the directory: {os.path.abspath(os.getcwd())}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
