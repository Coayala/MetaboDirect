#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2021-04-10
Purpose: Program to run MetaboDirect scripts
"""

import argparse


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

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Main body to run all MetaboDirect scripts"""

    args = get_args()
    data_file = args.data_file
    metadata_file = args.sampleinfo_file

    print(f'Data file is "{data_file}"')
    print(f'Metadata file is "{metadata_file}"')


# --------------------------------------------------
if __name__ == '__main__':
    main()
