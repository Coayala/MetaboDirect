#!/usr/bin/env python3
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2023-01-18
Purpose: Program to run MetaboDirect scripts
"""

import argparse
import os
import sys
import pandas as pd
import py4cytoscape as p4c

from metabodirect import preprocessing, r_control, transformations


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Program for creating molecular transformation networks, '
                    'based on previously calculated transformations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('outdir',
                        help='Output directory used to create networks with '
                        'metabodirect and the -t option',
                        metavar='OUTDIR',
                        type=str)

    parser.add_argument('metadata_file',
                        help='Metadata file used in the analysis, if a '
                        'filtered metadata was generated please enter '
                        'that one',
                        metavar='METADATA',
                        type=str)

    parser.add_argument('group',
                        help='Grouping variables for coloring and faceting '
                        'figures (Max 2)',
                        metavar='GROUP',
                        type=str,
                        nargs='+')

    args = parser.parse_args()

    if len(args.group) > 2:
        parser.error(
            'Incorrect number of variables for grouping. Please enter exactly'
            'two values')

    return args


# --------------------------------------------------
def main():
    """Create networks with Cytoscape"""

    args = get_args()
    list_dir = preprocessing.make_directories(args.outdir)
    node_table = pd.read_csv(os.path.join(list_dir[5], 'node_table.csv'))

    check = ''
    while check != 'You are connected to Cytoscape!':
        cytoscape = input(
            'Please open Cytoscape and press the ENTER key [q for quit].')
        if cytoscape == '':
            try:
                check = p4c.cytoscape_ping()
            except (ValueError, Exception):
                print('Cytoscape is not open.')
        elif cytoscape == 'q':
            sys.exit()
        else:
            check = 'Cytoscape not open'

    print('Starting network construction on Cytoscape')

    transformations.create_cytoscape_network(node_table, path=list_dir[5])

    network_stats_script = r_control.write_r_script(
        'network_stats_template.R',
        outdir=list_dir[5],
        metadata_file=args.metadata_file,
        groups=args.group,
        norm_method='max'
    )

    print(f'Running R script: {network_stats_script}')
    r_control.run_r(network_stats_script)
    print(
        f'Find results and R script in'
        f'the directory: {os.path.abspath(list_dir[5])}')

    print('------------------------\n'
          'Network creation and analysis finished'
          '\n------------------------\n')


# --------------------------------------------------
if __name__ == '__main__':
    main()
