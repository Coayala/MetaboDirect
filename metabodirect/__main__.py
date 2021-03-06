# __main__.py
"""
Author : Christian Ayala <cayalaortiz@email.arizona.edu>
Date   : 2023-01-18
Purpose: Program to run MetaboDirect scripts

To get help, use metabodirect -h or visit the website https://github.com/Coayala/MetaboDirect
"""

import os
import time
import datetime
import sys
import pandas as pd
import py4cytoscape as p4c
import numpy as np
from metabodirect import get_args, preprocessing, diagnostics, r_control, transformations


# --------------------------------------------------
def main():
    """Main body to run all MetaboDirect scripts"""

    start_time = time.perf_counter()
    args = get_args.get_args()

    print('========================\nWelcome to MetaboDirect\n========================\n')
    print('Analysis starting on {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")))
    print('Results will be saved in directory: {}\n'.format(os.path.abspath(args.outdir)))

    print(f'Data file is {args.data_file}')
    print(f'Metadata file is {args.metadata_file}')

    print('\n------------------------\nStarting data pre-processing\n------------------------\n')

    df = pd.read_csv(args.data_file)
    metadata = pd.read_csv(args.metadata_file)
    list_dir = preprocessing.make_directories(args.outdir)

    if args.filter_by:
        print(f'Filtering samples based on {args.filter_by[0]} = {args.filter_by[1]}')
        df, metadata = preprocessing.sample_filtering(df, metadata, args.filter_by, path=list_dir[0])
        print(f'Filtered tables can be found in the directory: {os.path.abspath(list_dir[0])}')
    else:
        print('Option -f not detected, all samples will be used')

    df = preprocessing.data_filtering(df, mass_filter=args.mass_filter,
                                      peak_filter=args.peak_filter,
                                      error_filter=args.error_filter)
    df = preprocessing.thermo_idx_and_classes(df)
    df, df_nonorm = preprocessing.data_normalization(df, args.norm_method, args.norm_subset, args.subset_parameter,
                                                     args.log_transform)

    df_nonorm = df_nonorm[df_nonorm['C'] > 0]
    filename = os.path.join(list_dir[0], 'Report_processed_noNorm.csv')
    df_nonorm.to_csv(filename, index=False)

    filename = os.path.join(list_dir[0], 'Report_processed.csv')
    df.to_csv(filename, index=False)

    df_formulas = df[df['C'] > 0]
    filename = os.path.join(list_dir[0], 'Report_processed_MolecFormulas.csv')
    df_formulas.to_csv(filename, index=False)
    print(f'Report saved as: {filename}')

    preprocessing.calculate_summaries(df_formulas, path=list_dir[0])
    matrix_features = preprocessing.get_matrix(df_formulas)
    matrix_features.to_csv(os.path.join(list_dir[0], 'matrix_features.csv'))

    colnames = [col for col in list(df_formulas.columns) if col not in list(metadata['SampleID'])]

    # Pivot and merge dataframe with all peaks
    df_all = df.melt(id_vars=colnames, 
                     var_name='SampleID', value_name='NormIntensity')
    df_all = df_all[df_all['NormIntensity'] != 0].reset_index(drop=True)
    df_all = df_all.merge(metadata, on='SampleID')
    
    # Pivot and merge dataframe of peaks with molecular formula assignment
    df = df_formulas.melt(id_vars=colnames,
                          var_name='SampleID', value_name='NormIntensity')
    df = df[df['NormIntensity'] != 0].reset_index(drop=True)
    df = df.merge(metadata, on='SampleID')

    print('\n------------------------\nData pre-processing finished\n------------------------\n')

    print('------------------------\nStarting data diagnostics\n------------------------\n')

    print('Calculating number of peaks detected per sample')
    diagnostics.peaks_per_sample(df_all, metadata, args.group, path=list_dir[1])
    print('Calculating number of assigned molecular formulas per sample')
    diagnostics.formula_per_sample(df, metadata, args.group, path=list_dir[1])
    print("Calculating error distribution per group/s: {}".format(args.group))
    diagnostics.error_per_group(df, args.group, path=list_dir[1])

    print('\n------------------------\nData diagnostics finished\n------------------------\n')

    print('------------------------\nStarting data exploration\n------------------------\n')

    data_exploration_script = r_control.write_r_script('data_exploration_template.R', outdir=list_dir[2],
                                                       metadata_file=args.metadata_file if not args.filter_by
                                                       else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                       groups=args.group)
    print(f'Running R script: {data_exploration_script}')
    r_control.run_r(data_exploration_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[2])}')

    if args.kegg_annotation:
        print('Starting annotation of molecular formulas using the KEGG database')
        kegg_annotation_script = r_control.write_r_script('KEGG_annotation_template.R', outdir=list_dir[2],
                                                          metadata_file=args.metadata_file if not args.filter_by
                                                          else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                          groups=args.group)
        r_control.run_r(kegg_annotation_script)
        print(f'\nFind results and R script in the directory: {os.path.abspath(list_dir[2])}')
    else:
        print('\nKEGG annotation not selected. If you wish to perform a KEGG annotation run the script again using'
              'the -k/--kegg_annotation option')

    print('\n------------------------\nData Exploration finished\n------------------------\n')

    print('------------------------\nStarting chemodiversity analysis\n------------------------\n')

    data_chemodiversity_script = r_control.write_r_script('data_chemodiversity_template.R', outdir=list_dir[3],
                                                          metadata_file=args.metadata_file if not args.filter_by
                                                          else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                          groups=args.group)
    print(f'Running R script: {data_chemodiversity_script}')
    r_control.run_r(data_chemodiversity_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[3])}')

    print('\n------------------------\nChemodiversity analysis finished\n------------------------\n')

    print('------------------------\nStarting statistical analysis\n------------------------\n')

    data_statistics_script = r_control.write_r_script('data_statistics_template.R', outdir=list_dir[4],
                                                      metadata_file=args.metadata_file if not args.filter_by
                                                      else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                      groups=args.group, norm_method=args.norm_method)
    print(f'Running R script: {os.path.abspath(data_statistics_script)}')
    r_control.run_r(data_statistics_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[4])}')

    print('------------------------\nStarting transformation network analysis\n------------------------\n')

    if not args.calculate_transformations:
        print('Calculate transformations not selected. '
              'If you wish to do calculate transformations based on biochemical key please set the option "-t"')
        end_time = time.perf_counter()
        print(f'MetaboDirect finished running in {(end_time - start_time) / 60:0.4f} minutes')
        sys.exit()

    keys = transformations.get_keys(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                                 'data',
                                                 'transf_key.csv')) if args.biochem_key == 'Default key' else \
        transformations.get_keys(args.biochem_key)
    transformations.calculate_transformations(df, keys, path=list_dir[6])
    transformations.summarize_transformations(path=list_dir[5])
    node_table = transformations.get_node_table(df, path=list_dir[5])

    print(f'Finished to calculate transformatios, please find transformation files in the directory:'
          f'{os.path.abspath(list_dir[6])}')

    if not args.create_networks:
        print(f'Create networks not selected.\n'
              f'If you wish to create networks automatically please set the option -c.\n'
              f'Otherwise to create networks from previously calculated transformations using the "create_networks"'
              f'companion script')
        end_time = time.perf_counter()
        print(f'MetaboDirect finished running in {(end_time - start_time) / 60:0.4f} minutes')
        sys.exit()

    check = ''
    while check != 'You are connected to Cytoscape!':
        cytoscape = input(f'Please open Cytoscape and press the ENTER key [q for quit].')
        if cytoscape == '':
            try:
                check = p4c.cytoscape_ping()
            except (ValueError, Exception):
                print('Cytoscape is not open.')
        elif cytoscape == 'q':
            end_time = time.perf_counter()
            print(f'MetaboDirect finished running in {(end_time - start_time) / 60:0.4f} minutes')
            sys.exit()
        else:
            check = 'Cytoscape not open'

    print('Starting network construction on Cytoscape')

    transformations.create_cytoscape_network(node_table, path=list_dir[5])

    network_stats_script = r_control.write_r_script('network_stats_template.R', outdir=list_dir[5],
                                                    metadata_file=args.metadata_file if not args.filter_by
                                                    else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                    groups=args.group, norm_method=args.norm_method)

    print(f'Running R script: {network_stats_script}')
    r_control.run_r(network_stats_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[5])}')

    print('------------------------\nTransformation network analysis finished\n------------------------\n')
    print('========================\nThanks for using MetaboDirect\n========================\n')

    end_time = time.perf_counter()
    print(f'MetaboDirect finished running in {(end_time - start_time)/60:0.4f} minutes')


# --------------------------------------------------
if __name__ == '__main__':
    main()
