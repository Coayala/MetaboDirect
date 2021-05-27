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
import glob
import sys
import time
from itertools import combinations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import py4cytoscape as p4c
from shutil import which
from platform import system
from functions_metabodirect import *


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

    parser.add_argument('-t',
                        '--transformation_analysis',
                        help='Set this option to perform a transformation netwokr analysis of the samples',
                        action='store_true',
                        default=False)

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
                      help='Available methods to normalize data',
                      metavar='STR',
                      type=str,
                      choices=['mean', 'median', 'zscore', 'sum', 'max', 'minmax', 'binary', 'none'],
                      default='max')

    norm.add_argument('--norm_subset',
                      help='Subset of the data to use for normalization purpouses. LOS uses peaks in the top L '
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
                      metavar='STR',
                      type=str
                      )

    norm.add_argument('--log_transform',
                      help='Set this option to log transform the data. (Program will fail if there are peaks with '
                           'intensities of 0. Consider tranforming this values into 1 if log transformation is desired',
                      action='store_true',
                      default=False
                      )

    args = parser.parse_args()

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

    return args


# --------------------------------------------------
def make_directories(outdir):
    """Create and return a list of directories for the outputs of each step of the pipeline"""
    project_name = outdir
    preprocess_dir = '1_preprocessing_output'
    diagnostics_dir = '2_diagnostics'
    exploratory_dir = '3_exploratory'
    stats_dir = '4_statistics'
    transf_dir = '5_transformations'
    transf_sub_dir = os.path.join(transf_dir, 'transf_by_sample')

    list_dir = [
        preprocess_dir, diagnostics_dir, exploratory_dir, stats_dir, transf_dir, transf_sub_dir
    ]

    list_dir = [os.path.join(project_name, x) for x in list_dir]

    if not os.path.exists(project_name):
        os.makedirs(project_name)

    for d in list_dir:
        if not os.path.exists(d):
            os.makedirs(d)

    return list_dir


# --------------------------------------------------
def sample_filtering(df, metadata, filter_by, path):
    """Filter samples based on selected features and values"""
    filter_col = filter_by[0]
    filter_values = filter_by[1].split(sep=',')

    filt_metadata = pd.DataFrame()
    for i in filter_values:
        filt_metadata = filt_metadata.append(metadata[metadata[filter_col] == i])
    filt_metadata.to_csv(os.path.join(path, 'filtered_metadata.csv'), index=False)

    from_formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]
    col_df = filt_metadata['SampleID'].to_list()
    from_formularity.extend(col_df)
    filt_df = df[from_formularity]
    filt_df.to_csv(os.path.join(path, 'filtered_input.csv'), index=False)

    return filt_df, filt_metadata


# --------------------------------------------------
def data_filtering(df, filter_values):
    """Filter data based on a specified m/z range, presence of isotopes and quality"""

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
def thermo_idx_and_classes(df, path):
    """Calculate thermodynamic indices and the associated class of compounds based on the molecular formula"""

    df = calculate_ratios(df)
    df = calculate_classes(df)
    df_formulas = df[df['C'] > 0]

    filename = os.path.join(path, 'Report_processed.csv')
    df.to_csv(filename, index=False)
    filename = os.path.join(path, 'Report_processed_MolecFormulas.csv')
    df_formulas.to_csv(filename, index=False)

    print(f'Report saved as: {filename}')

    return df_formulas


# --------------------------------------------------
def data_normalization(df, norm_method, norm_subset, subset_parameter=1, log=False):
    """Normalize direct injection data based on the specified parameters."""

    samples = get_list_samples(df)
    input_data = df[samples]

    # Pivoting the data, so it removes the masses that are not present in any of the samples
    input_data.insert(0, 'Mass', list(df['Mass']))
    input_data = pd.melt(input_data, id_vars='Mass', value_vars=samples, var_name='SampleID', value_name='Intensity')
    input_data = input_data[input_data['Intensity'] > 0]
    input_data = input_data.pivot(index='Mass', columns='SampleID', values='Intensity')
    npeaks, nsamples = np.shape(input_data)

    if log:
        assert 0 not in input_data, 'Cannot calculate log 0, please consider changing this values to 1'
        input_data = np.log(input_data)

    # Perform sample subset to calculate normalization factors
    if norm_subset == 'PPP':
        min_observations = int(subset_parameter * nsamples)
        selected_peaks = []
        for peak in range(0, npeaks):
            if sum(input_data.iloc[peak, :] > 0) > min_observations:
                keep = input_data.iloc[peak, :].name
                selected_peaks.append(keep)
    elif norm_subset == 'LOS':
        ntop_peaks = int(subset_parameter * npeaks)
        keep = []
        for sample in range(0, nsamples):
            ord_peaks = input_data.iloc[:, sample].sort_values(ascending=False)
            keep.extend(list(ord_peaks.iloc[0:ntop_peaks, ].index))
        selected_peaks = list(np.unique(keep))
    else:
        selected_peaks = list(input_data.index)

    # Calculate normalization factors
    sample_mean = input_data.loc[selected_peaks, :].mean(axis=0, skipna=True)
    sample_min = input_data.loc[selected_peaks, :].min(axis=0, skipna=True)
    sample_max = input_data.loc[selected_peaks, :].max(axis=0, skipna=True)
    sample_median = input_data.loc[selected_peaks, :].median(axis=0, skipna=True)
    sample_std = input_data.loc[selected_peaks, :].std(axis=0, skipna=True)
    sample_sum = input_data.loc[selected_peaks, :].sum(axis=0, skipna=True)

    # Normalize data based on chosen method
    if norm_method == 'mean':
        norm_data = (input_data - sample_mean) / (sample_max - sample_min)
    elif norm_method == 'median':
        norm_data = (input_data - sample_median) / (sample_max - sample_min)
    elif norm_method == 'zscore':
        norm_data = (input_data - sample_mean) / sample_std
    elif norm_method == 'sum':
        norm_data = input_data / sample_sum
    elif norm_method == 'max':
        norm_data = input_data / sample_max
    elif norm_method == 'minmax':
        norm_data = (input_data - sample_min) / (sample_max - sample_min)
    elif norm_method == 'binary':
        norm_data = input_data.copy()
        norm_data[norm_data > 0] = 1
    else:
        norm_data = input_data.copy()

    norm_data['Mass'] = norm_data.index
    norm_data = norm_data.reset_index(drop=True)
    norm_data = norm_data.replace(np.nan, 0)
    temp = df[df.columns[~df.columns.isin(samples)]]
    norm_data = temp.merge(norm_data, on='Mass')

    return norm_data


# --------------------------------------------------
def calculate_summaries(df, path):
    """Get summaries for class composition, elemental composition and thermodynamic indices"""

    class_comp = get_summary(df, on='Class')
    el_comp = get_summary(df, on='El_comp')

    idx_stats = pd.DataFrame(get_list_samples(df), columns=['SampleID'])
    for i in ['NOSC', 'GFE', 'DBE', 'AI']:
        temp = get_summary_indices(df, i)
        idx_stats = idx_stats.merge(temp, on='SampleID')

    class_comp.to_csv(os.path.join(path, 'class_composition.csv'), index=False)
    el_comp.to_csv(os.path.join(path, 'elemental_composition.csv'), index=False)
    idx_stats.to_csv(os.path.join(path, 'indices_statistics.csv'), index=False)

    return


# --------------------------------------------------
def formula_per_sample(df, metadata, path):
    """Calculate and plot the number of samples that have molecular formula"""

    stats_per_sample = pd.DataFrame(df.groupby(['SampleID'])[['Mass']].size()).rename(columns={0: 'Counts'})
    stats_per_sample = stats_per_sample.reset_index()
    stats_per_sample = stats_per_sample.merge(metadata, on='SampleID')
    stats_per_sample.to_csv(os.path.join(path, 'stats_formula_per_sample.csv'), index=False)

    sns.set_theme(style='white')
    sns.set_style('ticks')
    plt.rcParams['figure.figsize'] = (15, 5)

    p = sns.barplot(x="SampleID",
                    y="Counts",
                    hue=metadata.columns[2],
                    data=stats_per_sample,
                    order=stats_per_sample.sort_values(by='Counts')['SampleID'],
                    saturation=0.7,
                    errcolor='.2',
                    ci="sd",
                    capsize=0.2,
                    errwidth=1.25,
                    dodge=False)

    plt.setp(p.get_xticklabels(), rotation=90, size=3)
    plt.savefig(os.path.join(path, 'stats_formula_per_sample.png'), dpi=300, bbox_inches="tight")
    print('The average number of masses assigned a molecular formula per sample is: {:.0f}'.format(
        np.mean(stats_per_sample['Counts']), ))

    return


# --------------------------------------------------
def error_per_group(df, group, path):
    """Calculate and plot error distribution per sample"""
    error = df.groupby([group[0], 'Mass'] if len(group) == 1
                       else [group[0], group[1], 'Mass'])['Error_ppm'].agg(['mean']).reset_index()
    error.to_csv(os.path.join(path, 'error_distribution_per_group.csv'), index=False)

    sns.set_theme(style="white")
    sns.set_style("ticks")

    grid = sns.FacetGrid(error,
                         col=group[0],
                         row=None if len(group) == 1 else group[1],
                         hue=group[0],
                         height=5)

    grid.map(plt.scatter, "Mass", "mean", s=1)
    grid.add_legend()
    grid.set_axis_labels(y_var="Error (ppm)")
    plt.savefig(os.path.join(path, 'error_distribution_per_group.png'), dpi=300, bbox_inches="tight")

    return


# --------------------------------------------------
def write_r_script(rscript, outdir, metadata_file, groups, norm_method='max'):
    """Function to write R scripts based on the provided Rscripts templates"""
    current_dir = os.getcwd().replace('\\', '/')
    metadata_file = metadata_file.replace('\\', '/')
    metabo_home = os.path.split(os.path.realpath(__file__))[0].replace('\\', '/')
    r_in = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'R_scripts_templates', rscript)
    r_in = open(r_in)
    r_file = os.path.join(outdir, rscript.replace('_template', ''))
    r_fh = open(r_file, 'w')
    for line in r_in:
        r_fh.write(line.replace('%currentdir%', current_dir)
                   .replace('%metadata%', metadata_file)
                   .replace('%Metabo_HOME%', metabo_home)
                   .replace('%outdir%', os.path.split(outdir)[0])
                   .replace('%group1%', groups[0])
                   .replace('%group2%', groups[1] if len(groups) == 2 else 'NULL')
                   .replace('%norm_method%', norm_method)
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
        # if out:
        #     print(out)

        if err:
            print(err)
    else:
        print('Rscript not accesible')


# ---------------------------------------------------
def get_keys(keyfile):
    """Get Biochemical Transformation Keys"""

    key = pd.read_csv(keyfile)
    key['mf'] = key['mf'].astype(float).apply(lambda x: '%.6f' % x)
    key = key.sort_values(by=['mf'])
    key_tuples = list(zip(key.Group, key.Transformation, key.Formula, key.mf))

    return key_tuples


# --------------------------------------------------
def calculate_transformations(df, keys, path):
    """Function to calculate transformations for transformation networks"""
    df_transf = pd.pivot_table(df, values='NormIntensity', index=['Mass'],
                               columns=['SampleID']).reset_index()
    df_transf['Mass'] = df_transf['Mass'].astype(float).apply(lambda x: '%.6f' % x)
    df_transf.replace([0, np.nan], ['X', 'X'], inplace=True)
    for col in df_transf.columns:
        df_transf[col] = np.where((df_transf[col] != 'X'), df_transf['Mass'], df_transf[col])
    df_transf = df_transf.drop('Mass', axis=1)
    df_transf = df_transf.replace('X', 0)
    print('Calculating m/z differences per sample column can take a little while...\n')

    for sample in sorted(df_transf.columns):
        print(sample)

        mz_list = set(float(x) for x in df_transf[sample] if float(x) > 0)
        print('   Total m/z values', len(mz_list))

        # make m/z substractions in all versus all fashion
        # doing all vs all the subtractions and filter
        result_tuples = [(x, y, round(abs(x - y), 6)) for x, y in combinations(mz_list, 2) if 1 < x - y < 766]

        result_tuples = [(
            r[0], r[1], r[2],
            k[0], k[1], k[2], k[3]
        ) for r in result_tuples for k in keys
            if r[2] - 0.001 <= float(k[3]) <= r[2] + 0.001]

        # make np.array from list of lists
        result_tuples = np.vstack(result_tuples)

        # make pd df
        result_df = pd.DataFrame(result_tuples, columns=[
            'Feature_X', 'Feature_Y', 'Difference',
            'Group', 'Transformation', 'Formula', 'mf'])

        result_df['SampleID'] = sample

        print('   Saving results')
        filename = os.path.join(path, 'transformations_' + sample + '.csv')
        result_df.to_csv(filename, index=False)

        # Compile counts
        result_counts = pd.DataFrame(
            result_df.groupby(['SampleID', 'Group', 'Transformation']).size().reset_index(name='Counts'))

        result_counts = result_counts.sort_values(by="Counts")

        # Save final_counts
        filename = os.path.join(path, 'counts_' + sample + '.csv')
        result_counts.to_csv(filename, index=False)

    print("\u2713 Done!")
    return


# --------------------------------------------------
def summarize_transformations(path):
    """Create a table summarizing the number of transformations of each sample"""

    files_path = os.path.join(path, 'transf_by_sample', 'counts_*.csv')
    files = glob.glob(files_path)
    summary_counts = pd.DataFrame()

    for file in files:
        df = pd.read_csv(file)
        summary_counts = pd.concat([summary_counts, df], axis=0)

    filename = os.path.join(path, 'Transformations_summary_counts.csv')
    summary_counts.to_csv(filename, index=False)
    return


# --------------------------------------------------
def get_node_table(df, path):
    """Create a node table for the transformation networks"""
    node_table = df[['Mass', 'C', 'H', 'O', 'N', 'S', 'P', 'OC', 'HC', 'NOSC',
                     'GFE', 'Class', 'MolecularFormula', 'El_comp']].drop_duplicates('Mass')

    filename = os.path.join(path, 'node_table.csv')
    node_table.to_csv(filename, index=False)

    return node_table


# --------------------------------------------------
def create_cytoscape_network(node_table, path):
    """Create a cytoscape network using the node table"""
    node_table['Mass'] = round(node_table['Mass'], 4)
    node_table['Mass'] = node_table['Mass'].astype(str)
    node_table = node_table.rename(columns={'Mass': 'id'})

    mol_classes = list(np.unique(node_table['Class']))
    node_colors = sns.color_palette('Set3', len(mol_classes)).as_hex()

    network_stats = []
    files_path = os.path.join(path, 'transf_by_sample', 'transformations_*.csv')
    edge_files = glob.glob(files_path)

    for file in edge_files:
        edge_table = pd.read_csv(file)
        edge_table['Feature_X'] = round(edge_table['Feature_X'], 4)
        edge_table['Feature_X'] = edge_table['Feature_X'].astype(str)
        edge_table['Feature_Y'] = round(edge_table['Feature_Y'], 4)
        edge_table['Feature_Y'] = edge_table['Feature_Y'].astype(str)
        edge_table = edge_table.rename(columns={'Feature_X': 'source', 'Feature_Y': 'target'})
        print(f'Creating transformations network for sample: {edge_table.SampleID.unique()[0]}')
        with pd.option_context('mode.chained_assignment', None):
            p4c.create_network_from_data_frames(None, edge_table,
                                                source_id_list='source',
                                                target_id_list='target',
                                                interaction_type_list='Transformation',
                                                title=edge_table.SampleID.unique()[0],
                                                collection='Transformation networks')
            p4c.load_table_data(node_table, data_key_column='id')

        transformation_group = list(np.unique(edge_table['Group']))
        edge_colors = sns.color_palette('dark', len(transformation_group)).as_hex()

        p4c.set_node_shape_default('ELLIPSE')
        p4c.set_node_color_mapping('Class', mol_classes, node_colors, mapping_type='d')
        p4c.set_edge_color_mapping('Group', transformation_group, edge_colors, mapping_type='d')

        network_analyzer = p4c.commands_post('analyzer analyze')
        network_stats.append(
            {'SampleID': edge_table.SampleID.unique()[0],
             'netHE': network_analyzer['heterogeneity'],
             'nNodes': network_analyzer['nodeCount'],
             'nEdges': network_analyzer['edgeCount'],
             'meanNeighbors': network_analyzer['avNeighbors'],
             'clusterCoeff': network_analyzer['cc'],
             'netDensity': network_analyzer['density'],
             'netCentralization': network_analyzer['centralization'],
             'connectComponents': network_analyzer['ncc'],
             'pathLength': network_analyzer['avSpl'],
             'netDiameter': network_analyzer['diameter'],
             'netRadius': network_analyzer['radius']}
        )

        time.sleep(10)

    network_stats = pd.DataFrame(network_stats)
    filename = os.path.join(path, 'network_summary_statistics.csv')
    network_stats.to_csv(filename, index=False)

    return


# --------------------------------------------------
def main():
    """Main body to run all MetaboDirect scripts"""

    args = get_args()

    print('========================\nWelcome to MetaboDirect\n========================\n')
    print('Analysis starting on {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")))
    print('Results will be saved in directory: {}\n'.format(os.path.abspath(args.outdir)))

    print(f'Data file is {args.data_file}')
    print(f'Metadata file is {args.metadata_file}')

    print('\n------------------------\nStarting data pre-processing\n------------------------\n')

    df = pd.read_csv(args.data_file)
    metadata = pd.read_csv(args.metadata_file)
    list_dir = make_directories(args.outdir)

    if args.filter_by:
        print(f'Filtering samples based on {args.filter_by[0]} = {args.filter_by[1]}')
        df, metadata = sample_filtering(df, metadata, args.filter_by, path=list_dir[0])
        print(f'Filtered tables can be found in the directory: {os.path.abspath(list_dir[0])}')
    else:
        print('Option -f not detected, all samples will be used')

    df = data_filtering(df, args.mass_filter)
    df = thermo_idx_and_classes(df, path=list_dir[0])
    df = data_normalization(df, args.norm_method, args.norm_subset, args.subset_parameter, args.log_transform)
    calculate_summaries(df, path=list_dir[0])
    matrix_features = get_matrix(df)
    matrix_features.to_csv(os.path.join(list_dir[0], 'matrix_features.csv'))

    colnames = [col for col in list(df.columns) if col not in list(metadata['SampleID'])]
    df = df.melt(id_vars=colnames,
                 var_name='SampleID', value_name='NormIntensity')
    df = df[df['NormIntensity'] > 0].reset_index(drop=True)
    df = df.merge(metadata, on='SampleID')

    print('\n------------------------\nData pre-processing finished\n------------------------\n')

    print('------------------------\nStarting data diagnostics\n------------------------\n')

    print('Calculating number of assigned molecular formulas per sample')
    formula_per_sample(df, metadata, path=list_dir[1])
    print("Calculating error distribution per group/s: {}".format(args.group))
    error_per_group(df, args.group, path=list_dir[1])

    print('\n------------------------\nData diagnostics finished\n------------------------\n')

    print('------------------------\nStarting data exploration\n------------------------\n')

    data_exploration_script = write_r_script('data_exploration_template.R', outdir=list_dir[2],
                                             metadata_file=args.metadata_file if not args.filter_by
                                             else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                             groups=args.group)
    print(f'Running R script: {data_exploration_script}')
    run_r(data_exploration_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[2])}')

    if args.kegg_annotation:
        print('Starting annotation of molecular formulas using the KEGG database')
        kegg_annotation_script = write_r_script('KEGG_annotation_template.R', outdir=list_dir[2],
                                                metadata_file=args.metadata_file if not args.filter_by
                                                else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                groups=args.group)
        run_r(kegg_annotation_script)
        print(f'Find results and R script in the directory: {os.path.abspath(list_dir[2])}')
    else:
        print('KEGG annotation not selected. If you wish to perform a KEGG annotation run the script again using'
              'the -k/--kegg_annotation option')

    print('\n------------------------\nData Exploration finished\n------------------------\n')

    print('------------------------\nStarting statistical analysis\n------------------------\n')

    data_statistics_script = write_r_script('data_statistics_template.R', outdir=list_dir[3],
                                            metadata_file=args.metadata_file if not args.filter_by
                                            else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                            groups=args.group)
    print(f'Running R script: {os.path.abspath(data_statistics_script)}')
    run_r(data_statistics_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[3])}')

    print('------------------------\nStarting transformation network analysis\n------------------------\n')

    df.to_csv('df.csv', index=False)

    if not args.transformation_analysis:
        print('Transformation analysis not selected. '
              'If you wish to do a transformation analysis please set the option "-t"')
        sys.exit()

    keys = get_keys(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                 'data',
                                 'transf_key.csv')) if args.biochem_key == 'Default key' else get_keys(args.biochem_key)
    calculate_transformations(df, keys, path=list_dir[5])
    summarize_transformations(path=list_dir[4])
    node_table = get_node_table(df, path=list_dir[4])

    check = ''
    while check != 'You are connected to Cytoscape!':
        cytoscape = input(f'Please open Cytoscape and press the ENTER key [q for quit].')
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

    create_cytoscape_network(node_table, path=list_dir[4])

    network_stats_script = write_r_script('network_stats_template.R', outdir=list_dir[4],
                                          metadata_file=args.metadata_file if not args.filter_by
                                          else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                          groups=args.group, norm_method=args.norm_method)

    print(f'Running R script: {network_stats_script}')
    run_r(network_stats_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[4])}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
