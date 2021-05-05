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
from itertools import combinations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from shutil import which
from platform import system
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

    parser.add_argument('metadata_file',
                        help='Name of the fie with the sample information (metadata) in tabular format',
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
                        default=['group1', 'group2'])

    parser.add_argument('-f',
                        '--filter_by',
                        help='Filter samples based on metadata. First enter the name of the feature,'
                             'followed by the values associated with the samples you want to keep in the analysis.'
                             '(Example -f Habitat Bog,Palsa)',
                        metavar='STR',
                        type=str,
                        nargs=2)

    parser.add_argument('-k',
                        '--biochem_key',
                        help='File with the biochemical key to use for the transformation network',
                        metavar='STR',
                        type=str,
                        default='Default key')

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error(f'File {args.data_file} does not exist')

    if not os.path.isfile(args.metadata_file):
        parser.error(f'File {args.metadata_file} does not exist')

    if len(args.group) > 2:
        parser.error(f'Incorrect number of variables for grouping. Please enter exactly two values')

    if args.mass_filter:
        if args.mass_filter[1] < args.mass_filter[0]:
            parser.error(f'Incorrect order of values for filtering. Please enter min_mz first (e.g. -f 200 900)')

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
    df = normalize_intensities(df)
    df_formulas = df[df['C'] > 0]

    filename = os.path.join(path, 'Report_processed.csv')
    df.to_csv(filename, index=False)
    filename = os.path.join(path, 'Report_processed_MolecFormulas.csv')
    df_formulas.to_csv(filename, index=False)

    print(f'Report saved as: {filename}')

    return df_formulas


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
                         hue=group[1],
                         height=5)

    grid.map(plt.scatter, "Mass", "mean", s=1)
    grid.add_legend()
    grid.set_axis_labels(y_var="Error (ppm)")
    plt.savefig(os.path.join(path, 'error_distribution_per_group.png'), dpi=300, bbox_inches="tight")

    return


# --------------------------------------------------
def write_r_script(rscript, outdir, metadata_file, groups):
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
                   .replace('%groups%',
                            '"{}","{}"'.format(groups[0], groups[1]) if len(groups) == 2 else '"{}"'.format(groups[0])
                            )
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
                      'GFE', 'Class', 'Formula', 'Element']].drop_duplicates('Mass')

    node_table['Mass'] = node_table['Mass'].astype(float).apply(lambda x: '%.6f' % x).astype(str)

    filename = os.path.join(path, 'node_table.csv')
    node_table.to_csv(filename, index=False)

    return node_table


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

    keys = get_keys(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                 'data',
                                 'transf_key.csv')) if args.biochem_key == 'Default key' else get_keys(args.biochem_key)
    calculate_transformations(df, keys, path=list_dir[5])
    summarize_transformations(path=list_dir[4])
    node_table = get_node_table(df, path=list_dir[4])



# --------------------------------------------------
if __name__ == '__main__':
    main()
