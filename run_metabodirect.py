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
import subprocess
from statsmodels.stats.weightstats import DescrStatsW
from shutil import which
from platform import system
import matplotlib.pyplot as plt
import seaborn as sns
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

    parser.add_argument('-f',
                        '--filter',
                        help='Range to filter m/z data (min_mz, max_mz). The pipeline will not filter '
                             'm/z values by default',
                        metavar='mz',
                        type=float,
                        nargs=2,
                        default=None)

    parser.add_argument('-g',
                        '--group',
                        help='Grouping variables for coloring and faceting figures (Max 2)',
                        metavar='mz',
                        type=str,
                        nargs='+',
                        default=None)

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error(f'File {args.data_file} does not exist')

    if not os.path.isfile(args.metadata_file):
        parser.error(f'File {args.metadata_file} does not exist')

    if len(args.group) > 2:
        parser.error(f'Incorrect number of variables for grouping. Please enter exactly two values')

    if args.filter:
        if args.filter[1] < args.filter[0]:
            parser.error(f'Incorrect order of values for filtering. Please enter min_mz first (e.g. -f 200 900')

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

    list_dir = [
        preprocess_dir, diagnostics_dir, exploratory_dir, stats_dir, transf_dir
    ]

    list_dir = [os.path.join(project_name, x) for x in list_dir]

    if not os.path.exists(project_name):
        os.makedirs(project_name)

    for d in list_dir:
        if not os.path.exists(d):
            os.makedirs(d)

    return list_dir


# --------------------------------------------------
def data_filtering(df_file, filter_values):
    """Filter data based on a specified m/z range, presence of isotopes and quality"""

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

    samples = get_list_samples(df)
    samples.append('Mass')
    samples.append('Error_ppm')
    df = df[samples]
    df = df.melt(id_vars=['Mass', 'Error_ppm'], var_name='SampleID', value_name='NormIntensity')
    df = df[df['NormIntensity'] > 0].reset_index(drop=True)
    df = df.merge(metadata, on='SampleID')

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

    return df


# --------------------------------------------------
def error_per_sample(df, metadata, path):
    """Calculate and plot error distribution per sample"""
    error = df.groupby(['SampleID', 'Mass'])['Error_ppm'].agg(['mean']).reset_index()
    error = error.reset_index()
    error = error.merge(metadata, on='SampleID')
    error.to_csv(os.path.join(path, 'error_distribution_per_sample.csv'), index=False)

    sns.set_theme(style="white")
    sns.set_style("ticks")

    grid = sns.FacetGrid(error,
                         col="SampleID",
                         col_wrap=10,
                         hue=metadata.columns[2],
                         height=5)

    grid.map(plt.scatter, "Mass", "mean", s=1)
    grid.set_axis_labels(y_var="Error (ppm)")
    plt.savefig(os.path.join(path, 'error_distribution_per_sample.png'), dpi=300, bbox_inches="tight")

    return


# --------------------------------------------------
def write_r_script(rscript, outdir, metadata_file, groups):
    """Function to write R scripts based on the provided Rscripts templates"""
    r_in = os.path.join('R_scripts_templates', rscript)
    r_in = open(r_in)
    r_file = os.path.join(outdir, rscript.replace('_template', ''))
    r_fh = open(r_file, 'w')
    for line in r_in:
        r_fh.write(line.replace('%metadata%', metadata_file)
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

    list_dir = make_directories(args.outdir)
    df = data_filtering(args.data_file, args.filter)
    df = thermo_idx_and_classes(df, path=list_dir[0])
    calculate_summaries(df, path=list_dir[0])
    matrix_features = get_matrix(df)
    matrix_features.to_csv(os.path.join(list_dir[0], 'matrix_features.csv'))

    print('\n------------------------\nData pre-processing finished\n------------------------\n')

    print('------------------------\nStarting data diagnostics\n------------------------\n')

    metadata = pd.read_csv(args.metadata_file)
    print('Calculating number of assigned molecular formulas per sample')
    df_meta = formula_per_sample(df, metadata, path=list_dir[1])
    print("Calculating error distribution per sample")
    error_per_sample(df_meta, metadata, path=list_dir[1])

    print('\n------------------------\nData diagnostics finished\n------------------------\n')

    print('------------------------\nStarting data exploration\n------------------------\n')

    data_exploration_script = write_r_script('data_exploration_template.R', outdir=list_dir[2],
                                             metadata_file=args.metadata_file, groups=args.group)
    print(f'Running R script: {data_exploration_script}')
    run_r(data_exploration_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[2])}')

    print('\n------------------------\nData Exploration finished\n------------------------\n')

    print('------------------------\nStarting statistical analysis\n------------------------\n')

    data_statistics_script = write_r_script('data_statistics_template.R', outdir=list_dir[3],
                                            metadata_file=args.metadata_file, groups=args.group)
    print(f'Running R script: {os.path.abspath(data_statistics_script)}')
    run_r(data_statistics_script)
    print(f'Find results and R script in the directory: {os.path.abspath(list_dir[3])}')



# --------------------------------------------------
if __name__ == '__main__':
    main()
