# diagnostics.py

import os
import json
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from loguru import logger


# --------------------------------------------------
def peaks_per_sample(df, metadata, group, path):
    """Calculate and plot the number of samples that have molecular formula."""

    # Calculating the number of samples that were assigned a molecular formula
    stats_per_sample = pd.DataFrame(df.groupby(['SampleID'])[['Mass']].size()).rename(columns={0: 'Counts'})
    stats_per_sample = stats_per_sample.reset_index()
    stats_per_sample = stats_per_sample.merge(metadata, on='SampleID')
    stats_per_sample.to_csv(os.path.join(path[1], '2.1.1_stats_peaks_per_sample.csv'), index=False)

    sns.set_theme(style='white')
    sns.set_style('ticks')
    plt.rcParams['figure.figsize'] = (15, 5)

    # Plotting the number of molecular formulas per each of the samples
    p = sns.barplot(x="SampleID",
                    y="Counts",
                    hue=group[0],
                    data=stats_per_sample,
                    order=stats_per_sample.sort_values(by=[group[0], 'Counts'])['SampleID'],
                    saturation=0.7,
                    errcolor='.2',
                    errorbar="sd",
                    capsize=0.2,
                    errwidth=1.25,
                    dodge=False)

    plt.setp(p.get_xticklabels(), rotation=45, size=8, ha='right')
    plt.savefig(os.path.join(path[1], '2.1.2_stats_peaks_per_sample.png'), dpi=300, bbox_inches="tight")
    logger.info('The average number of peaks detected per sample is: {}',
                round(np.mean(stats_per_sample['Counts'])))
    plt.clf()
    
    # Saving metadata JSON
    
    figure_metadata = {
        'figure_id': '2.1.2_stats_peaks_per_sample',
        'analysis_module': '2. Diagnostics',
        'figure_type': 'Bar plot',
        'figure_title': 'Peaks per sample',
        'caption': 'Bar plot of detected metabolites',
        'description': f'This barplot shows how many peaks were detected in each sample. Bar heigths represent the amount of peaks in a sample. Bars are colored by {group[0]}',
        'insights': f'The average number of peaks detected per sample is: {round(np.mean(stats_per_sample['Counts']))}',
        'plot_info': {
            'x_axis_label': 'SampleID',
            'y_axis_label': 'Counts',
            'grouping_variables': {'hue' : group[0]},
            'modifiable_aesthetics': {'hue' : group[0]},
            'units': 'none',
            'legend': {'hue': group[0]}
        },
        'data_source ': {
            'data': os.path.join(path[0], 'Report_processed.csv')
            },
        'script_path': {
            'functions': 'diagnostics.py'},
        'functions_used ': {'plot': 'peaks_per_sample()'},
        'figure_file_info': {
            'figure_file': os.path.join(path[1], '2.1.2_stats_peaks_per_sample.png'),
            'last_modified': '',
            'resolution': '300',
            'width': '',
            'height': ''
        }   
    }
    
    with open(os.path.join(path[6], '2.1.2_stats_peaks_per_sample.json'), "w") as json_file:
        json.dump(figure_metadata, json_file)
        
    figure_index = {
        'figure_id': '2.1.2_stats_peaks_per_sample',
        'figure_title': 'Peaks per sample',
        'metadata_file': os.path.join(path[6], '2.1.2_stats_peaks_per_sample.json')
    }

    return figure_index


# --------------------------------------------------
def formula_per_sample(df, metadata, group, path):
    """Calculate and plot the number of samples that have molecular formula."""

    # Calculating the number of samples that were assigned a molecular formula
    stats_per_sample = pd.DataFrame(df.groupby(['SampleID'])[['Mass']].size()).rename(columns={0: 'Counts'})
    stats_per_sample = stats_per_sample.reset_index()
    stats_per_sample = stats_per_sample.merge(metadata, on='SampleID')
    stats_per_sample.to_csv(os.path.join(path[1], '2.2.1_stats_formula_per_sample.csv'), index=False)

    sns.set_theme(style='white')
    sns.set_style('ticks')
    plt.rcParams['figure.figsize'] = (15, 5)

    # Plotting the number of molecular formulas per each of the samples
    f = sns.barplot(x="SampleID",
                    y="Counts",
                    hue=group[0],
                    data=stats_per_sample,
                    order=stats_per_sample.sort_values(by=[group[0], 'Counts'])['SampleID'],
                    saturation=0.7,
                    errcolor='.2',
                    errorbar="sd",
                    capsize=0.2,
                    errwidth=1.25,
                    dodge=False)

    plt.setp(f.get_xticklabels(), rotation=45, size=8, ha='right')
    plt.savefig(os.path.join(path[1], '2.2.2_stats_formula_per_sample.png'), dpi=300, bbox_inches="tight")
    logger.info('The average number of masses assigned a molecular formula per sample is: {}',
                round(np.mean(stats_per_sample['Counts'])))
    plt.clf()
    
    # Saving metadata JSON

    figure_metadata = {
        'figure_id': '2.2.2_stats_formula_per_sample',
        'analysis_module': '2. Diagnostics',
        'figure_type': 'Bar plot',
        'figure_title': 'Formulas per sample',
        'caption': 'Bar plot of detected metabolites that have molecular formula',
        'description': f'This barplot shows how many metabolites were assigned a molecular formula in each sample. Bar heigths represent the amount of peaks in a sample. Bars are colored by {group[0]}',
        'insights': f'The average number of masses assigned a molecular formula per sample is: {round(np.mean(stats_per_sample['Counts']))}',
        'plot_info': {
            'x_axis_label': 'SampleID',
            'y_axis_label': 'Counts',
            'grouping_variables': {'hue': group[0]},
            'modifiable_aesthetics': {'hue': group[0]},
            'units': 'none',
            'legend': {'hue': group[0]}
        },
        'data_source ': {
            'data': os.path.join(path[0], 'Report_processed_MolecFormulas.csv')
        },
        'script_path': {
            'functions': 'diagnostics.py'},
        'functions_used ': {'plot': 'formula_per_sample()'},
        'figure_file_info': {
            'figure_file': os.path.join(path[1], '2.2.2_stats_formula_per_sample.png'),
            'last_modified': '',
            'resolution': '300',
            'width': '',
            'height': ''
        }
    }

    with open(os.path.join(path[6], '2.2.2_stats_formula_per_sample.json'), "w") as json_file:
        json.dump(figure_metadata, json_file)

    figure_index = {
        'figure_id': '2.2.2_stats_formula_per_sample',
        'figure_title': 'Peaks per sample',
        'metadata_file': os.path.join(path[6], '2.2.2_stats_formula_per_sample.json')
    }

    return figure_index


# --------------------------------------------------
def error_per_group(df, group, path):
    """Calculate and plot error distribution per sample."""
    error = df.groupby([group[0], 'Mass'] if len(group) == 1
                       else [group[0], group[1], 'Mass'])['Error_ppm'].agg(['mean']).reset_index()
    error.to_csv(os.path.join(path[1], '2.3.1_error_distribution_per_group.csv'), index=False)

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
    plt.savefig(os.path.join(path[1], '2.3.2_error_distribution_per_group.png'), dpi=300, bbox_inches="tight")
    
    figure_metadata = {
        'figure_id': '2.3.2_error_distribution_per_group',
        'analysis_module': '2. Diagnostics',
        'figure_type': 'Scatter plot',
        'figure_title': 'Error distribution per group',
        'caption': 'Scatter plot showing the distribution of errors when assigning molecular formula',
        'description': f'This barplot shows the error in ppm between the detected mass and the predicted molecular formula. Points are colored by {group[0]}',
        'insights': f'This barplot shows the error in ppm between the detected mass and the predicted molecular formula. Points are colored by {group[0]}',
        'plot_info': {
            'x_axis_label': 'Mass',
            'y_axis_label': 'Error (ppm)',
            'grouping_variables': {'hue': group[0],
                                   'facet': group[0]},
            'modifiable_aesthetics': {'hue': group[0],
                                      'facet': group[0]},
            'units': 'none',
            'legend': {'hue': group[0]}
        },
        'data_source ': {
            'data': os.path.join(path[0], 'Report_processed.csv')
        },
        'script_path': {
            'functions': 'diagnostics.py'},
        'functions_used ': {'plot': 'error_per_group()'},
        'figure_file_info': {
            'figure_file': os.path.join(path[1], '2.3.2_error_distribution_per_group.png'),
            'last_modified': '',
            'resolution': '300',
            'width': '',
            'height': ''
        }
    }

    with open(os.path.join(path[6], '2.3.2_error_distribution_per_group.json'), "w") as json_file:
        json.dump(figure_metadata, json_file)

    figure_index = {
        'figure_id': '2.3.2_error_distribution_per_group',
        'figure_title': 'Peaks per sample',
        'metadata_file': os.path.join(path[6], '2.3.2_error_distribution_per_group.json')
    }

    return figure_index
