# diagnostics.py

import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
def peaks_per_sample(df, metadata, group, path):
    """Calculate and plot the number of samples that have molecular formula."""

    # Calculating the number of samples that were assigned a molecular formula
    stats_per_sample = pd.DataFrame(df.groupby(['SampleID'])[['Mass']].size()).rename(columns={0: 'Counts'})
    stats_per_sample = stats_per_sample.reset_index()
    stats_per_sample = stats_per_sample.merge(metadata, on='SampleID')
    stats_per_sample.to_csv(os.path.join(path, 'stats_peaks_per_sample.csv'), index=False)

    sns.set_theme(style='white')
    sns.set_style('ticks')
    plt.rcParams['figure.figsize'] = (15, 5)

    # Plotting the number of molecular formulas per each of the samples
    p = sns.barplot(x="SampleID",
                    y="Counts",
                    hue=group[0],
                    data=stats_per_sample,
                    order=stats_per_sample.sort_values(by='Counts')['SampleID'],
                    saturation=0.7,
                    errcolor='.2',
                    ci="sd",
                    capsize=0.2,
                    errwidth=1.25,
                    dodge=False)

    plt.setp(p.get_xticklabels(), rotation=45, size=8, ha='right')
    plt.savefig(os.path.join(path, 'stats_peaks_per_sample.png'), dpi=300, bbox_inches="tight")
    print('The average number of peaks detected per sample is: {:.0f}'.format(
        np.mean(stats_per_sample['Counts']), ))
    plt.clf()

    return


# --------------------------------------------------
def formula_per_sample(df, metadata, group, path):
    """Calculate and plot the number of samples that have molecular formula."""

    # Calculating the number of samples that were assigned a molecular formula
    stats_per_sample = pd.DataFrame(df.groupby(['SampleID'])[['Mass']].size()).rename(columns={0: 'Counts'})
    stats_per_sample = stats_per_sample.reset_index()
    stats_per_sample = stats_per_sample.merge(metadata, on='SampleID')
    stats_per_sample.to_csv(os.path.join(path, 'stats_formula_per_sample.csv'), index=False)

    sns.set_theme(style='white')
    sns.set_style('ticks')
    plt.rcParams['figure.figsize'] = (15, 5)

    # Plotting the number of molecular formulas per each of the samples
    f = sns.barplot(x="SampleID",
                    y="Counts",
                    hue=group[0],
                    data=stats_per_sample,
                    order=stats_per_sample.sort_values(by='Counts')['SampleID'],
                    saturation=0.7,
                    errcolor='.2',
                    ci="sd",
                    capsize=0.2,
                    errwidth=1.25,
                    dodge=False)

    plt.setp(f.get_xticklabels(), rotation=45, size=8, ha='right')
    plt.savefig(os.path.join(path, 'stats_formula_per_sample.png'), dpi=300, bbox_inches="tight")
    print('The average number of masses assigned a molecular formula per sample is: {:.0f}'.format(
        np.mean(stats_per_sample['Counts']), ))
    plt.clf()

    return


# --------------------------------------------------
def error_per_group(df, group, path):
    """Calculate and plot error distribution per sample."""
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
