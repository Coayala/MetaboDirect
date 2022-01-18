# get_args.py

import os
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW


# --------------------------------------------------------------
def filter_c13(df):
    """ Filter predicted formulas with 13C.
    Returns filtered df and n excluded """

    shape_i = df.shape[0]

    df = df[df['C13'] == 0]
    df = df.reset_index(drop=True)

    shape_f = df.shape[0]

    n_excluded = shape_i - shape_f

    return df, n_excluded


# --------------------------------------------------------------
def filter_mz(df, min_mz=None, max_mz=None):
    """ Optional filtering. Filter between a range of m/z """

    df = df[df['Mass'].between(min_mz, max_mz)]
    df = df.reset_index(drop=True)

    return df


# --------------------------------------------------------------
def filter_error_ppm(df, err_range=0.5):
    """ Filter based on error range of - err_range to + err_range """

    df = df[df['Error_ppm'].between(-err_range, err_range)]

    return df


# --------------------------------------------------------------
def calculate_ratios(df):
    """ Calculate ratios of O:C and H:C and indices of NOSC, GFE, DBE, AI, et al. """

    df[['C', 'H', 'O', 'N', 'S', 'P']] = df[['C', 'H', 'O', 'N', 'S',
                                             'P']].replace(np.nan, 0)

    df[['C', 'H', 'O', 'N', 'S', 'P']] = df[['C', 'H', 'O', 'N', 'S',
                                             'P']].astype(float)

    df['OC'] = df['O'] / df['C']
    df['HC'] = df['H'] / df['C']

    # Thermodynamics
    df['NOSC'] = -((4 * df['C'] + df['H'] - 3 * df['N'] - 2 * df['O'] +
                    5 * df['P'] - 2 * df['S']) / (df['C'])) + 4

    df['GFE'] = -(28.5 * df['NOSC']) + 60.3

    # Koch, B. P. and Dittmar, T. doi:10.1002/rcm.2386, 2006.
    df['DBE'] = 1 + 0.5 * (2 * df['C'] - df['H'] + df['N'] + df['P'])

    df['DBE_O'] = (1 + 0.5 *
                   (2 * df['C'] - df['H'] + df['N'] + df['P'])) - df['O']

    df['AI'] = (1 + df['C'] - df['O'] - df['S'] -
                ((df['H'] + df['P'] + df['N']) * 0.5)) / (
                       df['C'] - df['O'] - df['S'] - df['N'] - df['P'])

    df['AI_mod'] = (1 + df['C'] - (df['O'] * 0.5) - df['S'] - (
            (df['H'] + df['P'] + df['N']) * 0.5)) / (df['C'] - (df['O'] * 0.5) -
                                                     df['S'] - df['N'] - df['P'])

    df['DBE_AI'] = 1 + df['C'] - df['O'] - df['S'] - (
            0.5 * (df['H'] + df['N'] + df['P']))

    list_new_columns = [
        'OC', 'HC', 'NOSC', 'GFE', 'DBE', 'DBE_O', 'AI', 'AI_mod', 'DBE_AI'
    ]

    df[list_new_columns] = df[list_new_columns].replace([-np.inf, +np.inf],
                                                        [0, 0])

    list_dbes = ['DBE', 'DBE_O', 'DBE_AI']

    df[list_dbes] = df[list_dbes].replace(1, np.nan)

    df = molecular_formula(df)

    return df


# --------------------------------------------------------------
def calculate_classes(df):
    """ Calculate compound classes based on M. Tfaily boundaries """

    boundaries = [
        (df['OC'].between(0, 0.3)) & (df['HC'].between(1.5, 2.5)),
        (df['OC'].between(0, 0.125)) & (df['HC'].between(0.8, 1.5)),
        (df['OC'].between(0, 0.95)) & (df['HC'].between(0.2, 0.8)),
        (df['OC'].between(0.3, 0.55)) & (df['HC'].between(1.5, 2.3)),
        (df['OC'].between(0.55, 0.7)) & (df['HC'].between(1.5, 2.2)),
        (df['OC'].between(0.7, 1.5)) & (df['HC'].between(1.5, 2.5)),
        (df['OC'].between(0.125, 0.65)) & (df['HC'].between(0.8, 1.5)),
        (df['OC'].between(0.65, 1.1)) & (df['HC'].between(0.8, 1.5))
    ]

    choices = [
        'Lipid', 'Unsaturated hydrocarbon', 'Condensed hydrocarbon', 'Protein',
        'Amino sugar', 'Carbohydrate', 'Lignin', 'Tannin'
    ]

    df['Class'] = np.select(boundaries, choices, default='Other')

    df = reorder_columns(df)

    return df


# --------------------------------------------------------------
def get_list_samples(df):
    """ Returns a list of sample names and a list of Formularity columns """

    from_formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]

    calculated_indices = [
        'MolecularFormula', 'OC', 'HC', 'NOSC', 'GFE', 'DBE', 'DBE_O', 'AI', 'AI_mod', 'DBE_AI'
    ]

    from_formularity.extend(calculated_indices)

    extended_list = from_formularity

    list_samples = [x for x in df.columns.to_list() if x not in extended_list]

    return list_samples


# --------------------------------------------------------------
def reorder_columns(df):
    """ Reorder columns so indices are not in the end of the """

    from_formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]

    calculated_indices = [
        'MolecularFormula', 'OC', 'HC', 'NOSC', 'GFE', 'DBE', 'DBE_O', 'AI', 'AI_mod', 'DBE_AI'
    ]

    from_formularity.extend(calculated_indices)

    samples_list = get_list_samples(df)

    correct_order = from_formularity

    correct_order.extend(samples_list)

    df = df[correct_order]

    return df


# --------------------------------------------------------------
def molecular_formula(df):
    """ Get a molecular formula """

    elements = df[['C', 'H', 'O', 'N', 'S', 'P']]
    elements = elements.applymap(int)
    elements = elements.applymap(str)

    for col in elements.columns.to_list():
        elements[col] = col + elements[col]

    for col in elements.columns.to_list():
        elements[col] = elements[col].replace(f"{col}1", col)
        elements[col] = elements[col].replace(f"{col}0", '')

    df['MolecularFormula'] = elements['C'] + elements['H'] + elements[
        'O'] + elements['N'] + elements['S'] + elements['P']

    return df


# --------------------------------------------------------------
def get_summary(df, on='Class'):
    """ Get summary of class or elemental composition """

    samples = get_list_samples(df)

    samples.append(on)

    t = df[samples]

    t = t.melt(id_vars=[on], var_name='SampleID', value_name='NormIntensity')

    t = t[t['NormIntensity'] != 0].reset_index(drop=True)

    comp = t.groupby(['SampleID', on
                      ]).size().reset_index().rename(columns={0: 'Count'})

    comp = pd.pivot_table(comp, index='SampleID', values='Count', columns=on)

    comp_p = round(comp.div(comp.sum(axis=1), axis=0) * 100, 2).reset_index()

    return comp_p


# --------------------------------------------------------------
def get_summary_indices(df, on='NOSC'):
    """ Get the summary stats for the indices: median, mean, std, weighted mean and weighted std """

    samples = get_list_samples(df)

    samples.append(on)

    t = df[samples]

    t = t.melt(id_vars=[on], var_name='SampleID', value_name='NormIntensity')

    t = t[t['NormIntensity'] > 0].reset_index(drop=True)

    t_agg = t.groupby(['SampleID']).agg({on: ['median', 'mean', 'std']})

    t_agg.columns = t_agg.columns.map('_'.join)

    t_agg = t_agg.reset_index()

    t_agg[[on + '_w_mean', on + '_w_std']] = ''

    for sample in t['SampleID'].unique():
        # print(sample)
        temp = t[t['SampleID'] == sample]
        wdf = DescrStatsW(temp[on], weights=temp['NormIntensity'])
        t_agg.loc[t_agg['SampleID'] == sample, on + '_w_mean'] = wdf.mean
        t_agg.loc[t_agg['SampleID'] == sample, on + '_w_std'] = wdf.std

    return t_agg


# --------------------------------------------------------------
def get_matrix(df):
    """ Make a matrix for multivariate analysis like NMDS, PCA.
    Matrix based on m/z and not molec formula """

    samples = get_list_samples(df)

    samples.append('Mass')

    t = df[samples].set_index('Mass').rename(columns={'index': 'SampleID'})

    t = t.replace(0, np.nan)

    n_samples = np.ceil((len(samples) - 1) / 5)

    # print(t.shape)

    t = t.dropna(axis=0, thresh=n_samples)

    # print(t.shape)

    t = t.replace(np.nan, 0)

    matrix = t.copy()

    return matrix


# --------------------------------------------------
def make_directories(outdir):
    """Create and return a list of directories for the outputs of each step of the pipeline."""
    project_name = outdir
    preprocess_dir = '1_preprocessing_output'
    diagnostics_dir = '2_diagnostics'
    exploratory_dir = '3_exploratory'
    diversity_dir = '4_chemodiversity'
    stats_dir = '5_statistics'
    transf_dir = '6_transformations'
    transf_sub_dir = os.path.join(transf_dir, 'transf_by_sample')

    list_dir = [
        preprocess_dir, diagnostics_dir, exploratory_dir, diversity_dir, stats_dir, transf_dir, transf_sub_dir
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
    """Filter samples based on selected features and values."""

    # Get the variable a values specified for sample filtering
    filter_col = filter_by[0]
    filter_values = filter_by[1].split(sep=',')

    # Saving a new metadata file containing only the samples remaining after filtering
    filt_metadata = pd.DataFrame()
    for i in filter_values:
        filt_metadata = filt_metadata.append(metadata[metadata[filter_col] == i])
    filt_metadata.to_csv(os.path.join(path, 'filtered_metadata.csv'), index=False)

    # Saving a new input file containing only the samples remaining after filtering
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
    """Filter data based on a specified m/z range, presence of isotopes and quality."""

    print(f'Number of m/z in provided file: {df.shape[0]}')
    # Filtering masses determined to be isotopic
    filt_df, n_excluded = filter_c13(df)
    print(f'Number of masses excluded because of C13 isotope: {n_excluded}')
    print(f'Number of m/z remaining after isotope-filtering: {filt_df.shape[0]}')

    # Filter peaks based in m/z range
    filt_df = filter_mz(filt_df, filter_values[0], filter_values[1]) if filter_values else filt_df
    print(f'Number of m/z remaining after m/z-filtering: {filt_df.shape[0]}')

    # Filter samples with error higher than 0.5 ppm
    filt_df = filter_error_ppm(filt_df, err_range=0.5)
    print(f'Number of m/z after error filtering (0.5 ppm): {filt_df.shape[0]}')

    return filt_df


# --------------------------------------------------
def thermo_idx_and_classes(df
                           ):
    """Calculate thermodynamic indices and the associated class of compounds based on the molecular formula."""

    df = calculate_ratios(df)
    df = calculate_classes(df)

    return df


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
        norm_data = (input_data-sample_mean) / (sample_max-sample_min)
    elif norm_method == 'median':
        norm_data = (input_data-sample_median) / (sample_max-sample_min)
    elif norm_method == 'zscore':
        norm_data = (input_data-sample_mean) / sample_std
    elif norm_method == 'sum':
        norm_data = input_data / sample_sum
    elif norm_method == 'max':
        norm_data = input_data / sample_max
    elif norm_method == 'minmax':
        norm_data = (input_data-sample_min) / (sample_max-sample_min)
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

    # Non-normalized data for chemodiversity analysis
    input_data['Mass'] = input_data.index
    nonorm_data = input_data.reset_index(drop=True)
    nonorm_data = nonorm_data.replace(np.nan, 0)
    nonorm_data = temp.merge(nonorm_data, on='Mass')

    return norm_data, nonorm_data


# --------------------------------------------------
def calculate_summaries(df, path):
    """Get summaries for class composition, elemental composition and thermodynamic indices."""

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
