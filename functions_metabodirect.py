#
# Functions for processing of raw output from Formularity
# modified from functions_metabodirect.py by Nathalia Graf Grachet (https://github.com/nathaliagg/fusarium_wilt_lettuce_di)
#

import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW


# --------------------------------------------------------------
def filter_C13(df):
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
def filter_error_ppm(df, err_range = 0.5):
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

    list_DBEs = ['DBE', 'DBE_O', 'DBE_AI']

    df[list_DBEs] = df[list_DBEs].replace(1, np.nan)
    
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

    from_Formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]

    calculated_indices = [
        'MolecularFormula', 'OC', 'HC', 'NOSC', 'GFE', 'DBE', 'DBE_O', 'AI', 'AI_mod', 'DBE_AI'
    ]

    from_Formularity.extend(calculated_indices)

    extended_list = from_Formularity

    list_samples = [x for x in df.columns.to_list() if x not in extended_list]

    return list_samples


# --------------------------------------------------------------
def reorder_columns(df):
    """ Reorder columns so indices are not in the end of the """

    from_Formularity = [
        'Mass', 'C', 'H', 'O', 'N', 'C13', 'S', 'P', 'Na', 'El_comp', 'Class',
        'NeutralMass', 'Error_ppm', 'Candidates'
    ]

    calculated_indices = [
        'MolecularFormula', 'OC', 'HC', 'NOSC', 'GFE', 'DBE', 'DBE_O', 'AI', 'AI_mod', 'DBE_AI'
    ]

    from_Formularity.extend(calculated_indices)

    samples_list = get_list_samples(df)

    correct_order = from_Formularity

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
def normalize_intensities(df):
    """ Normalize intensities (divide by the max value) """
    
    samples = get_list_samples(df)

    for col in samples:
        max_value = df[col].max()
        df[col] = df[col]/max_value
        
    return df


# --------------------------------------------------------------
def get_summary(df, on = 'Class'):
    """ Get summary of class or elemental composition """
    
    samples = get_list_samples(df)
    
    samples.append(on)

    t = df[samples]
    
    t = t.melt(id_vars = [on], var_name = 'SampleID', value_name = 'NormIntensity')

    t = t[t['NormIntensity']>0].reset_index(drop=True)
    
    comp = t.groupby(['SampleID', on
                        ]).size().reset_index().rename(columns={ 0 : 'Count'})

    comp = pd.pivot_table(comp, index = 'SampleID', values='Count', columns=on)

    # comp['Total_Peaks_with_Formula'] = comp.sum(axis=1)

    comp_p = round(comp.div(comp.sum(axis=1), axis=0)*100,2).reset_index()

    # comp_p = comp_p.melt(id_vars = ['SampleID'], var_name = on, value_name = 'Percent')
    
    return comp_p


# --------------------------------------------------------------
def get_summary_indices(df, on='NOSC'):
    """ Get the summary stats for the indices: median, mean, std, weighted mean and weighted std """

    samples = get_list_samples(df)

    samples.append(on)

    t = df[samples]

    t = t.melt(id_vars = [on], var_name = 'SampleID', value_name = 'NormIntensity')

    t = t[t['NormIntensity']>0].reset_index(drop=True)

    t_agg = t.groupby(['SampleID']).agg({on : ['median', 'mean', 'std']})

    t_agg.columns = t_agg.columns.map('_'.join)

    t_agg = t_agg.reset_index()

    t_agg[[on+'_w_mean', on+'_w_std']] = ''


    for sample in t['SampleID'].unique():
        # print(sample)
        temp = t[t['SampleID']==sample]
        wdf = DescrStatsW(temp[on], weights=temp['NormIntensity'])
        t_agg.loc[t_agg['SampleID']==sample , on+'_w_mean'] = wdf.mean
        t_agg.loc[t_agg['SampleID']==sample , on+'_w_std'] = wdf.std

    return t_agg


# --------------------------------------------------------------
def get_matrix(df):
    """ Make a matrix for multivariate analysis like NMDS, PCA. 
    Matrix based on m/z and not molec formula """

    samples = get_list_samples(df)

    samples.append('Mass')

    t = df[samples].set_index('Mass').rename(columns={'index':'SampleID'})

    t = t.replace(0, np.nan)

    n_samples = np.ceil((len(samples)-1)/5)

    # print(t.shape)

    t = t.dropna(axis=0, thresh = n_samples)

    # print(t.shape)

    t = t.replace(np.nan, 0)
    
    matrix = t.copy()

    return matrix


# --------------------------------------------------------------
