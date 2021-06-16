# transformations.py
import os
import glob
import time
import datetime
import seaborn as sns
import pandas as pd
import py4cytoscape as p4c
import numpy as np
from itertools import combinations


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

    # Create a dataframe that has the masses of the peaks that are present in each sample, and 0 in the peaks that are
    # not in that sample
    df_transf = pd.pivot_table(df, values='NormIntensity', index=['Mass'],
                               columns=['SampleID']).reset_index()
    df_transf['Mass'] = df_transf['Mass'].astype(float).apply(lambda x: '%.6f' % x)
    df_transf.replace([0, np.nan], ['X', 'X'], inplace=True)
    for col in df_transf.columns:
        df_transf[col] = np.where((df_transf[col] != 'X'), df_transf['Mass'], df_transf[col])
    df_transf = df_transf.drop('Mass', axis=1)
    df_transf = df_transf.replace('X', 0)
    print('Calculating m/z differences per sample column can take a little while...\n')
    i = 1

    for sample in sorted(df_transf.columns):
        print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")}]\t{i}\\{len(df_transf.columns)}\t{sample}')

        mz_list = set(float(x) for x in df_transf[sample] if float(x) > 0)
        print('   Total m/z values', len(mz_list))

        # make m/z substractions in all versus all fashion
        # doing all vs all the subtractions and filter
        result_tuples = [(x, y, round(abs(x - y), 6)) for x, y in combinations(mz_list, 2) if 1 < abs(x - y) < 766]

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
            result_df.groupby(['SampleID', 'Group', 'Transformation', 'Formula']).size().reset_index(name='Counts'))

        total_transformations = sum(result_counts['Counts'])

        result_counts['Perc_Counts'] = result_counts['Counts'] / total_transformations
        result_counts = result_counts.sort_values(by="Counts")

        # Save final_counts
        filename = os.path.join(path, 'counts_' + sample + '.csv')
        result_counts.to_csv(filename, index=False)
        i = i + 1

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

    files_path = os.path.join(path, 'transf_by_sample', 'transformations_*.csv')
    files = glob.glob(files_path)
    summary_transf = pd.DataFrame()

    for file in files:
        df = pd.read_csv(file)
        summary_transf = pd.concat([summary_transf, df], axis=0)

    filename = os.path.join(path, 'Transformations_summary_all.csv')
    summary_transf.to_csv(filename, index=False)
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

    # Create a vector of colors for the compound classes
    mol_classes = list(np.unique(node_table['Class']))
    node_colors = sns.color_palette('Set3', len(mol_classes)).as_hex()

    network_stats = []

    # Create a list of the transformation files to be used as the edge tables for the networks
    files_path = os.path.join(path, 'transf_by_sample', 'transformations_*.csv')
    edge_files = glob.glob(files_path)
    i = 1

    for file in edge_files:
        edge_table = pd.read_csv(file)
        edge_table['Feature_X'] = round(edge_table['Feature_X'], 4)
        edge_table['Feature_X'] = edge_table['Feature_X'].astype(str)
        edge_table['Feature_Y'] = round(edge_table['Feature_Y'], 4)
        edge_table['Feature_Y'] = edge_table['Feature_Y'].astype(str)
        # Renaming columns for clarity
        edge_table = edge_table.rename(columns={'Feature_X': 'source', 'Feature_Y': 'target'})
        print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")}]\t{i}\\{len(edge_files)}'
              f'\tCreating network for sample: {edge_table.SampleID.unique()[0]}')
        with pd.option_context('mode.chained_assignment', None):
            p4c.create_network_from_data_frames(None, edge_table,
                                                source_id_list='source',
                                                target_id_list='target',
                                                interaction_type_list='Transformation',
                                                title=edge_table.SampleID.unique()[0],
                                                collection='Transformation networks')
            # Load node table with data information
            p4c.load_table_data(node_table, data_key_column='id')

        # Create a vector of colors for the transformation groups
        transformation_group = list(np.unique(edge_table['Group']))
        edge_colors = sns.color_palette('dark', len(transformation_group)).as_hex()

        # Coloring the table based on compound classes and transformation groups
        p4c.set_node_shape_default('ELLIPSE')
        p4c.set_node_color_mapping('Class', mol_classes, node_colors, mapping_type='d')
        p4c.set_edge_color_mapping('Group', transformation_group, edge_colors, mapping_type='d')

        # Run the Network Analyzer command for cytoscape and capture its output in a list of lists
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

        time.sleep(5)
        i = i + 1

    network_stats = pd.DataFrame(network_stats)
    filename = os.path.join(path, 'network_summary_statistics.csv')
    network_stats.to_csv(filename, index=False)

    return
