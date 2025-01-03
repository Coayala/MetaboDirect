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
    key['mf'] = key['mf'].astype(float)
    key = key.sort_values(by=['mf'])
    key_tuples = list(zip(key.Group, key.Transformation, key.Formula, key.mf))

    return key_tuples


# --------------------------------------------------
def calculate_transformations(df, keys, path, err_thresh=0.001):
    import bisect
    import operator
    """Function to calculate transformations for transformation networks"""

    # Create a dataframe with m/z and per-sample intensity columns.

    df_transf = pd.pivot_table(df, values='NormIntensity', index=['Mass'],
                               columns=['SampleID']).reset_index()
    i = 1

    for sample in df_transf.columns[1:]:  # Skip mass column
        print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")}]'
              f'{i}\\{len(df_transf.columns)}\t{sample}')

        # For each sample, extract masses for all compounds present
        # (i.e. where intensity > 0)
        mz_list = df_transf.Mass[df_transf[sample] > 0]        
        print('\t\tTotal m/z values', len(mz_list))

        # Create a distance matrix for the m/z values, where:
        #   diff[i,j] = mz_list[j] - mz_list[i]
        # It will be antisymmetric, with the upper triangle positive.

        diffs = (mz_list.values - mz_list.values[:, None])

        # Create a set of X,Y coordinates for diffs. We can use these
        # to find the original masses.
        X, Y = np.indices((len(mz_list), len(mz_list)))

        # Now, for all mass differences (in the range of 1-766 Da),
        # create a list of (delta, src_idx, target_idx) and sort by mass delta.
        valid_locs = np.where((diffs > 1) & (diffs < 766), True, False)
        candidates = zip(diffs[valid_locs], X[valid_locs], Y[valid_locs])

        # Sort by mass delta (first element in tuple).
        # Using operator.itemgetter is a skosh faster than
        # using key=lambda x: x[0], but only if we create it beforehand.
        getter = operator.itemgetter(0)
        candidates = sorted(candidates, key=getter)

        c_min = 0
        c_max = len(candidates)

        # Now, since we have two sorted lists (candidates and keys), we can
        # compare them. Timing tests show that the fastest method is to
        # use binary search ('bisect')

        results = []
        for grp, xform, formula, mf in keys:
            # For each transformation, find all of our candidates within the
            # range of delta +/- err_thresh.
            # We need to convert the 'mf +/- err) to a tuple, as bisect can
            # only compare like types:            
            lo = bisect.bisect_left(candidates, (mf - err_thresh,), lo=c_min)
            hi = bisect.bisect_left(candidates, (mf + err_thresh,), lo=c_min)
            if lo == hi:
                # No candidates fall in range for this transform
                continue
            if lo == c_max:
                # We've exhausted all our candidates
                break
            # Every candiate between lo and hi has a mass that matches our
            # transform.
            results.append([(mz_list.iloc[cand[1]],  # Source m/z
                             mz_list.iloc[cand[2]],  # Target m/z
                             cand[0],  # Actual delta m/z
                             grp, xform, formula, mf)
                            for cand in candidates[lo:hi]])
            # We can speed things up a bit by only searching masses greater than
            # or equal to what we've searched so far by setting c_min.
            # If we set this to 'hi', it will exclude all previous matches, so compound
            # pairs will only match a single transformation. By setting it to 'lo',
            # all possible transformations are found.
            c_min = lo

        # Note: in this algorithm, if two formulas have the same mass delta
        # (such as leucine and isoleucine), or if two mass deltas are within
        # err_thresh of each other, only a single transformation will be
        # discovered (whichever comes first).
        # In the previous version, *all* matching transforms are returned;
        # in particular, all leucine transforms are double counted as
        # isoleucine.
        # We could easily revert to the previous behavior without loss of speed
        # if that is desired (just set 'c_min=lo' at the end of above loop).
        
        # Convert our list of lists to a flat list
        result_tuples = sum(results, start=[])
        
        if len(result_tuples) == 0:
            print('No transformations were found for this sample, moving to the next one')

        else:
            # make pd df
            result_df = pd.DataFrame.from_records(result_tuples, columns=[
                'Feature_X', 'Feature_Y', 'Difference',
                'Group', 'Transformation', 'Formula', 'mf'])

            result_df['SampleID'] = sample

            print('   Saving results')
            filename = os.path.join(path, 'transformations_' + sample + '.csv')
            result_df.to_csv(filename, index=False, float_format='%0.6f')

            # Compile counts
            result_counts = pd.DataFrame(
                result_df.groupby(
                    ['SampleID', 'Group', 'Transformation', 'Formula']
                ).size().reset_index(name='Counts'))

            total_transformations = sum(result_counts['Counts'])

            result_counts['Perc_Counts'] = (result_counts['Counts']
                                            / total_transformations)
            
            result_counts = result_counts.sort_values(by="Counts")

            # Save final_counts
            filename = os.path.join(path, 'counts_' + sample + '.csv')
            result_counts.to_csv(filename, index=False, float_format='%0.6g')
        i = i + 1

    print("Done!")
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
    summary_transf.to_csv(filename, index=False, float_format='%.6f')
    return


# --------------------------------------------------
def get_node_table(df, path):
    """Create a node table for the transformation networks"""

    node_table = df[['Mass', 'C', 'H', 'O', 'N', 'S', 'P', 'OC', 'HC', 'NOSC',
                     'GFE', 'Class', 'MolecularFormula', 'El_comp']
                    ].drop_duplicates('Mass')

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

    for ii, file in enumerate(edge_files):
        edge_table = pd.read_csv(file)
        edge_table['Feature_X'] = round(edge_table['Feature_X'], 4)
        edge_table['Feature_X'] = edge_table['Feature_X'].astype(str)
        edge_table['Feature_Y'] = round(edge_table['Feature_Y'], 4)
        edge_table['Feature_Y'] = edge_table['Feature_Y'].astype(str)
        # Renaming columns for clarity
        edge_table = edge_table.rename(columns={'Feature_X': 'source', 'Feature_Y': 'target'})
        print(f'[{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p")}]'
              f'{ii}\\{len(edge_files)}\tCreating network for sample: '
              f'{edge_table.SampleID.unique()[0]}'
              f'{edge_table.SampleID.iloc[0]}')

        p4c.create_network_from_data_frames(None, edge_table,
                                            source_id_list='source',
                                            target_id_list='target',
                                            interaction_type_list='Transformation',
                                            title=edge_table.SampleID.unique()[0],
                                            collection='Transformation networks')
        with pd.option_context('mode.chained_assignment', None):
            # Load node table with data information
            p4c.load_table_data(node_table, data_key_column='id')

        if (ii == 0):
            # These mappings only need to be set the first time; they slow things
            # down quite a bit if we do them for every network.
            # Coloring the table based on compound classes and transformation groups
            p4c.set_node_shape_default('ELLIPSE')
            p4c.set_node_color_mapping('Class', mol_classes, node_colors, mapping_type='d')


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

        # Saving a partial file with network stats in case the job did not finish
        temp_network_stats = pd.DataFrame(network_stats)
        filename = os.path.join(path, 'temp_network_summary_statistics.csv')
        temp_network_stats.to_csv(filename, index=False)

        time.sleep(5) # Ideally, this isn't necessary, but on slower machines we
        # have reports of cytoscape crashing if not here.
        

    network_stats = pd.DataFrame(network_stats)
    filename = os.path.join(path, 'network_summary_statistics.csv')
    network_stats.to_csv(filename, index=False)

    # Removing partial network summary statistics
    os.remove(os.path.join(path, 'temp_network_summary_statistics.csv'))

    return
