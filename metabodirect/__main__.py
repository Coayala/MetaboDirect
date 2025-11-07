# __main__.py
"""
Author : Christian Ayala <cayalaortiz@arizona.edu>
Date   : 2023-08-23
Purpose: Program to run MetaboDirect scripts

To get help, use metabodirect -h or visit the website:
https://github.com/Coayala/MetaboDirect
"""

import os
import time
import datetime
import sys
import shutil
import json
import pandas as pd
import py4cytoscape as p4c
from loguru import logger
from math import ceil
from metabodirect import get_args, preprocessing, diagnostics, r_control, transformations


# --------------------------------------------------
def main():
    """Main body to run all MetaboDirect scripts"""
    start_time = time.perf_counter()
    args = get_args.get_args()

    logger.remove()
    logger.add('log_{time:MMDDYY-HH_mm_ss}.log',
               format='[{time:MM/DD/YYYY HH:mm}] {level: <8}| <lvl>{message}</lvl>')

    logger.add(sys.stdout,
               format='[{time:MM/DD/YYYY HH:mm}] {level: <8}| <lvl>{message}</lvl>')

    # Starting the pipeline
    logger.opt(colors=True)
    logger.level("PROCESS", no=38, color="<magenta>")
    logger.log('PROCESS', 'WELCOME TO METABODIRECT\n')
    logger.info('Analysis starting on {}',
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"))
    logger.success('Run was initiated with the command: {}', sys.argv)
    logger.info('Results will be saved in directory: {}',
                os.path.abspath(args.outdir))

    logger.info('Data file is {}', args.data_file)
    logger.info('Metadata file is {}\n', args.metadata_file)

    # Starting the data pre-processing step
    logger.log('PROCESS', 'Start data pre-processing\n')

    df = pd.read_csv(args.data_file)
    metadata = pd.read_csv(args.metadata_file)
    list_dir = preprocessing.make_directories(args.outdir)

    # Filtering samples based on metadata
    if args.filter_by:
        logger.info('Filtering samples based on {} = {}',
                    args.filter_by[0], args.filter_by[1])
        df, metadata = preprocessing.sample_filtering(df, metadata,
                                                      args.filter_by,
                                                      path=list_dir[0],
                                                      args=args)
        logger.info('Filtered tables can be found in the directory: {}',
                    os.path.abspath(list_dir[0]))
    else:
        logger.warning('Option -f not detected, all samples will be used')

    # Filtering peaks
    df = preprocessing.data_filtering(df,
                                      mass_filter=args.mass_filter,
                                      peak_filter=args.peak_filter,
                                      error_filter=args.error_filter,
                                      args=args)

    # Calculating indices
    logger.info('Calculating thermodynamic indexes')
    df = preprocessing.thermo_idx_and_classes(df, args)

    # Normalizing data
    logger.info('Normalizing data using the {} method', args.norm_method)

    df, df_nonorm = preprocessing.data_normalization(df,
                                                     args,
                                                     args.norm_method,
                                                     args.norm_subset,
                                                     args.subset_parameter,
                                                     args.log_transform)

    # Remove peaks without formulas
    filename = os.path.join(list_dir[0], 'Report_processed_noNorm.csv')
    df_nonorm.to_csv(filename, index=False)

    df_nonorm = df_nonorm[df_nonorm['C'] > 0]
    filename = os.path.join(list_dir[0], 'Report_processed_noNorm_MolecFormulas.csv')
    df_nonorm.to_csv(filename, index=False)

    filename = os.path.join(list_dir[0], 'Report_processed.csv')
    df.to_csv(filename, index=False)

    df_formulas = df[df['C'] > 0]
    filename = os.path.join(list_dir[0], 'Report_processed_MolecFormulas.csv')
    df_formulas.to_csv(filename, index=False)
    logger.info('Report saved as: {}\n', filename)

    preprocessing.calculate_summaries(df_formulas, args, path=list_dir[0])
    matrix_features = preprocessing.get_matrix(df_formulas, args)
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

    logger.log('PROCESS', 'Data pre-processing finished\n')

    
    # Starting Data diagnostics ste[]
    logger.log('PROCESS', 'Starting data diagnostics\n')

    logger.info('Calculating number of peaks detected per sample')
    index_1 = diagnostics.peaks_per_sample(df_all, metadata, args.group, path=list_dir)
    logger.info('Calculating number of assigned molecular formulas per sample')
    index_2 = diagnostics.formula_per_sample(df, metadata, args.group, path=list_dir)
    logger.info("Calculating error distribution per group/s: {}", args.group)
    index_3 = diagnostics.error_per_group(df, args.group, path=list_dir)
    
    figure_index = {
        'analysis_module': '2. Diagnostics',
        'script_name': 'diagnostics.py',
        'last_run': str(datetime.datetime.now()),
        'figures': [
            index_1,
            index_2,
            index_3
        ]
    }
    
    with open(os.path.join(args.outdir, 'index.json'), "w") as json_file:
        json.dump(figure_index, json_file)

    logger.log('PROCESS', 'Data diagnostics finished\n')

    if not args.skip_analyses:

        # Starting data exploration step
        logger.log('PROCESS', 'Starting data exploration\n')

        # Copying script with functions
        functions_file = os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                      'R_scripts_templates/custom_functions.R')
        shutil.copy(functions_file, args.outdir)

        e_script_name = 'data_exploration_template.R' if len(args.group) == 1 else 'data_exploration_template_2_groups.R'

        d_explor_script = r_control.write_r_script(e_script_name,
                                                   outdir=list_dir[2],
                                                   metadata_file=args.metadata_file if not args.filter_by
                                                   else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                   groups=args.group)
        logger.info('Running R script: {}', d_explor_script)

        r_control.run_r(d_explor_script)
        logger.info('Find results and R script in the directory: {}',
                    os.path.abspath(list_dir[2]))

        # Annotating with KEGG if option was selected
        if args.kegg_annotation:
            logger.info('Starting annotation of molecular formulas using the KEGG database')
            kegg_annot_script = r_control.write_r_script('KEGG_annotation_template.R',
                                                         outdir=list_dir[2],
                                                         metadata_file=args.metadata_file if not args.filter_by
                                                         else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                         groups=args.group)
            r_control.run_r(kegg_annot_script)
            logger.info('Find results and R script in the directory: {}',
                        os.path.abspath(list_dir[2]))
        else:
            logger.warning(
                'KEGG annotation not selected. If you wish to perform a KEGG '
                'annotation run the script again using the '
                '"-k/--kegg_annotation" option')

        logger.log('PROCESS', 'Data exploration finished\n')

        # Starting chemodiversity analysis step
        logger.log('PROCESS', 'Starting chemodiversity analysis\n')

        data_chemodiversity_script = r_control.write_r_script('data_chemodiversity_template.R', outdir=list_dir[3],
                                                              metadata_file=args.metadata_file if not args.filter_by
                                                              else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                              groups=args.group)
        logger.info('Running R script: {}', data_chemodiversity_script)
        r_control.run_r(data_chemodiversity_script)
        logger.info('Find results and R script in the directory: {}\n',
                    os.path.abspath(list_dir[3]))

        logger.log('PROCESS', 'Chemodiversity analysis finished\n')

        # Starting Statistical analysis step
        logger.log('PROCESS', 'Starting statistical analysis\n')

        s_script_name = 'data_statistics_template.R' if len(args.group) == 1 else 'data_statistics_template_2_groups.R'
        data_statistics_script = r_control.write_r_script(s_script_name,
                                                          outdir=list_dir[4],
                                                          metadata_file=args.metadata_file if not args.filter_by
                                                          else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                          groups=args.group, norm_method=args.norm_method)
        logger.info('Running R script: {}', os.path.abspath(data_statistics_script))
        r_control.run_r(data_statistics_script)
        logger.info('Find results and R script in the directory: {}\n',
                    os.path.abspath(list_dir[4]))

        logger.log('PROCESS', 'Statistical analysis finished\n')
    else:
        logger.log('PROCESS', 'Skipping analysis steps.\n')

        
    # Starting data transformation networks step
    logger.log('PROCESS', 'Starting transformation networks analysis\n')

    if not args.calculate_transformations:
        logger.warning('Calculate transformations not selected.\n'
                       'If you wish to do calculate transformations based on biochemical key please set the option "-t"\n')
        end_time = time.perf_counter()
        logger.success('MetaboDirect finished running in {} minutes', ceil((end_time - start_time) / 60))
        sys.exit()

    # Calculating potential transformations
    keys = transformations.get_keys(os.path.join(os.path.split(os.path.realpath(__file__))[0],
                                                 'data',
                                                 'transf_key.csv')) if args.biochem_key == 'Default key' else \
        transformations.get_keys(args.biochem_key)
    transformations.calculate_transformations(df,
                                              keys,
                                              path=list_dir[6],
                                              err_thresh=args.transformation_threshold)
    transformations.summarize_transformations(path=list_dir[5])
    node_table = transformations.get_node_table(df, path=list_dir[5])

    logger.info('Finished to calculate transformatios, please find transformation files in the directory:{}',
                os.path.abspath(list_dir[6]))

    # Build transformation networks
    if not args.create_networks:
        logger.warning('Create networks not selected.\n'
                       'If you wish to create networks automatically please set the option -c.\n'
                       'Otherwise to create networks from previously calculated transformations using the "create_networks"'
                       'companion script\n')
        end_time = time.perf_counter()
        logger.success('MetaboDirect finished running in {} minutes', ceil((end_time - start_time) / 60))
        sys.exit()

    # Opening Cytoscape
    check = ''
    while check != 'You are connected to Cytoscape!':
        cytoscape = input('Please open Cytoscape and press the ENTER key [q for quit].')
        if cytoscape == '':
            try:
                check = p4c.cytoscape_ping()
            except (ValueError, Exception):
                print('Cytoscape is not open.')
        elif cytoscape == 'q':
            end_time = time.perf_counter()
            logger.success('MetaboDirect finished running in {} minutes', ceil((end_time - start_time) / 60))
            sys.exit()
        else:
            check = 'Cytoscape not open'

    # Start trasnformation networks
    logger.info('Starting network construction on Cytoscape')

    transformations.create_cytoscape_network(node_table, path=list_dir[5])

    network_stats_script = r_control.write_r_script('network_stats_template.R', outdir=list_dir[5],
                                                    metadata_file=args.metadata_file if not args.filter_by
                                                    else os.path.join(list_dir[0], 'filtered_metadata.csv'),
                                                    groups=args.group, norm_method=args.norm_method)

    logger.info('Running R script: {}', network_stats_script)
    r_control.run_r(network_stats_script)
    logger.info('Find results and R script in the directory: {}',
                os.path.abspath(list_dir[5]))

    logger.log('PROCESS', 'Transformation network analysis finished\n')

    # Finishing
    logger.success('Thanks for using MetaboDirect\n')

    end_time = time.perf_counter()
    logger.success('MetaboDirect finished running in {} minutes',
                   ceil((end_time - start_time) / 60))


# --------------------------------------------------
if __name__ == '__main__':
    main()
