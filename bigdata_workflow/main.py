import logging
import re
import os
import subprocess
import shutil
import pandas as pd
import numpy as np
from . import utils
logger = logging.getLogger(__name__)


def process_folder(args):
    logger.info('Starting data analysis...')
    infolder = os.path.abspath(args['indir'])
    path_to_cfg = args['config']
    modes = utils.parse_mode(args['mode'])
    # fmods_val = args['fixed_mods']
    path_to_output_wild_reversed_fasta = args['wdatabase']
    path_to_outout_wild_and_target_peptides_fasta = args['cdatabase']
    path_to_genemap = args['genemap']
    path_to_cormap = args['cormap']
    if 1 in modes:
        total_cnt = 0
        for root, _, files in os.walk(infolder):
            for infilename in files:
                infile = os.path.join(root, infilename)
                if infile.lower().endswith('.raw'):
                    total_cnt += 1
        logger.info(
            '%s raw files are found, starting convertion...',
            total_cnt
            )
        for root, _, files in os.walk(infolder):
            for infilename in files:
                infile = os.path.join(root, infilename)
                if infile.lower().endswith('.raw') and utils.file_exist_and_nonempty(infile):
                    # infile_mzml = os.path.splitext(infile)[0] + '.mzML'
                    infile_mzml = os.path.join(infolder, os.path.splitext(infilename)[0] + '.mzML')
                    if args['overwrite'] or \
                            not utils.file_exist_and_nonempty(infile_mzml):
                        subprocess.run(
                            [
                                "msconvert.exe",
                                infile,
                                '--mzML',
                                '--filter',
                                '"peakPicking true 1-"',
                                '--filter',
                                '"MS2Deisotope"',
                                '--filter',
                                '"zeroSamples removeExtra"',
                                '--filter',
                                'threshold absolute 1 most-intense',
                                '-o',
                                infolder
                            ]
                        )
                    # logger.info('%d%d')

    if 2 in modes:
        # Wild search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.mzml'):

                pepxml_tmp = os.path.splitext(infile)[0] + '_identipy.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'

                if args['overwrite'] or \
                        not utils.file_exist_and_nonempty(pepxml_wild):
                    enz = "'[RK]|{P}'"
                    mc = 1
                    utils.run_identipy(
                        infile,
                        path_to_cfg,
                        path_to_output_wild_reversed_fasta,
                        enz,
                        # fmods_val,
                        mc,
                        dino=True)
                    shutil.move(pepxml_tmp, pepxml_wild)
                    utils.run_scavager(
                        pepxml_wild,
                        path_to_output_wild_reversed_fasta)

    if 3 in modes:
        # Variant search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            # print(infile)
            if infile.lower().endswith('_identipy.mgf'):

                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_variant = pepxml_tmp.split('.pep.xml')[0] + \
                    '_variant.pep.xml'
                # print(infile, utils.file_exist_and_nonempty(pepxml_variant))

                if args['overwrite'] or \
                        not utils.file_exist_and_nonempty(pepxml_variant):

                    enz = "'{X}|{X}'"
                    mc = 0
                    utils.run_identipy(
                        infile,
                        path_to_cfg,
                        path_to_outout_wild_and_target_peptides_fasta,
                        enz,
                        # fmods_val,
                        mc,
                        dino=False)
                    shutil.move(pepxml_tmp, pepxml_variant)
                    utils.run_scavager(pepxml_variant)

    if 4 in modes:
        # Brute force search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('_identipy.mgf'):

                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'
                pepxml_bruteforce = pepxml_tmp.split('.pep.xml')[0] + \
                    '_bruteforce.pep.xml'

                if args['overwrite'] or \
                        not utils.file_exist_and_nonempty(pepxml_bruteforce):
                    path_to_top100_wild_fasta = utils.make_top100_fasta(
                        pepxml_wild,
                        path_to_output_wild_reversed_fasta)
                    if path_to_top100_wild_fasta:
                        enz = "'[RK]|{P}'"
                        mc = 0
                        utils.run_identipy(
                            infile,
                            path_to_cfg,
                            path_to_top100_wild_fasta,
                            enz,
                            # fmods_val,
                            mc,
                            dino=False,
                            snp=True)
                        shutil.move(pepxml_tmp, pepxml_bruteforce)
                        utils.run_scavager(pepxml_bruteforce)
    if 5 in modes:
        # Prepare output variant tables
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('_identipy.mgf'):
                try:
                    table_variant = os.path.splitext(infile)[0] + \
                        '_variant_PSMs_full.tsv'
                    table_bruteforce = os.path.splitext(infile)[0] + \
                        '_bruteforce_PSMs_full.tsv'
                    table_output_variant = os.path.splitext(infile)[0] + \
                        '_variant_final.tsv'

                    brute_counter_percent = utils.get_brute_counts(
                        table_bruteforce)

                    df_variants = utils.get_filtered_variants(
                        table_variant,
                        brute_counter_percent)

                    df_variants.to_csv(
                        path_or_buf=table_output_variant,
                        sep='\t',
                        index=False)
                except:
                    print('Missing tsv files')

    if 6 in modes:
        # Prepare output variant table for whole folder
        flag = 1
        folder_name = os.path.basename(os.path.normpath(infolder))
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.mzml'):
                fn = os.path.splitext(infile)[0]
                table_output_variant = os.path.splitext(infile)[0] + \
                    '_identipy_variant_final.tsv'
                try:
                    df1 = pd.read_csv(table_output_variant, sep='\t')
                    df1['foldername'] = folder_name
                    df1['filename'] = os.path.basename(fn)
                    if flag:
                        dfc = df1.copy()
                        flag = 0
                    else:
                        dfc = dfc.append(df1)
                        dfc.reset_index(inplace=True, drop=True)
                except:
                    print('Missing tsv files')

        if not flag:
            dfc = utils.get_final_table(dfc)
            dfc.to_csv(path_or_buf=table_final, sep='\t', index=False)

    if 7 in modes:
        
        if path_to_genemap:
            cos_map = dict()
            for l in open(path_to_genemap, 'r'):
                cosmname, genename = l.strip().split('\t')
                cos_map[cosmname] = genename
        else:
            cos_map = False

        if path_to_cormap:
            df0 = pd.read_csv(path_to_cormap)
            df0 = df0.sort_values(by='correlation', ascending=False)
            df0 = df0.drop_duplicates(subset='peptides')
            cor_dict = {}
            for pep, cor in df0[['peptides', 'correlation']].values:
                cor_dict[pep] = float(cor)
        else:
            cor_dict = False


        # Process output variant table
        folder_name = os.path.basename(os.path.normpath(infolder))
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        table_final_output = os.path.join(infolder, folder_name + '_final.tsv')


        df1 = pd.read_table(table_final)
        df1['group'] = 'common_group'#df1['filename'].apply(lambda x: x.split('_')[0].lower())
        dfo = df1.copy()
        dfo = dfo[['peptide', 'filename', 'database', 'gene', 'aach', 'group', 'MS1Intensity', 'brute_count', 'length', 'PEP']]
        dfo['MS1Intensity'] = np.log10(dfo['MS1Intensity'])
        dfo['total_psms'] = dfo.groupby('peptide')['peptide'].transform('count')
        # dfo['total_files'] = dfo.groupby(('peptide', 'filename'))['peptide'].transform('count')
        dfo['PSM count'] = dfo.groupby(['peptide', 'filename'])['peptide'].transform('count')
        dfo['total_intensity_max'] = dfo.groupby('peptide')['MS1Intensity'].transform('max')
        df1 = dfo.reset_index(level=0).rename(columns={'level_0': 'group'})
        if cos_map:
            df1['gene'] = df1['gene'].apply(utils.remap_gene, cos_map=cos_map)
        df = df1.groupby(['group', 'peptide']).agg({'aach': 'first', 'gene': 'first', 'database': 'first',
                                            'filename': lambda x: len(set(x)),
                                            'PSM count': 'sum',
                                            'MS1Intensity': 'max',
                                            'brute_count': 'max',
                                            'length': 'first',
                                            'PEP': 'min'})
        df = df.rename(columns={"filename": "filecount"})
        dfd = df.reset_index(level=0)
        cols_to_save = ['aach', 'gene', 'database', 'brute_count', 'length', 'PEP']
        info = dfd.loc[~dfd.index.duplicated(), cols_to_save]
        udf = df.unstack(level=0)
        gdf = udf.swaplevel(axis=1).sort_index(axis=1).loc[:, (slice(None), ['PSM count', 'MS1Intensity', 'filecount'])].fillna(0)#.astype(int)
        gdf.columns.names = (None, None)
        gdf[cols_to_save] = info[cols_to_save]
        # gdf['PSM count'] = gdf.loc[:, (slice(None), 'PSM count')]
        gdf['PSM count'] = gdf.loc[:, (slice(None), 'PSM count')].sum(axis=1)
        gdf['filecount'] = gdf.loc[:, (slice(None), 'filecount')].sum(axis=1)
        gdf = gdf.reset_index()
        if cor_dict:
            gdf['correlation'] = gdf['peptide'].apply(lambda x: cor_dict.get(x, 0))
        else:
            gdf['correlation'] = 0
        c = list(gdf.columns)
        order = {'peptide': 0, 'PSM count': 1, 'filecount': 2, 'aach': 3, 'database': 4, 'gene': 5, 'correlation': 6,
                'length': 7, 'brute_count': 8, 'PEP': 9}
        gdf = gdf[sorted(c, key=lambda x: order.get(x[0], 1e5))].copy()
        gdf = gdf.sort_values(by='filecount', ascending = False)
        gdf = gdf.replace([np.inf, -np.inf], 0)
        cols_to_int = {}
        for coln in gdf.columns:
            if 'PSM count' in coln or 'filecount' in coln:
                cols_to_int[coln] = 'int32'
        cols_to_int
        gdf = gdf.astype(cols_to_int)
        gdf.to_csv(path_or_buf=table_final_output, sep='\t', index=False, float_format='%.2f')

