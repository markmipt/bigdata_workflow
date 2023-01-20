import logging
import re
import os
import subprocess
import shutil
import pandas as pd
import numpy as np
import tempfile
from pyteomics import auxiliary as aux
from scipy.stats import scoreatpercentile
from scipy.optimize import curve_fit
from scipy import exp

from .prositCorrelationUtils.PrositPipeline import PrositPipeline
from . import utils
logger = logging.getLogger(__name__)



def noisygaus(x, a, x0, sigma, b):
    return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

def calibrate_RT_gaus(bwidth, mass_left, mass_right, true_md):

    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]


    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), bwidth * 5, 1])
    mass_shift, mass_sigma = popt[1], abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]

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

                # pepxml_tmp = os.path.splitext(infile)[0] + '_identipy.pep.xml'
                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'

                if args['overwrite'] or \
                        not utils.file_exist_and_nonempty(pepxml_wild):

                    try:

                        enz = "'[RK]|{P}'"
                        mc = 1
                        utils.run_identipy(
                            infile,
                            path_to_cfg,
                            path_to_output_wild_reversed_fasta,
                            enz,
                            # fmods_val,
                            mc,
                            dino=False,
                            # dino=args['featurefinder'],
                            path_to_identipy=args['identipy'])
                        shutil.move(pepxml_tmp, pepxml_wild)
                        utils.run_scavager(
                            pepxml_wild,
                            path_to_output_wild_reversed_fasta, path_to_scavager=args['scavager'])


                    except:
                        logger.info('Cannot process file %s...', infile)
                        trash_dir = os.path.join(infolder, 'unprocessed/')

                        if not os.path.exists(trash_dir):
                            os.makedirs(trash_dir)

                        infile_to_trash = os.path.join(trash_dir, infilename)

                        shutil.move(infile, infile_to_trash)

    if 3 in modes:
        # Variant search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            # print(infile)
            # if infile.lower().endswith('_identipy.mgf'):
            if infile.lower().endswith('.mzml'):

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
                        dino=False,
                        path_to_identipy=args['identipy'])
                    shutil.move(pepxml_tmp, pepxml_variant)
                    utils.run_scavager(pepxml_variant, path_to_scavager=args['scavager'])

                        

    if 4 in modes:
        # Brute force search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            # if infile.lower().endswith('_identipy.mgf'):
            if infile.lower().endswith('.mzml'):

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
                            snp=True,
                            path_to_identipy=args['identipy'])
                        shutil.move(pepxml_tmp, pepxml_bruteforce)
                        utils.run_scavager(pepxml_bruteforce, path_to_scavager=args['scavager'])
    if 5 in modes:
        # Prepare output variant tables
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            # if infile.lower().endswith('_identipy.mgf'):
            if infile.lower().endswith('.mzml'):
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
                    print('Missing tsv files for mode #5')

    if 6 in modes:
        # Prepare output variant table for whole folder
        flag = 1
        folder_name = os.path.basename(os.path.normpath(infolder))
        print(folder_name)
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.mzml'):
                fn = os.path.splitext(infile)[0]
                # table_output_variant = os.path.splitext(infile)[0] + \
                #     '_identipy_variant_final.tsv'
                table_output_variant = os.path.splitext(infile)[0] + \
                    '_variant_final.tsv'
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
                    print('Missing tsv files for mode #6')

        # print('HERE', flag)
        if not flag:
            dfc = utils.get_filtered_variants_fdr_only(dfc, 0.05)
            # print(dfc)
            dfc = utils.get_final_table(dfc)
            dfc.to_csv(path_or_buf=table_final, sep='\t', index=False)

    if 7 in modes:
        folder_name = os.path.basename(os.path.normpath(infolder))

        prosit_path = args['prosit']
        MODEL_SPECTRA = args['prosit_model']
        MODEL_IRT = args['prosit_irt_model']
        prosit_pipeline = PrositPipeline()
        prosit_pipeline.main_prosit(infolder, folder_name, prosit_path, MODEL_SPECTRA, MODEL_IRT)

    if 8 in modes:
        deeplc_path = args['deeplc']
        deeplc_path = deeplc_path.strip()

        folder_name = os.path.basename(os.path.normpath(infolder))
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        table_wilds = os.path.join(infolder, folder_name + '_wilds.tsv')
        df1 = pd.read_table(table_final)
        df2 = pd.read_table(table_wilds)

        df1['mods_for_deepLC'] = df1.apply(utils.mods_for_deepLC, axis=1)
        df2['mods_for_deepLC'] = df2.apply(utils.mods_for_deepLC, axis=1)

        df3 = df1.append(df2, ignore_index=True)


        outtrain = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outcalib = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outres = tempfile.NamedTemporaryFile(suffix='.txt', mode='w')
        outres_name = outres.name
        outres.close()
        ns = df3['peptide'].values
        nr = df3['RT exp'].values
        nm = df3['mods_for_deepLC'].values
        print(df3['mods_for_deepLC'])
        print('Peptides used for RT prediction: %d' % (len(ns), ))
        ns2 = df2['peptide'].values
        nr2 = df2['RT exp'].values
        nm2 = df2['mods_for_deepLC'].values
        calibset = set(ns2)

        outtrain.write('seq,modifications,tr\n')
        for seq, RT, mods_tmp in zip(ns2, nr2, nm2):
            # mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outtrain.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
        outtrain.flush()

        outcalib.write('seq,modifications,tr\n')
        for seq, RT, mods_tmp in zip(ns, nr, nm):
            # mods_tmp = '|'.join([str(idx+1)+'|Carbamidomethyl' for idx, aa in enumerate(seq) if aa == 'C'])
            outcalib.write(seq + ',' + str(mods_tmp) + ',' + str(RT) + '\n')
        outcalib.flush()

        subprocess.call([deeplc_path, '--file_pred', outcalib.name, '--file_cal', outtrain.name, '--file_pred_out', outres_name])
        pepdict = dict()
        train_RT = []
        train_seq = []
        for x in open(outres_name).readlines()[1:]:
            _, seq, _, RTexp, RT = x.strip().split(',')
            pepdict[seq] = float(RT)

            if seq in calibset:
                train_seq.append(seq)
                train_RT.append(float(RTexp))


        train_RT = np.array(train_RT)
        RT_pred = np.array([pepdict[s] for s in train_seq])

        rt_diff_tmp = RT_pred - train_RT
        RT_left = -min(rt_diff_tmp)
        RT_right = max(rt_diff_tmp)

        try:
            start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 100
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
        except:
            start_width = (scoreatpercentile(rt_diff_tmp, 95) - scoreatpercentile(rt_diff_tmp, 5)) / 50
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(start_width, RT_left, RT_right, rt_diff_tmp)
        if np.isinf(covvalue):
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(0.1, RT_left, RT_right, rt_diff_tmp)
        if np.isinf(covvalue):
            XRT_shift, XRT_sigma, covvalue = calibrate_RT_gaus(1.0, RT_left, RT_right, rt_diff_tmp)
        print('Calibrated RT shift: ', XRT_shift)
        print('Calibrated RT sigma: ', XRT_sigma)

        aa, bb, RR, ss = aux.linear_regression(RT_pred, train_RT)

        df1['RT pred DeepLC'] = df1['peptide'].apply(lambda x: pepdict[x])
        df1['RT diff DeepLC'] = df1['RT exp'] - df1['RT pred DeepLC']
        df1['RT Z abs diff DeepLC'] = df1['RT diff DeepLC'].abs() / XRT_sigma
        df1.to_csv(path_or_buf=table_final, sep='\t', index=False)

    if 9 in modes:
        
        if path_to_genemap:
            cos_map = dict()
            for l in open(path_to_genemap, 'r'):
                cosmname, genename = l.strip().split('\t')
                cos_map[cosmname] = genename
        else:
            cos_map = False

        # if path_to_cormap:
        #     df0 = pd.read_csv(path_to_cormap)
        #     df0 = df0.sort_values(by='correlation', ascending=False)
        #     df0 = df0.drop_duplicates(subset='peptides')
        #     cor_dict = {}
        #     for pep, cor in df0[['peptides', 'correlation']].values:
        #         cor_dict[pep] = float(cor)
        # else:
        #     cor_dict = False


        # Process output variant table
        folder_name = os.path.basename(os.path.normpath(infolder))
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        table_final_output = os.path.join(infolder, folder_name + '_final.tsv')

        try:
            table_wilds = os.path.join(infolder, folder_name + '_wilds.tsv')
            df2 = pd.read_table(table_wilds)
            wild_threshold = scoreatpercentile(df2['correlation'], 5)
        except:
            wild_threshold = -100

        df1 = pd.read_table(table_final)
        # df1 = df1.rename(columns={"modified_sequence": "peptide"})
        df1['group'] = 'common_group'#df1['filename'].apply(lambda x: x.split('_')[0].lower())
        dfo = df1.copy()

        if 'correlation' not in dfo.columns:
            dfo['correlation'] = -1
        if 'RT Z abs diff DeepLC' not in dfo.columns:
            dfo['RT Z abs diff DeepLC'] = -1

        dfo = dfo[['peptide', 'filename', 'database', 'gene', 'aach', 'group', 'MS1Intensity', 'brute_count', 'length', 'PEP', 'correlation', 'comment', 'RT Z abs diff DeepLC']]
        dfo['MS1Intensity'] = np.log10(dfo['MS1Intensity'])
        dfo['PSM count'] = 1#dfo.groupby('peptide')['peptide'].transform('count')
        # dfo['total_files'] = dfo.groupby(('peptide', 'filename'))['peptide'].transform('count')
        # dfo['PSM count'] = dfo.groupby(['peptide', 'filename'])['peptide'].transform('count')
        # print(dfo[dfo['peptide'] == 'GEGEPCGGGGAGGGYCAPGMECVK'])
        dfo['total_intensity_max'] = dfo.groupby('peptide')['MS1Intensity'].transform('max')


        dfo['comment'] = ''
        w_95 = wild_threshold
        dfo.loc[dfo['correlation'] < w_95, 'comment'] = 'unreliable'
        dfo.loc[dfo['RT Z abs diff DeepLC'] > 3.0, 'comment'] = 'unreliable'
        # dfo = dfo[dfo['comment'] != 'unreliable']
        if len(dfo):

            df1 = dfo.reset_index(level=0).rename(columns={'level_0': 'group'})
            # # df1 = 
            # print(df1[pd.isna(df1['gene'])])
            df1.loc[pd.isna(df1['gene']), 'gene'] = 'NAgene'
            print(df1[pd.isna(df1['gene'])])
            if cos_map:
                df1['gene'] = df1['gene'].apply(utils.remap_gene, cos_map=cos_map)
            df = df1.groupby(['group', 'peptide']).agg({'aach': 'first', 'gene': 'first', 'database': 'first',
                                                'filename': lambda x: len(set(x)),
                                                'PSM count': 'sum',
                                                'MS1Intensity': 'max',
                                                'brute_count': 'max',
                                                'length': 'first',
                                                'PEP': 'min',
                                                'comment': 'min',
                                                'correlation': 'max',
                                                'RT Z abs diff DeepLC': 'min'})
            df = df.rename(columns={"filename": "filecount"})
            dfd = df.reset_index(level=0)
            cols_to_save = ['aach', 'gene', 'database', 'brute_count', 'length', 'PEP', 'correlation', 'comment', 'RT Z abs diff DeepLC']
            info = dfd.loc[~dfd.index.duplicated(), cols_to_save]
            udf = df.unstack(level=0)
            gdf = udf.swaplevel(axis=1).sort_index(axis=1).loc[:, (slice(None), ['PSM count', 'MS1Intensity', 'filecount'])].fillna(0)#.astype(int)
            gdf.columns.names = (None, None)
            gdf[cols_to_save] = info[cols_to_save]
            # gdf['PSM count'] = gdf.loc[:, (slice(None), 'PSM count')]
            gdf['PSM count'] = gdf.loc[:, (slice(None), 'PSM count')].sum(axis=1)
            gdf['filecount'] = gdf.loc[:, (slice(None), 'filecount')].sum(axis=1)
            gdf = gdf.reset_index()
            # if cor_dict:
            #     gdf['correlation'] = gdf['peptide'].apply(lambda x: cor_dict.get(x, 0))
            # else:
            #     gdf['correlation'] = 0
            c = list(gdf.columns)

            # w_95 = scoreatpercentile(df2['correlation'], 5)
            # gdf.loc[gdf['correlation'] < w_95, 'comment'] = 'unreliable'
            # gdf.loc[gdf['RT Z abs diff DeepLC'] > 3.0, 'comment'] = 'unreliable'
            # gdf.loc[gdf['brute_count'] > 1.0, 'comment'] = 'unreliable'

            order = {'peptide': 0, 'PSM count': 1, 'filecount': 2, 'aach': 3, 'database': 4, 'gene': 5, 'correlation': 6,
                    'length': 7, 'brute_count': 8, 'PEP': 9, 'comment': 10, 'RT Z abs diff DeepLC': 11}
            gdf = gdf[sorted(c, key=lambda x: order.get(x[0], 1e5))].copy()
            gdf = gdf.replace([np.inf, -np.inf], 0)
            gdf['comment'] = gdf['comment'].fillna('')
            gdf = gdf.sort_values(by=['comment', 'filecount'], ascending = [True, False])
            cols_to_int = {}
            for coln in gdf.columns:
                if 'PSM count' in coln or 'filecount' in coln:
                    cols_to_int[coln] = 'int32'
            gdf = gdf.astype(cols_to_int)
            gdf.to_csv(path_or_buf=table_final_output, sep='\t', index=False, float_format='%.2f')
        else:
            pd.DataFrame(columns=['peptide', 'PSM count', 'filecount', 'aach', 'database', 'gene',
            'correlation', 'length', 'brute_count', 'PEP', 'comment',
            'RT Z abs diff DeepLC', 'common_group', 'common_group.1',
            'common_group.2']).to_csv(path_or_buf=table_final_output, sep='\t', index=False, float_format='%.2f')


    if 10 in modes:
        # Wild search union

        list_of_pepxml = []
        folder_name = os.path.basename(os.path.normpath(infolder))

        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.mzml'):

                # pepxml_tmp = os.path.splitext(infile)[0] + '_identipy.pep.xml'
                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'

                if utils.file_exist_and_nonempty(pepxml_wild):
                    list_of_pepxml.append(pepxml_wild)

        if len(list_of_pepxml):
            utils.run_scavager_union(list_of_pepxml, folder_name, path_to_output_wild_reversed_fasta, path_to_scavager=args['scavager'])

                # if args['overwrite'] or \
                #         not utils.file_exist_and_nonempty(pepxml_wild):
                    # enz = "'[RK]|{P}'"
                    # mc = 1
                    # utils.run_identipy(
                    #     infile,
                    #     path_to_cfg,
                    #     path_to_output_wild_reversed_fasta,
                    #     enz,
                    #     # fmods_val,
                    #     mc,
                    #     dino=True)
                    # shutil.move(pepxml_tmp, pepxml_wild)
                    # utils.run_scavager(
                    #     pepxml_wild,
                    #     path_to_output_wild_reversed_fasta)

