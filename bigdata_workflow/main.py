import logging
import re
import os
import subprocess
import shutil
import pandas as pd
from . import utils
logger = logging.getLogger(__name__)

def process_folder(args):
    infolder = os.path.abspath(args['indir'])
    modes = utils.parse_mode(args['mode'])
    fmods_val = args['fixed_mods']
    path_to_output_wild_reversed_fasta = args['wdatabase']
    path_to_outout_wild_and_target_peptides_fasta = args['cdatabase']
    if 1 in modes:
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.raw'):
                infile_mzml = os.path.splitext(infile)[0] + '.mzML'
                if args['overwrite'] or not utils.file_exist_and_nonempty(infile_mzml):
                    subprocess.run(["msconvert.exe",
                    infile,
                    '--singleThreaded',
                    '--mzML',
                    '--filter',
                    '"peakPicking true 1-"',
                    '--filter',
                    '"MS2Deisotope"',
                    '--filter',
                    '"zeroSamples removeExtra"',
                    '-o',
                    infolder             
                    ])

    if 2 in modes:
        # Wild search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.mzml'):

                pepxml_tmp = os.path.splitext(infile)[0] + '_identipy.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'

                if args['overwrite'] or not utils.file_exist_and_nonempty(pepxml_wild):


                    enz = "'[RK]|{P}'"
                    mc = 1
                    utils.run_identipy(infile, path_to_output_wild_reversed_fasta, enz, fmods_val, mc, dino=True)
                    shutil.move(pepxml_tmp, pepxml_wild)
                    utils.run_scavager(pepxml_wild, path_to_output_wild_reversed_fasta)

    if 3 in modes:
        # Variant search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('_identipy.mgf'):

                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_variant = pepxml_tmp.split('.pep.xml')[0] + '_variant.pep.xml'

                if args['overwrite'] or not utils.file_exist_and_nonempty(pepxml_variant):

                    enz = "'{X}|{X}'"
                    mc = 0
                    utils.run_identipy(infile, path_to_outout_wild_and_target_peptides_fasta, enz, fmods_val, mc, dino=False)
                    shutil.move(pepxml_tmp, pepxml_variant)
                    utils.run_scavager(pepxml_variant)

    if 4 in modes:
        # Brute force search
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('_identipy.mgf'):

                pepxml_tmp = os.path.splitext(infile)[0] + '.pep.xml'
                pepxml_wild = pepxml_tmp.split('.pep.xml')[0] + '_wild.pep.xml'
                pepxml_bruteforce = pepxml_tmp.split('.pep.xml')[0] + '_bruteforce.pep.xml'

                if args['overwrite'] or not utils.file_exist_and_nonempty(pepxml_bruteforce):
                    path_to_top100_wild_fasta = utils.make_top100_fasta(pepxml_wild, path_to_output_wild_reversed_fasta)
                    enz = "'[RK]|{P}'"
                    mc = 0
                    utils.run_identipy(infile, path_to_top100_wild_fasta, enz, fmods_val, mc, dino=False, snp=True)
                    shutil.move(pepxml_tmp, pepxml_bruteforce)
                    utils.run_scavager(pepxml_bruteforce)


    if 5 in modes:
        # Prepare output variant tables
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('_identipy.mgf'):
                table_variant = os.path.splitext(infile)[0] + '_variant_PSMs_full.tsv'
                table_bruteforce = os.path.splitext(infile)[0] + '_bruteforce_PSMs_full.tsv'
                table_output_variant = os.path.splitext(infile)[0] + '_variant_final.tsv'

                brute_counter_percent = utils.get_brute_counts(table_bruteforce)

                df_variants = utils.get_filtered_variants(table_variant, brute_counter_percent)

                df_variants.to_csv(path_or_buf = table_output_variant, sep='\t', index=False)

    if 6 in modes:
        # Prepare output variant table for whole folder
        flag = 1
        folder_name = os.path.basename(os.path.normpath(infolder))
        table_final = os.path.join(infolder, folder_name + '_variants.tsv')
        for infilename in os.listdir(infolder):
            infile = os.path.join(infolder, infilename)
            if infile.lower().endswith('.raw'):
                fn = os.path.splitext(infile)[0]
                table_output_variant = os.path.splitext(infile)[0] + '_identipy_variant_final.tsv'
                df1 = pd.read_csv(table_output_variant, sep='\t')
                df1['foldername'] = folder_name
                df1['filename'] = os.path.basename(fn)
                if flag:
                    dfc = df1.copy()
                    flag = 0
                else:
                    dfc = dfc.append(df1)
                    dfc.reset_index(inplace=True, drop=True)

        if not flag:
            cols = dfc.columns.tolist()
            cols.remove('foldername')
            cols.insert(1, 'foldername')
            cols.remove('filename')
            cols.insert(2, 'filename')
            cols.remove('gene')
            cols.insert(3, 'gene')
            cols.remove('aach')
            cols.insert(4, 'aach')
            dfc = dfc[cols]
            dfc.to_csv(path_or_buf = table_final, sep='\t', index=False)
