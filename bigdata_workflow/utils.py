import os
import subprocess
import pandas as pd
from collections import Counter
from pyteomics import fasta, parser, auxiliary as aux
import ast


def rule_for_mutant(x):
    return all('mut' in z for z in x)

def rule_for_aa_change(x):
    return x[0].split(':')[2].split(',', 3)[-1].split(';')[0]

def get_aachange(zz):
    aachange = zz[0].split('_')[-1]
    aachange = aachange[0] + '>' + aachange.split('to')[-1][0]
    return aachange

def get_brute_counts(table_bruteforce):
    df3 = pd.read_csv(table_bruteforce, sep='\t')

    df3 = df3[~df3['decoy1']]
    df3['protein'] = df3['protein'].apply(lambda x: ast.literal_eval(x))
    df3 = df3[df3['protein'].apply(lambda z: not any(zz.endswith('_wild') for zz in z))]
    df3['aachange'] = df3['protein'].apply(get_aachange)
    df3_v_f = aux.filter(df3, key='ML score', is_decoy='decoy2', fdr=0.1)

    brute_counter = Counter()
    for z in df3_v_f.drop_duplicates(subset='peptide')[['aachange']].values:
        brute_counter[z[0]] += 1

    brute_counter_percent = {}
    koeff = sum(brute_counter.values())
    for k, v in brute_counter.items():
        brute_counter_percent[k] = float(v) / koeff * 100
    return brute_counter_percent

def get_filtered_variants(table_variant, brute_counter_percent):
    # Read all PSMs
    df2_v = pd.read_csv(table_variant, sep='\t')

    df2_v['protein'] = df2_v['protein'].apply(lambda x: ast.literal_eval(x))
    df2_v = df2_v[df2_v['protein'].apply(rule_for_mutant)]
    df2_v = df2_v[~df2_v['decoy1']]

    df2_v['len'] = df2_v['peptide'].apply(lambda z: len(z))
    df2_v['MS1Intensity'] = df2_v['MS1Intensity'].astype(int)
    df2_v['seq_IL'] = df2_v['peptide'].apply(lambda z: z.replace('L', 'I'))
    df2_v['AAchange'] = df2_v['protein'].apply(lambda z: z[0].split(':')[2].split(',', 3)[-1].split(';')[0])
    df2_v['database'] = df2_v['protein'].apply(lambda z: z[0].split(':')[0].split('_')[-1])
    df2_v['aach'] = df2_v['AAchange'].apply(lambda z: z.split(',')[0])
    df2_v['brute_count'] = df2_v['aach'].apply(lambda z: brute_counter_percent.get(z, 0))
    df2_v['mc'] = df2_v['peptide'].apply(lambda z: parser.num_sites(z, parser.expasy_rules['trypsin']))
    df2_v = df2_v[df2_v['mc'] == 0] ### May be not to filter to 0 missed cleavages ???
    df2_v = df2_v[df2_v['brute_count'] <= 5.0] ### Another rule instead ???
    df2_v = df2_v[df2_v['AAchange'].apply(lambda z: 'L>I' not in z and 'I>L' not in z)]
    df2_v_f = aux.filter(df2_v, key='PEP', is_decoy='decoy2', fdr=0.05)
    df2_v_f['gene'] = df2_v_f['protein'].apply(lambda z: z[0].split(':')[1].split(',')[0])
    return df2_v_f


#Reversing of variant peptides
def smart_reverse(ts):
    ns = fasta.reverse(ts, keep_nterm=True, keep_cterm=True)
    return ns

def make_top100_fasta(pepxml_wild, path_to_fasta):
    proteins_wild = pepxml_wild.split('.pep.xml')[0] + '_protein_groups.tsv'
    path_to_top100_wild_fasta = pepxml_wild.split('.pep.xml')[0] + '_top100wild_reversed.fasta'
    df0 = pd.read_csv(proteins_wild, sep='\t')
    cnt_prots = Counter()
    for z in df0[['dbname']].values:
        cnt_prots[z[0].split()[0]] += 1

    top_wild_prots = []
    abund_prots = set([z[0] for z in cnt_prots.most_common()[:100]])
    for p in fasta.read(path_to_fasta):
        if p[0].split()[0] in abund_prots:
            rev_seq = smart_reverse(p[1])
            top_wild_prots.append(p)
            top_wild_prots.append(('DECOY_' + p[0], rev_seq))

    fasta.write(top_wild_prots, output=open(path_to_top100_wild_fasta, 'w')).close()
    return path_to_top100_wild_fasta


def parse_mode(mode_str):
    if '-' in mode_str:
        if mode_str.startswith('-'):
            mode_str = '1' + mode_str
        elif mode_str.endswith('-'):
            mode_str = mode_str + '9'
        start_idx, end_idx = (int(md) for md in mode_str.split('-'))
        return set(range(start_idx, end_idx+1, 1))
    else:
        return set(int(md) for md in mode_str.strip().split(','))

def file_exist_and_nonempty(infile):
    return os.path.isfile(infile) and os.stat(infile).st_size

def run_identipy(infile, path_to_fasta, enz, fmods_val, mc, dino=False, snp=False):
    array_for_subprocess = ["identipy",
    infile,
    '-db',
    path_to_fasta,
    '-at',
    '-ptol',
    '1',
    '-ftol',
    '0.05',
    '-e',
    enz,
    '-mc',
    str(mc),
    '-fmods',
    fmods_val,
    '-dyn',
    '100',
    '-maxp',
    '50',
    '-cmin',
    '1',
    '-cmax',
    '6',
    '-ime',
    '0',
    '-lmin',
    '6',
    '-deis',
    'no',
    ]
    if dino:
        array_for_subprocess.extend(['-dino', '/home/lab006/Dinosaur/Dinosaur-1.1.3.free.jar'])
    if snp:
        array_for_subprocess.extend(['-snp', '1'])
    subprocess.run(array_for_subprocess)
    return 

def run_scavager(pepxml_wild, path_to_fasta=False):
    array_for_subprocess = ["scavager",
    pepxml_wild,
    '-fdr',
    '5.0',
    ]
    if path_to_fasta:
        array_for_subprocess.extend(['-db', path_to_fasta])
    subprocess.run(array_for_subprocess)
    return 

def get_final_table(dfc):
    columns = [
        'peptide',
        'foldername',
        'filename',
        'database',
        'gene',
        'aach',
        'brute_count',
        'length',
        'PEP',
        'massdiff_ppm',
        'modifications',
        'MS1Intensity',
        'RT exp',
        'RT pred',
        'assumed_charge',
        'num_missed_cleavages',
        'ISOWIDTHDIFF',
        'protein',
        'protein_descr',
        'calc_neutral_pep_mass',
        'num_tol_term',
        'q',
        'ML score',
    ]
    return dfc[columns]