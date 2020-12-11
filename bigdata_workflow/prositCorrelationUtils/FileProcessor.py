import pandas as pd
import os
from pyteomics import mzml
import ast

def oxidation_positions(raw):
    tmp = ast.literal_eval(raw)
    tmp = [z.split('@') for z in tmp]
    return set([int(z[1]) for z in tmp if z[0].startswith('147.')])

def new_sequence(raw):
    tmp = ''
    for idx, aa in enumerate(raw['peptide']):
        if aa == 'M' and idx+1 in raw['oxpos']:
            tmp += 'M(ox)'
        else:
            tmp += aa
    return tmp
    

class FileProcessor:

    @staticmethod
    def get_collision_energy(path):
        files = os.listdir(path)
        for file in files:
            if file.endswith('.mzML'):
                file_to_compute = file
                break
        m = mzml.read(path + '/' + file_to_compute)
        tmp_mzml = m.next()
        while not tmp_mzml['ms level'] == 2:
            tmp_mzml = m.next()
        collision_energy = float(tmp_mzml['precursorList']['precursor'][0]['activation']['collision energy'])
        return collision_energy

    @staticmethod
    def make_prosit_file(name, wild, collision_energy):
        if wild:
            wild_or_variant = "Wild"
        else:
            wild_or_variant = "Variant"
        pre_prosit_file = pd.read_csv(name)
        prosit_file = pre_prosit_file[['modified_sequence', 'assumed_charge']]
        prosit_file['collision_energy'] = [collision_energy for i in range(len(pre_prosit_file))]
        prosit_file.columns = ['modified_sequence', 'precursor_charge', 'collision_energy']
        return_name = "tmpPrositDir/" + wild_or_variant + "PrositFile.csv"
        prosit_file.to_csv(return_name, sep=',', index=False)
        return return_name

    @staticmethod
    def change_and_save_identipy_copy(name, wild):
        if wild:
            wild_or_variant = "Wild"
            identipy_file = pd.read_table(name)
            identipy_file = identipy_file[:1000]
        else:
            wild_or_variant = "Variant"
            identipy_file = pd.read_table(name + '_variants.tsv')
        # identipy_file = identipy_file.rename(columns={"peptide": "modified_sequence"})       
        identipy_file['length'] = identipy_file['peptide'].str.len()
        identipy_file = identipy_file.loc[identipy_file['length'] <= 30]
        identipy_file['oxpos'] = identipy_file['modifications'].apply(oxidation_positions)
        identipy_file['modified_sequence'] = identipy_file.apply(new_sequence, axis=1)
        
        return_name = 'tmpPrositDir/identipy' + wild_or_variant + 'TmpFile.csv'
        identipy_file.to_csv(return_name)
        return return_name

    @staticmethod
    def get_most_frequent_file_of_wild(file):
        files = os.listdir(file)
        most_len = 0
        best_wild_file = ''
        for i in files:
            if i.endswith("_wild_peptides.tsv"):
                tmp_df = pd.read_table(file + "/" + i)
                tmp_len = len(tmp_df)
                if tmp_len > most_len:
                    most_len = tmp_len
                    best_wild_file = i
        return best_wild_file
