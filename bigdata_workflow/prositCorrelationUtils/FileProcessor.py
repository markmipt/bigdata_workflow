import pandas as pd
import os


class FileProcessor:

    @staticmethod
    def make_prosit_file(name, wild):
        if wild:
            wild_or_variant = "Wild"
        else:
            wild_or_variant = "Variant"
        pre_prosit_file = pd.read_csv(name)
        prosit_file = pd.DataFrame(pre_prosit_file['modified_sequence'])
        prosit_file['precursor_charge'] = [1 for i in range(len(pre_prosit_file))]
        prosit_file['collision_energy'] = [30 for i in range(len(pre_prosit_file))]
        prosit_file.columns = ['modified_sequence', 'precursor_charge', 'collision_energy']
        return_name = "tmpPrositDir/" + wild_or_variant + "PrositFile.csv"
        prosit_file.to_csv(return_name, sep=',', index=False)
        return return_name

    @staticmethod
    def change_and_save_identipy_copy(name, wild):
        if wild:
            wild_or_variant = "Wild"
            identipy_file = pd.read_table(name)
        else:
            wild_or_variant = "Variant"
            identipy_file = pd.read_table(name + '_variants.tsv')
        identipy_file = identipy_file.rename(columns={"peptide": "modified_sequence"})
        identipy_file['length'] = identipy_file['modified_sequence'].str.len()
        identipy_file = identipy_file.loc[identipy_file['length'] <= 30]
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
