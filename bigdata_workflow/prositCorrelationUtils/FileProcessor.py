import pandas as pd
import os


class FileProcessor:

    @staticmethod
    def make_prosit_file(file):
        name = file.split('/')[-1]
        # FIXME Change file to process
        pre_prosit_file = pd.read_csv('tmpPrositDir/identipyTmpFile.csv')
        print(pre_prosit_file.columns)
        prosit_file = pd.DataFrame(pre_prosit_file['modified_sequence'])
        prosit_file['precursor_charge'] = [1 for i in range(len(pre_prosit_file))]
        prosit_file['collision_energy'] = [30 for i in range(len(pre_prosit_file))]
        prosit_file.columns = ['modified_sequence', 'precursor_charge', 'collision_energy']
        prosit_file.to_csv('tmpPrositDir/prositFile.csv', sep=',', index=False)
        return 'tmpPrositDir/prositFile.csv'

    @staticmethod
    def change_and_save_identipy_copy(file):
        # FIXME Change file to process

        identipy_file = pd.read_table('/home/results/identipy/PXD005921/PXD005921_variants.tsv')

        identipy_file = identipy_file.rename(columns={"peptide": "modified_sequence"})
        identipy_file['length'] = identipy_file['modified_sequence'].str.len()
        identipy_file = identipy_file.loc[identipy_file['length'] <= 30]
        identipy_file.to_csv('tmpPrositDir/identipyTmpFile.csv')
        return 'tmpPrositDir/identipyTmpFile.csv'

    @staticmethod
    def get_most_frequent_file_of_wild(file):
        files = os.listdir(file)
        wilds = filter(lambda x: x.endswith('_wild_peptides.csv '), files)
        most_len = 0
        best_wild_file = ''
        for i in wilds:
            tmp_df = pd.read_csv(i)
            tmp_len = len(tmp_df)
            if tmp_len > most_len:
                most_len = tmp_len
                best_wild_file = i
        return best_wild_file
