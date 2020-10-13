import subprocess
import os
import sys
import time
import pandas as pd
import numpy as np
from pyteomics import mgf
import math


class FileProcessor:

    @staticmethod
    def make_prosit_file(file):
        name = file.split('/')[-1]
        # FIXME Change file to process
        prePrositFile = pd.read_csv('tmpPrositDir/identipyTmpFile.csv')
        print(prePrositFile.columns)
        prositFile = pd.DataFrame(prePrositFile['modified_sequence'])
        prositFile['precursor_charge'] = [1 for i in range(len(prePrositFile))]
        prositFile['collision_energy'] = [30 for i in range(len(prePrositFile))]
        prositFile.columns = ['modified_sequence', 'precursor_charge', 'collision_energy']
        prositFile.to_csv('tmpPrositDir/prositFile.csv', sep=',', index=False)

    def change_and_save_identipy_copy(self, file):
        # FIXME Change file to process

        identipy_file = pd.read_table('/home/results/identipy/PXD005921/PXD005921_variants.tsv')

        identipy_file = identipy_file.rename(columns={"peptide": "modified_sequence"})
        identipy_file['length'] = identipy_file['modified_sequence'].str.len()
        identipy_file = identipy_file.loc[identipy_file['length'] <= 30]
        identipy_file.to_csv('tmpPrositDir/identipyTmpFile.csv')

    @staticmethod
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    @staticmethod
    def cos_correlation(theoretical_list, experimental_list):

        suit_len = min(len(theoretical_list), len(experimental_list))

        theoretical_list = theoretical_list[:suit_len]
        experimental_list = experimental_list[:suit_len]
        top = 0
        bottom = math.sqrt(sum([numb * numb for numb in theoretical_list])) * \
                 math.sqrt(sum([numb * numb for numb in experimental_list]))
        for i1, i2 in zip(theoretical_list, experimental_list):
            top += i1 * i2

        return top / bottom

    def make_correlations(self, old_file, new_file):

        test = pd.read_csv(old_file)
        ready = pd.read_table(new_file)

        full_list = []
        second_list = []
        new_list = []
        corr_list = []
        name_list = []
        i = 0
        test_len = len(test)
        for fn, spn, pep_name in test[['filename', 'spectrum', 'modified_sequence']].values:
            print(str(i) + '/' + str(test_len))
            i += 1
            reader = mgf.read('/home/results/identipy/PXD005921/' + fn + '_identipy.mgf')
            mz_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Masses'])
            int_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Intensities'])
            full_list.append([mz_tmp, int_tmp])
            mz_spec = list(reader[spn]['m/z array'])
            int_spec = list(reader[spn]['intensity array'])
            second_list.append([mz_spec, int_spec, pep_name])

        for idx, i in enumerate(full_list):
            new_mz = []
            new_int = []
            for each in np.array(i[0][0].split(';'), dtype=float):
                each = np.float(each)
                new_id = self.find_nearest(second_list[idx][0], each)
                if abs(second_list[idx][0][new_id] - each) <= 0.05:
                    new_mz.append(second_list[idx][0][new_id])
                    new_int.append(second_list[idx][1][new_id])
                else:
                    new_mz.append(0)
                    new_int.append(0)
            new_list.append([new_mz, new_int, second_list[idx][2]])

        for idx, i in enumerate(new_list):
            corr_list.append(self.cos_correlation(i[1], np.array(full_list[idx][1][0].split(';'), dtype=float)))
            name_list.append(i[2])

        test['correlation'] = corr_list
        test.to_csv(old_file+'_variants.tsv')
