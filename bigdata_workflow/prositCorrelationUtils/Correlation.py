import math
import numpy as np
import pandas as pd
from pyteomics import mgf
from progress.bar import IncrementalBar


class Correlation:

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

    def make_correlations(self, old_file, new_file, directory, folder_name, wild_or_variant, wild_name):

        test = pd.read_csv(old_file)
        ready = pd.read_table(new_file)
        print(ready.columns)
        full_list = []
        second_list = []
        new_list = []
        corr_list = []
        name_list = []
        i = 0
        test_len = len(test)
        print(test.columns)
        if wild_or_variant == "variant":
            bar = IncrementalBar('Countdown', max = test_len)
            for fn, spn, pep_name in test[['filename', 'spectrum', 'modified_sequence']].values:
                bar.next()
                i += 1
                reader = mgf.read(directory + '/' + fn + '_identipy.mgf')
                mz_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Masses'])
                int_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Intensities'])
                full_list.append([mz_tmp, int_tmp])
                mz_spec = list(reader[spn]['m/z array'])
                int_spec = list(reader[spn]['intensity array'])
                second_list.append([mz_spec, int_spec, pep_name])
            bar.finish()
        else:
            bar = IncrementalBar('Countdown', max = test_len)
            for spn, pep_name in test[[ 'spectrum', 'modified_sequence']].values:
                bar.next()
                fn = wild_name[:-len('_identipy_wild_peptides.tsv')]
                i += 1
                reader = mgf.read(directory + '/' + fn + '_identipy.mgf')
                mz_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Masses'])
                int_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Intensities'])
                full_list.append([mz_tmp, int_tmp])
                mz_spec = list(reader[spn]['m/z array'])
                int_spec = list(reader[spn]['intensity array'])
                second_list.append([mz_spec, int_spec, pep_name])
            bar.finish()
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

        if wild_or_variant == 'variant':
            test.to_csv(directory + '/' + folder_name + '_variants.tsv', sep='\t', index=False)
        else:
            test.to_csv(directory + '/' + folder_name + '_wilds.tsv', sep='\t', index=False)
