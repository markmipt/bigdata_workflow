import pandas as pd


class FileProcessor:

    @staticmethod
    def make_prosit_file(file):
        name = file.split('/')[-1]
        pre_prosit_file = pd.read_table(file + '/' + name + '_variants.tsv')
        prosit_file = pd.DataFrame(pre_prosit_file['peptide'])
        prosit_file['precursor_charge'] = [1 for i in range(len(pre_prosit_file))]
        prosit_file['collision_energy'] = [30 for i in range(len(pre_prosit_file))]
        prosit_file.columns = ['modified_sequence', 'precursor_charge', 'collision_energy']
        prosit_file.to_csv('prositFile.csv', sep=',', index=False)

    def makeCorrelations(self, oldFile, newFile):
        test = pd.read_csv(oldFile)
        ready = pd.read_csv(newFile)

        full_list = []
        second_list = []
        new_list = []
        corr_list = []
        name_list = []
        i = 0
        test_len = len(test)
        for fn, spn, pep_name in test[['file_name', 'spectrum', 'modified_sequence']].values:
            print(str(i) + '/' + str(test_len))
            i += 1
            reader = mgf.read(fn.replace('_variant_final.tsv', '.mgf'))
            #     reader = mgf.read('test_data/e0590_lt_08_LyscTrypsin_identipy.mgf')
            mz_tmp = list(ready.loc[ready['LabeledPeptide'] == pep_name]['FragmentMz'])
            int_tmp = list(ready.loc[ready['LabeledPeptide'] == pep_name]['RelativeIntensity'])
            full_list.append([mz_tmp, int_tmp])
            mz_spec = list(reader[spn]['m/z array'])
            int_spec = list(reader[spn]['intensity array'])
            second_list.append([mz_spec, int_spec, pep_name])