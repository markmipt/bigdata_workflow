import pandas as pd
import os
from pyteomics import mzml
import ast
from ms2pip.ms2pipC import MS2PIP
import numpy as np
import math
from . import utils


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def cos_correlation(theoretical_list, experimental_list):
    suit_len = min(len(theoretical_list), len(experimental_list))

    theoretical_list = theoretical_list[:suit_len]
    experimental_list = experimental_list[:suit_len]
    top = 0
    bottom = math.sqrt(sum([numb * numb for numb in theoretical_list])) * \
                math.sqrt(sum([numb * numb for numb in experimental_list]))
    for i1, i2 in zip(theoretical_list, experimental_list):
        top += i1 * i2
    if top == 0:
        return 0

    return top / bottom

# def make_correlations(old_file, new_file, directory, folder_name, wild_or_variant, wild_name):



#     test = pd.read_csv(old_file)
#     ready = pd.read_table(new_file)
#     print(ready.columns)
#     full_list = []
#     second_list = []
#     new_list = []
#     corr_list = []
#     name_list = []
#     i = 0
#     test_len = len(test)
#     print(test.columns)
#     if wild_or_variant == "variant":
#         for fn, spn, pep_name in test[['filename', 'spectrum', 'modified_sequence']].values:
#             i += 1
#             # reader = mgf.read(directory + '/' + fn + '_identipy.mgf')
#             reader = mzml.read(directory + '/' + fn + '.mzML')
#             mz_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Masses'])
#             int_tmp = list(ready.loc[ready['Modified Sequence'] == pep_name]['Intensities'])
#             full_list.append([mz_tmp, int_tmp])
#             mz_spec = list(reader[spn]['m/z array'])
#             int_spec = list(reader[spn]['intensity array'])
#             second_list.append([mz_spec, int_spec, pep_name])


#     for idx, i in enumerate(full_list):
#         new_mz = []
#         new_int = []
#         for each in np.array(i[0][0].split(';'), dtype=float):
#             each = np.float(each)
#             new_id = self.find_nearest(second_list[idx][0], each)
#             if abs(second_list[idx][0][new_id] - each) <= 0.05:
#                 new_mz.append(second_list[idx][0][new_id])
#                 new_int.append(second_list[idx][1][new_id])
#             else:
#                 new_mz.append(0)
#                 new_int.append(0)
#         new_list.append([new_mz, new_int, second_list[idx][2]])

#     for idx, i in enumerate(new_list):
#         corr_list.append(self.cos_correlation(i[1], np.array(full_list[idx][1][0].split(';'), dtype=float)))
#         name_list.append(i[2])

#     test['correlation'] = corr_list

#     if wild_or_variant == 'variant':
#         test.to_csv(directory + '/' + folder_name + '_variants.tsv', sep='\t', index=False)
#     else:
#         test.to_csv(directory + '/' + folder_name + '_wilds.tsv', sep='\t', index=False)

def make_correlations(predictions, directory):

    # corr_list = []
    corr_dict = {}
    # name_list = []

    all_mzml_files = set(predictions['filename'].values)
    all_mzml_files = list(z + '.mzML' for z in all_mzml_files)
    for mzml_file in all_mzml_files:


        full_list = []
        second_list = []
        new_list = []

        print(mzml_file)

        mzml_file_short = mzml_file.replace('.mzML', '')

        idx_filename = predictions['filename'] == mzml_file_short
        prediction_tmp = predictions[idx_filename]
        all_sp_ids = set(prediction_tmp['spec_id'])
        

        # reader = mzml.read(os.path.join(directory, mzml_file))
        reader = mzml.PreIndexedMzML(os.path.join(directory, mzml_file))
        # for spectrum in reader:
        for spectrum_id in all_sp_ids:

            spectrum = reader.get_by_id(spectrum_id)

            # spectrum_id = spectrum['id']

            # if spectrum_id in all_sp_ids:

            idx_local = prediction_tmp['spec_id'] == spectrum_id

            mz_tmp = list(prediction_tmp.loc[idx_local]['mz'])
            int_tmp = list(prediction_tmp.loc[idx_local]['prediction'].apply(lambda x: ((2 ** x) - 0.001)))
            full_list.append([mz_tmp, int_tmp])


            mz_spec = list(spectrum['m/z array'])
            int_spec = list(spectrum['intensity array'])
            second_list.append([mz_spec, int_spec, spectrum_id])


        for idx, i in enumerate(full_list):
            new_mz = []
            new_int = []
            for mz_1, intensity_1 in zip(*i):
                new_id = find_nearest(second_list[idx][0], mz_1)
                if abs(second_list[idx][0][new_id] - mz_1) <= 0.05:
                    new_mz.append(second_list[idx][0][new_id])
                    new_int.append(second_list[idx][1][new_id])
                else:
                    new_mz.append(0)
                    new_int.append(0)
            new_list.append([new_mz, new_int, second_list[idx][2]])

        for idx, i in enumerate(new_list):

            corr_value = cos_correlation(i[1], full_list[idx][1])
            corr_dict[(mzml_file_short, i[2])] = corr_value

            # print(corr_value)

            # corr_list.append(cos_correlation(i[1], np.array(full_list[idx][1][0].split(';'), dtype=float)))
            # name_list.append(i[2])

    return corr_dict

    # test['correlation'] = corr_list

    # if wild_or_variant == 'variant':
    #     test.to_csv(directory + '/' + folder_name + '_variants.tsv', sep='\t', index=False)
    # else:
    #     test.to_csv(directory + '/' + folder_name + '_wilds.tsv', sep='\t', index=False)

def main_ms2pip(directory, folder_name):

    directory_tmp = os.path.join(directory, folder_name)

    df_variants = pd.read_table(os.path.join(directory, folder_name + '_variants.tsv'))

    df_wilds = pd.read_table(os.path.join(directory, folder_name + '_wilds.tsv'))

    # dfx = df_variants.drop_duplicates(subset=['peptide', 'assumed_charge']).copy()
    dfx = df_variants.copy()

    dfx = pd.concat([dfx, df_wilds])
    dfx.reset_index(inplace=True, drop=True)

    dfx['spec_id'] = dfx.apply(lambda x: x['filename'] + '_AND_' + x['spectrum'], axis=1)
    dfx['modifications'] = dfx.apply(utils.mods_for_deepLC, axis=1)
    dfx['charge'] = dfx['assumed_charge']

    params = {
        "ms2pip": {
            "ptm": [
                "Oxidation,15.994915,opt,M",
                "Carbamidomethyl,57.021464,opt,C",
                "Acetyl,42.010565,opt,N-term",
                "TMT6plex,229.162932,opt,N-term",
                "TMT6plex,229.162932,opt,K",
            ],
    #         "frag_method": "HCD",
            "frag_method": "TMT",
            "frag_error": 0.02,
            "out": "csv",
            "sptm": [], "gptm": [],
        }
    }
    ms2pip = MS2PIP(dfx, params=params, num_cpu=8, return_results=True)
    predictions = ms2pip.run()
    predictions[['filename', 'spec_id']] = predictions['spec_id'].str.split('_AND_', n=1, expand=True)

    corr_dict = make_correlations(predictions, directory)

    df_variants['correlation'] = df_variants.apply(lambda x: corr_dict[(x['filename'], x['spectrum'])], axis=1)
    df_wilds['correlation'] = df_wilds.apply(lambda x: corr_dict[(x['filename'], x['spectrum'])], axis=1)

    df_variants.to_csv(os.path.join(directory, folder_name + '_variants.tsv'), sep='\t', index=False)
    df_wilds.to_csv(os.path.join(directory, folder_name + '_wilds.tsv'), sep='\t', index=False)

    # print(folder_name)
    # # template_dir = TemplateDir()
    # # template_dir.create_prosit_template_dir()

    # print(directory_tmp)
    # print(directory)
    # # file_processor = FileProcessor()
    # collision_energy = get_collision_energy(directory)#file_processor.get_collision_energy('/'.join(directory_tmp.split('/')[:-1]))
    # identipy_file = file_processor.change_and_save_identipy_copy(directory_tmp, False)
    # prosit_file_variant = file_processor.make_prosit_file(identipy_file, False, collision_energy)


    # most_freq_wild_file = file_processor.get_most_frequent_file_of_wild(directory_tmp[:-len(folder_name)])
    # saved_wild_name = most_freq_wild_file
    # most_freq_wild_file = file_processor.change_and_save_identipy_copy(directory_tmp[:-len(folder_name)] + most_freq_wild_file, True)
    # prosit_file_wild = file_processor.make_prosit_file(most_freq_wild_file, True, collision_energy)

    # # # start Prosit server
    # # prosit_server = PrositServer()
    # # prosit_server.run_prosit_server(prosit_path, MODEL_SPECTRA, MODEL_IRT)
    # # after_prosit_variant = prosit_server.send_file(prosit_file_variant, False)
    # # after_prosit_wild = prosit_server.send_file(prosit_file_wild, True)
    # # prosit_server.stop_prosit_server()

    # # count correlations and changing file
    # correlation = Correlation()
    # correlation.make_correlations(identipy_file, after_prosit_variant, directory, folder_name, 'variant', saved_wild_name)
    # correlation.make_correlations(most_freq_wild_file, after_prosit_wild, directory, folder_name, 'wild', saved_wild_name)
    # template_dir.delete_prosit_template_dir()

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