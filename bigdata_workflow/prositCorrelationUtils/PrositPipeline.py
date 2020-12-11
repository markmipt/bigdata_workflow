from .Correlation import Correlation
from .FileProcessor import FileProcessor
from .PrositServer import PrositServer
from .TemplateDir import TemplateDir
from os import path


class PrositPipeline:
    @staticmethod
    def main_prosit(directory, folder_name, prosit_path, MODEL_SPECTRA, MODEL_IRT):

        directory_tmp = path.join(directory, folder_name)


        print(folder_name)
        template_dir = TemplateDir()
        template_dir.create_prosit_template_dir()

        print(directory_tmp)
        print(directory)
        file_processor = FileProcessor()
        collision_energy = file_processor.get_collision_energy('/'.join(directory_tmp.split('/')[:-1]))
        identipy_file = file_processor.change_and_save_identipy_copy(directory_tmp, False)
        prosit_file_variant = file_processor.make_prosit_file(identipy_file, False, collision_energy)
        most_freq_wild_file = file_processor.get_most_frequent_file_of_wild(directory_tmp[:-len(folder_name)])
        saved_wild_name = most_freq_wild_file
        most_freq_wild_file = file_processor.change_and_save_identipy_copy(directory_tmp[:-len(folder_name)] + most_freq_wild_file, True)
        prosit_file_wild = file_processor.make_prosit_file(most_freq_wild_file, True, collision_energy)

        # start Prosit server
        prosit_server = PrositServer()
        prosit_server.run_prosit_server(prosit_path, MODEL_SPECTRA, MODEL_IRT)
        after_prosit_variant = prosit_server.send_file(prosit_file_variant, False)
        after_prosit_wild = prosit_server.send_file(prosit_file_wild, True)
        prosit_server.stop_prosit_server()

        # count correlations and changing file
        correlation = Correlation()
        correlation.make_correlations(identipy_file, after_prosit_variant, directory, folder_name, 'variant', saved_wild_name)
        correlation.make_correlations(most_freq_wild_file, after_prosit_wild, directory, folder_name, 'wild', saved_wild_name)
        template_dir.delete_prosit_template_dir()
