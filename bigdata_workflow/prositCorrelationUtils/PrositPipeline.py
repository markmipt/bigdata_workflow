from .Correlation import Correlation
from .FileProcessor import FileProcessor
from .PrositServer import PrositServer
from .TemplateDir import TemplateDir


class PrositPipeline:
    @staticmethod
    def main_prosit(directory, prosit_path, MODEL_SPECTRA, MODEL_IRT):

        name = directory.split('/')[-2]

        print(name)
        template_dir = TemplateDir()
        template_dir.create_prosit_template_dir()

        print(directory)
        file_processor = FileProcessor()
        collision_energy = file_processor.get_collision_energy(directory)
        identipy_file = file_processor.change_and_save_identipy_copy(directory, False)
        prosit_file_variant = file_processor.make_prosit_file(identipy_file, False, collision_energy)
        most_freq_wild_file = file_processor.get_most_frequent_file_of_wild(directory[:-len(name)])
        saved_wild_name = most_freq_wild_file
        most_freq_wild_file = file_processor.change_and_save_identipy_copy(directory[:-len(name)] + most_freq_wild_file, True)
        prosit_file_wild = file_processor.make_prosit_file(most_freq_wild_file, True, collision_energy)

        # start Prosit server
        prosit_server = PrositServer()
        prosit_server.run_prosit_server(prosit_path, MODEL_SPECTRA, MODEL_IRT)
        after_prosit_variant = prosit_server.send_file(prosit_file_variant, False)
        after_prosit_wild = prosit_server.send_file(prosit_file_wild, True)
        prosit_server.stop_prosit_server()

        # count correlations and changing file
        correlation = Correlation()
        correlation.make_correlations(identipy_file, after_prosit_variant, name, 'variant', saved_wild_name)
        correlation.make_correlations(most_freq_wild_file, after_prosit_wild, name, 'wild', saved_wild_name)
        template_dir.delete_prosit_template_dir()
