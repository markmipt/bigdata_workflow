from prositCorrelationUtils.Correlation import Correlation
from prositCorrelationUtils.FileProcessor import FileProcessor
from prositCorrelationUtils.PrositServer import PrositServer
from prositCorrelationUtils.TemplateDir import TemplateDir


class PrositPipeline:
    @staticmethod
    def main_prosit(file):

        name = file.split('/')[-1]

        # make template dir
        template_dir = TemplateDir()
        template_dir.create_prosit_template_dir()

        # process prepare prosit file
        file_processor = FileProcessor()
        identipy_file = file_processor.change_and_save_identipy_copy(file)
        prosit_file_variant = file_processor.make_prosit_file(identipy_file)
        most_freq_wild_file = file_processor.get_most_frequent_file_of_wild(file)
        most_freq_wild_file = file_processor.change_and_save_identipy_copy(most_freq_wild_file)
        prosit_file_wild = file_processor.make_prosit_file(most_freq_wild_file)

        # start Prosit server
        prosit_server = PrositServer()
        prosit_server.run_prosit_server()
        after_prosit_variant = prosit_server.send_file(prosit_file_variant)
        after_prosit_wild = prosit_server.send_file(prosit_file_wild)
        prosit_server.stop_prosit_server()

        # count correlations and changing file
        correlation = Correlation()
        correlation.make_correlations(identipy_file, after_prosit_variant, name, 'variant')
        correlation.make_correlations(most_freq_wild_file, after_prosit_wild, name, 'wild')
