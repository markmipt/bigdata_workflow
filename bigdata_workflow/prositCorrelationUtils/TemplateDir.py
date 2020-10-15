import os


class TemplateDir:

    @staticmethod
    def create_prosit_template_dir():
        os.system('mkdir tmpPrositDir')

    @staticmethod
    def delete_prosit_template_dir():
        os.system('rm -r tmpPrositDir')