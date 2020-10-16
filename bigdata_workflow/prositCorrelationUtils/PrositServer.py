import os
import time


class PrositServer:

    @staticmethod
    def run_prosit_server(prosit_path, MODEL_SPECTRA, MODEL_IRT):
        orig_wd = os.getcwd()
        os.chdir(prosit_path)
        stream = os.popen(
            "make server MODEL_SPECTRA=" + MODEL_SPECTRA + " "
            "MODEL_IRT=" + MODEL_IRT + " HOSTPORT=3000")
        time.sleep(40)
        os.chdir(orig_wd)

    @staticmethod
    def stop_prosit_server():
        os.system('docker stop $(docker ps -q  --filter ancestor=prosit)')
        time.sleep(5)

    @staticmethod
    def send_file(file, wild):
        if wild:
            print(file)
            os.system('curl -o tmpPrositDir/WildAfterProsit.csv -F "peptides=@' + file + '" http://127.0.0.1:3000/predict/msms')
            time.sleep(15)
            return 'tmpPrositDir/WildAfterProsit.csv'
        else:
            print(file)
            os.system('curl -o tmpPrositDir/VariantAfterProsit.csv -F "peptides=@' + file + '" http://127.0.0.1:3000/predict/msms')
            time.sleep(15)
            return 'tmpPrositDir/VariantAfterProsit.csv'
