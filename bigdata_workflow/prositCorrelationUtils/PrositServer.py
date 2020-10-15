import os
import time


class PrositServer:

    @staticmethod
    def run_prosit_server():
        orig_wd = os.getcwd()
        os.chdir('/home/zloydanny/test_data/prosit/')
        stream = os.popen(
            "make server MODEL_SPECTRA=/home/zloydanny/test_data/prosit/model_spectra/ "
            "MODEL_IRT=/home/zloydanny/test_data/prosit/model_irt_prediction HOSTPORT=3000")
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
            return 'tmpPrositDir/WildAfterProsit.csv'
        else:
            print(file)
            os.system('curl -o tmpPrositDir/VariantAfterProsit.csv -F "peptides=@' + file + '" http://127.0.0.1:3000/predict/msms')
            return 'tmpPrositDir/VariantAfterProsit.csv'
