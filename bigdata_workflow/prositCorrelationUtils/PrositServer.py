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
        time.sleep(20)
        os.chdir(orig_wd)

    @staticmethod
    def stop_prosit_server():
        os.system('docker stop $(docker ps -q  --filter ancestor=prosit)')
        time.sleep(5)

    @staticmethod
    def send_file(file):
        os.system('curl -o tmp.csv -F "peptides=@mega_new.csv" http://127.0.0.1:3000/predict/msms')
