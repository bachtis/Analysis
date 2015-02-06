import os
def getFullPath(path):
    return os.environ['CMSSW_BASE']+"/src/CMGTools/KalmanCalibration/data/"+path
