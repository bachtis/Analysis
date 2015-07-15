import ROOT
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()

def fetch(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    hh.SetName(name)
    f.cd()
    hh.Write()


def fetchCovariance(filename,histo,f):
    ff=ROOT.TFile(filename)
    cov=ff.Get(histo)
    N=cov.GetNbinsX()
    matrix = ROOT.TMatrixDSym(N)

    for i in range(1,N+1):
        for j in range(1,N+1):
            bin=cov.GetBin(i,j)
            matrix[i-1][j-1] = cov.GetBinContent(bin)

    choleskyD = ROOT.TDecompChol(matrix)
    choleskyD.Decompose()
    cholesky=choleskyD.GetU()

    #Also Diagonalize the matrix
    eigen = ROOT.TMatrixDSymEigen(matrix)
    eigenvalues = eigen.GetEigenValues()
    eigenvectors = eigen.GetEigenVectors()

    #create maps
    binMap = ROOT.TH1I("covBinMap","covBinMap",N,0,N)
    histoMap = ROOT.TH1I("covHistoMap","covBinMap",N,0,N)
    for i in range(0,N):
        label = cov.GetXaxis().GetBinLabel(i+1)
        obj = label.split('_')
        binMap.SetBinContent(i+1,int(obj[1]))

        if obj[0] =='A':
            histoMap.SetBinContent(i+1,1)
        if obj[0] =='K':
            histoMap.SetBinContent(i+1,2)
        if obj[0] =='M':
            histoMap.SetBinContent(i+1,3)
        if obj[0] =='B':
            histoMap.SetBinContent(i+1,4)
    


            


    f.cd()
    cholesky.Write("cholesky")
    eigenvalues.Write("eigenvalues")
    eigenvectors.Write("eigenvectors")

    binMap.Write()
    histoMap.Write()




f=ROOT.TFile("dataInputs.root","RECREATE")
fetch("step0_magneticMap/mapCalibration.root","mapCorrection","magnetic",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","A_1","A1",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","K_1","A2",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","B_1","B",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","M_1","e",f)
fetchCovariance("KalmanFilter_Ext/kalmanScale_data.root","COVARIANCE_0",f)
fetch("ClosureSyst/closure.root","closure","closure",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_A_target","sigma_A_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_B_target","sigma_B_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_C_target","sigma_C_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_A_target","sigma_A_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_B_target","sigma_B_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_C_target","sigma_C_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_A","ebe_A",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_B","ebe_B",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_C","ebe_C",f)


f.Close()

f=ROOT.TFile("mcInputs.root","RECREATE")
fetch("KalmanFilter_Ext/kalmanScale_mc.root","A_0","A1",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","K_0","A2",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","B_0","B",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","M_0","e",f)
fetchCovariance("KalmanFilter_Ext/kalmanScale_mc.root","COVARIANCE_0",f)
fetch("ClosureSyst/closure.root","closure","closure",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_A_target","sigma_A_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_B_target","sigma_B_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_C_target","sigma_C_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_A_target","sigma_A_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_B_target","sigma_B_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_C_target","sigma_C_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_A","ebe_A",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_B","ebe_B",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_C","ebe_C",f)


f.Close()










