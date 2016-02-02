import ROOT
import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()
import math

def fetch(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    hh.SetName(name)
    f.cd()
    hh.Write()

def fetchEtaProj(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    h=hh.ProjectionY(name,1,1)
    f.cd()
    h.Write()



def fetchDiff(filename1,histo1,filename2,histo2,name,f):
    ff=ROOT.TFile(filename1)
    hh=ff.Get(histo1).Clone()
    fff=ROOT.TFile(filename2)
    hhh=fff.Get(histo2).Clone()
    new =ROOT.TH2D(name,name,hh.GetNbinsX(),hh.GetXaxis().GetXbins().GetArray(),hh.GetNbinsY(),hh.GetYaxis().GetXbins().GetArray())
    for i in range(1,hh.GetNbinsX()+1):
        for j in range(1,hh.GetNbinsY()+1):
                bin=hh.GetBin(i,j)
                print 'content',hh.GetBinContent(bin),hhh.GetBinContent(bin),hh.GetBinContent(bin)-hhh.GetBinContent(bin)
                new.SetBinContent(bin,hh.GetBinContent(bin)-hhh.GetBinContent(bin))
                new.SetBinError(bin,math.sqrt(hh.GetBinError(bin)*hh.GetBinError(bin)+hhh.GetBinError(bin)*hhh.GetBinError(bin)))
    new.SetName(name)
    f.cd()
    new.Write()


def fetchCovariance(filename,histo,f):
    ff=ROOT.TFile(filename)
    cov=ff.Get(histo)
    N=cov.GetNbinsX()
    matrix = ROOT.TMatrixDSym(N)

    for i in range(1,N+1):
        for j in range(1,N+1):
            bin=cov.GetBin(i,j)
            matrix[i-1,j-1] = cov.GetBinContent(bin)

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
        if obj[0] =='B0':
            histoMap.SetBinContent(i+1,4)
        if obj[0] =='B1':
            histoMap.SetBinContent(i+1,5)
        if obj[0] =='B2':
            histoMap.SetBinContent(i+1,6)
        if obj[0] =='C1':
            histoMap.SetBinContent(i+1,7)
        if obj[0] =='C2':
            histoMap.SetBinContent(i+1,8)

    f.cd()
    eigenvalues.Write("eigenvalues")
    eigenvectors.Write("eigenvectors")
    binMap.Write()
    histoMap.Write()


def fetchCovarianceBinned(filename,histo,f):
    ff=ROOT.TFile(filename)
    cov=ff.Get(histo)
    N=cov.GetNbinsX()
    matrix = ROOT.TMatrixDSym(N)

    for i in range(1,N+1):
        for j in range(1,N+1):
            bin=cov.GetBin(i,j)
            matrix[i-1,j-1] = cov.GetBinContent(bin)

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
        if obj[0] =='B0':
            histoMap.SetBinContent(i+1,4)
    f.cd()
    eigenvalues.Write("eigenvalues")
    eigenvectors.Write("eigenvectors")
    binMap.Write()
    histoMap.Write()



def fake(name,f):
    h=ROOT.TH3F(name,name,1,0,1,1,0,1,1,0,1)
    f.cd()
    h.Write()





f=ROOT.TFile("DATA_76X_13TeV.root","RECREATE")
fetch("MagneticMap/mapCorrections.root","magnetic","magnetic",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","A_0","A1",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","K_0","A2",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","M_0","e",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","B0_0","B0",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","B1_0","B1",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","B2_0","B2",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","C1_0","C1",f)
fetch("KalmanFilter/calibrationDATA13TeV.root","C2_0","C2",f)
fetchCovariance("KalmanFilter/calibrationDATA13TeV.root","COVARIANCE_0",f)

fetchEtaProj("Resolution/resultsDATA.root","a_0","aSRC",f)
fetchEtaProj("Resolution/resultsDATA.root","b_0","bSRC",f)
fetchEtaProj("Resolution/resultsDATA.root","c_0","cSRC",f)
fetchEtaProj("Resolution/resultsDATA.root","d_0","dSRC",f)

fetch("EbE/mcResolution.root","a_0","aTARGET",f)
fetch("EbE/mcResolution.root","b_0","bTARGET",f)
fetch("EbE/mcResolution.root","c_0","cTARGET",f)
fetch("EbE/mcResolution.root","d_0","dTARGET",f)

fetch("EbE/dataEbE.root","a_0","aEBE",f)
fetch("EbE/dataEbE.root","b_0","bEBE",f)
fetch("EbE/dataEbE.root","c_0","cEBE",f)
fetch("EbE/dataEbE.root","d_0","dEBE",f)



fetch("ClosureSyst/closure.root","closure","closure",f)



#fake("closure",f)

#fake("A1",f)
#fake("A2",f)
#fake("e",f)
#fake("B0",f)
#fake("B1",f)
#fake("B2",f)
#fake("C1",f)
#fake("C2",f)
#fake("b2",f)
f.Close()



f=ROOT.TFile("MC_76X_13TeV.root","RECREATE")
fetch("MagneticMap/mapCorrections.root","magnetic","magnetic",f)
fetch("KalmanFilter/calibrationMC13TeV.root","A_0","A1",f)
fetch("KalmanFilter/calibrationMC13TeV.root","K_0","A2",f)
fetch("KalmanFilter/calibrationMC13TeV.root","M_0","e",f)
fetch("KalmanFilter/calibrationMC13TeV.root","B0_0","B0",f)
fetch("KalmanFilter/calibrationMC13TeV.root","B1_0","B1",f)
fetch("KalmanFilter/calibrationMC13TeV.root","B2_0","B2",f)
fetch("KalmanFilter/calibrationMC13TeV.root","C1_0","C1",f)
fetch("KalmanFilter/calibrationMC13TeV.root","C2_0","C2",f)
fetchCovariance("KalmanFilter/calibrationMC13TeV.root","COVARIANCE_0",f)

fetchEtaProj("Resolution/resultsDATA.root","a_0","aTARGET",f)
fetchEtaProj("Resolution/resultsDATA.root","b_0","bTARGET",f)
fetchEtaProj("Resolution/resultsDATA.root","c_0","cTARGET",f)
fetchEtaProj("Resolution/resultsDATA.root","d_0","dTARGET",f)

fetch("EbE/mcResolution.root","a_0","aSRC",f)
fetch("EbE/mcResolution.root","b_0","bSRC",f)
fetch("EbE/mcResolution.root","c_0","cSRC",f)
fetch("EbE/mcResolution.root","d_0","dSRC",f)

fetch("EbE/mcEbE.root","a_0","aEBE",f)
fetch("EbE/mcEbE.root","b_0","bEBE",f)
fetch("EbE/mcEbE.root","c_0","cEBE",f)
fetch("EbE/mcEbE.root","d_0","dEBE",f)
fake("closure",f)

fetch("ClosureSyst/closure.root","closure","closure",f)



#fake("A1",f)
#fake("A2",f)
#fake("e",f)
#fake("B0",f)
#fake("B1",f)
#fake("B2",f)
#fake("C1",f)
#fake("C2",f)
#fake("b2",f)
f.Close()











