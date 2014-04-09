import ROOT
ROOT.gSystem.Load("libCMGToolsMuonCalibration")


def findOutliers(filename,hist='Cpos_0',threshold=0.01):
    f=ROOT.TFile(filename)
    h=f.Get(hist)
    max=0
    maxBin=-1
    for i in range(1,h.GetNbinsX()+1):
        for j in range(1,h.GetNbinsY()+1):
            for k in range(1,h.GetNbinsZ()+1):
                bin=h.GetBin(i,j,k)
                content=abs(h.GetBinContent(bin))
                if content>threshold:
                    print 'bin above threshold',bin
                if content>max:
                    max=content
                    maxBin=bin
    print 'Maximum at',maxBin,'with content =',max



def checkFile(filename):
    f=ROOT.TFile(filename)
    w=f.Get("w")
    w.var('scale').Print()
    frame=w.var("massRaw").frame()
    w.data("data").plotOn(frame)
    w.pdf("model").plotOn(frame)
    c=ROOT.TCanvas("c")
    c.cd()
    frame.Draw()
    c.Draw()
    c.Flush()
    return c,frame




def refit(filename):
    f=ROOT.TFile(filename,"UPDATE")
    f.cd()
    w=f.Get("w")
    w.pdf("model").fitTo(w.data("data"))
    w.saveSnapshot('result',w.allVars())
    w.Write("",ROOT.TObject.kOverwrite)
    f.Close()






