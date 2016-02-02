import ROOT
def makePlot(filename,dataset,var,mini,maxi):
    f=ROOT.TFile(filename)
    d=f.Get(dataset)
    d.get().find(var).setMin(mini)
    d.get().find(var).setMax(maxi)
    frame=d.get().find(var).frame()
    d.plotOn(frame)
    return f,d,frame
