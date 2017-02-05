import ROOT
from math import sqrt
def quadrature(histo1,histo2):
    h=ROOT.TH2D(histo1)
    h.SetName("h")
    for i in range(1,histo1.GetNbinsX()+1):
        for j in range(1,histo1.GetNbinsY()+1):
            bin=histo1.GetBin(i,j)
            v1=histo1.GetBinContent(bin)
            e1=histo1.GetBinError(bin)
            v2=histo2.GetBinContent(bin)
            e2=histo2.GetBinError(bin)
            val = v1*v1-v2*v2
            error = 2*sqrt(v1*v1*e1*e1+v2*v2*e2*e2)
            h.SetBinContent(bin,val)
            h.SetBinError(bin,error)
    return h
 

def fit(histo,suffix=""):
    a=ROOT.TH1D(histo.ProjectionY())
    a.SetName("a"+suffix)
    b=ROOT.TH1D(histo.ProjectionY())
    b.SetName("b"+suffix)
    c=ROOT.TH1D(histo.ProjectionY())
    c.SetName("c"+suffix)
    d=ROOT.TH1D(histo.ProjectionY())
    d.SetName("d"+suffix)
    func=ROOT.TF1('func','sqrt([0]+[1]*x*x+[2]/(1+[3]/(x*x)))',0,200)

    func.SetParLimits(0,0,1)
    func.SetParLimits(1,0,1)
    func.SetParLimits(2,0,1)
    func.SetParLimits(3,0,10000)

    for i in range(1,histo.GetNbinsY()+1):
        proje=histo.ProjectionX("q",i,i)
        func.SetParameters(0,0,0,0)
        proje.Fit(func)
        proje.Fit(func)
        proje.Fit(func)
        a.SetBinContent(i,func.GetParameter(0))
        a.SetBinError(i,func.GetParError(0))
        b.SetBinContent(i,func.GetParameter(1))
        b.SetBinError(i,func.GetParError(1))
        c.SetBinContent(i,func.GetParameter(2))
        c.SetBinError(i,func.GetParError(2))
        d.SetBinContent(i,func.GetParameter(3))
        d.SetBinError(i,func.GetParError(3))
    return a,b,c,d 


def fitEBE(histo,suffix,dH):
    a=ROOT.TH1D(histo.ProjectionY())
    a.SetName("a"+suffix)
    b=ROOT.TH1D(histo.ProjectionY())
    b.SetName("b"+suffix)
    c=ROOT.TH1D(histo.ProjectionY())
    c.SetName("c"+suffix)

    for i in range(1,histo.GetNbinsY()+1):
        proje=histo.ProjectionX("q",i,i)
        func=ROOT.TF1('func','[0]+[1]*x*x+[2]/(1+{d}/(x*x))'.format(d=dH.GetBinContent(i)),0,200)
        func.SetParameters(0,0,0,0)
        func.SetParLimits(0,-1e-3,1e-3)
        func.SetParLimits(1,-1e-5,1e-5)
        func.SetParLimits(2,-1e-3,1e-3)
        proje.Fit(func)
        proje.Fit(func)
        proje.Fit(func)
        a.SetBinContent(i,func.GetParameter(0))
        a.SetBinError(i,func.GetParError(0))
        b.SetBinContent(i,func.GetParameter(1))
        b.SetBinError(i,func.GetParError(1))
        c.SetBinContent(i,func.GetParameter(2))
        c.SetBinError(i,func.GetParError(2))
    return a,b,c 

        

f=ROOT.TFile("results.root")
res=f.Get("resolution")
resRef=f.Get("resolutionRef")
res.FitSlicesZ()
resRef.FitSlicesZ()


sigma=ROOT.gDirectory.Get("resolution_2")
sigmaRef=ROOT.gDirectory.Get("resolutionRef_2")
h=quadrature(sigma,sigmaRef)

fOUT=ROOT.TFile("resolutionMC.root","RECREATE")
fOUT.cd()
a,b,c,d=fit(sigma)
a.Write("aRes")
b.Write("bRes")
c.Write("cRes")
d.Write("dRes")
d.Write("dEbE") # we use the same d and we change c only

a,b,c=fitEBE(h,"EBE",d)
a.Write("aEbE")
b.Write("bEbE")
c.Write("cEbE")


h.Write("diff2")

fOUT.Close()

