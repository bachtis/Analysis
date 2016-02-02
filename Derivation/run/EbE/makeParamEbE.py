import ROOT
f=ROOT.TFile("outputDATA13TeV.root")
histo=f.Get("histoEbE")
histo.FitSlicesZ()
sigma = ROOT.gDirectory.Get("histoEbE_2")

func=ROOT.TF1("func","sqrt([0]+[1]/(1+[3]/x)+[2]*x)",0,100000)
func.SetParLimits(0,0,1000e-6)
func.SetParLimits(1,0,1000e-6)
func.SetParLimits(2,0,100e-7)
func.SetParLimits(3,0,100000)

histo_a = ROOT.TH1D("a_0","a",sigma.GetYaxis().GetNbins(),sigma.GetYaxis().GetXbins().GetArray())
histo_b = ROOT.TH1D("b_0","a",sigma.GetYaxis().GetNbins(),sigma.GetYaxis().GetXbins().GetArray())
histo_c = ROOT.TH1D("c_0","a",sigma.GetYaxis().GetNbins(),sigma.GetYaxis().GetXbins().GetArray())
histo_d = ROOT.TH1D("d_0","a",sigma.GetYaxis().GetNbins(),sigma.GetYaxis().GetXbins().GetArray())

for i in range(1,sigma.GetNbinsY()+1):
    func.SetParameters(1e-6,1e-6,1e-8,900)
    proje=sigma.ProjectionX("q",i,i)
    proje.Fit(func,"","",25,1000000)
    proje.Fit(func,"","",25,1000000)
    histo_a.SetBinContent(i,func.GetParameter(0))
    histo_a.SetBinError(i,func.GetParError(0))

    histo_b.SetBinContent(i,func.GetParameter(1))
    histo_b.SetBinError(i,func.GetParError(1))

    histo_c.SetBinContent(i,func.GetParameter(2))
    histo_c.SetBinError(i,func.GetParError(2))

    histo_d.SetBinContent(i,func.GetParameter(3))
    histo_d.SetBinError(i,func.GetParError(3))



f2=ROOT.TFile("dataEbE.root","RECREATE")
f2.cd()
histo_a.Write()
histo_b.Write()
histo_c.Write()
histo_d.Write()
f2.Close()
