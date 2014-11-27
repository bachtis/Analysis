from defs import *

def recalibrateMap(data):
    f=ROOT.TFile('../step0_magneticMap/mapCalibration.root')
    h = f.Get('mapCorrection')
    newData = ROOT.RooDataSet(data.GetName()+'cal','',data.get())

    for i in range(0,data.numEntries()):
        line =data.get(i)
        bin1 = h.GetBin(h.GetXaxis().FindBin(line.find('phiRaw1').getVal()),h.GetYaxis().FindBin(line.find('etaRaw1').getVal()))
        bin2 = h.GetBin(h.GetXaxis().FindBin(line.find('phiRaw2').getVal()),h.GetYaxis().FindBin(line.find('etaRaw2').getVal()))
        
            
        factor1 = h.GetBinContent(bin1)
        factor2 = h.GetBinContent(bin2)
        valNew1 = line.find('curvRaw1').getVal()*factor1
        valNew2 = line.find('curvRaw2').getVal()*factor2
        line.find('curvRaw1').setVal(valNew1)
        line.find('curvRaw2').setVal(valNew2)

        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                        line.find('etaRaw1').getVal(),
                        line.find('phiRaw1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                        line.find('etaRaw2').getVal(),
                        line.find('phiRaw2').getVal(),
                        0.1056583715)
        line.find('massRaw').setVal((v1+v2).M())
        
        newData.add(line)
        
    return newData



def loop(data,h1,h2,H):
    h1.Sumw2()
    h2.Sumw2()
    H.Sumw2()
    for i in range(0,data.numEntries()):
        line =data.get(i)

        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                        line.find('etaRaw1').getVal(),
                        line.find('phiRaw1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                        line.find('etaRaw2').getVal(),
                        line.find('phiRaw2').getVal(),
                        0.1056583715)
        
        h1.Fill(v1.Pt())
        h2.Fill(v2.Pt())
        H.Fill((v1+v2).Pt())


builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = recalibrateMap(builder.tree)
dataPt1 = ROOT.TH1F("dataPt1","mcPt1",40,0,40)
dataPt2 = ROOT.TH1F("dataPt2","mcPt",40,0,40)
dataPt  = ROOT.TH1F("dataPt","mcPt",40,0,40)
loop(builder.tree,dataPt1,dataPt2,dataPt)


builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)

mcPt1 = ROOT.TH1F("mcPt1","mcPt1",40,0,40)
mcPt2 = ROOT.TH1F("mcPt2","mcPt",40,0,40)
mcPt  = ROOT.TH1F("mcPt","mcPt",40,0,40)
loop(builder.tree,mcPt1,mcPt2,mcPt)


mcPt1.SetLineWidth(2)
mcPt2.SetLineWidth(2)
mcPt.SetLineWidth(2)

dataPt1.SetMarkerStyle(20)
dataPt2.SetMarkerStyle(20)
dataPt.SetMarkerStyle(20)


c1=ROOT.TCanvas("c1","c1")
c1.cd()
mcPt1.DrawNormalized("HIST")
dataPt1.DrawNormalized("Psame")

c2=ROOT.TCanvas("c2","c2")
c2.cd()
mcPt2.DrawNormalized("HIST")
dataPt2.DrawNormalized("Psame")

c3=ROOT.TCanvas("c3","c2")
c3.cd()
mcPt.DrawNormalized("HIST")
dataPt.DrawNormalized("Psame")



