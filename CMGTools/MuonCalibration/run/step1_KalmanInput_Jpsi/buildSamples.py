from defs import *

print 'This step is needed only for data to apply the bfield map corrections'

def recalibrateMap(data):
    f=ROOT.TFile('../0_magneticfieldmap/mapCalibration.root')
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


builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = recalibrateMap(builder.tree)
builder.tree = builder.tree.reduce('1.0/curvRaw1>4.&& 1.0/curvRaw2>4.') 

builder.buildVsDiLepton()
builder.save("JDATA_Input.root")

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000,False,True)
builder.tree = builder.tree.reduce('1.0/curvRaw1>4.&& 1.0/curvRaw2>4.') 
builder.buildVsDiLepton()
builder.save("JGEN_Input.root")

