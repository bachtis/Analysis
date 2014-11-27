from defs import *

print 'This step is needed only for data to apply the bfield map corrections'

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



print 'J /psi data'
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = recalibrateMap(builder.tree)
f2=ROOT.TFile('JData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'J /psi MC'

builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
f2=ROOT.TFile('JMC_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()



print 'Z data'

builder = DataSetBuilder(pmap,w,'../../data/ZDATA.root','data',10000000)
#builder.tree = recalibrateMap(builder.tree)
builder.tree = builder.tree.reduce('abs(eta)<5') 
f2=ROOT.TFile('ZData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'Z MC1'

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5') 
f2=ROOT.TFile('ZMC1_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'Z MC2'

builder = DataSetBuilder(pmap,w,'../../data/ZMC2.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5') 
f2=ROOT.TFile('ZMC2_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'Z GEN'

builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5') 
f2=ROOT.TFile('ZGEN_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()

print 'J GEN'

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5') 
f2=ROOT.TFile('JGEN_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()
