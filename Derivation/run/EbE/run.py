import ROOT
from array import array
from KaMuCa.Derivation.tools.Calibrator import Calibrator


ptArr=[5*5,7*7,9*9,10*10,12*12,15*15,20*20,30*30,40*40,50*50,60*60,80*80,100*100]
etaArr=[-2.5,-2.0,-1.7,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.7,2.0,2.5]
resArr=[]
for i in range(0,51):
    resArr.append(-0.2+i*0.4/50.0)

histo = ROOT.TH3D("histo","histo",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr),len(resArr)-1,array('d',resArr))
histoEbE = ROOT.TH3D("histoEbE","histo",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr),len(resArr)-1,array('d',resArr))



pullArr=[]
for i in range(0,51):
    pullArr.append(-4+i*8/50.0)
pull = ROOT.TH3D("pull","histo",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr),len(pullArr)-1,array('d',pullArr))


random=ROOT.TRandom(101082)

def fill(filename,histo,histoEbE,pull):
    calibrator = Calibrator("MC_76X_13TeV",True,False,False)
    f=ROOT.TFile(filename)
    data=f.Get("data")
    for i in range(0,data.numEntries()):
        line=data.get(i)
        calibrator.correct(line)
        m=line.find("mass").getVal()        

        c1=line.find("c1").getVal()        
        eta1=line.find("eta1").getVal()        
        gc1=line.find("gc1").getVal()
        err1=line.find("cErr1").getVal()
        cSmear=random.Gaus(c1,err1*c1)
        histoEbE.Fill(1.0/(c1*c1),eta1,(cSmear-c1)/c1)


        if gc1>0:
            histo.Fill(1.0/(c1*c1),eta1,(c1-gc1)/gc1)
            pull.Fill(1.0/(c1*c1),eta1,(c1-gc1)/(err1*c1))


        c1=line.find("c2").getVal()        
        eta1=line.find("eta2").getVal()        
        gc1=line.find("gc2").getVal()
        err1=line.find("cErr2").getVal()
        cSmear=random.Gaus(c1,err1*c1)
        histoEbE.Fill(1.0/(c1*c1),eta1,(cSmear-c1)/c1)

        if gc1>0:
            histo.Fill(1.0/(c1*c1),eta1,(c1-gc1)/gc1)
            pull.Fill(1.0/(c1*c1),eta1,(c1-gc1)/(err1*c1))

def fillData(filename,histo,histoEbE,pull):
    f=ROOT.TFile(filename)
    data=f.Get("data")
    for i in range(0,data.numEntries()):
        line=data.get(i)
        m=line.find("mass").getVal()        

        c1=line.find("c1").getVal()        
        eta1=line.find("eta1").getVal()        
        err1=line.find("cErr1").getVal()
        cSmear=random.Gaus(c1,err1*c1)
        histoEbE.Fill(1.0/(c1*c1),eta1,(cSmear-c1)/c1)

        c1=line.find("c2").getVal()        
        eta1=line.find("eta2").getVal()        
        err1=line.find("cErr2").getVal()
        cSmear=random.Gaus(c1,err1*c1)
        histoEbE.Fill(1.0/(c1*c1),eta1,(cSmear-c1)/c1)



#fill("/data/bachtis/CALIB/13TeV_76X/JMC.root",histo,histoEbE,pull)            
#fill("/data/bachtis/CALIB/13TeV_76X/ZMC.root",histo,histoEbE,pull)            
#fill("/data/bachtis/CALIB/13TeV_76X/YMC.root",histo,histoEbE,pull)            


fillData("/data/bachtis/CALIB/13TeV_76X/JDATA.root",histo,histoEbE,pull)            
fillData("/data/bachtis/CALIB/13TeV_76X/ZDATA.root",histo,histoEbE,pull)            
fillData("/data/bachtis/CALIB/13TeV_76X/YDATA.root",histo,histoEbE,pull)            

f=ROOT.TFile("outputDATA13TeV.root","RECREATE")
f.cd()
histo.Write()
histoEbE.Write()
pull.Write()
f.Close()
