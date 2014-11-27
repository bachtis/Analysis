import ROOT
from array import array
import copy
import math
class PartitionMap(object):
    def __init__(self,ptInv,eta,phi,log = ""):
        self.map =ROOT.TH3F('map','map',len(ptInv)-1,array('f',ptInv),len(eta)-1,array('f',eta),len(phi)-1,array('f',phi))

        self.loader = None
        self.data={}
        if log != "":
            self.log = ROOT.TFile(log,"RECREATE")
        else:    
            self.log = None

        

    def limits(self):
        xmin = self.map.GetXaxis().GetXmin()
        xmax = self.map.GetXaxis().GetXmax()
        ymin = self.map.GetYaxis().GetXmin()
        ymax = self.map.GetYaxis().GetXmax()
        zmin = self.map.GetZaxis().GetXmin()
        zmax = self.map.GetZaxis().GetXmax()
        return {'curv':[xmin,xmax],'eta':[ymin,ymax],'phi':[zmin,zmax]}

    def load(self,filename,histo,data):
        if self.loader is None:
            self.loader = ROOT.TFile(filename)
        h=self.loader.Get(histo)
        self.data[data] = h.Clone()
        

    def save(self,prefix = ""):
        if self.log is not None:
            for name,histo in self.data.iteritems():
                self.log.cd()
                histo.Write(name+"_"+prefix)
            self.log.Flush()    


            

    def average(self,name,datasets):
        sumw=0.0
        sum=0.0
        for binx in range(1,self.map.GetXaxis().GetNbins()+1):
            for biny in range(1,self.map.GetYaxis().GetNbins()+1):
                for binz in range(1,self.map.GetYaxis().GetNbins()+1):
                    bin = self.bin(binx,biny,binz)
                    sum=sum+datasets[bin].sumEntries()*self.data[name].GetBinContent(bin)
                    sumw=sumw+datasets[bin].sumEntries()
        return sum/sumw             
        
            
    def declareData(self,name,defaultVal = 0.0):
        newMap = copy.deepcopy(self.map)
        newMap.SetName(name)
        
        for i in range(0,self.bins_curv()+2):
            for j in range(0,self.bins_eta()+2):
                for k in range(0,self.bins_phi()+2):
                    bin = self.bin(i,j,k)
                    newMap.SetBinContent(bin,defaultVal)
        self.data[name] = newMap
                    


    def setData(self,name,bin,data,error = 0.0):
        self.data[name].SetBinContent(bin,data)
        self.data[name].SetBinError(bin,error)

    def getData(self,name,bin):
        return self.data[name].GetBinContent(bin)

    def addData(self,name,bin,data,error = 0.0):
        self.data[name].SetBinContent(bin,self.data[name].GetBinContent(bin)+data)
        self.data[name].SetBinError(bin,error)

    def multiplyData(self,name,bin,data,error=0.0):
        self.data[name].SetBinContent(bin,self.data[name].GetBinContent(bin)*data)
        self.data[name].SetBinError(bin,error)

    def multiplyDataProj(self,name,binx,biny,binz,data):
        if binx == 'all':
            binsx = range(1,self.map.GetXaxis().GetNbins()+1)
        else:
            binsx = [binx]
        if biny == 'all':
            binsy = range(1,self.map.GetYaxis().GetNbins()+1)
        else:
            binsy = [biny]
        if binz == 'all':
            binsz = range(1,self.map.GetZaxis().GetNbins()+1)
        else:
            binsz = [binz]

        for bx  in binsx:
            for by  in binsy:
                for bz  in binsz:
                    self.data[name].SetBinContent(self.bin(bx,by,bz),self.data[name].GetBinContent(self.bin(bx,by,bz))*data)

    def addDataProj(self,name,binx,biny,binz,data):
        if binx == 'all':
            binsx = range(1,self.map.GetXaxis().GetNbins()+1)
        else:
            binsx = [binx]
        if biny == 'all':
            binsy = range(1,self.map.GetYaxis().GetNbins()+1)
        else:
            binsy = [biny]
        if binz == 'all':
            binsz = range(1,self.map.GetZaxis().GetNbins()+1)
        else:
            binsz = [binz]

        for bx  in binsx:
            for by  in binsy:
                for bz  in binsz:
                    self.data[name].SetBinContent(self.bin(bx,by,bz),self.data[name].GetBinContent(self.bin(bx,by,bz))+data)
            

        


    def getMap(self):
        return self.map

    def bin(self,x,y,z):
        return self.map.GetBin(x,y,z)

    def binXYZ(self,bin):
        binx = ROOT.Long(0)
        biny = ROOT.Long(0)
        binz = ROOT.Long(0)
        self.map.GetBinXYZ(bin,binx,biny,binz)

        
        return binx,biny,binz 

    def binVals(self,bin):
        binx,biny,binz = self.binXYZ(bin)
        return self.map.GetXaxis().GetBinCenter(binx), \
               self.map.GetYaxis().GetBinCenter(biny), \
               self.map.GetZaxis().GetBinCenter(binz)
    



    def binFromVals (self,x,y,z):
        binx = self.map.GetXaxis().FindBin(x)
        biny = self.map.GetYaxis().FindBin(y)
        binz = self.map.GetZaxis().FindBin(z)
        return self.bin(binx,biny,binz)


    def bins_curv(self):
        return self.map.GetNbinsX()

    def bins_eta(self):
        return self.map.GetNbinsY()

    def bins_phi(self):
        return self.map.GetNbinsZ()

    def boundaries_curv(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)

        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetXaxis().GetBinLowEdge(binC),self.map.GetXaxis().GetBinUpEdge(binC))

    def boundaries_eta(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)

        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetYaxis().GetBinLowEdge(binEta),self.map.GetYaxis().GetBinUpEdge(binEta))

    def boundaries_phi(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)

        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetZaxis().GetBinLowEdge(binPhi),self.map.GetZaxis().GetBinUpEdge(binPhi))
    
    def recalibrate(self,data,what = 'C',calib = 'sum',leg=2):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            bin1 = self.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
            
            factor1 = self.getData(what+'pos',bin1)
            bin2 = self.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())
            
            factor2 = self.getData(what+'neg',bin2)


            if calib =="sum":
                valNew1 = line.find('curvRaw1').getVal()+factor1
                valNew2 = line.find('curvRaw2').getVal()+factor2
            elif calib =="prod":
                valNew1 = line.find('curvRaw1').getVal()*factor1
                valNew2 = line.find('curvRaw2').getVal()*factor2
            elif calib =="sum2":
                valNew1 = line.find('curvRaw1').getVal()*(1+factor1*line.find('curvRaw1').getVal())
                valNew2 = line.find('curvRaw2').getVal()*(1+factor2*line.find('curvRaw2').getVal())
            if leg==0 or leg ==2:
                line.find('curvRaw1').setVal(valNew1)
            if leg==1 or leg ==2:
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


    def recalibrateGlobal(self,factor):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            valNew1 = line.find('curvRaw1').getVal()*factor
            valNew2 = line.find('curvRaw2').getVal()*factor
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


    def recalibrateBoth(self,data,leg=2):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            bin1 = self.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
            
            Cpos = self.getData('Cpos',bin1)
            Upos = self.getData('Upos',bin1)
            
            bin2 = self.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())
            
            Cneg = self.getData('Cneg',bin2)
            Uneg = self.getData('Uneg',bin2)
            if leg ==0 or leg==2:
                valNew1 = line.find('curvRaw1').getVal()*Upos+Cpos
                line.find('curvRaw1').setVal(valNew1)
            if leg ==1 or leg==2:
                valNew2 = line.find('curvRaw2').getVal()*Uneg+Cneg
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


    def recalibrateTriple(self,data,leg=2):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            bin1 = self.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
            
            Apos = self.getData('A',bin1)
            Bpos = self.getData('B',bin1)
            Mpos = self.getData('M',bin1)
            
            bin2 = self.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())

            Aneg = self.getData('A',bin2)
            Bneg = self.getData('B',bin2)
            Mneg = self.getData('M',bin2)
            
            c1 = line.find('curvRaw1').getVal()
            c2 = line.find('curvRaw2').getVal()

            if leg ==0 or leg==2:
                valNew1 = (Apos-1.0)*c1 + c1/(1+Mpos*c1)+Bpos
                line.find('curvRaw1').setVal(valNew1)
            if leg ==1 or leg==2:
                valNew2 = (Aneg-1.0)*c2 + c2/(1+Mneg*c2)-Bneg
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


    def recalibrateBothMEV(self,data,leg=2):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            bin1 = self.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
            
            Cpos = self.getData('Cpos',bin1)
            Upos = self.getData('Upos',bin1)
            
            bin2 = self.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())
            
            Cneg = self.getData('Cneg',bin2)
            Uneg = self.getData('Uneg',bin2)
            if leg ==0 or leg==2:
                p = 1./line.find('curvRaw1').getVal()+Upos
                valNew1 = 1./p+Cpos
                line.find('curvRaw1').setVal(valNew1)
            if leg ==1 or leg==2:
                p = 1./line.find('curvRaw2').getVal()+Uneg
                valNew2 = 1/p+Cneg
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


    def recalibrateEbE(self,data):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            bin1 = self.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
            
            bin2 = self.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())

            factor1 = self.getData('Rpos',bin1)
            factor2 = self.getData('Rneg',bin2)
            line.find('massErrRaw1').setVal(line.find('massErrRaw1').getVal()*factor1)
            line.find('massErrRaw2').setVal(line.find('massErrRaw2').getVal()*factor2)
            line.find('massErrRaw').setVal(math.sqrt(line.find('massErrRaw1').getVal()*line.find('massErrRaw1').getVal()+\
                                                  line.find('massErrRaw2').getVal()*line.find('massErrRaw2').getVal()))
            newData.add(line)

        return newData




    def smearEbE2D(self,data,w,shift = 0.0,errorScale=1.0):
        random = ROOT.TRandom3(101082)
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            errx = errorScale*2*line.find('massErrRaw1').getVal()*line.find('curvRaw1').getVal()/line.find('massRaw').getVal()
            erry = errorScale*2*line.find('massErrRaw2').getVal()*line.find('curvRaw2').getVal()/line.find('massRaw').getVal()
            
#            print 'before',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal()
            c1 = shift*line.find('curvRaw1').getVal()+random.Gaus(0.0,errx)
            c2 = shift*line.find('curvRaw2').getVal()+random.Gaus(0.0,erry)

            if c1<=0:
                c1=1e-19
            if c2<=0:
                c2=1e-19
            
            line.find('curvRaw1').setVal(c1)
            line.find('curvRaw2').setVal(c2)
#            print 'after',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal(),errx,erry

            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                                   line.find('etaRaw1').getVal(),
                                   line.find('phiRaw1').getVal(),
                                   w.var('muMass').getVal())
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                                   line.find('etaRaw2').getVal(),
                                   line.find('phiRaw2').getVal(),
                                   w.var('muMass').getVal())

            line.find('massRaw').setVal((v1+v2).M())

            
            newData.add(line)

        return newData


    def smear2D(self,data,w,resolution=0.02):
        random = ROOT.TRandom3(101082)
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            errx = resolution*line.find('curvRaw1').getVal()
            erry = resolution*line.find('curvRaw2').getVal()


#            print 'before',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal(),line.find('massRaw').getVal()

            line.find('curvRaw1').setVal(line.find('curvRaw1').getVal()+random.Gaus(0.0,errx))
            line.find('curvRaw2').setVal(line.find('curvRaw2').getVal()+random.Gaus(0.0,erry))

            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                                   line.find('etaRaw1').getVal(),
                                   line.find('phiRaw1').getVal(),
                                   w.var('muMass').getVal())
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                                   line.find('etaRaw2').getVal(),
                                   line.find('phiRaw2').getVal(),
                                   w.var('muMass').getVal())

            line.find('massRaw').setVal((v1+v2).M())

#            print 'after',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal(),line.find('massRaw').getVal()
            
            newData.add(line)

        return newData





