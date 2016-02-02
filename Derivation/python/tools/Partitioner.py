import ROOT

from array import array
import copy
import math

class Partitioner(object):
    def __init__(self,x,y,z,log = ""):

        self.dimensions={'x':x,'y':y,'z':z}
        self.map =ROOT.TH3F('map','map',len(x)-1,array('f',x),len(y)-1,array('f',y),len(z)-1,array('f',z))
        self.loader = None
        self.data={}
        if log != "":
            self.log = ROOT.TFile(log,"RECREATE")
        else:    
            self.log = None

        

    def addLog(self,log):
            self.log = ROOT.TFile(log,"RECREATE")
            self.log.cd()
        
    def limits(self):
        xmin = self.map.GetXaxis().GetXmin()
        xmax = self.map.GetXaxis().GetXmax()
        ymin = self.map.GetYaxis().GetXmin()
        ymax = self.map.GetYaxis().GetXmax()
        zmin = self.map.GetZaxis().GetXmin()
        zmax = self.map.GetZaxis().GetXmax()
        return {'x':[xmin,xmax],'y':[ymin,ymax],'z':[zmin,zmax]}


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
        
        for i in range(0,self.bins_x()+2):
            for j in range(0,self.bins_y()+2):
                for k in range(0,self.bins_z()+2):
                    bin = self.bin(i,j,k)
                    newMap.SetBinContent(bin,defaultVal)
        self.data[name] = newMap

    def setData(self,name,bin,data,error = 0.0):
        
        self.data[name].SetBinContent(bin,data)
        self.data[name].SetBinError(bin,error)

    def getData(self,name,bin):
        return self.data[name].GetBinContent(bin)

    def getDataError(self,name,bin):
        return self.data[name].GetBinError(bin)

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
        if binx==0 or biny==0 or binz==0:
            return -1
        if binx== self.map.GetNbinsX()+1 or biny==self.map.GetNbinsY()+1 or binz==self.map.GetNbinsZ()+1:
            return -1

        return self.bin(binx,biny,binz)


    def bins_x(self):
        return self.map.GetNbinsX()

    def bins_y(self):
        return self.map.GetNbinsY()

    def bins_z(self):
        return self.map.GetNbinsZ()

    def boundaries_x(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)
        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetXaxis().GetBinLowEdge(binC),self.map.GetXaxis().GetBinUpEdge(binC))

    def boundaries_y(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)

        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetYaxis().GetBinLowEdge(binEta),self.map.GetYaxis().GetBinUpEdge(binEta))

    def boundaries_z(self,bin):
        binC=ROOT.Long(0)
        binEta=ROOT.Long(0)
        binPhi=ROOT.Long(0)

        self.map.GetBinXYZ(bin,binC,binEta,binPhi)
        return (self.map.GetZaxis().GetBinLowEdge(binPhi),self.map.GetZaxis().GetBinUpEdge(binPhi))
    
