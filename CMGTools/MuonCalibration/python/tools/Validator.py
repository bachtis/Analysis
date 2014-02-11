import ROOT
import os


class Validator:
    def __init__(self,fname):
        self.file = ROOT.TFile(fname)
#        ROOT.gROOT.ProcessLine(".x "+os.environ['CMSSW_BASE']+"/src/CMGTools/MuonCalibration/python/tools/tdrstyle.C")

        
    def getTrend(self,var,bin,mini = 0 , maxi=5):
        g = ROOT.TGraph(maxi-mini)
        for i in range(mini,maxi):
            h = self.file.Get(var+"_"+str(i))
            g.SetPoint(i,i,h.GetBinContent(bin))
        return g    

    def getTrendXYZ(self,var,binx,biny,binz,mini = 0 , maxi=5):
            h = self.file.Get(var+"_"+str(mini))
            bin = h.GetBin(binx,biny,binz)
            return self.getTrend(var,bin,mini,maxi)

    def getTrendPtEtaPhi(self,var,pt,eta,phi,mini = 0 , maxi=5):
            h = self.file.Get(var+"_"+str(mini))
            binx = h.GetXaxis().FindBin(1./pt)
            biny = h.GetYaxis().FindBin(eta)
            binz = h.GetZaxis().FindBin(phi)
            return self.getTrendXYZ(var,binx,biny,binz,mini,maxi)
        

    def getMultiTrend(self,var,mini,maxi):
            h = self.file.Get(var+"_"+str(mini))
            c = ROOT.TCanvas('c','c')
            obj=[]
            for binx in range(1,h.GetNbinsX()+1):
                for biny in range(1,h.GetNbinsY()+1):
                    for binz in range(1,h.GetNbinsZ()+1):
                        bin = h.GetBin(binx,biny,binz)
                        obj.append(self.getTrend(var,bin,mini,maxi))
                        if len(obj)==1:
                            obj[-1].Draw("AC")
                            obj[-1].GetYaxis().SetRangeUser(0.97,1.03)
                        else:    
                            obj[-1].Draw("Csame")
            return c,obj                

    def getAverageProjection(self,histo,axis = 'X'):
        h = self.file.Get(histo)
        axes = ['X','Y','Z']
        axes.remove(axis)
        bins=1
        for ax in axes:
            bins=bins*getattr(h,'Get'+ax+'axis')().GetNbins()
        hp = getattr(h,'Projection'+axis)()
        hp.SetName(histo+'proj'+axis)
        hp.Scale(1./bins)
        return hp

    def getAverageProjection3D(self,histo,axisp = 'xy',axisbins = 'Z'):
        h = self.file.Get(histo)
        hp = h.Project3D(axisp)
        bins=getattr(h,'Get'+axisbins+'axis')().GetNbins()
        hp.Scale(1./bins)
        return hp



    def findBinsAboveLimit(self,var,iter,mini,maxi):
            h = self.file.Get(var+"_"+str(iter))
            for binx in range(1,h.GetNbinsX()+1):
                for biny in range(1,h.GetNbinsY()+1):
                    for binz in range(1,h.GetNbinsZ()+1):
                        bin = h.GetBin(binx,biny,binz)
                        if h.GetBinContent(bin)>maxi or h.GetBinContent(bin)<mini:
                            print 'Extreme Values for (pt,eta,phi) = ',1/h.GetXaxis().GetBinCenter(binx),h.GetYaxis().GetBinCenter(biny),h.GetZaxis().GetBinCenter(binz)
