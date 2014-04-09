import ROOT
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.MCClosureTools import MCClosureTools

class Unfolding(object):
    def __init__(self,pmap,samples,bins,min,max,sign='pos'):
        w = ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        self.unfoldMap={}
        if sign=='pos':
            curv='curvRaw1'
            gencurv='genCurvRaw1'
        else:    
            curv='curvRaw2'
            gencurv='genCurvRaw2'

        closureTools = MCClosureTools()    
        random = ROOT.TRandom3(101082)    
        offset=(max-min)/bins
        for bin,data in samples.iteritems():
            p = ROOT.TProfile("map","map",bins,min,max,min/1.2,max*1.2)
            for i in range(0,bins):
                print 'bin:',i
                shift = min+bin*offset
                smearedData = pmap.smearEbE2D(data,w,shift)
                mean,meanErr,spread,spreadErr = closureTools.getMeanAndSpreadFit(sign,smearedData)
                p.Fill(mean,shift)
            self.unfoldMap[bin] =p
            
                    
    def unfold3D(self,hist):                 
        for i in range(1,hist.GetNbinsX()+1):
            for j in range(1,hist.GetNbinsY()+1):
                for k in range(1,hist.GetNbinsZ()+1):
                    bin = hist.GetBin(i,j,k)
                    content = hist.GetBinContent(bin)
                    binP = self.unfoldMap[bin].FindBin(content)
                    newContent = self.unfoldMap[bin].GetBinContent(binP)
                    hist.SetBinContent(bin,content-newContent)

                
            
            
        
        
