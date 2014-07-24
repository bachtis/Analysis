import ROOT
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.MCClosureTools import MCClosureTools

class Unfolding(object):
    def __init__(self,pmap,tree,bins,min,max,sign='pos'):
        w = ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        self.unfoldProf={}
        self.unfoldMap={}
        if sign=='pos':
            curv='curvRaw1'
            eta='etaRaw1'
            phi='phiRaw1'
            massErr='massErrRaw1'
            gencurv='curvGenRaw1'
        else:    
            curv='curvRaw2'
            eta='etaRaw2'
            phi='phiRaw2'
            massErr='massErrRaw2'
            gencurv='curvGenRaw2'

        closureTools = MCClosureTools()    
        random = ROOT.TRandom3(101082)    
        offset=(max-min)/bins
        for i in range(0,pmap.bins_curv()+2):
            for j in range(0,pmap.bins_eta()+2):
                for k in range(0,pmap.bins_phi()+2):
                    bin = pmap.bin(i,j,k)
                    self.unfoldProf[bin]=ROOT.TProfile("p"+str(bin),"p",bins,min,max,min*1.2,max*1.2)
        for i in range(0,tree.numEntries()):
            line =tree.get(i)
            etaVal = line.find(eta).getVal()    
            phiVal = line.find(phi).getVal()    
            genCurvVal = line.find(gencurv).getVal()
            resolution = 2*line.find(massErr).getVal()*line.find(curv).getVal()/line.find('massRaw').getVal()

            for b in range(0,bins):
                x=min+b*offset
                curvVal=genCurvVal+random.Gaus(x,resolution)
                bin=pmap.binFromVals(curvVal,etaVal,phiVal)
                self.unfoldProf[bin].Fill(x,genCurvVal-curvVal)

        out=ROOT.TFile("unfoldMaps.root","RECREATE")
        out.cd()
        for bin,profile in self.unfoldProf.iteritems():
            self.unfoldMap[bin]=ROOT.TGraph("u_"+str(bin))
            for i in range(1,profile.GetNbinsX()+1):
                self.unfoldMap[bin].SetPoint(i-1,profile.GetBinContent(i),profile.GetXaxis().GetBinCenter(i))
            self.unfoldMap[bin].Write()
        out.Close()
            
                    
    def unfold3D(self,hist):                 
        for i in range(1,hist.GetNbinsX()+1):
            for j in range(1,hist.GetNbinsY()+1):
                for k in range(1,hist.GetNbinsZ()+1):
                    bin = hist.GetBin(i,j,k)
                    content = hist.GetBinContent(bin)
                    newContent = self.unfoldMap[bin].Eval(content)
                    hist.SetBinContent(bin,newContent)
                    
                    
            
            
        
        
