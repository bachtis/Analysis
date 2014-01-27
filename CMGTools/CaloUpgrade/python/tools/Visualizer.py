import ROOT
from array import array
class Visualizer(object):
    def __init__(self,filename):
        self.f = ROOT.TFile(filename)
        self.hcal3D = self.createHCALHistogram()
        self.legend,self.hcalLayers=self.createHCALLayers()

        
    def exit(self):
        self.f.Close()


    def createHCALHistogram(self):
        h = ROOT.TH3D("hcal","hcal",16,-8,8,16,-8,8,5,0,5)
        h.GetXaxis().SetTitle("i #eta")
        h.GetYaxis().SetTitle("i #phi")
        h.GetZaxis().SetTitle("depth")
        return h

    def createHCALLayers(self):
        layers=[]
        legend =ROOT.TLegend(0.8,0.8,0.9,0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        for i in range(1,5+1):
            layers.append(ROOT.TH2D("hcal_"+str(i),"hcal",16,-8,8,16,-8,8))
            layers[-1].GetXaxis().SetTitle("i #eta")
            layers[-1].GetYaxis().SetTitle("i #phi")
            legend.AddEntry(layers[-1],"Layer "+str(i),"l")
        layers[0].SetLineColor(ROOT.kRed)    
        layers[1].SetLineColor(ROOT.kGreen)    
        layers[2].SetLineColor(ROOT.kMagenta)    
        layers[3].SetLineColor(ROOT.kYellow)    
        layers[3].SetLineColor(ROOT.kOrange)    
        layers[3].SetLineColor(ROOT.kAzure)    
        return legend,layers
        

    def process(self,event,shower):
        self.hcal3D.Reset()
        for layer in self.hcalLayers:
            layer.Reset()
            
        tree = self.f.Get("t_"+str(event)+"_"+str(shower))
        iEta=0
        iPhi=0
        for event in tree:
            if event.type>-0.5 and event.type<0.5:
                iEta=event.ieta
                iPhi=event.iphi
                print 'Track energy',event.energy
                
            if event.type>1.5 and event.type<2.5:
                deltaPhi = event.iphi-iPhi
                if deltaPhi>40:
                    deltaPhi = deltaPhi-72
                elif deltaPhi<-40:    
                    deltaPhi = deltaPhi+72

                self.hcal3D.Fill(event.ieta-iEta,deltaPhi,event.depth,event.energy)
                self.hcalLayers[int(event.depth)-1].Fill(event.ieta-iEta,deltaPhi,event.energy)
                print event.ieta,event.iphi,event.depth,event.energy
                
        self.canvases=[(ROOT.TCanvas("c3d",""))]
        self.canvases[-1].cd()
        self.hcal3D.Draw("box")

        self.canvases.append(ROOT.TCanvas("c32d",""))
        self.canvases[-1].cd()
        sortedLayers=sorted(self.hcalLayers,key=lambda x : x.GetMaximum(),reverse=True)
        for i,layer in enumerate(sortedLayers):
            if i==0:
                layer.Draw("box")
            else:    
                layer.Draw("box,same")
                
        self.legend.Draw()        
        

