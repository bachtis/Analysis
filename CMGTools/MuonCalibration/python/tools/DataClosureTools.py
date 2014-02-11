import ROOT
import numpy
import math
from array import array
def divideGraphs(g1,g2):
    gOut = ROOT.TGraphErrors()
    x1 = ROOT.Double(0)
    x2 = ROOT.Double(0)
    y1 = ROOT.Double(0)
    y2 = ROOT.Double(0)

    for i in range(0,g1.GetN()):
        g1.GetPoint(i,x1,y1)
        g2.GetPoint(i,x2,y2)
        gOut.SetPoint(i, x1,y1/y2)
    return gOut    

def multiplyGraphs(g1,g2):
    gOut = ROOT.TGraphErrors()
    x1 = ROOT.Double(0)
    x2 = ROOT.Double(0)
    y1 = ROOT.Double(0)
    y2 = ROOT.Double(0)

    for i in range(0,g1.GetN()):
        g1.GetPoint(i,x1,y1)
        g2.GetPoint(i,x2,y2)
        gOut.SetPoint(i, x1,y1*y2)
    return gOut    

def subtractGraphs(g1,g2):
    gOut = ROOT.TGraphErrors()
    x1 = ROOT.Double(0)
    x2 = ROOT.Double(0)
    y1 = ROOT.Double(0)
    y2 = ROOT.Double(0)

    for i in range(0,g1.GetN()):
        g1.GetPoint(i,x1,y1)
        g2.GetPoint(i,x2,y2)
        gOut.SetPoint(i, x1,y1-y2)
    return gOut    

class DataClosureTools (object):

    def __init__(self,filename ='ZMC.root' ,tree = "data"):
        self.f = ROOT.TFile(filename)
        self.cache = ROOT.TFile("__cache__.root","RECREATE")
        
        self.tree = self.f.Get(tree)
        

    def getMeanAndSpread(self,vary,varx,bins,min,max,miny,maxy,cut = ""):
        offset = (max-min)/bins
        gmean = ROOT.TGraphErrors()
        gres  = ROOT.TGraphErrors()

        for i in range(1,bins+1):
            low  = min+(i-1)*offset
            high = min+i*offset
            c ="("+cut+")*("+ varx+'>'+str(low)+'&&'+varx+'<'+str(high)+")"
            mean = self.tree.mean(self.tree.get().find(vary),c)
            sigma = self.tree.sigma(self.tree.get().find(vary),c)
            
            gmean.SetPoint(i-1,low+(high-low)/2.,mean)
            gres.SetPointError(i-1,(high-low)/2.,sigma)

        return gmean,gres    


    def getDerivativeOverContent(self,h,x):
        b = h.GetXaxis().FindBin(x)
        content = h.GetBinContent(b)
        if content==0:
            return 0.0

        if b==1:
            return (h.GetBinContent(b+1)-h.GetBinContent(b))/(content*(h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b)))
        elif b==h.GetNbinsX():
            return (h.GetBinContent(b)-h.GetBinContent(b-1))/(content*(h.GetXaxis().GetBinCenter(b)-h.GetXaxis().GetBinCenter(b-1)))
        else:
            return (h.GetBinContent(b+1)-h.GetBinContent(b-1))/(content*(h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b-1)))

    def getDerivative(self,h,x):
        b = h.GetXaxis().FindBin(x)
        content = h.GetBinContent(b)

        if b==1:
            return (h.GetBinContent(b+1)-h.GetBinContent(b))/((h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b)))
        elif b==h.GetNbinsX():
            return (h.GetBinContent(b)-h.GetBinContent(b-1))/((h.GetXaxis().GetBinCenter(b)-h.GetXaxis().GetBinCenter(b-1)))
        else:
            return (h.GetBinContent(b+1)-h.GetBinContent(b-1))/((h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b-1)))


    def getDerivativeOverContent2D(self,h,x,y,direction = 'X'):
        b={}
        bnext={}
        bprev={}
            
        otherdirection='X'
        if direction =='X':
            otherdirection='Y'


        b['X'] = h.GetXaxis().FindBin(x)
        b['Y'] = h.GetYaxis().FindBin(y)
        content = h.GetBinContent(h.GetBin(b['X'],b['Y']))
        
        if b[direction]==getattr(h,'GetNbins'+direction)():
            bnext[direction]=b[direction]
            bprev[direction]=b[direction]-1
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]
        elif b[direction]==1:
            bnext[direction]=b[direction]+1
            bprev[direction]=b[direction]
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]
        else:
            bnext[direction]=b[direction]+1
            bprev[direction]=b[direction]-1
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]

        bnext2D=h.GetBin(bnext['X'],bnext['Y'])
        bprev2D=h.GetBin(bprev['X'],bprev['Y'])
       
        deriv=(h.GetBinContent(bnext2D)-h.GetBinContent(bprev2D))/(getattr(h,'Get'+direction+'axis')().GetBinCenter(bnext[direction])-getattr(h,'Get'+direction+'axis')().GetBinCenter(bprev[direction]))
        if content>0:
            return deriv/content
        else:
            return 0.0


    def makeHistogram(self,dataset,var,bins,min,max,trim=1):
        dataset.get().find(var).setBins(bins)
        dataset.get().find(var).setMin(min)
        dataset.get().find(var).setMax(max)

        datahist =  ROOT.RooDataHist(dataset.GetName(),'',ROOT.RooArgSet(dataset.get().find(var)),dataset)
        histogram = datahist.createHistogram(var,bins) 
        for i in range(1,trim+1):
            histogram.SetBinContent(i,0)
            histogram.SetBinContent(bins-(i-1),0)
        return histogram
    


    def getBiasFromSpectrum(self,bins,min,max,cut = "",points = 20,step=0.000005,minStep=0.6e-4,trim=1):
        dataPlus = self.tree.reduce(cut+'&&abs(phiRaw1)<0.26')
        dataMinus = self.tree.reduce(cut+'&&abs(phiRaw2)<0.26')
        w = ROOT.RooWorkspace('w','w')
        getattr(w,'import')(dataPlus.get().find('curvRaw1'))
        getattr(w,'import')(dataMinus.get().find('curvRaw2'))
        chi2 = ROOT.TGraph()
        for xpoint in range(0,points):
            offset = minStep+xpoint*step
            w.factory('expr::shift_1_'+str(xpoint)+"('curvRaw1+"+str(offset)+"',curvRaw1)")
            w.factory('expr::shift_2_'+str(xpoint)+"('curvRaw2-"+str(offset)+"',curvRaw2)")

            dataPlus.addColumn(w.function('shift_1_'+str(xpoint)))
            dataMinus.addColumn(w.function('shift_2_'+str(xpoint)))
            histoPlus  = self.makeHistogram(dataPlus,'shift_1_'+str(xpoint),bins,min,max,trim)
            histoMinus = self.makeHistogram(dataMinus,'shift_2_'+str(xpoint),bins,min,max,trim)
            if xpoint==0:
                histoRP = histoPlus.Clone()
                histoRN = histoMinus.Clone()
                histoRP.SetName("pos")
                histoRN.SetName("neg")
                
            histoPlus.Scale(1./histoPlus.Integral())
            histoMinus.Scale(1./histoMinus.Integral())
            
            chi = histoPlus.Chi2Test(histoMinus,"WW,P,CHI2")
            chi2.SetPoint(xpoint,offset,chi)
        return chi2,histoRP,histoRN



    def getBiasFromSpectrumKOL(self,cut = "",points = 20,step=0.000005,minStep=0.6e-4):
        dataPlus = self.tree.reduce(cut+'&&abs(phiRaw1)<0.26')
        dataMinus = self.tree.reduce(cut+'&&abs(phiRaw2)<0.26')
        w = ROOT.RooWorkspace('w','w')
        chi2 = ROOT.TGraph()
        for xpoint in range(0,points):
            print 'testing',xpoint 
            offset = minStep+xpoint*step

            plus=[]
            minus=[]
            for i in range(0,dataPlus.numEntries()):
                line=dataPlus.get(i)
                plus.append(line.find('curvRaw1').getVal()+offset)
            for i in range(0,dataMinus.numEntries()):
                line=dataMinus.get(i)
                minus.append(line.find('curvRaw2').getVal()-offset)

            plus.sort()
            minus.sort()
            prob = ROOT.TMath.KolmogorovTest(len(plus),array('d',plus),len(minus),array('d',minus),"D")    
            chi2.SetPoint(xpoint,offset,prob)
        return chi2
    
