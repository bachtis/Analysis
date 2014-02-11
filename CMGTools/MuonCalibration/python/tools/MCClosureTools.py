import ROOT
import numpy
import math
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
        error1 = g1.GetErrorY(i)
        error2 = g2.GetErrorY(i)
        
        gOut.SetPoint(i, x1,y1-y2)
        gOut.SetPointError(i, 0.0,math.sqrt(error1*error1+error2*error2))
    return gOut    

class MCClosureTools (object):

    def __init__(self,filename ='/data/bachtis/sandbox/W/Z.root' ,tree = "ZTreeProducer",flatten = False,cut = "nMuons==2&&(MuNeg_pt>20&&MuPos_pt>20)&&MuPosIsTightAndIso&&MuNegIsTightAndIso&&MuPos_pt<200&&MuNeg_pt<200&&Z_mass>70&&Z_mass<120"):
        self.f = ROOT.TFile(filename)
        self.cache = ROOT.TFile("/tmp/bachtis__cacheG__.root","RECREATE")
        treeIn = self.f.Get(tree)
        self.tree = treeIn.CopyTree(cut)

        self.cut=""
        
        if flatten:
            histoPlus = ROOT.TH1F("histoPlus","",100,1/100., 1/25.)
            histoMinus = ROOT.TH1F("histoMinus","",100,1/100., 1/25.)
            self.tree.Draw("1./MuPosGenStatus1_pt>>histoPlus")
            self.tree.Draw("1./MuNegGenStatus1_pt>>histoMinus")
            histoPlus.Scale(1./histoPlus.Integral())
            histoMinus.Scale(1./histoMinus.Integral())
            n = numpy.zeros(1,float)
            n2 = numpy.zeros(1,float)
            branch = self.tree.Branch("weightPos",n,'weightPos/D')
            branch2 = self.tree.Branch("weightNeg",n2,'weightPos/D')

            counter=0
            for event in self.tree:
                if counter>10000:
                    print "Done 10000 events"
                    counter = 0
                counter=counter+1
                bin = histoPlus.GetXaxis().FindBin(1./event.MuPosGenStatus1_pt)
                n[0] = 1./histoPlus.GetBinContent(bin)
                bin = histoMinus.GetXaxis().FindBin(1./event.MuNegGenStatus1_pt)
                n2[0] = 1./histoMinus.GetBinContent(bin)
                branch.Fill()
                branch2.Fill()
                
                
                
        

    def getMeanAndSpread(self,vary,varx,bins,min,max,miny,maxy,cut = ""):
        offset = (max-min)/bins
        gmean = ROOT.TGraphErrors()
        gres  = ROOT.TGraphErrors()
        for i in range(1,bins+1):
            low  = min+(i-1)*offset
            high = min+i*offset
            h=ROOT.TH1F("h","h",500,miny,maxy)
            self.tree.Draw(vary+'>>h',cut+'*('+varx+'>'+str(low)+'&&'+varx+'<'+str(high)+")")

            gaus = ROOT.TF1("g","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",miny,maxy)
            gaus.SetParameter(0,1)
            gaus.SetParameter(1,0)
            gaus.SetParameter(2,1)
            h.Fit(gaus)
            mean = gaus.GetParameter(1)
            meanError =gaus.GetParError(1) 
            sigma = gaus.GetParameter(2)
            sigmaError =gaus.GetParError(2) 
            
            print i,h.GetMean(),h.GetRMS(),mean,sigma

#            gmean.SetPoint(i-1,low+(high-low)/2.,mean)
#            gmean.SetPointError(i-1,(high-low)/2.,meanError)
#            gres.SetPoint(i-1,low+(high-low)/2.,sigma)
#            gres.SetPointError(i-1,(high-low)/2.,sigmaError)

            gmean.SetPoint(i-1,low+(high-low)/2.,h.GetMean())
            gmean.SetPointError(i-1,(high-low)/2.,h.GetMeanError())
            gres.SetPoint(i-1,low+(high-low)/2.,h.GetRMS())
            gres.SetPointError(i-1,(high-low)/2.,h.GetRMSError())

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


    def getMeanAndSpreadUnfolded(self,vary,varx,bins,min,max,miny,maxy,cut = ""):
        print 'will make the plot 1/', vary, '- 1/',varx,'vs',varx
        offset = (max-min)/bins
        gmean = ROOT.TGraphErrors()
        gres  = ROOT.TGraphErrors()

        #first draw the spectrum
        spectrum = ROOT.TH1F("spectrum","spectrum",bins,min,max)
        self.tree.Draw('1./'+varx+'>>spectrum',cut+'*(1./'+varx+'>'+str(min)+'&&1./'+varx+'<'+str(max)+")")
        
        for i in range(1,bins+1):
            low  = min+(i-1)*offset
            high = min+i*offset
            #for each bin loop  draw the resolution

            hh=ROOT.TH1F("hh","h",500,-0.01,0.01)
            self.tree.Draw('1./'+varx+'-'+'1./'+vary+'>>hh',cut+'*(1./'+varx+'>'+str(low)+'&&1./'+varx+'<'+str(high)+")")
            resolution = hh.GetRMS()

            #now unfold and fill
            cache2 = ROOT.TFile("__cacheGU__.root","RECREATE")
            cache2.cd()
            reduced = self.tree.CopyTree(cut+'*(1./'+varx+'>'+str(low)+'&&1./'+varx+'<'+str(high)+")")
            xrec = low+(high-low)/2.
            hhh=ROOT.TH1F("hhh","h",500,miny,maxy)

            for event in reduced:
                x=1./getattr(event,varx)
                y=1./getattr(event,vary)

                derivative = self.getDerivativeOverContent(spectrum,x)
                unfoldedx = x+resolution*resolution*derivative
                hhh.Fill(y-unfoldedx)
            cache2.Close()    
                
#            if xrec>1./50:    
#                import pdb;pdb.set_trace()
#            print xrec
            gmean.SetPoint(i-1,xrec,hhh.GetMean())
            gmean.SetPointError(i-1,0.0,hhh.GetMeanError())
            gres.SetPoint(i-1,xrec,hhh.GetRMS())
            gres.SetPointError(i-1,0.0,hhh.GetRMSError())

        return gmean,gres    
            


