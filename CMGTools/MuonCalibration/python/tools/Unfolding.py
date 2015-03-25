import ROOT
import math

def dfx(histo,x,y):
    binx=histo.GetXaxis().FindBin(x)
    biny=histo.GetYaxis().FindBin(y)

    f=histo.GetBinContent(histo.GetBin(binx,biny))

    if f==0:
        return 0

    x=histo.GetXaxis().GetBinCenter(binx)
    if binx+1<=histo.GetNbinsX():
        fp=histo.GetBinContent(histo.GetBin(binx+1,biny))
        xp=histo.GetXaxis().GetBinCenter(binx+1)
    else:
        fp=f
        xp=x
    if binx-1>=1:
        fm=histo.GetBinContent(histo.GetBin(binx-1,biny))
        xm=histo.GetXaxis().GetBinCenter(binx-1)
    else:
        fm=f
        xm=x
    return (1.0/f)*(fp-fm)/(xp-xm)


def dfy(histo,x,y):
    binx=histo.GetXaxis().FindBin(x)
    biny=histo.GetYaxis().FindBin(y)

    f=histo.GetBinContent(histo.GetBin(binx,biny))

    if f==0:
        return 0

    x=histo.GetYaxis().GetBinCenter(biny)
    if biny+1<=histo.GetNbinsY():
        fp=histo.GetBinContent(histo.GetBin(binx,biny+1))
        xp=histo.GetYaxis().GetBinCenter(biny+1)
    else:
        fp=f
        xp=x
    if biny-1>=1:
        fm=histo.GetBinContent(histo.GetBin(binx,biny-1))
        xm=histo.GetYaxis().GetBinCenter(biny-1)
    else:
        fm=f
        xm=x
    return (1.0/f)*(fp-fm)/(xp-xm)
    

def getDerivativeHisto(histo,f):
    histo2 = histo.Clone()
    histo2.SetName("df")
    for i in range(1,histo.GetNbinsX()+1):
        for j in range(1,histo.GetNbinsY()+1):
            bin=histo.GetBin(i,j)
            histo2.SetBinContent(bin,f(histo,histo.GetXaxis().GetBinCenter(i),histo.GetYaxis().GetBinCenter(j)))
    return histo2 


class Unfolding(object):
    def __init__(self,filename):
        self.f=ROOT.TFile(filename)
        self.h = self.f.Get('unfolding')

    def getVal(self,c1,c2,sigma):
        bin = self.h.GetBin(self.h.GetXaxis().FindBin(c1),self.h.GetYaxis().FindBin(c2),self.h.GetZaxis().FindBin(sigma/c1))
        return c1+self.h.GetBinContent(bin)


    def unfoldDataSet(self,data):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            c1=line.find('curvRaw1').getVal()
            c2=line.find('curvRaw2').getVal()
            e1=line.find('massErrRaw1').getVal()
            e2=line.find('massErrRaw2').getVal()
            m=line.find('massRaw').getVal()
            

            bin1 = self.h.GetBin(self.h.GetXaxis().FindBin(c1),self.h.GetYaxis().FindBin(c2),self.h.GetZaxis().FindBin(2*e1/m))
            c1=c1+self.h.GetBinContent(bin1)
            bin2 = self.h.GetBin(self.h.GetXaxis().FindBin(c2),self.h.GetYaxis().FindBin(c1),self.h.GetZaxis().FindBin(2*e2/m))
            c2=c2+self.h.GetBinContent(bin2)

            line.find('curvRaw1').setVal(c1)
            line.find('curvRaw2').setVal(c2)
            newData.add(line)

        return newData




class AnalyticalUnfolding(object):
    def __init__(self,data,binscurv,mincurv,maxcurv):
        self.spectrum = ROOT.TH2F("spectrum","spectrum",binscurv,mincurv,maxcurv,binscurv,mincurv,maxcurv)
        self.ebeSpectrum = ROOT.TProfile("ebeSpectrum","spectrum",2*binscurv,mincurv,maxcurv,0,10)
        self.x = ROOT.RooRealVar("x","x",1/1000.,1./5.)
        self.y = ROOT.RooRealVar("y","y",1./1000.,1./5.)


        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            c1=line.find('curvRaw1').getVal()
            c2=line.find('curvRaw2').getVal()
            e1=2*line.find('massErrRaw1').getVal()*c1/line.find('massRaw').getVal()
            e2=2*line.find('massErrRaw1').getVal()*c1/line.find('massRaw').getVal()

            self.spectrum.Fill(c1,c2)
            self.ebeSpectrum.Fill(c1,e1*e1)
            self.ebeSpectrum.Fill(c2,e2*e2)
            

        self.datahist = ROOT.RooDataHist("datahist","datahist",ROOT.RooArgList(self.x,self.y),self.spectrum)
        self.histpdf = ROOT.RooHistPdf("histpdf","histpdf",ROOT.RooArgSet(self.x,self.y),self.datahist,2)




        self.dx = self.histpdf.derivative(self.x)
        self.dy = self.histpdf.derivative(self.y)



    def ebeDX(self,x):

        binx=self.ebeSpectrum.GetXaxis().FindBin(x)
        f=self.ebeSpectrum.GetBinContent(self.ebeSpectrum.GetBin(binx))
            
        if f==0:
            return 0

        x=self.ebeSpectrum.GetXaxis().GetBinCenter(binx)
        if binx+1<=self.ebeSpectrum.GetNbinsX():
            fp=self.ebeSpectrum.GetBinContent(self.ebeSpectrum.GetBin(binx+1))
            xp=self.ebeSpectrum.GetXaxis().GetBinCenter(binx+1)
        else:
            fp=f
            xp=x
        if binx-1>=1:
            fm=self.ebeSpectrum.GetBinContent(self.ebeSpectrum.GetBin(binx-1))
            xm=self.ebeSpectrum.GetXaxis().GetBinCenter(binx-1)
        else:
            fm=f
            xm=x

        return (fp-fm)/(xp-xm)



    def getVal(self,c1,c2,sigma):
        self.x.setVal(c1)
        self.y.setVal(c2)
        fv = self.histpdf.getVal()
        dfx = self.dx.getVal()
#        import pdb;pdb.set_trace()


        term1 = -sigma*sigma*dfx/fv
        term2 = self.ebeDX(c1)

        print 'terms',term1,term2

        return c1+term1+term2

    def unfoldDataSet(self,data,updateMass = False):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())

        bad= 0
        n =  data.numEntries()

        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            c1=line.find('curvRaw1').getVal()
            c2=line.find('curvRaw2').getVal()
            e1=line.find('massErrRaw1').getVal()
            e2=line.find('massErrRaw2').getVal()
            m=line.find('massRaw').getVal()

            sigma1 = 2*e1*c1/m
            sigma2 = 2*e2*c2/m

            line.find("curvRaw1").setVal(self.getVal(c1,c2,sigma1))
            line.find("curvRaw2").setVal(self.getVal(c2,c1,sigma2))
            print 'before',c1,c2,'after',line.find("curvRaw1").getVal(),line.find("curvRaw2").getVal()

            if updateMass:

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



class AnalyticalUnfolding1D(object):
    def __init__(self,data,binscurv,mincurv,maxcurv):
        self.spectrum = ROOT.TH1F("spectrum","spectrum",binscurv,mincurv,maxcurv)
        self.x = ROOT.RooRealVar("x","x",1/1000.,1./5.)

        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            c1=line.find('curvRaw1').getVal()
            self.spectrum.Fill(c1)

            

        self.datahist = ROOT.RooDataHist("datahist","datahist",ROOT.RooArgList(self.x),self.spectrum)
        self.histpdf = ROOT.RooHistPdf("histpdf","histpdf",ROOT.RooArgSet(self.x),self.datahist,2)
        self.dx = self.histpdf.derivative(self.x)



    def getVal(self,c,sigma):
        self.x.setVal(c)
        fvx = self.histpdf.getVal()
        dfx = self.dx.getVal()

        if not (dfx==0 or fvx==0 ):
            fodf = fvx/dfx
            sqr = (1-4*sigma*sigma/(fodf*fodf)) 
            if sqr<0:
                print 'negative square root'
                return 0
                    
        return c-0.5*fodf*(1-math.sqrt(sqr))



        

    def unfoldDataSet(self,data,updateMass=False):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())

        bad= 0
        n =  data.numEntries()

        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            c1=line.find('curvRaw1').getVal()
            c2=line.find('curvRaw2').getVal()
            e1=line.find('massErrRaw1').getVal()
            e2=line.find('massErrRaw2').getVal()
            m=line.find('massRaw').getVal()

            sigma1 = 2*e1*c1/m
            sigma2 = 2*e2*c2/m

            self.x.setVal(c1)
            fvx = self.histpdf.getVal()
            dfx = self.dx.getVal()

            self.x.setVal(c2)
            fvy = self.histpdf.getVal()
            dfy = self.dx.getVal()

            if not (dfx==0 or fvx==0 ):
                fodf = fvx/dfx
                sqr = (1-4*sigma1*sigma1/(fodf*fodf)) 
                if sqr<0:
                    bad=bad+1
                    continue
                    
                uc = c1-0.5*fodf*(1-math.sqrt(sqr))
                line.find('curvRaw1').setVal(uc)

            if not (dfy==0 or fvy==0 ):
                fodf = fvy/dfy
                sqr = (1-4*sigma2*sigma2/(fodf*fodf)) 
                if sqr<0:
                    bad=bad+1
                    continue
                    
                uc = c2-0.5*fodf*(1-math.sqrt(sqr))
                line.find('curvRaw2').setVal(uc)

            if updateMass:

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
        print 'bad=',bad,'all',2*n
        return newData




class AnalyticalUnfoldingMass(object):
    def __init__(self,data,binsmass,minmass,maxmass):
        self.spectrum = ROOT.TH1F("spectrum","spectrum",binsmass,minmass,maxmass)
        self.x = ROOT.RooRealVar("x","x",1.,120)

        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            m=line.find('massRaw').getVal()
            self.spectrum.Fill(m)

            

        self.datahist = ROOT.RooDataHist("datahist","datahist",ROOT.RooArgList(self.x),self.spectrum)
        self.histpdf = ROOT.RooHistPdf("histpdf","histpdf",ROOT.RooArgSet(self.x),self.datahist,2)
        self.dx = self.histpdf.derivative(self.x)



    def getVal(self,m,sigma):
        self.x.setVal(m)
        fvx = self.histpdf.getVal()
        dfx = self.dx.getVal()

        if not (dfx==0 or fvx==0 ):
            fodf = fvx/dfx
            sqr = (1-4*sigma*sigma/(fodf*fodf)) 
            if sqr<0:
                print 'negative square root'
                return 0
                    
        return m-0.5*fodf*(1-math.sqrt(sqr))



        

    def unfoldDataSet(self,data):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())

        bad= 0
        n =  data.numEntries()

        for evt in range(0,data.numEntries()):
            line = data.get(evt)

            m=line.find('massRaw').getVal()
            e=line.find('massErrRaw').getVal()
            line.find('massRaw').setVal(self.getVal(m,e))
            newData.add(line)
        print 'bad=',bad,'all',2*n
        return newData
