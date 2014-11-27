import ROOT
import os
import math


class LineshapeFitter(object):
    def __init__(self,workspace):
        self.w = workspace
        self.constraints=[]
        self.poi=[]
        self.observables=[]
        self.nuisances=[]



    def save(self):
        self.w.defineSet('nuisances',','.join(self.nuisances))
        self.w.defineSet('constraints',','.join(self.constraints))
        self.w.defineSet('poi',','.join(self.poi))
        self.w.defineSet('observables',','.join(self.observables))


    def load(self):
        for setName in ['constraints','poi','observables','nuisances']:
            set = self.w.set(setName)
            if set.getSize()==0:
                continue
            iter=set.createIterator()
            for i in range(0,set.getSize()):
                arg=iter.Next()
                getattr(self,setName).append(arg.GetName())


            
    def workspace(self):
        return self.w

    def addObservable(self,obs):
        self.observables.append(obs)

    def addPOI(self,name,mini,maxi):
        self.w.factory(name+'['+str(mini)+','+str(maxi)+']')
        self.poi.append(name)


    def addUniformConstraint(self,poi):
        self.w.factory('RooUniform::'+poi+'Constraint('+poi+')')
        self.constraints.append(poi+'Constraint')

    def addGaussianConstraint(self,poi):
        self.w.factory('RooGaussian::'+poi+'Constraint('+poi+','+str(self.w.var(poi).getVal())+','+str(6*self.w.var(poi).getError())+')')
        self.constraints.append(poi+'Constraint')
        
        
    def formula(self,name,formula,dependents):
        self.w.factory('expr::'+name+"('"+formula+"',"+dependents+')')


    def importData(self,data):
        nom = data.GetName()
        getattr(self.w,'import')(data)
        self.w.data(nom).SetName("data")


    def buildZModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[1,0.8,1.5]')
        self.w.factory('error2[1.0,0.8,1.5]')
        self.poi.append('scale')
        
#        pdf = ROOT.RooGaussianSumPdfDouble(name,name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma2.Class(),1)

#        pdf = ROOT.RooGaussianSumPdfDouble(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')


        pdf = ROOT.RooGaussianSumPdfWithSigma2(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw1','massErrRaw2')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0.0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])




    def buildJDualModelMAT(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        hPP = ROOT.TH1F("hPP","hP",100,0.005,0.2)
        hPP.Sumw2()
        hPM = ROOT.TH1F("hPM","hM",100,0.005,0.2)
        hPM.Sumw2()

        hMP = ROOT.TH1F("hMP","hP",100,0.005,0.2)
        hMP.Sumw2()
        hMM = ROOT.TH1F("hMM","hM",100,0.005,0.2)
        hMM.Sumw2()

        hGPP = ROOT.TH1F("hGPP","hP",100,0.005,0.2)
        hGPP.Sumw2()
        hGPM = ROOT.TH1F("hGPM","hM",100,0.005,0.2)
        hGPM.Sumw2()

        hGMP = ROOT.TH1F("hGMP","hP",100,0.005,0.2)
        hGMP.Sumw2()
        hGMM = ROOT.TH1F("hGMM","hM",100,0.005,0.2)
        hGMM.Sumw2()


        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hPP.Fill(line.find("curvRaw1").getVal())
            hPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hMM.Fill(line.find("curvRaw2").getVal())
            hMP.Fill(line.find("curvRaw1").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hGPP.Fill(line.find("curvRaw1").getVal())
            hGPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hGMM.Fill(line.find("curvRaw2").getVal())
            hGMP.Fill(line.find("curvRaw1").getVal())


        pp = hPP.GetMean()
        ppErr = hPP.GetMeanError()
        pm = hPM.GetMean()
        pmErr = hPM.GetMeanError()
        mp = hMP.GetMean()
        mpErr = hMP.GetMeanError()
        mm = hMM.GetMean()
        mmErr = hMM.GetMeanError()

        gpp = hGPP.GetMean()
        gppErr = hGPP.GetMeanError()
        gpm = hGPM.GetMean()
        gpmErr = hGPM.GetMeanError()
        gmp = hGMP.GetMean()
        gmpErr = hGMP.GetMeanError()
        gmm = hGMM.GetMean()
        gmmErr = hGMM.GetMeanError()

        

        ppMmm = pp-mm

        pmMmp = mp-pm

        gppMmm = gpp-gmm

        gpmMmp = gmp-gpm

        ppMmmErr = math.sqrt(ppErr*ppErr+mmErr*mmErr+gppErr*gppErr+gmmErr*gmmErr)
        pmMmpErr = math.sqrt(pmErr*pmErr+mpErr*mpErr+gpmErr*gpmErr+gmpErr*gmpErr)


        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[0.,-200e-6,200e-6]')
        self.w.factory('a[1.,0.995,1.005]')
        self.w.factory('b[0.,-200e-6,200e-6]')


        self.w.factory("expr::plus('sqrt((a+b/curvRaw1+c/curvRaw1)*(a+b/curvRaw2-u/curvRaw2))',u,c,a,b,curvRaw1,curvRaw2)")
        self.w.factory("expr::minus('sqrt((a+b/curvRaw2-c/curvRaw2)*(a+b/curvRaw1+u/curvRaw1))',u,c,a,b,curvRaw1,curvRaw2)")

        self.w.factory('error[0.01,0.003,0.1]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')


        getattr(self.w,'import')(dataPU,ROOT.RooFit.Rename('dataPlus'))
        getattr(self.w,'import')(dataMU,ROOT.RooFit.Rename('dataMinus'))

                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')

        self.w.factory("expr::ppMmm('"+str(gppMmm)+"+2*c',c)")
        self.w.factory("expr::pmMmp('"+str(gpmMmp)+"+2*u',u)")


        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(ppMmm)+',ppMmm,'+str(ppMmmErr)+'),RooGaussian('+str(pmMmp)+',pmMmp,'+str(pmMmpErr)+'))')

        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('curvRaw1'),self.w.var('curvRaw2')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')))

        getattr(self.w,'import')(data)




    def buildJDualModelBinned(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        hPP = ROOT.TH1F("hPP","hP",100,0.005,0.2)
        hPP.Sumw2()
        hPM = ROOT.TH1F("hPM","hM",100,0.005,0.2)
        hPM.Sumw2()

        hMP = ROOT.TH1F("hMP","hP",100,0.005,0.2)
        hMP.Sumw2()
        hMM = ROOT.TH1F("hMM","hM",100,0.005,0.2)
        hMM.Sumw2()

        hGPP = ROOT.TH1F("hGPP","hP",100,0.005,0.2)
        hGPP.Sumw2()
        hGPM = ROOT.TH1F("hGPM","hM",100,0.005,0.2)
        hGPM.Sumw2()

        hGMP = ROOT.TH1F("hGMP","hP",100,0.005,0.2)
        hGMP.Sumw2()
        hGMM = ROOT.TH1F("hGMM","hM",100,0.005,0.2)
        hGMM.Sumw2()


        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hPP.Fill(line.find("curvRaw1").getVal())
            hPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hMM.Fill(line.find("curvRaw2").getVal())
            hMP.Fill(line.find("curvRaw1").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hGPP.Fill(line.find("curvRaw1").getVal())
            hGPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hGMM.Fill(line.find("curvRaw2").getVal())
            hGMP.Fill(line.find("curvRaw1").getVal())


        pp = hPP.GetMean()
        ppErr = hPP.GetMeanError()
        pm = hPM.GetMean()
        pmErr = hPM.GetMeanError()
        mp = hMP.GetMean()
        mpErr = hMP.GetMeanError()
        mm = hMM.GetMean()
        mmErr = hMM.GetMeanError()

        gpp = hGPP.GetMean()
        gppErr = hGPP.GetMeanError()
        gpm = hGPM.GetMean()
        gpmErr = hGPM.GetMeanError()
        gmp = hGMP.GetMean()
        gmpErr = hGMP.GetMeanError()
        gmm = hGMM.GetMean()
        gmmErr = hGMM.GetMeanError()

        

        ppMmm = pp-mm

        pmMmp = mp-pm

        gppMmm = gpp-gmm

        gpmMmp = gmp-gpm

        ppMmmErr = math.sqrt(ppErr*ppErr+mmErr*mmErr+gppErr*gppErr+gmmErr*gmmErr)
        pmMmpErr = math.sqrt(pmErr*pmErr+mpErr*mpErr+gpmErr*gpmErr+gmpErr*gmpErr)


        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[0.,-200e-6,200e-6]')
        self.w.factory('a[1.,0.995,1.005]')
        self.w.factory('b[0.,-0.05,0.05]')


        self.w.factory("expr::plus('sqrt((a-b*curvRaw1+c/curvRaw1)*(a-b*curvRaw2-u/curvRaw2))',u,c,a,b,curvRaw1,curvRaw2)")
        self.w.factory("expr::minus('sqrt((a-b*curvRaw2-c/curvRaw2)*(a-b*curvRaw1+u/curvRaw1))',u,c,a,b,curvRaw1,curvRaw2)")

        self.w.factory('error[0.01,0.003,0.1]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')



                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')

        self.w.factory("expr::ppMmm('"+str(gppMmm)+"+2*c',c)")
        self.w.factory("expr::pmMmp('"+str(gpmMmp)+"+2*u',u)")


        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(ppMmm)+',ppMmm,'+str(ppMmmErr)+'),RooGaussian('+str(pmMmp)+',pmMmp,'+str(pmMmpErr)+'))')


        ####
        ####
        ####


        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var('curvRaw1'),self.w.var('curvRaw2'),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var('curvRaw1'),self.w.var('curvRaw2'),self.w.var("weight")),"weight")

        dataPB = self.convertToBinned3D(dataPU,['massRaw','curvRaw1','curvRaw2'],50,20,20)
        dataMB = self.convertToBinned3D(dataMU,['massRaw','curvRaw1','curvRaw2'],50,20,20)
        for i in range(0,dataPB.numEntries()):
            dataP.add(dataPB.get(i),dataPB.weight())
        for i in range(0,dataMB.numEntries()):
            dataM.add(dataMB.get(i),dataMB.weight())
            

        getattr(self.w,'import')(dataP,ROOT.RooFit.Rename('dataPlus'))
        getattr(self.w,'import')(dataM,ROOT.RooFit.Rename('dataMinus'))


        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('curvRaw1'),self.w.var('curvRaw2'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))

        getattr(self.w,'import')(data)


    def buildJModelAL(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[1.0,0.5,1.5]')
        self.w.factory('error2[1,0.5,1.5]')
        self.poi.append('scale')
        
#        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma2.Class(),1)

        pdf = ROOT.RooGaussianSumPdfWithSigma2(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw1','massErrRaw2')
#        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildJModelCBParam2(self,name,isMC = False):
        self.w.factory("a[0.00001,0.0,0.1]")
        self.w.factory("b[0.00001,0.0,0.1]")
        self.poi.append('a')
        self.poi.append('b')

        self.w.factory("expr::sigma('sqrt(a*a/(curvRaw1*curvRaw1*curvRaw1*curvRaw1)+a*a/(curvRaw2*curvRaw2*curvRaw2*curvRaw2)+b*b/(curvRaw1*curvRaw1)+b*b/(curvRaw2*curvRaw2))',a,b,curvRaw1,curvRaw2)")
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory('RooCBShape::'+name+postfix+'(massRaw,mass[3.09,3.10],sigma,alpha[1,0.05,20],n[5])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NSIG','NBKG','mass'])
        else:    
            self.nuisances.extend(['mass'])

    def buildJModelCBMAT(self,name):
        self.poi.append('m')
        self.poi.append('u')

        self.w.factory("u[1.0,0.995,1.005]")
        self.w.factory("m[0.0,-0.02,0.02]")
        self.w.factory("expr::mass('3.09423*sqrt((u-1+1.0/(1.0+m*curvRaw1))*(u-1+1.0/(1.0+m*curvRaw2)))',m,u,curvRaw1,curvRaw2)")
        self.w.factory('RooCBShape::modelSig(massRaw,mass,error[0.01,0.001,0.2],alpha[1,0,10],n[5,0,30])')
        self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])

    def buildJModelCB(self,name):
        self.poi.append('scale')
        massGEN=3.093165684198312
        self.w.factory("scale[1.0,0.992,1.008]")
        self.w.factory("expr::mass('3.093165684198312*scale',scale)")
        self.w.factory('RooCBShape::modelSig(massRaw,mass,error[0.01,0.001,0.2],alpha[1,-10,10],n[5,0,10])')
        self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])



    def buildZDualModel(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        hPP = ROOT.TH1F("hPP","hP",100,0.005,0.2)
        hPP.Sumw2()
        hPM = ROOT.TH1F("hPM","hM",100,0.005,0.2)
        hPM.Sumw2()

        hMP = ROOT.TH1F("hMP","hP",100,0.005,0.2)
        hMP.Sumw2()
        hMM = ROOT.TH1F("hMM","hM",100,0.005,0.2)
        hMM.Sumw2()

        hGPP = ROOT.TH1F("hGPP","hP",100,0.005,0.2)
        hGPP.Sumw2()
        hGPM = ROOT.TH1F("hGPM","hM",100,0.005,0.2)
        hGPM.Sumw2()

        hGMP = ROOT.TH1F("hGMP","hP",100,0.005,0.2)
        hGMP.Sumw2()
        hGMM = ROOT.TH1F("hGMM","hM",100,0.005,0.2)
        hGMM.Sumw2()


        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hPP.Fill(line.find("curvRaw1").getVal())
            hPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hMM.Fill(line.find("curvRaw2").getVal())
            hMP.Fill(line.find("curvRaw1").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hGPP.Fill(line.find("curvRaw1").getVal())
            hGPM.Fill(line.find("curvRaw2").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hGMM.Fill(line.find("curvRaw2").getVal())
            hGMP.Fill(line.find("curvRaw1").getVal())


        pp = hPP.GetMean()
        ppErr = hPP.GetMeanError()
        pm = hPM.GetMean()
        pmErr = hPM.GetMeanError()
        mp = hMP.GetMean()
        mpErr = hMP.GetMeanError()
        mm = hMM.GetMean()
        mmErr = hMM.GetMeanError()

        gpp = hGPP.GetMean()
        gppErr = hGPP.GetMeanError()
        gpm = hGPM.GetMean()
        gpmErr = hGPM.GetMeanError()
        gmp = hGMP.GetMean()
        gmpErr = hGMP.GetMeanError()
        gmm = hGMM.GetMean()
        gmmErr = hGMM.GetMeanError()

        

        ppMmm = pp-mm

        pmMmp = mp-pm

        gppMmm = gpp-gmm

        gpmMmp = gmp-gpm

        ppMmmErr = math.sqrt(ppErr*ppErr+mmErr*mmErr+gppErr*gppErr+gmmErr*gmmErr)
        pmMmpErr = math.sqrt(pmErr*pmErr+mpErr*mpErr+gpmErr*gpmErr+gmpErr*gmpErr)


        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[0.,-200e-6,200e-6]')
        self.w.factory('a[1.,0.995,1.005]')
        self.w.factory('b[0.,-200e-6,200e-6]')
        self.w.factory("expr::plus('sqrt((a+b/curvRaw1+c/curvRaw1)*(a+b/curvRaw2-u/curvRaw2))',u,c,a,b,curvRaw1,curvRaw2)")
        self.w.factory("expr::minus('sqrt((a+b/curvRaw2-c/curvRaw2)*(a+b/curvRaw1+u/curvRaw1))',u,c,a,b,curvRaw1,curvRaw2)")
        self.w.factory('error[1.0,0.7,2.0]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')


        getattr(self.w,'import')(dataPU,ROOT.RooFit.Rename('dataPlus'))
        getattr(self.w,'import')(dataMU,ROOT.RooFit.Rename('dataMinus'))

                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')

        self.w.factory("expr::ppMmm('"+str(gppMmm)+"+2*c',c)")
        self.w.factory("expr::pmMmp('"+str(gpmMmp)+"+2*u',u)")


        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(ppMmm)+',ppMmm,'+str(ppMmmErr)+'),RooGaussian('+str(pmMmp)+',pmMmp,'+str(pmMmpErr)+'))')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('curvRaw1'),self.w.var('curvRaw2')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')))
        getattr(self.w,'import')(data)




    def convertToBinned(self,data,variable,bins=100):
        self.w.var(variable).setBins(bins)
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variable)),data)

    def convertToBinned2D(self,data,variables,binsx,binsy):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1])),data)

    def buildZDualModel2(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        ###FIND THE MEANS first#########################
        hP = ROOT.TH1F("hP","hP",100,0.005,0.05)
        hP.Sumw2()
        hM = ROOT.TH1F("hM","hM",100,0.005,0.05)
        hM.Sumw2()

        hPG = ROOT.TH1F("hPG","hP",100,0.005,0.05)
        hPG.Sumw2()

        hMG = ROOT.TH1F("hMG","hM",100,0.005,0.05)
        hMG.Sumw2()

        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hP.Fill(line.find("curvRaw1").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hM.Fill(line.find("curvRaw2").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hPG.Fill(line.find("curvRaw1").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hMG.Fill(line.find("curvRaw2").getVal())

        plus = hP.GetMean()
        plusErr = hP.GetMeanError()
        minus = hM.GetMean()
        minusErr = hM.GetMeanError()
        plusGEN = hPG.GetMean()
        plusGENErr = hPG.GetMeanError()
        minusGEN = hMG.GetMean()
        minusGENErr = hMG.GetMeanError()

        PminusM = plus-minus
        PminusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr)
        PplusM = plus+minus
        PplusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr)
        PminusMG = plusGEN-minusGEN
        PminusMGErr = math.sqrt(plusGENErr*plusGENErr+minusGENErr*minusGENErr)
        PplusMG = plusGEN+minusGEN
        PplusMGErr = math.sqrt(plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        print 'plus - minus = ',PminusM,' +-', PminusMErr
        print 'plus + minus = ',PplusM,' +-', PplusMErr
        print 'plus - minus [GEN] = ',PminusMG,' +-', PminusMGErr
        print 'plus + minus [GEN] = ',PplusMG,' +-', PplusMGErr


        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[0.,-15e-6,15e-6]')
        self.w.factory("expr::plus('(1+(u+c)/"+str(plus)+")',u,c)")
        self.w.factory("expr::minus('(1+(u-c)/"+str(minus)+")',u,c)")

        self.w.factory('error[1,0.7,1.3]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')

        dataPlus = self.convertToBinned(dataPU,'massRaw',100)
        dataMinus = self.convertToBinned(dataMU,'massRaw',100)
        
        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)




        self.w.factory("expr::plusMinus1('"+str(PminusMG)+"+2*c',c)")
        self.w.factory("expr::plusMinus2('"+str(PplusMG)+"+2*u',u)")

        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(PminusM)+',plusMinus1,'+str(PminusMErr)+'),RooGaussian('+str(PplusM)+',plusMinus2,'+str(PplusMErr)+'))')
        
            



    def buildZDualModel3(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        ###FIND THE MEANS first#########################
        hP = ROOT.TH1F("hP","hP",100,0.005,0.05)
        hP.Sumw2()
        hM = ROOT.TH1F("hM","hM",100,0.005,0.05)
        hM.Sumw2()

        hPG = ROOT.TH1F("hPG","hP",100,0.005,0.05)
        hPG.Sumw2()

        hMG = ROOT.TH1F("hMG","hM",100,0.005,0.05)
        hMG.Sumw2()

        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hP.Fill(line.find("curvRaw1").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hM.Fill(line.find("curvRaw2").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hPG.Fill(line.find("curvRaw1").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hMG.Fill(line.find("curvRaw2").getVal())

        plus = hP.GetMean()
        plusErr = hP.GetMeanError()
        minus = hM.GetMean()
        minusErr = hM.GetMeanError()
        plusGEN = hPG.GetMean()
        plusGENErr = hPG.GetMeanError()
        minusGEN = hMG.GetMean()
        minusGENErr = hMG.GetMeanError()

        PminusM = plus-minus
        PminusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr)
        PplusM = plus+minus
        PplusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr)
        PminusMG = plusGEN-minusGEN
        PminusMGErr = math.sqrt(plusGENErr*plusGENErr+minusGENErr*minusGENErr)
        PplusMG = plusGEN+minusGEN
        PplusMGErr = math.sqrt(plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        print 'plus - minus = ',PminusM,' +-', PminusMErr
        print 'plus + minus = ',PplusM,' +-', PplusMErr
        print 'plus - minus [GEN] = ',PminusMG,' +-', PminusMGErr
        print 'plus + minus [GEN] = ',PplusMG,' +-', PplusMGErr


        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[0.,-400e-6,400e-6]')
        self.w.factory("expr::plus('sqrt(1+(u+c)/"+str(plus)+")',u,c)")
        self.w.factory("expr::minus('sqrt(1+(u-c)/"+str(minus)+")',u,c)")

        self.w.factory('error[1,0.7,2.0]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')

        dataPlus = self.convertToBinned(dataPU,'massRaw',100)
        dataMinus = self.convertToBinned(dataMU,'massRaw',100)
        
        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)




        self.w.factory("expr::plusMinus1('"+str(PminusMG)+"+2*c',c)")
        self.w.factory("expr::plusMinus2('"+str(PplusMG)+"+2*u',u)")


            
        
        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(PminusM)+',plusMinus1,'+str(PminusMErr)+'),RooGaussian('+str(PplusM)+',plusMinus2,'+str(PplusMErr)+'))')


    def buildZDualModel4(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        ###FIND THE MEANS first#########################
        hP = ROOT.TH1F("hP","hP",100,0.005,0.05)
        hP.Sumw2()
        hM = ROOT.TH1F("hM","hM",100,0.005,0.05)
        hM.Sumw2()

        hPG = ROOT.TH1F("hPG","hP",100,0.005,0.05)
        hPG.Sumw2()

        hMG = ROOT.TH1F("hMG","hM",100,0.005,0.05)
        hMG.Sumw2()

        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hP.Fill(line.find("curvRaw1").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hM.Fill(line.find("curvRaw2").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hPG.Fill(line.find("curvRaw1").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hMG.Fill(line.find("curvRaw2").getVal())

        plus = hP.GetMean()
        plusErr = hP.GetMeanError()
        minus = hM.GetMean()
        minusErr = hM.GetMeanError()
        plusGEN = hPG.GetMean()
        plusGENErr = hPG.GetMeanError()
        minusGEN = hMG.GetMean()
        minusGENErr = hMG.GetMeanError()

        PminusM = plus-minus
        PminusMG = plusGEN-minusGEN
        PminusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        PplusM=plus*minus
        PplusMG = plusGEN*minusGEN
        PplusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)



 

        print 'plus - minus = ',PminusM,PminusMG,' +-', PminusMErr
        print 'plus + minus = ',PplusM,PplusMG,' +-', PplusMErr



        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[1.,0.995,1.005]')
        self.w.factory("expr::plus('(u+c/"+str(plus)+")',u,c)")
        self.w.factory("expr::minus('(u-c/"+str(minus)+")',u,c)")

        self.w.factory('error[1,0.7,2.0]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')

        dataPlus = self.convertToBinned(dataPU,'massRaw',100)
        dataMinus = self.convertToBinned(dataMU,'massRaw',100)
        
        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)




        self.w.factory("expr::plusMinus1('"+str(PminusM)+"-2*c',c)")
        self.w.factory("expr::plusMinus2('"+str(PplusM)+"*u',u)")


            
        
        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(PminusMG)+',plusMinus1,'+str(PminusMErr)+'),RooGaussian('+str(PplusMG)+',plusMinus2,'+str(PplusMErr)+'))')




    def buildJDualModel(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        ###FIND THE MEANS first#########################
        hP = ROOT.TH1F("hP","hP",100,1./20.,1./5.)
        hP.Sumw2()
        hM = ROOT.TH1F("hM","hM",100,1./20.,1./5.)
        hM.Sumw2()

        hPG = ROOT.TH1F("hPG","hP",100,1./20.,1./5.)
        hPG.Sumw2()

        hMG = ROOT.TH1F("hMG","hM",100,1./20.,1./5.)
        hMG.Sumw2()

        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hP.Fill(line.find("curvRaw1").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hM.Fill(line.find("curvRaw2").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hPG.Fill(line.find("curvRaw1").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hMG.Fill(line.find("curvRaw2").getVal())

        plus = hP.GetMean()
        plusErr = hP.GetMeanError()
        minus = hM.GetMean()
        minusErr = hM.GetMeanError()
        plusGEN = hPG.GetMean()
        plusGENErr = hPG.GetMeanError()
        minusGEN = hMG.GetMean()
        minusGENErr = hMG.GetMeanError()

        PminusM = plus-minus
        PminusMG = plusGEN-minusGEN
        PminusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        PplusM=plus*minus
        PplusMG = plusGEN*minusGEN
        PplusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        print 'plus - minus = ',PminusM,PminusMG,' +-', PminusMErr
        print 'plus + minus = ',PplusM,PplusMG,' +-', PplusMErr

        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[1.,0.995,1.005]')
#        self.w.factory("expr::plus('3.09173281695*(u+c/"+str(plus)+")',u,c)")
#        self.w.factory("expr::minus('3.09173281695*(u-c/"+str(minus)+")',u,c)")
        self.w.factory("expr::plus('3.09423*(u+c/"+str(plus)+")',u,c)")
        self.w.factory("expr::minus('3.09423*(u-c/"+str(minus)+")',u,c)")

        self.w.factory('error[0.01,0.001,0.5]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')

        dataPlus = self.convertToBinned(dataPU,'massRaw',100)
        dataMinus = self.convertToBinned(dataMU,'massRaw',100)
        
        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



        self.w.factory('RooCBShape::'+name+'Sig1(massRaw,plus,error,alpha[1,0,30],n[3,0,20])')
        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,10])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
        self.w.factory('RooCBShape::'+name+'Sig2(massRaw,minus,error,alpha2[1,0,30],n2[3,0,20])')
        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope2[-1,-8.,10.])')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')
        self.nuisances.extend(['NSIG2','NBKG2'])
        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)
        self.w.factory("expr::plusMinus1('"+str(PminusM)+"-2*c',c)")
        self.w.factory("expr::plusMinus2('"+str(PplusM)+"*2*u',u)")
        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(PminusMG)+',plusMinus1,'+str(PminusMErr)+'),RooGaussian('+str(PplusMG)+',plusMinus2,'+str(PplusMErr)+'))')




    def buildJDualModelKERNEL(self,name,lineshapePlus,lineshapeMinus,dataPU,dataMU):

        ####OK NOW THE PLUS MINUS CONSTRAINS############!
        ###FIND THE MEANS first#########################
        hP = ROOT.TH1F("hP","hP",100,1./20.,1./5.)
        hP.Sumw2()
        hM = ROOT.TH1F("hM","hM",100,1./20.,1./5.)
        hM.Sumw2()

        hPG = ROOT.TH1F("hPG","hP",100,1./20.,1./5.)
        hPG.Sumw2()

        hMG = ROOT.TH1F("hMG","hM",100,1./20.,1./5.)
        hMG.Sumw2()

        for i in range(0,dataPU.numEntries()):
            line = dataPU.get(i)
            hP.Fill(line.find("curvRaw1").getVal())
        for i in range(0,dataMU.numEntries()):
            line = dataMU.get(i)
            hM.Fill(line.find("curvRaw2").getVal())

        for i in range(0,lineshapePlus.numEntries()):
            line = lineshapePlus.get(i)
            hPG.Fill(line.find("curvRaw1").getVal())
        for i in range(0,lineshapeMinus.numEntries()):
            line = lineshapeMinus.get(i)
            hMG.Fill(line.find("curvRaw2").getVal())

        plus = hP.GetMean()
        plusErr = hP.GetMeanError()
        minus = hM.GetMean()
        minusErr = hM.GetMeanError()
        plusGEN = hPG.GetMean()
        plusGENErr = hPG.GetMeanError()
        minusGEN = hMG.GetMean()
        minusGENErr = hMG.GetMeanError()

        PminusM = plus-minus
        PminusMG = plusGEN-minusGEN
        PminusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)

        PplusM=plus*minus
        PplusMG = plusGEN*minusGEN
        PplusMErr = math.sqrt(plusErr*plusErr+minusErr*minusErr+plusGENErr*plusGENErr+minusGENErr*minusGENErr)



 

        print 'plus - minus = ',PminusM,PminusMG,' +-', PminusMErr
        print 'plus + minus = ',PplusM,PplusMG,' +-', PplusMErr



        self.w.factory('c[0,-200e-6,200e-6]')
        self.w.factory('u[1.,0.995,1.005]')
#        self.w.factory("expr::plus('((u+1./(1.+m*"+str(plus)+")+c/"+str(plus)+")',u,m,c)")
#        self.w.factory("expr::minus('((u-1)+1./(1.+m*"+str(minus)+")-c/"+str(minus)+")',u,m,c)")
        self.w.factory("expr::plus('(u+c/"+str(plus)+")',u,c)")
        self.w.factory("expr::minus('(u-c/"+str(minus)+")',u,c)")


        self.w.factory('error[0.01,0.001,2.0]')
        self.w.factory('error2[0.0]')

        self.poi.append('c')

        dataPlus = self.convertToBinned(dataPU,'massRaw',100)
        dataMinus = self.convertToBinned(dataMU,'massRaw',100)
        
        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



                
        pdf1 = ROOT.RooGaussianSumPdf(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,10])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdf(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('error'),self.w.var('error2'),lineshapeMinus,'massRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope)')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'Shape(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')
        
        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)




        self.w.factory("expr::plusMinus1('"+str(PminusM)+"-2*c',c)")
        self.w.factory("expr::plusMinus2('"+str(PplusM)+"*u',u)")


            
        
        self.w.factory('PROD::'+name+'('+name+'Shape,RooGaussian('+str(PminusMG)+',plusMinus1,'+str(PminusMErr)+'),RooGaussian('+str(PplusMG)+',plusMinus2,'+str(PplusMErr)+'))')



    def buildZDualModelUnc(self,name,lineshapePlus,lineshapeMinus,dataPlus,dataMinus):
        self.w.factory('scalePlus[1,0.99,1.01]')
        self.w.factory('scaleMinus[1,0.99,1.01]')


        self.w.factory('error[1,0.8,1.5]')
        self.w.factory('error2[0.0]')



        self.poi.append('scalePlus')
        self.poi.append('scaleMinus')

        ###CONVERT the hists in weighted datasets
        self.w.factory('weight[1,0,100000000]')
        dataP = ROOT.RooDataSet("dataPlus","dataP",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        dataM = ROOT.RooDataSet("dataMinus","dataM",ROOT.RooArgSet(self.w.var("massRaw"),self.w.var("weight")),"weight")
        for i in range(0,dataPlus.numEntries()):
            dataP.add(dataPlus.get(i),dataPlus.weight())
        for i in range(0,dataMinus.numEntries()):
            dataM.add(dataMinus.get(i),dataMinus.weight())
            
        getattr(self.w,'import')(dataP)
        getattr(self.w,'import')(dataM)



                
        pdf1 = ROOT.RooGaussianSumPdfWithSigma(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.var('scalePlus'),self.w.var('error'),self.w.var('error2'),lineshapePlus,'massRaw','massErrRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope','NSIG1','NBKG1'])




    def buildZModelEbE(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[0,-0.5,0.5]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
#        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)

        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildZModelNobkg(self,name, var,dataset):
        self.w.factory('scale[1.0,0.5,1.5]')
        self.w.factory('error1[1,0.5,5.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        self.w.var("scale").setError(0.5)
        
        pdf = ROOT.RooGaussianSumPdf(name,name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf)

    def buildZModelParametric(self,name, var,dataset,polyTerms=1,fourierTerms=1):
        fourierPlus='1'
        fourierMinus='1'
        fourierList=['k']
        self.w.factory('k[1.0,0.0,2.]')
        
        for i in range(0,fourierTerms):
            self.w.factory('A_'+str(i+1)+'[0.0,-1,1]')
            self.w.factory('B_'+str(i+1)+'[0.0,-1,1]')
            fourierPlus=fourierPlus+'+A_'+str(i+1)+"*cos(k*"+str(i+1)+"*phiRaw1+B_"+str(i+1)+")"
            fourierMinus=fourierMinus+'+A_'+str(i+1)+"*cos(k*"+str(i+1)+"*phiRaw2+B_"+str(i+1)+"+"+str(math.pi)+")"

            fourierList.append('A_'+str(i+1))
            fourierList.append('B_'+str(i+1))


        self.w.factory('C_0[1.0,0.99,1.01]')

        polyPlus='C_0'
        polyMinus='C_0'
        polyList=['C_0']
        for i in range(0,polyTerms):
            self.w.factory('C_'+str(i+1)+'[0.0,-0.005,0.005]')
            polyPlus=polyPlus+'+C_'+str(i+1)+"*"+'*'.join(['curvRaw1']*(i+1))
            polyMinus=polyMinus+'+C_'+str(i+1)+"*"+'*'.join(['curvRaw2']*(i+1))
            polyList.append('C_'+str(i+1))


        self.w.factory("expr::fourierPlus('"+fourierPlus+"',"+','.join(fourierList+['phiRaw1'])+")")
        self.w.factory("expr::fourierMinus('"+fourierMinus+"',"+','.join(fourierList+['phiRaw2'])+")")
        self.w.factory("expr::polyPlus('"+polyPlus+"',"+','.join(polyList+['curvRaw1'])+")")
        self.w.factory("expr::polyMinus('"+polyMinus+"',"+','.join(polyList+['curvRaw2'])+")")


        self.w.factory("expr::scaleFactor('sqrt(polyPlus*polyMinus*fourierPlus*fourierMinus)',polyPlus,polyMinus,fourierPlus,fourierMinus)")
        self.w.factory('errorScale[1,0.1,3.]')
        self.w.factory("expr::error1('errorScale*massErrRaw',errorScale,massErrRaw)")
        self.w.factory('error2[0.0]')

        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.function('scaleFactor'),self.w.function('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')

        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])




    def buildZModelParametricKernel(self,name, var,dataset,polyTerms=1,fourierTerms=1):
        fourierPlus='0'
        fourierMinus='0'
        fourierList=['k']
        self.w.factory('k[1.0,0.0,2.]')
        
        for i in range(0,fourierTerms):
            self.w.factory('B_'+str(i+1)+'[0.0,0,4]')

            fourierPlus=fourierPlus+"+cos(k*"+str(i+1)+"*phiRaw1+B_"+str(i+1)+")"
            fourierMinus=fourierMinus+"+cos(k*"+str(i+1)+"*phiRaw2+B_"+str(i+1)+"+"+str(math.pi)+")"
            fourierList.append('B_'+str(i+1))


        self.w.factory('C_0[1.0,0.99,1.01]')

        polyPlus='C_0'
        polyMinus='C_0'
        polyList=['C_0']
        for i in range(0,polyTerms):
            self.w.factory('C_'+str(i+1)+'[0.0,-0.005,0.005]')
            polyPlus=polyPlus+'+C_'+str(i+1)+"*"+'*'.join(['curvRaw1']*(i+1))
            polyMinus=polyMinus+'+C_'+str(i+1)+"*"+'*'.join(['curvRaw2']*(i+1))
            polyList.append('C_'+str(i+1))


        self.w.factory("expr::fourierPlus('"+fourierPlus+"',"+','.join(fourierList+['phiRaw1'])+")")
        self.w.factory("expr::fourierMinus('"+fourierMinus+"',"+','.join(fourierList+['phiRaw2'])+")")
        self.w.factory("expr::polyPlus('"+polyPlus+"',"+','.join(polyList+['curvRaw1'])+")")
        self.w.factory("expr::polyMinus('"+polyMinus+"',"+','.join(polyList+['curvRaw2'])+")")


        self.w.factory("expr::scaleFactor('sqrt(polyPlus*polyMinus*fourierPlus*fourierMinus)',polyPlus,polyMinus,fourierPlus,fourierMinus)")
        self.w.factory('errorScale[1,0.1,3.]')
        self.w.factory("expr::error1('errorScale*massErrRaw',errorScale,massErrRaw)")
        self.w.factory('error2[0.0]')

        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.function('scaleFactor'),self.w.function('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')

        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])



        
    def extend(self,name,pdf,normName):
        self.w.factory(normName+'[0.,1000000]')
        self.w.factory('RooExtendPdf::'+name+'('+pdf+','+normName+')')
        


    def fit(self,model,data, hint = False, verbose =0,snapshot = "result"):
        if verbose==0 :
            ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

        hintResult=None
        fitResult=None

        if self.w.allCats().getSize()>0:
            print 'Weighted data input'
#            for poi in self.poi:
#                self.w.var(poi).setConstant(1)
#            self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.SumW2Error(1))
#            if self.w.var('error').getError()>0.01:
#                self.w.var('sigmaErr0').setVal(self.w.var('error').getError())
#                self.w.var('sigma0').setVal(self.w.var('error').getVal())
#            else:    
#                self.w.var('sigmaErr0').setVal(1.0)
#                self.w.var('sigma0').setVal(0.3)
#
#            for poi in self.poi:
#                self.w.var(poi).setConstant(0)

            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2),ROOT.RooFit.SumW2Error(1))
#            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2),ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(self.w.var("curvRaw1"),self.w.var("curvRaw2"))))

        else:
#            self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(6,0))
#            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(6,0),ROOT.RooFit.Save(1))

            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2),ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(self.w.var("curvRaw1"),self.w.var("curvRaw2"),self.w.var("massErrRaw1"),self.w.var("massErrRaw2"))))
        self.w.saveSnapshot(snapshot,','.join(self.poi))    
        if hintResult is not None:
            print '-------HINT------'
            hintResult.Print()
            
        print '-------FINAL RESULT------'
        if fitResult is not None:
            fitResult.Print()
            print '-COVARIANCE-'
            fitResult.correlationMatrix().Print()
        return (hintResult,fitResult)        



    def fitConditional(self,model,data, hint = False, verbose =0,snapshot = "result"):
        if verbose==0 :
            ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

        hintResult=None
        fitResult=None
        

        fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(1),ROOT.RooFit.NumCPU(8,0),ROOT.RooFit.Offset(1),ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(self.w.var('massErrRaw'),self.w.var("curvRaw1"),self.w.var("curvRaw2"),self.w.var("phiRaw1"),self.w.var("phiRaw2"))))
            



    def convertToBinned3D(self,data,variables,binsx,binsy,binsz):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        self.w.var(variables[2]).setBins(binsz);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1]),self.w.var(variables[2])),data)
