
import ROOT
import os
import math
from KaMuCa.Derivation.tools.DataSetSplitter import *

class LineshapeFitter(object):
    def __init__(self,workspace):
        self.w = workspace
        self.constraints=[]
        self.poi=[]
        self.observables=[]
        self.nuisances=[]
        self.conditionals=[]


    def save(self):
        self.w.defineSet('nuisances',','.join(self.nuisances))
        self.w.defineSet('constraints',','.join(self.constraints))
        self.w.defineSet('poi',','.join(self.poi))
        self.w.defineSet('observables',','.join(self.observables))
        self.w.defineSet('conditionals',','.join(self.conditionals))


    def load(self):
        for setName in ['constraints','poi','observables','nuisances','conditionals']:
            seti = self.w.set(setName)
            if seti.getSize()==0:
                continue
            iter=seti.createIterator()
            for i in range(0,seti.getSize()):
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
        getattr(self.w,'import')(data,ROOT.RooFit.Rename("data"))





    def buildJModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[0.0,-2e-4,5e-4]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass','massErr')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildYModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[1e-9,1e-10,0.5]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass','massErr')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildJModelSimple(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[0.02,0.001,0.1]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')

        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass')
        getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(name+'Sig'))
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildPullModel(self,name, var,dataset):
        self.w.factory('scale[0,-1,1]')
        self.w.factory('error1[1,0.5,1.5]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooPullSumPdf.Class(),1)
        pdf = ROOT.RooPullSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass')
        getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(name+'Sig'))
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildYModelSimple(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[0.1,0.04,0.3]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')

        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])



    def buildJModelLazy(self,name, var,dataset):

        #create lineshape histogram"
        #--------------------------
        
        lineHisto=ROOT.TH1F("lineshape","lineshape",var.getBinning().numBins(),var.getMin(),var.getMax())
        for i in range(0,dataset.numEntries()):
            line=dataset.get(i)
            lineHisto.Fill(line.find("mass").getVal())


        self.w.factory('scale[0,-0.1,0.1]')
        self.w.factory('error1[0.02,0.001,0.1]')
        getattr(self.w,'importClassCode')(ROOT.DynamicBinnedSmearingPdf.Class(),1)
        pdf = ROOT.DynamicBinnedSmearingPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),lineHisto)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])

    def buildResolutionUnbinned(self,name, var,lineHisto):

        #create lineshape histogram"
        #--------------------------
        self.w.factory('scale[0,-0.2,0.2]')
        self.w.factory('b[0.0000,-2.5e-5,1e-4]')
        self.w.factory("expr::sigma('sqrt(massErr*massErr/(mass*mass)+b)*mass',massErr,mass,b)")

        getattr(self.w,'importClassCode')(ROOT.DynamicBinnedSmearingPdf.Class(),1)
        pdf = ROOT.DynamicBinnedSmearingPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.function('sigma'),lineHisto)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])
        self.conditionals.extend(['massErr'])     




    def buildZModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.99,1.01]')

        self.w.factory('error1[0.0,-1e-4,5e-4]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass','massErr')
        getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(name+'Sig'))
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildZModelNoBKG(self,name, var,dataset):
        self.w.factory('scale[1.0,0.99,1.01]')

        self.w.factory('error1[0.000001,0.0,2e-4]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name,name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass','massErr')
        getattr(self.w,'import')(pdf)



    def buildZModelSimple(self,name, var,dataset):
        self.w.factory('scale[1.0,0.99,1.01]')
        self.w.factory('error1[1,0.6,4.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass')
        getattr(self.w,'import')(pdf,ROOT.RooFit.Rename(name+'Sig'))
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[0,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[0.000000001,0,1000000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])

    def buildZModelSimpleNoBKG(self,name, var,dataset):
        self.w.factory('scale[1.0,0.99,1.01]')
        self.w.factory('error1[1,0.6,4.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        
        pdf = ROOT.RooGaussianSumPdf(name,name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'mass')
        getattr(self.w,'import')(pdf)


    def buildJModelCBParam(self,name,isMC = False):
        self.w.factory("a[50e-6,1e-6,0.001]")
        self.w.factory("b[1e-7,0,1e-6]")
        self.w.factory("c[1e-4,0,1e-3]")
        self.w.factory("d[200,10,1000]")

        self.w.factory("expr::sigma('sqrt(0.5*a+0.25*b/(1+d*c1*c1)+0.25*b/(1+d*c2*c2)+(0.25/(c1*c1))*c+(0.25/(c2*c2))*c)*mass',a,b,c,d,c1,c2,mass)")
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory('RooCBShape::'+name+postfix+'(mass,mean[3.091,3.08,3.12],sigma,alpha[3,0.5,20],n[4.77])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(mass,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NSIG','NBKG','mean','a','b','c','d'])
        else:    
            self.nuisances.extend(['mean','a','b','c','d'])
        self.conditionals.extend(['c1','c2'])     





    def buildSimModelParam(self,dataZ,dataJ,isMC=False):
        self.w.factory("a[50e-6,1e-6,1000e-6]")
        self.w.factory("b[40e-6,0,1000e-6]")
        self.w.factory("c[30e-9,12e-9,5e-7]")
        self.w.factory("d[500,50,10000]")
        self.w.factory("expr::sigmaZ('sqrt(0.5*a+0.25*b/(1+d*c1*c1)+0.25*b/(1+d*c2*c2)+(0.25/(c1*c1))*c+(0.25/(c2*c2))*c)*massZ',a,b,c,d,c1,c2,massZ)")
        self.w.factory("expr::sigmaJ('sqrt(0.5*a+0.25*b/(1+d*c1*c1)+0.25*b/(1+d*c2*c2)+(0.25/(c1*c1))*c+(0.25/(c2*c2))*c)*massJ',a,b,c,d,c1,c2,massJ)")
        self.w.factory('RooVoigtian::modelZSig(massZ,meanZ[91,80,100],zwidth[2.4952],sigmaZ)')
        self.w.factory('RooExponential::modelZBkg(massZ,bkgSlopeZ[-1,-8.,5])')
        self.w.factory('SUM::modelZ(NZSIG[0,100000000]*modelZSig,NZBKG[1,0,1000000]*modelZBkg)')
        self.nuisances.extend(['bkgSlopeZ','NZSIG','NZBKG','meanZ','a','b','c','d'])


        
        if isMC:
            self.w.factory('RooCBShape::modelJ(massJ,meanJ[3.091,3.08,3.12],sigmaJ,alpha[3,0.5,20],n[4.77])')
            self.nuisances.extend(['meanJ'])
        else:
            self.w.factory('RooCBShape::modelJSig(massJ,meanJ[3.091,3.08,3.12],sigmaJ,alpha[3,0.5,20],n[4.77])')
            self.w.factory('RooExponential::modelJBkg(massJ,bkgSlopeJ[-1,-8.,5])')
            self.w.factory('SUM::modelJ(NJSIG[0,100000000]*modelJSig,NJBKG[1,0,1000000]*modelJBkg)')
            self.nuisances.extend(['bkgSlopeJ','NJSIG','NJBKG','meanJ'])

        getattr(self.w,'import')(dataZ,ROOT.RooFit.Rename('dataZ'),ROOT.RooFit.RenameVariable('mass','massZ'))
        getattr(self.w,'import')(dataJ,ROOT.RooFit.Rename('dataJ'),ROOT.RooFit.RenameVariable('mass','massJ'))

        self.w.factory('pois{a,b,c,d}')

        self.w.factory('resonance[Z,J]')
        self.w.factory('SIMUL::model(resonance,Z=modelZ,J=modelJ)')
        
        varset=self.w.data('dataZ').get()
        varset.add(self.w.data('dataJ').get().find('massJ'))

        data=ROOT.RooDataSet('data','data',varset,ROOT.RooFit.Index(self.w.cat('resonance')),ROOT.RooFit.Import("Z",self.w.data('dataZ')),ROOT.RooFit.Import("J",self.w.data('dataJ')))
        getattr(self.w,'import')(data,ROOT.RooFit.Rename('data'))
        self.conditionals.extend(['c1','c2'])     








    def buildJModelCBParamEbE(self,name,isMC = False):
        
        self.w.factory("a[1e-6,-20e-6,300e-6]")
        self.w.factory("b[1e-7,-1e-7,1e-7]")
        self.w.factory("c[10e-6,-300e-6,300e-6]")
        self.w.factory("d[300,20,500000]")

        self.w.factory("expr::sigma('sqrt(massErr*massErr/(mass*mass)+0.5*a+0.25*b/(1+d*c1*c1)+0.25*b/(1+d*c2*c2)+(0.25/(c1*c1))*c+(0.25/(c2*c2))*c)*mass',a,b,c,d,c1,c2,massErr,mass)")
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory('RooCBShape::'+name+postfix+'(mass,mean[3.091,3.08,3.12],sigma,alpha[3,0.5,20],n[4.77])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(mass,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NSIG','NBKG','mean','a','b','c','d'])
        else:    
            self.nuisances.extend(['mean','a','b','c','d'])
        self.conditionals.extend(['c1','c2','massErr'])     


    def buildJModelCBSimple(self,name,isMC = False):
        self.w.factory("error1[0.01,0.005,0.1]")
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory("scale[0.995,1.005]")
        self.w.factory("expr::mass('3.097*scale',scale)")
        self.poi.append('scale')

                    
        self.w.factory('RooCBShape::'+name+postfix+'(mass,mass,error1,alpha[3,0.5,20],n[4.77,0,20])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(mass,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NBKG','error1','alpha','n'])
        else:    
            self.nuisances.extend(['error1','alpha','n'])






    def convertToBinned(self,data,variable,bins=100):
        self.w.var(variable).setBins(bins)
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variable)),data)

    def convertToBinned2D(self,data,variables,binsx,binsy):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1])),data)

        
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
            self.w.var('alpha').setConstant(0)
            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.ConditionalObservables(self.w.set('conditionals')))
            self.w.var('alpha').setConstant(1)
            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Minos(self.w.set('pois')),ROOT.RooFit.NumCPU(1,0),ROOT.RooFit.Save(1),ROOT.RooFit.ConditionalObservables(self.w.set('conditionals')))          
        else:
            if len(self.conditionals)>0:
                fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(8,0),ROOT.RooFit.Save(1),ROOT.RooFit.ConditionalObservables(self.w.set('conditionals')))

            else:
                fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2))
#                for n in self.nuisances:
#                    self.w.var(n).setConstant(1)
#                fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2))



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




    def convertToBinned(self,data,variable,bins=100):
        self.w.var(variable).setBins(bins)
        data.get().find(variable).setBins(bins)

        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(data.get().find(variable)),data)

    def convertToBinned2D(self,data,variables,binsx,binsy):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1])),data)



    def convertToBinned3D(self,data,variables,binsx,binsy,binsz):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        self.w.var(variables[2]).setBins(binsz);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1]),self.w.var(variables[2])),data)
