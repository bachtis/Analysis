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
        getattr(self.w,'import')(data)
        self.w.data(nom).SetName("data")




    def buildJModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')
        self.w.factory('error1[1e-9,0,0.02]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
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
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])



    def buildZModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.995,1.005]')

        self.w.factory('error1[1e-9,0,0.2]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildZModelSimple(self,name, var,dataset):
        self.w.factory('scale[1.0,0.996,1.005]')
        self.w.factory('error1[1,0.01,10.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        self.w.factory("expr::relativeError('error1*massErrRaw',error1,massErrRaw")
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])



    def buildJModelCBParam(self,name,isMC = False):
        self.w.factory("b[3e-6,0,14e-6]")
        self.poi.append('b')
        self.w.factory("expr::sigma('1e-12+sqrt(massErrRaw1*massErrRaw1/(massRaw*massRaw)+massErrRaw2*massErrRaw2/(massRaw*massRaw)+2*b)*massRaw',massErrRaw1,massErrRaw2,b,curvRaw1,curvRaw2,massRaw)")
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory('RooCBShape::'+name+postfix+'(massRaw,mass[3.091,3.08,3.12],sigma,alpha[3,0.5,20],n[4.77,0.1,100])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NSIG','NBKG','mass'])
        else:    
            self.nuisances.extend(['mass'])

        self.conditionals.extend(['massErrRaw1','massErrRaw2','curvRaw1','curvRaw2'])     


    def buildJModelCBSimple(self,name,isMC = False):
        self.w.factory("error1[0.01,0.0001,0.1]")
        self.poi.append('error1')
        if isMC:
            postfix=''
        else:    
            postfix='Sig'

        self.w.factory('RooCBShape::'+name+postfix+'(massRaw,mass[3.091,3.08,3.12],error1,alpha[3,0.5,20],n[4.77,0.1,100])')
        if not isMC:
            self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
            self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
            self.nuisances.extend(['bkgSlope','NSIG','NBKG','mass'])
        else:    
            self.nuisances.extend(['mass'])




    def buildZModelCBParam(self,name):
        self.w.factory("a[1e-6,1e-12,1e-4]")
        self.poi.append('a')
#        self.w.factory("expr::sigma('1e-12+sqrt(massErrRaw1*massErrRaw1/(massRaw2*massRaw2)+massErrRaw2*massErrRaw2/(massRaw2*massRaw2)+a/(curvRaw1*curvRaw1)+a/(curvRaw2*curvRaw2))*massRaw2',massErrRaw1,massErrRaw2,a,curvRaw1,curvRaw2,massRaw2)")
        self.w.factory("expr::sigma('1e-12+sqrt(massErrRaw1*massErrRaw1/(massRaw2*massRaw2)+massErrRaw2*massErrRaw2/(massRaw2*massRaw2)+2*a)*massRaw2',massErrRaw1,massErrRaw2,a,curvRaw1,curvRaw2,massRaw2)")
        self.w.factory('RooBreitWigner::lineshape(massRaw,90.86,2.4952)')
        self.w.factory('RooCBShape::resolution(massRaw,shift[0.0,-5.,5.],sigma,alpha[3],n[4.77,0.1,100])')
        convolved = ROOT.RooFFTConvPdf(name+'Sig','',self.w.var('massRaw'),self.w.pdf('lineshape'),self.w.pdf('resolution'),2)
        convolved.setBufferFraction(0.2)
        getattr(self.w,'import')(convolved)

        self.w.factory('RooExponential::'+name+'Bkg(massRaw,bkgSlope[-1,-8.,5])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,1000000]*'+name+'Bkg)')
        self.nuisances.extend(['bkgSlope','NSIG','NBKG','shift','n'])
        self.conditionals.extend(['massErrRaw1','massErrRaw2','curvRaw1','curvRaw2','massRaw2'])     



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
            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2),ROOT.RooFit.SumW2Error(1))
        else:
            if len(self.conditionals)>0:
                fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(0),ROOT.RooFit.Strategy(2),ROOT.RooFit.ConditionalObservables(self.w.set('conditionals')))
            else:
                fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2))




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
