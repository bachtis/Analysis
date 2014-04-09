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

        self.w.defineSet('constraints',','.join(self.constraints))
        self.w.defineSet('poi',','.join(self.poi))
        self.w.defineSet('observables',','.join(self.observables))
        self.w.defineSet('nuisances',','.join(self.nuisances))


    def load(self):
        for setName in ['constraints','poi','observables','nuisances']:
            set = self.w.set(setName)
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
        self.w.factory('scale[1.0,0.5,1.5]')
        self.w.factory('error1[1,0.1,5.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        self.w.var("scale").setError(0.5)
        
        pdf = ROOT.RooGaussianSumPdfWithSigma(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw','massErrRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')

        self.nuisances.extend(['bkgSlope','NSIG','NBKG'])


    def buildZDualModel(self,name,lineshapePlus,lineshapeMinus,dataPlus,dataMinus,averagePlus,averageMinus):
        self.w.factory('c[-0.0005,0.0005]')
        self.w.factory('u[-0.0005,0.0005]')

        self.w.factory("expr::plus('sqrt((1+c/"+str(averagePlus)+")*(1-u))',c,u)")
        self.w.factory("expr::minus('sqrt((1-c/"+str(averageMinus)+")*(1+u))',c,u)")

        self.w.factory('errorPlus[1,0.1,5.]')
        self.w.factory('errorMinus[1,0.1,5.]')
        self.w.factory('error2[0.0]')



        self.poi.append('c')

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



                
        pdf1 = ROOT.RooGaussianSumPdfWithSigma(name+'Sig1',name+'Sig1',self.w.var("massRaw"),self.w.function('plus'),self.w.var('errorPlus'),self.w.var('error2'),lineshapePlus,'massRaw','massErrRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdfWithSigma.Class(),1)
        getattr(self.w,'import')(pdf1)

        self.w.factory('RooExponential::'+name+'Bkg1(massRaw,bkgSlope1[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Plus(NSIG1[0,100000000]*'+name+'Sig1,NBKG1[1,0,100000]*'+name+'Bkg1)')

        self.nuisances.extend(['bkgSlope1','NSIG1','NBKG1'])


        
                
        pdf2 = ROOT.RooGaussianSumPdfWithSigma(name+'Sig2',name+'Sig2',self.w.var("massRaw"),self.w.function('minus'),self.w.var('errorMinus'),self.w.var('error2'),lineshapeMinus,'massRaw','massErrRaw')
        getattr(self.w,'import')(pdf2)

        self.w.factory('RooExponential::'+name+'Bkg2(massRaw,bkgSlope2[-1,-8.,0])')
        self.w.factory('SUM::'+name+'Minus(NSIG2[0,100000000]*'+name+'Sig2,NBKG2[1,0,100000]*'+name+'Bkg2)')

        self.nuisances.extend(['bkgSlope2','NSIG2','NBKG2'])

        self.w.factory('sign[PLUS,MINUS]')
        self.w.factory('SIMUL::'+name+'(sign,PLUS='+name+'Plus,MINUS='+name+'Minus)')

        data = ROOT.RooDataSet('data','data',ROOT.RooArgSet(self.w.var('massRaw'),self.w.var('weight')),ROOT.RooFit.Index(self.w.cat('sign')),ROOT.RooFit.Import('PLUS',self.w.data('dataPlus')),ROOT.RooFit.Import('MINUS',self.w.data('dataMinus')),ROOT.RooFit.WeightVar('weight'))
        getattr(self.w,'import')(data)

        

    def buildZModelEbE(self,name, var,dataset):
        self.w.factory('scale[1.0,0.5,1.5]')
        self.w.factory('errorScale[1,0.1,3.]')
        self.w.factory("expr::error1('errorScale*massErrRaw',errorScale,massErrRaw)")
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        self.w.var("scale").setError(0.5)

        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.function('error1'),self.w.var('error2'),dataset,'massRaw')
        getattr(self.w,'importClassCode')(ROOT.RooGaussianSumPdf.Class(),1)
        getattr(self.w,'import')(pdf)
        self.w.factory('RooExponential::'+name+'Bkg('+var.GetName()+',bkgSlope[-1,-8.,0])')
        self.w.factory('SUM::'+name+'(NSIG[0,100000000]*'+name+'Sig,NBKG[1,0,100000]*'+name+'Bkg)')

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
        
        self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Strategy(2),ROOT.RooFit.Minos(1),ROOT.RooFit.SumW2Error(1))
        fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Save(1),ROOT.RooFit.SumW2Error(1))
            
        self.w.saveSnapshot(snapshot,','.join(self.poi))    

        if hintResult is not None:
            print '-------HINT------'
            hintResult.Print()
            
        print '-------FINAL RESULT------'
        if fitResult is not None:
            fitResult.Print()

        return (hintResult,fitResult)        



    def fitConditional(self,model,data, hint = False, verbose =0,snapshot = "result"):
        if verbose==0 :
            ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

        hintResult=None
        fitResult=None
        

        fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(1),ROOT.RooFit.NumCPU(8,0),ROOT.RooFit.Offset(1),ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(self.w.var('massErrRaw'),self.w.var("curvRaw1"),self.w.var("curvRaw2"),self.w.var("phiRaw1"),self.w.var("phiRaw2"))))
            
