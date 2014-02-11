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

    def buildMassFromVectors(self,name = 'mass'):
        self.w.factory("expr::px1('pt1*cos(phi1)',pt1,phi1)")
        self.w.factory("expr::py1('pt1*sin(phi1)',pt1,phi1)")
        self.w.factory("expr::theta1('2*atan(exp(-eta1))',eta1)")
        self.w.factory("expr::pz1('pt1/tan(theta1)',pt1,theta1)")
        self.w.factory("expr::E1('sqrt(pz1*pz1+pt1*pt1+muMass*muMass)',pt1,pz1,muMass)")

        self.w.factory("expr::px2('pt2*cos(phi2)',pt2,phi2)")
        self.w.factory("expr::py2('pt2*sin(phi2)',pt2,phi2)")
        self.w.factory("expr::theta2('2*atan(exp(-eta2))',eta2)")
        self.w.factory("expr::pz2('pt2/tan(theta2)',pt2,theta2)")
        self.w.factory("expr::E2('sqrt(pz2*pz2+pt2*pt2+muMass*muMass)',pt2,pz2,muMass)")

        self.w.factory("expr::crossProduct('px1*px2+py1*py2+pz1*pz2',px1,px2,py1,py2,pz1,pz2)")

        self.w.factory("expr::"+name+"('sqrt(2*muMass*muMass+2*(E1*E2-crossProduct))',muMass,E1,E2,crossProduct)")


    def buildZModel(self,name, var,dataset):
        self.w.factory('scale[1.0,0.5,1.5]')
        self.w.factory('error1[1,0.5,5.]')
        self.w.factory('error2[0.0]')
        self.poi.append('scale')
        self.w.var("scale").setError(0.5)
        
        pdf = ROOT.RooGaussianSumPdf(name+'Sig',name+'Sig',var,self.w.var('scale'),self.w.var('error1'),self.w.var('error2'),dataset,'massRaw')
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


        
    def extend(self,name,pdf,normName):
        self.w.factory(normName+'[0.,1000000]')
        self.w.factory('RooExtendPdf::'+name+'('+pdf+','+normName+')')
        


    def fitConditional(self,verbose = 1):
        nll = self.w.pdf('model').createNLL(self.w.data('data'),ROOT.RooFit.NumCPU(5,1),ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(self.w.var('curvRaw1'),self.w.var('curvRaw2'),self.w.var('etaRaw1'),self.w.var('etaRaw2'),self.w.var('phiRaw1'),self.w.var('phiRaw2'))),ROOT.RooFit.Verbose(verbose),ROOT.RooFit.Offset(1))

        minuit = ROOT.RooMinuit(nll)
        minuit.setVerbose(verbose)
        minuit.setStrategy(0)
        minuit.setEps(1e-8)
#        minuit.simplex()
        minuit.migrad()
        minuit.hesse()
        minuit.cleanup()

        
    def fit(self,model,data, hint = False, verbose =0,snapshot = "result"):
        minimum = ROOT.Double(0.0)
        maximum = ROOT.Double(0.0)
        if verbose==0 :
            ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

        hintResult=None
        fitResult=None
        
        #first extract a  guess
        if not hint:
            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Warnings(verbose),ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Strategy(2),ROOT.RooFit.Minos(1),ROOT.RooFit.Save(1))
            fitResult=self.w.pdf(model).fitTo(data,ROOT.RooFit.Warnings(verbose),ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Strategy(2),ROOT.RooFit.Minos(1),ROOT.RooFit.Save(1))

        if hint:
            for poi in self.poi:
                self.w.var(poi).setConstant(1)
            hintResult = self.w.pdf(model).fitTo(data,ROOT.RooFit.Warnings(0),ROOT.RooFit.Verbose(0),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(4,0),ROOT.RooFit.Minos(1),ROOT.RooFit.Save(1))
        #now constrain the nuisances within 5 sigma of this guess
            for nuis in self.nuisances:
                self.addGaussianConstraint(nuis)

        #now release the POIs
            for poi in self.poi:
                self.w.var(poi).setConstant(0)
                
            #now build the new model
            pdfs = ROOT.RooArgList()
            for pdf in self.constraints+[model]:
                pdfs.add(self.w.pdf(pdf))

            constrained = ROOT.RooProdPdf('constrained','constrained',pdfs)

            fitResult = constrained.fitTo(data,ROOT.RooFit.Warnings(verbose),ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(6,0),ROOT.RooFit.Minos(1),ROOT.RooFit.Save(1))
            fitResult = constrained.fitTo(data,ROOT.RooFit.Warnings(verbose),ROOT.RooFit.Verbose(verbose),ROOT.RooFit.PrintLevel(verbose),ROOT.RooFit.Timer(0),ROOT.RooFit.NumCPU(6,0),ROOT.RooFit.Minos(1),ROOT.RooFit.Save(1))

            self.w.saveSnapshot(snapshot,','.join(self.poi))    

        if hintResult is not None:
            print '-------HINT------'
            hintResult.Print()
            
        print '-------FINAL RESULT------'
        if fitResult is not None:
            fitResult.Print()

        return (hintResult,fitResult)        
