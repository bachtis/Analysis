
import ROOT
import time

class MegaFitter(object):
    def __init__(self,inputData):
        self.data = inputData
        self.w = ROOT.RooWorkspace("w")

        binsPosArr=[]
        binsNegArr=[]

        for N,inp in inputData.iteritems():
            binsPosArr.append(inp['bin1'])
            binsNegArr.append(inp['bin2'])
            

        self.binsPos = set(binsPosArr)
        self.binsNeg = set(binsNegArr)

        #create nuisances
        for bin in self.binsPos:
            self.w.factory('pos_'+str(bin)+'[0.0,-0.001,0.001]')
        for bin in self.binsNeg:
            self.w.factory('neg_'+str(bin)+'[0.0,-0.001,0.001]')
            
            
        self.pdfs=[]    
        self.obs=[]
        self.obsErr=[]

    def checkCoverage(self,partitionMap):
        print 'HARVESTING RESULTS-------------------------------'
        for i in range(1,partitionMap.bins_curv()+1):
            for j in range(1,partitionMap.bins_eta()+1):
                for k in range(1,partitionMap.bins_phi()+1):
                    bin = partitionMap.bin(i,j,k)
                    if not (bin in self.binsPos):
                        print 'Positive constrain for bin:',bin,i,j,k,'missing'
                    if not (bin in self.binsNeg):
                        print 'Negative constrain for bin:',bin,i,j,k,'missing'
        print '-------------------------------------------------'
        
        
        

    def buildModel(self,maxErr =0.2):
        for N,inp in self.data.iteritems():
            var1 = 'pos_'+str(inp['bin1'])+'/'+str(inp['curv1'])
            var2 = 'neg_'+str(inp['bin2'])+'/'+str(inp['curv2'])
            value = 1./(inp['value']*inp['value'])-1
            error = 2*abs(inp['error']/(inp['value']*inp['value']*inp['value']))
            if abs(error/value)>maxErr:
                continue
            expression = 'expr::Constraint_'+str(N)+"('"+var1+"+"+var2+"',pos_"+str(inp['bin1'])+",neg_"+str(inp['bin2'])+")"
            self.w.factory(expression)
            self.w.factory("Scale_"+str(N)+"["+str(value)+"]")
            self.w.factory("ScaleErr_"+str(N)+"["+str(error)+"]")
            self.obs.append("Scale_"+str(N))
            self.obsErr.append("ScaleErr_"+str(N))
            self.w.factory("RooGaussian::Constraint_"+str(N)+"_PDF(Scale_"+str(N)+",Constraint_"+str(N)+",ScaleErr_"+str(N)+")")
            self.pdfs.append('Constraint_'+str(N)+'_PDF')
            print 'Added Gaussian with mean',value ,' and error',error , 'with inputerr',abs(inp['error']/inp['value'])

    def buildModelProd(self):
        for N,inp in self.data.iteritems():
            var1 = 'pos_'+str(inp['bin1'])
            var2 = 'neg_'+str(inp['bin2'])
            value = 1./(inp['value']*inp['value'])
            error = inp['error']/(inp['value']*inp['value']*inp['value'])
            expression = 'expr::Constraint_'+str(N)+"('"+var1+'*'+var2+"',pos_"+str(inp['bin1'])+",neg_"+str(inp['bin2'])+")"
            self.w.factory(expression)
            self.w.factory("Scale_"+str(N)+"["+str(value)+"]")
            self.w.factory("ScaleErr_"+str(N)+"["+str(error)+"]")
            self.obs.append("Scale_"+str(N))
            self.obsErr.append("ScaleErr_"+str(N))
            self.w.factory("RooGaussian::Constraint_"+str(N)+"_PDF(Scale_"+str(N)+",Constraint_"+str(N)+",ScaleErr_"+str(N)+")")
            self.pdfs.append('Constraint_'+str(N)+'_PDF')
            


        


    def fit(self,errorScale = 1.0):
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
        list=ROOT.RooArgList("pdfList")
        for pdf in self.pdfs:
            list.add(self.w.pdf(pdf))
        product = ROOT.RooProdPdf('model','model',list)
        getattr(self.w,'import')(product)
        argset=ROOT.RooArgSet('set')
        for obs in self.obs:
            argset.add(self.w.var(obs))
        data = ROOT.RooDataSet('data','data',argset)
        data.add(argset)
        getattr(self.w,'import')(data)
        #first save a snapshot of the error
        errSet = ROOT.RooArgSet('errSet')
        for err in self.obsErr:
            errSet.add(self.w.var(err))


        self.w.saveSnapshot('NOMINAL',errSet)
        for err in self.obsErr:
            self.w.var(err).setVal(self.w.var(err).getVal()*errorScale)
            self.w.var(err).setConstant(1)


        self.w.pdf('model').fitTo(data)

        self.w.loadSnapshot('NOMINAL')
        



#        self.w.saveSnapshot('NOMINAL',errSet)
#        self.w.pdf('model').fitTo(data)
        
        

    def writeToMap(self,pmap):
        for bin in self.binsPos:
            value = self.w.var('pos_'+str(bin)).getVal()
            error = self.w.var('pos_'+str(bin)).getError()
            pmap.setData('Cpos',bin,value,error)
        for bin in self.binsNeg:
            value = self.w.var('neg_'+str(bin)).getVal()
            error = self.w.var('neg_'+str(bin)).getError()
            pmap.setData('Cneg',bin,value,error)
