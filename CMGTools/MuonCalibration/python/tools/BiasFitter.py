import ROOT
from array import array

class BiasFitter (object):
    def __init__(self,workspace):
        self.posData=[]
        self.negData=[]
        self.w=workspace
        self.w.factory('bias[1e-6,-0.001,0.001]')

    def prepareModel(self):
        posData = self.w.data('pos')
        negData = self.w.data('neg')
        
        for i in range(0,posData.numEntries()):
            line = posData.get(i)
            self.posData.append(line.find('curvRaw1').getVal())

        for i in range(0,negData.numEntries()):
            line = negData.get(i)
            self.negData.append(line.find('curvRaw2').getVal())

            
        self.pdf = ROOT.RooPosNegBiasEstimator('model','biasE',self.w.var('bias'),len(self.posData),array('d',self.posData),len(self.negData),array('d',self.negData))
        

    def importData(self,dataPos,dataNeg):
        getattr(self.w,'import')(dataPos)
        getattr(self.w,'import')(dataNeg)
        self.w.data(dataPos.GetName()).SetName("pos")
        self.w.data(dataNeg.GetName()).SetName("neg")

    def optimize(self):
        print 'optimizing range'
        offset=5e-6

        mini = self.w.var("bias").getMin()
        maxi = self.w.var("bias").getMax()
        bins =int((maxi-mini)/offset)
        
        inputs=[]

        for i in range(0,bins):
            x = mini+i*offset
            self.w.var("bias").setVal(x)
            prob =self.pdf.getProbability() 
            if prob!=0:
                inputs.append(x)

        mini =  min(inputs)       
        maxi =  max(inputs)       
        self.w.var("bias").setVal(mini+(maxi-mini)/2.)        
        self.w.var("bias").setMin(mini)        
        self.w.var("bias").setMax(maxi)
                
        

        
    def fit(self,verbose=1):
        minuit = ROOT.RooMinuit(self.pdf)
        minuit.setVerbose(verbose)
        minuit.setStrategy(2)
        minuit.setEps(1e-7)
#        minuit.seek()
        self.optimize()
        minuit.migrad()
        minuit.improve()
        minuit.hesse()
        minuit.minos()
        result = minuit.save()
        minuit.cleanup()
        return result
        

            
            
