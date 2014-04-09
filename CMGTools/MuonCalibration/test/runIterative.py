import ROOT
import math
ROOT.gSystem.Load("libCMGToolsMuonCalibration")


from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler
from CMGTools.MuonCalibration.tools.MegaFitter import MegaFitter


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(70)
w.var('massRaw').setMax(120)




curvArr= [0.005,0.050]

phiArr=[]
for i in range(0,13):
    phiArr.append(-math.pi+2*math.pi*i/12.)

etaArr = [-2.4,-2.3,-1.9,-1.5,-1.1,-0.9,0.0,0.9,1.1,1.5,1.9,2.3,2.4]

#etaArr = [-1.1,-0.9,0.0,0.9,1.1]
etaArr = [-2.4,-1.1,1.1,2.4]


pmap = PartitionMap(curvArr,etaArr,phiArr,"results_IR_1_broadEta.root")

#pmap = PartitionMap([0.005,0.05],[-1.3,1.3],[-math.pi,math.pi],"results_IR.root")


pmap.declareData('CposScale',1.00)
pmap.declareData('CnegScale',1.00)
pmap.declareData('CposRes',0.00)
pmap.declareData('CnegRes',0.00)
pmap.declareData('Cpos',0.00)
pmap.declareData('Cneg',0.00)
pmap.declareData('Rpos',1.00)
pmap.declareData('Rneg',1.00)


builder = DataSetBuilder(pmap,w,'ZMC.root','data',500000000,True)
#builder = DataSetBuilder(pmap,w,'ZMC.root','data',200000,False,False)
builderMC = DataSetBuilder(pmap,w,'ZGEN.root','data',500000000)
#builderMC = DataSetBuilder(pmap,w,'ZGEN.root','data',10000)


builder.build()
builderMC.build(2000)



for N in range(0,10):
#    for sign in ['pos','neg']:
    dataC={'pos':{},'neg':{}}
    for sign in ['pos','neg']:
        for bin,data in builder.data(sign).iteritems():
            dataC[sign][bin]= pmap.recalibrate(data,w)

    sched={}
    for sign in ['pos','neg']:
#    for sign in ['pos']:
        sched[sign] = Scheduler(['IR',sign,str(N)],2)
       ###create Scheduler here
        for bin,data in dataC[sign].iteritems():
            w=ROOT.RooWorkspace('w','w')
            prepareWorkspace(w)
            w.var('massRaw').setBins(100)
            w.var('massRaw').setMin(70)
            w.var('massRaw').setMax(120)
            fitter = LineshapeFitter(w)
            fitter.addObservable('massRaw')
            dataBinned=builder.convertToBinned2D(data,['massRaw','massErrRaw'],100,30)
            prelineShape = builderMC.data(sign)[bin]
            lineshape = prelineShape.reduce("massRaw>65&&massRaw<125")
            fitter.buildZModelEbE('model',w.var('massRaw'),lineshape)
            fitter.importData(dataBinned)
            fitter.save()
            w.writeToFile(sign+"_"+str(bin)+".root")
            sched[sign].declareJob(sign+"_"+str(bin)+".root",bin)
 
#        sched[sign].submitLOCAL()
        sched[sign].submit('8nh')
    sched['pos'].wait()
    for sign in ['pos','neg']:
        for bin,data in builder.data(sign).iteritems():
            if sign == "pos":
                othersign = "neg"
                othercurv = "curvRaw2"
                curv = "curvRaw1"
                massErr = 'massErrRaw1'
                otherMassErr = 'massErrRaw2'
            else:
                othersign = "pos"
                othercurv = "curvRaw1"
                curv = 'curvRaw2'
X                massErr = 'massErrRaw2'
                otherMassErr = 'massErrRaw1'
                
            value,error = sched[sign].harvest(bin)


            pmap.setData('C'+sign+'Scale',bin,value,error )

            average = data.mean(data.get().find(curv))
            otheraverage = data.mean(data.get().find(othercurv))

            weightOfBin =  (0.25)*average

            bias =weightOfBin*((value*value)-1.0)
            biasError = 2*abs(error*weightOfBin*value)
            pmap.addData('C'+sign,bin,bias,biasError)


            sigma,sigmaErr = sched[sign].harvest(bin,'errorScale')
            pmap.setData('C'+sign+'Res',bin,sigma,sigmaErr)

#            averageMassErr = data.mean(data.get().find('massErrRaw'))

 #           sigmaFactor = sigma/averageMassErr
 #           sigmaFactorErr = sigmaErr/averageMassErr
            
            avMassErr1 = data.mean(data.get().find(massErr))
            avMassErr2 = data.mean(data.get().find(otherMassErr))
            

            errorScale = 1+(sigma-1)*avMassErr1/(avMassErr1+avMassErr2)
            errorScaleError = sigmaErr*avMassErr1/((avMassErr1+avMassErr2))
            pmap.multiplyData('R'+sign,bin,errorScale,errorScaleError )

    
            

    pmap.save(str(N))    
