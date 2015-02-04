from defs import *

pmapZ = PartitionMap(curvArr,etaArr,phiArr,"zResolution2D.root")


pmapZ.declareData('scale',0.00)
pmapZ.declareData('sigmaCorr',0.00)
pmapZ.declareData('sigmaAbs',0.00)
pmapZ.declareData('ebeAbs',0.00)


builder = DataSetBuilder(pmapZ,w,'../../data/JMC.root','data',1000000000)
builder.load("ZData_Input.root")

builderG = DataSetBuilder(pmapZ,w,'../../data/JGEN.root','data',1000000000)
builderG.load("ZGEN_Input.root")


sched = Scheduler(['Z'],1)
for bin,data in builder.data('pos').iteritems():
    w=ROOT.RooWorkspace('w','w')
    prepareWorkspace(w)
    w.var('massRaw').setMin(80)
    w.var('massRaw').setMax(120)
    w.var('massRaw').setBins(100)
    fitter = LineshapeFitter(w)
    fitter.addObservable('massRaw')
    dataU.append(builder.data('neg')[bin] 
    dataU = data.reduce('massRaw>80&&massRaw<110')
    dataU.get().find("massRaw").setMin(80)
    dataU.get().find("massRaw").setMax(120)
    dataU.get().find("massRaw").setBins(100)
    dataBinned=builder.convertToBinned(dataU,'massRaw',100)
    lineshape=builderG.data('pos')[bin].reduce('massRaw>75&&massRaw<125').reduce(ROOT.RooFit.EventRange(0,8500))
    lineshape.append(builderG.data('neg')[bin].reduce('massRaw>75&&massRaw<125').reduce(ROOT.RooFit.EventRange(0,8500)))                 
    fitter.buildZModelSimple('model',w.var('massRaw'),lineshape)
    fitter.importData(dataBinned)
    fitter.save()
    w.writeToFile("fit_"+str(bin)+".root")
    sched.declareJob("fit_"+str(bin)+".root",bin)
#sched.submit('2nd')
#sched.wait()
sched.submitLOCAL()

for bin,data in builder.data('pos').iteritems():
    value,error = sched.harvest(bin,'scale')
    pmapZ.setData('scale',bin,value,error)
    value,error = sched.harvest(bin,'error1')
    massAvg = data.mean(data.get().find("massRaw"))
    errAvg = data.mean(data.get().find("massErrRaw"))
    errAvgPos = data.mean(data.get().find("massErrRaw1"))

    ratio =errAvgPos/errAvg 

    if (data.numEntries())<=0:
        continue
    pmapZ.setData('sigmaAbs',bin,(value*value*ratio*ratio)/(massAvg*massAvg) ,2*value*ratio*ratio*error/(massAvg*massAvg))
    pmapZ.setData('ebeAbs',bin,(errAvgPos*errAvgPos)/(massAvg*massAvg) ,0.0)
    pmapZ.setData('sigmaCorr',bin,(value*value*ratio*ratio -errAvgPos*errAvgPos)/(massAvg*massAvg) ,2*value*ratio*ratio*error/(massAvg*massAvg))

pmapZ.save('0')

    
