from defs import *

pmapJ = PartitionMap(curvArr,etaArr,phiArr,"JpsiResolutionMC2D.root")



pmapJ.declareData('sigmaCorr',0.00)
pmapJ.declareData('sigmaAbs',0.00)
pmapJ.declareData('ebeAbs',0.00)

builder = DataSetBuilder(pmapJ,w,'../../data/JDATA.root','data',1000000000)
builder.load("JMC_Input.root")

builderG = DataSetBuilder(pmapJ,w,'../../data/JGEN.root','data',1000000000)
builderG.load("JGEN_Input.root")


sched = Scheduler(['JPSIMC'],1)
for bin,data in builder.data('pos').iteritems():
    w=ROOT.RooWorkspace('w','w')
    prepareWorkspace(w)
    w.var('massRaw').setMin(2.9)
    w.var('massRaw').setMax(3.2)
    w.var('massRaw').setBins(50)
    fitter = LineshapeFitter(w)
    fitter.addObservable('massRaw')

    dataU = data.reduce('massRaw>2.9&&massRaw<3.2')
    data2 = builder.data('neg')[bin].reduce("massRaw>2.9&&massRaw<3.2")
    dataU.append(data2)

    dataU.get().find("massRaw").setMin(2.9)
    dataU.get().find("massRaw").setMax(3.2)
    dataU.get().find("massRaw").setBins(50)
    dataBinned=builder.convertToBinned(dataU,'massRaw',50)
#    lineshape=builderG.data('pos')[bin].reduce('massRaw>2.8&&massRaw<3.3').reduce(ROOT.RooFit.EventRange(0,1000))
#    fitter.buildJModelSimple('model',w.var('massRaw'),lineshape)
    fitter.buildJModelCBSimple('model',False)
    fitter.importData(dataBinned)
    fitter.save()
    w.writeToFile("fit_"+str(bin)+".root")
    sched.declareJob("fit_"+str(bin)+".root",bin)
#sched.submit('2nd')
#sched.wait()
sched.submitLOCAL()

for bin,data in builder.data('pos').iteritems():
    if (data.numEntries())<=0:
        continue

#    value,error = sched.harvest(bin,'scale')
#    pmapJ.setData('scale',bin,value,error)

    massAvg,error = sched.harvest(bin,'mass')

#    massAvg = data.mean(data.get().find("massRaw"))


    value,error = sched.harvest(bin,'error1')
    errAvg = data.mean(data.get().find("massErrRaw"))
    errAvgPos = data.mean(data.get().find("massErrRaw1"))

    ratio = errAvgPos/(errAvg)
    
    
    
    pmapJ.setData('ebeAbs',bin,errAvgPos*errAvgPos/(massAvg*massAvg),0.0)

    #the reall error square is err+^2 = err^2*()
    pmapJ.setData('sigmaAbs',bin,(value*value*ratio*ratio)/(massAvg*massAvg) ,2*value*ratio*ratio*error/(massAvg*massAvg))
    pmapJ.setData('sigmaCorr',bin,(value*value*ratio*ratio -errAvgPos*errAvgPos)/(massAvg*massAvg) ,2*value*ratio*ratio*error/(massAvg*massAvg))
pmapJ.save('0')

    
