from defs import *


pmap = PartitionMap(curvArr,etaArr,phiArr,"resolutionData.root")
pmap.declareData('b2',0.00)

builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',1000000000)
builder.load("JData_Input.root")
isMC=False


sched = Scheduler(['PARAMETRICJ'],1)
for bin,data in builder.data('pos').iteritems():
    w=ROOT.RooWorkspace('w','w')
    prepareWorkspace(w)
    w.var('massRaw').setMin(2.9)
    w.var('massRaw').setMax(3.2)
    w.var('curvRaw1').setVal(1./10.)
    w.var('curvRaw2').setVal(1./10.)
    w.var('curvRaw1').setMin(1./60.)
    w.var('curvRaw1').setMax(1./5.)
    w.var('curvRaw2').setMin(1./60.)
    w.var('curvRaw2').setMax(1./5.)
    w.var('massErrRaw1').setVal(0.03)
    w.var('massErrRaw2').setVal(0.03)


    data.get().find('massRaw').setMin(2.9)
    data.get().find('massRaw').setMax(3.2)
    data.get().find('curvRaw1').setMin(1./60.)
    data.get().find('curvRaw1').setMax(1./5.)
    data.get().find('curvRaw2').setMin(1./60.)
    data.get().find('curvRaw2').setMax(1./5.)



    fitter = LineshapeFitter(w)
    fitter.addObservable('massRaw')
    dataU = data.reduce('massRaw>2.9&&massRaw<3.2')

    fitter.buildJModelCBParam('model',isMC)
    fitter.importData(dataU)
    fitter.save()
    w.writeToFile("fit_"+str(bin)+".root")
    sched.declareJob("fit_"+str(bin)+".root",bin)
sched.submitLOCAL()
#sched.wait()

for bin,data in builder.data('pos').iteritems():
    value,error = sched.harvest(bin,'b')
    pmap.setData('b2',bin,value,error)
pmap.save('0')

    
