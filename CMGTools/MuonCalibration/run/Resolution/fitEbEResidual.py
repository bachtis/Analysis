from defs import *


pmap = PartitionMap(curvArr,etaArr,phiArr,"paramResult.root")
pmap.declareData('a2',0.00)
pmap.declareData('b2',0.00)
pmap.declareData('a',0.00)
pmap.declareData('b',0.00)



builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',1000000000)
builder.load("Data_Input.root")
isMC=False


sched = Scheduler(['PARAMETRIC'],1)
for bin,data in builder.data('pos').iteritems():
    w=ROOT.RooWorkspace('w','w')
    prepareWorkspace(w)
    w.var('massRaw').setMin(2.9)
    w.var('massRaw').setMax(3.3)
    w.var('curvRaw1').setVal(1./10.)
    w.var('curvRaw2').setVal(1./10.)
    w.var('curvRaw1').setMin(1./100.)
    w.var('curvRaw1').setMax(1./5.)
    w.var('curvRaw2').setMin(1./100.)
    w.var('curvRaw2').setMax(1./5.)

    fitter = LineshapeFitter(w)
    fitter.addObservable('massRaw')
    dataU = data.reduce('massRaw>2.9&&massRaw<3.3')
    fitter.buildJModelCBParam('model',isMC)
    fitter.importData(dataU)
    fitter.save()
    w.writeToFile("fit_"+str(bin)+".root")
    sched.declareJob("fit_"+str(bin)+".root",bin)
sched.submitLOCAL()

for bin,data in builder.data('pos').iteritems():
#    value,error = sched.harvest(bin,'a')
#    pmap.setData('a2',bin,value,error )
#    value=value*math.sqrt(abs(value))/abs(value)
#    error=0.5*error/(2*math.sqrt(abs(value)))

#    pmap.setData('a',bin,value,error )

    value,error = sched.harvest(bin,'b')
    pmap.setData('b2',bin,value,error)
    value=value*math.sqrt(abs(value))/abs(value)
    error=0.5*error/(math.sqrt(abs(value)))
    pmap.setData('b',bin,value,error )

pmap.save('0')

    
