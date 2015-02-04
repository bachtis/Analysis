from defs import *


pmap = PartitionMap(curvArr,etaArr,phiArr,"resolutionHighPtMC.root")
pmap.declareData('a2',0.00)

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',1000000000)
builder.load("ZData_Input.root")
isMC=False


sched = Scheduler(['PARAMETRICZDATA'],1)
for bin,data in builder.data('pos').iteritems():
    w=ROOT.RooWorkspace('w','w')
    prepareWorkspace(w)
    w.var('massRaw').setMin(80)
    w.var('massRaw').setMax(120)
    w.var('massRaw').setBins(1000,'cache')
    w.var('curvRaw1').setVal(1./10.)
    w.var('curvRaw2').setVal(1./10.)
    w.var('curvRaw1').setMin(1./100.)
    w.var('curvRaw1').setMax(1./5.)
    w.var('curvRaw2').setMin(1./100.)
    w.var('curvRaw2').setMax(1./5.)
    w.var('massErrRaw1').setVal(0.03)
    w.var('massErrRaw2').setVal(0.03)

    w.factory('massRaw2[80,120]')
    fitter = LineshapeFitter(w)
    fitter.addObservable('massRaw')
    dataU = data.reduce('massRaw>80&&massRaw<120')
    dataU.get().find("massRaw").setBins(1000,'cache')

    fitter.buildZModelCBParam('model')

    #aDD A COLUMN WITH THE COPY OF THE MASS
    m2=ROOT.RooFormulaVar("massRaw2","massRaw2",'massRaw',ROOT.RooArgList(w.var('massRaw')))
    dataU.addColumn(m2)

    fitter.importData(dataU)
    fitter.save()
    w.writeToFile("fit_"+str(bin)+".root")
    sched.declareJob("fit_"+str(bin)+".root",bin)
#sched.submitLOCAL()
sched.submit('8nh')
sched.wait()

for bin,data in builder.data('pos').iteritems():
    value,error = sched.harvest(bin,'a')
    pmap.setData('a2',bin,value,error)


pmap.save('0')

    
