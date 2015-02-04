from defs import *


pmap = PartitionMap(curvArr,etaArr,phiArr,"resolutionHighPtData.root")
pmap.declareData('a2',0.00)

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',1000000000)
builder.load("ZData_Input.root")



sched = Scheduler(['PARAMETRICZDATA'],1)
for bin,data in builder.data('pos').iteritems():
    sched.declareFile("fit_"+str(bin)+".root",bin)

for bin,data in builder.data('pos').iteritems():
    value,error = sched.harvest(bin,'a')
    pmap.setData('a2',bin,value,error)


pmap.save('0')

    
