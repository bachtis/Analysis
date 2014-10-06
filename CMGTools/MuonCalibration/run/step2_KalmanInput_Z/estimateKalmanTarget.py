from defs import *


#load the ratio of signal and background and the slopes
pmap = PartitionMap(curvArr,etaArr,phiArr,"kalmanTargetZ.root")
pmap.declareData('mass',0.0)
pmap.declareData('width',0.0)

builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',50000000)
builder.load("ZGEN_Input.root")

w.var('massRaw').setBins(50)
w.var('massRaw').setMin(70)
w.var('massRaw').setMax(120)

builderS = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',50000000)
builderS.load("ZGEN_InputSmeared.root")





def estimate(minMass,maxMass):
    for bin,data in builderS.data('pos').iteritems():
        dataset=data.reduce('massRaw>'+str(minMass)+'&&massRaw<'+str(maxMass))
        pmap.setData('mass',bin,dataset.mean(w.var('massRaw')),0.0)
    for bin,data in builder.data('pos').iteritems():
        dataset=data.reduce('massRaw>'+str(minMass)+'&&massRaw<'+str(maxMass))
        pmap.setData('width',bin,dataset.sigma(w.var('massRaw')),0.0)
          
    pmap.save('fit')      
          


