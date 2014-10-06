from defs import *


pmap = PartitionMap(curvArr,etaArr,phiArr,"JDataFits.root")
pmap.declareData('scale',0.00)
pmap.declareData('BKGoverSIGNAL',0.00)
pmap.declareData('slope',0.00)
pmap.declareData('error1',0.00)
pmap.declareData('error2',0.00)
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',1000000000)
builderMC = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',50000000)
builder.load("JDATA_Input.root")
#builder.load("samplesGENFull.root")
builderMC.load("JGEN_Input.root")

def harvest(sched):
    for bin,data in builder.data('pos').iteritems():
            value,error = sched.harvest(bin)
            pmap.setData('scale',bin,value,error )
            value,error = sched.harvest(bin,'error1')
            pmap.setData('error1',bin,value,error )
            value,error = sched.harvest(bin,'error2')
            pmap.setData('error2',bin,value,error )
            value,error = sched.harvest(bin,'bkgSlope')
            pmap.setData('slope',bin,value,error )
            NSIG,error = sched.harvest(bin,'NSIG')
            NBKG,error = sched.harvest(bin,'NBKG')
            pmap.setData('BKGoverSIGNAL',bin,NBKG/NSIG,0.0)
    pmap.save('fit')


def run():

    scheduler = Scheduler(['JMASSFits'],1)
       ###create Scheduler here
    for bin,data in builder.data('pos').iteritems():
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('massRaw').setBins(50)
        w.var('massRaw').setMin(2.9)
        w.var('massRaw').setMax(3.3)
        fitter = LineshapeFitter(w)
        fitter.addObservable('massRaw')
        dataU = data.reduce('massRaw>2.9&&massRaw<3.3')
        dataBinned=builder.convertToBinned(dataU,'massRaw',50)
        lineshape = builderMC.data('pos')[bin].reduce("massRaw>2.8&&massRaw<3.4")
        fitter.buildJModelAL('model',w.var('massRaw'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("dataFit_"+str(bin)+".root")
        scheduler.declareJob("dataFit_"+str(bin)+".root",bin)
 
#        scheduler.submit('8nh')
    scheduler.submitLOCAL()
#    sched['pos'].wait()
    harvest(scheduler)

run()
    
