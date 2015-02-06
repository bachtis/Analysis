from defs import *
import math

#load the ratio of signal and background and the slopes


def estimate(minMass,maxMass,bkg=False):

    filename="kalmanTargetJpsi_Data.root"
    if  bkg:
        filename="kalmanTargetJpsi_Data_Bkg.root"
        
    pmap = PartitionMap(curvArr,etaArr,phiArr,filename)

    pmap.load('JDataFits.root','slope_fit','slope')
    pmap.load('JDataFits.root','BKGoverSIGNAL_fit','ratio')
    pmap.declareData('mass',0.0)


    builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',50000000)
    builder.load("JGENDataSmear.root")

    w.var('massRaw').setBins(50)
    w.var('massRaw').setMin(2.9)
    w.var('massRaw').setMax(3.3)

    w.factory('RooExponential::pdf(massRaw,slope[-1,-8.,10])')

    for bin,data in builder.data('pos').iteritems():
        NSIG = data.numEntries()
        ratio = pmap.getData('ratio',bin)
        
        #for this bin generate a background sample
        slope = pmap.getData('slope',bin)
        w.var('slope').setVal(slope)
        dataset = w.pdf('pdf').generate(ROOT.RooArgSet(w.var('massRaw')),NSIG*ratio)
        dataset.append(data)
        dataset=dataset.reduce('massRaw>'+str(minMass)+'&&massRaw<'+str(maxMass))
        if not bkg:
            dataset = data.reduce('massRaw>'+str(minMass)+'&&massRaw<'+str(maxMass))
        pmap.setData('mass',bin,dataset.mean(w.var('massRaw')),dataset.sigma(w.var('massRaw'))/math.sqrt(dataset.numEntries()))
          
    pmap.save('fit')      


estimate(3.0,3.2,False)          


