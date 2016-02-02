
import ROOT
from KaMuCa.Derivation.tools.LineshapeFitter import LineshapeFitter
from KaMuCa.Derivation.tools.DataSetSplitter import DataSetSplitter
from KaMuCa.Derivation.tools.Scheduler import Scheduler
from KaMuCa.Derivation.tools.workspaceTools import prepareWorkspace
import math



def fitJSimple(pmap,outputFile,dataIn,genIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('sigma',0.00)
    pmap.declareData('sigma2',0.00)
    pmap.declareData('ebe',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)

    genFile = ROOT.TFile(genIn)
    genData = genFile.Get('data')
    lineshape=genData.reduce('mass>2.8&&mass<3.4').reduce(ROOT.RooFit.EventRange(0,10000))

    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        data.get().find("mass").setMin(2.9)
        data.get().find("mass").setMax(3.3)

        data.get().find("c1").setMax(1.0)
        data.get().find("c1").setMin(0)

        data.get().find("c2").setMax(1.0)
        data.get().find("c2").setMin(0)


        w.var('mass').setMin(data.get().find("mass").getMin())
        w.var('mass').setMax(data.get().find("mass").getMax())

        w.var('c1').setMin(data.get().find("c1").getMin())
        w.var('c1').setMax(data.get().find("c1").getMax())

        w.var('c2').setMin(data.get().find("c2").getMin())
        w.var('c2').setMax(data.get().find("c2").getMax())

        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>2.9&&mass<3.3')
        dataBinned=fitter.convertToBinned(dataU,'mass',50)
        fitter.buildJModelSimple('model',w.var('mass'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        ebeMean=data.mean(data.get().find("massErr"))
        ebeRMS=data.sigma(data.get().find("massErr"))

        pmap.setData('ebe',bin,ebeMean,ebeRMS/math.sqrt(data.numEntries()))


        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('sigma',bin,value,error)
        pmap.setData('sigma2',bin,value*value,2*value*error)
    pmap.save('0')


def fitJParam(pmap,outputFile,dataIn,grid=False,isMC=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('mean',0.00)
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(2.95)
        w.var('mass').setMax(3.25)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>2.95&&mass<3.25')
        dataU.get().find("mass").setMin(2.95)
        dataU.get().find("mass").setMax(3.25)
        fitter.buildJModelCBParam('model',isMC)
        fitter.importData(dataU)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'mean')
        pmap.setData('mean',bin,value,error)
        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')



def fitZParam(pmap,outputFile,dataIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('mean',0.00)
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(80)
        w.var('mass').setMax(110)
#        w.var('mass').setBins(1000,'cache')
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>80&&mass<110')
        dataU.get().find("mass").setMin(80)
        dataU.get().find("mass").setMax(110)
#        dataU.get().find("mass").setBins(1000,'cache')
        fitter.buildZModelCBParam('model')
#        w.var('mass').setMin(80)
#        w.var('mass').setMax(110)
#        dataU.get().find("mass").setMin(80)
#        dataU.get().find("mass").setMax(110)
        fitter.importData(dataU)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'mean')
        pmap.setData('mean',bin,value,error)
        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')


def simFitParam(pmap,outputFile,dataIn1,dataIn2,MC=False,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn1)

    dataSplitter2 = DataSetSplitter(pmap)
    dataSplitter2.load(dataIn2)


    sched = Scheduler([outputFile.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]
        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue


        w=ROOT.RooWorkspace('w','w')
        w.factory('c1[0,1]')
        w.factory('c2[0,1]')
        w.factory('N1[0,200]')
        w.factory('N2[0,200]')
        w.factory('eta1[-2.4,2.4]')
        w.factory('eta2[-2.4,2.4]')
        w.factory('massErrZ[0.01,0,3.]')
        w.factory('massErrJ[0.01,0,3.]')
        w.factory('cErr1[0.01,0.000001,3.]')
        w.factory('cErr2[0.01,0.000001,3.]')
        w.factory('weight[1,0.0,1000000.]')
        w.factory('phi1[-'+str(math.pi)+','+str(math.pi)+ ']')
        w.factory('phi2[-'+str(math.pi)+','+str(math.pi)+ ']')
        w.factory('muMass[0.1056583715]')
        w.factory('massZ[80,110]')
        w.factory('massJ[2.95,3.25]')

        dataZ = data2.reduce('mass>80&&mass<110')
        dataZ.get().find("mass").setMin(80)
        dataZ.get().find("mass").setMax(110)

        dataJ = data.reduce('mass>2.95&&mass<3.25')
        dataJ.get().find("mass").setMin(2.95)
        dataJ.get().find("mass").setMax(3.25)
#        w.var('mass').setBins(1000,'cache')
        fitter = LineshapeFitter(w)
        fitter.addObservable('massZ')
        fitter.addObservable('massJ')

        fitter.buildSimModelParam(dataZ,dataJ,MC)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('2nd')
        return
#        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]
        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')


def simFitParamRelative(pmap,outputFile,dataIn1,dataIn2,MC=False,grid=False,inputFile=''):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)

    f2=ROOT.TFile(inputFile)
    p_a = f2.Get('a_0')
    p_b = f2.Get('b_0')
    p_c = f2.Get('c_0')
    p_d = f2.Get('d_0')


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn1)

    dataSplitter2 = DataSetSplitter(pmap)
    dataSplitter2.load(dataIn2)


    sched = Scheduler([outputFile.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]
        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue


        w=ROOT.RooWorkspace('w','w')
        w.factory('c1[0,1]')
        w.factory('c2[0,1]')
        w.factory('N1[0,200]')
        w.factory('N2[0,200]')
        w.factory('eta1[-2.4,2.4]')
        w.factory('eta2[-2.4,2.4]')
        w.factory('massErrZ[0.01,0,3.]')
        w.factory('massErrJ[0.01,0,3.]')
        w.factory('cErr1[0.01,0.000001,3.]')
        w.factory('cErr2[0.01,0.000001,3.]')
        w.factory('weight[1,0.0,1000000.]')
        w.factory('phi1[-'+str(math.pi)+','+str(math.pi)+ ']')
        w.factory('phi2[-'+str(math.pi)+','+str(math.pi)+ ']')
        w.factory('muMass[0.1056583715]')
        w.factory('massZ[80,110]')
        w.factory('massJ[2.95,3.25]')

        dataZ = data2.reduce('mass>80&&mass<110')
        dataZ.get().find("mass").setMin(80)
        dataZ.get().find("mass").setMax(110)

        dataJ = data.reduce('mass>2.95&&mass<3.25')
        dataJ.get().find("mass").setMin(2.95)
        dataJ.get().find("mass").setMax(3.25)
#        w.var('mass').setBins(1000,'cache')
        fitter = LineshapeFitter(w)
        fitter.addObservable('massZ')
        fitter.addObservable('massJ')

        bx,by,bz=pmap.binXYZ(bin)
        a = p_a.GetBinContent(by)
        b = p_b.GetBinContent(by)
        c = p_c.GetBinContent(by)
        d = p_d.GetBinContent(by)
        
        fitter.w.factory('a[{val},{mini},{maxi}]'.format(val=a,mini=0.5*a,maxi=2*a))
        fitter.w.factory('b[{val},{mini},{maxi}]'.format(val=b,mini=0.5*b,maxi=2*b))
        fitter.w.factory('c[{val},{mini},{maxi}]'.format(val=c,mini=0.5*c,maxi=2*c))
        fitter.w.factory('d[{val},{mini},{maxi}]'.format(val=d,mini=0.5*d,maxi=2*d))

        print 'BUIULDING MODEL'
        fitter.buildSimModelParam(dataZ,dataJ,MC)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('2nd')
        return
#        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]
        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')


def fitJParamEbE(pmap,outputFile,dataIn,grid=False,isMC=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('mean',0.00)
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(2.95)
        w.var('mass').setMax(3.25)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>2.95&&mass<3.25')
        dataU.get().find("mass").setMin(2.95)
        dataU.get().find("mass").setMax(3.25)
        fitter.buildJModelCBParamEbE('model',isMC)
        fitter.importData(dataU)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'mean')
        pmap.setData('mean',bin,value,error)
        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')



def fitPull(pmap,outputFile,dataIn,genIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('sigma',0.00)
    pmap.declareData('sigma2',0.00)
    pmap.declareData('muMinus1',0.00)

    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)

    genFile = ROOT.TFile(genIn)
    genData = genFile.Get('data_13')
    genData.get().find("mass").setMin(-8)
    genData.get().find("mass").setMax(8)
    lineshape=genData.reduce('mass>-8&&mass<8').reduce(ROOT.RooFit.EventRange(0,10000))

    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(-5)
        w.var('mass').setMax(5)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>-5&&mass<5')
        dataU.get().find("mass").setMin(-5)
        dataU.get().find("mass").setMax(5)
        dataU.get().find("mass").setBins(50)
        dataBinned=fitter.convertToBinned(dataU,'mass',50)
        fitter.buildPullModel('model',w.var('mass'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        print 'calculating error'
        avg=0.0
        N=data.numEntries()
        for i in range(0,N):
            line=data.get(i)
            cErr1 = line.find("cErr1").getVal()
            cErr2 = line.find("cErr2").getVal()
            avg=avg+ 0.25*(cErr1*cErr1+cErr2*cErr2)
        massErr2=avg/float(N)    

        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('sigma',bin,value,error)
        pmap.setData('sigma2',bin,value*value,2*value*error)

        pmap.setData('muMinus1',bin,massErr2*(value*value-1),2*massErr2*value*error)
    pmap.save('0')


def fitJ(pmap,outputFile,dataIn,genIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('b',0.00)



    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    genSplitter = DataSetSplitter(pmap)
    genSplitter.load(genIn)



    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(2.9)
        w.var('mass').setMax(3.3)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>2.9&&mass<3.3')
        dataU.get().find("mass").setMin(2.9)
        dataU.get().find("mass").setMax(3.3)
        dataU.get().find("mass").setBins(50)
        lineshape = genSplitter.data[bin].reduce("mas>2.8&&mass<3.4&&massErr/mass>0.002").reduce(ROOT.RooFit.EventRange(0,4000))
        dataBinned=fitter.convertToBinned(dataU,'mass',50)
        fitter.buildJModel('model',w.var('mass'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)

    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('b',bin,value,error)

    pmap.save('0')


def fitZSimple(pmap,outputFile,dataIn,genIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('sigma',0.00)
    pmap.declareData('sigma2',0.00)
    pmap.declareData('ebe',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)

    genSplitter = DataSetSplitter(pmap)
    genSplitter.load(genIn)

    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(85)
        w.var('mass').setMax(100.)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>85.0&&mass<100.0')
        dataU.get().find("mass").setMin(85.0)
        dataU.get().find("mass").setMax(100.)
        dataU.get().find("mass").setBins(50)
        dataBinned=fitter.convertToBinned(dataU,'mass',50)
        lineshape = genSplitter.data[bin].reduce('mass>78&&mass<112').reduce(ROOT.RooFit.EventRange(0,10000))
        fitter.buildZModelSimple('model',w.var('mass'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)
#    sched.submitLOCAL()
    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue


        ebeMean=data.mean(data.get().find("massErr"))
        ebeRMS=data.sigma(data.get().find("massErr"))
        pmap.setData('ebe',bin,ebeMean,ebeRMS/math.sqrt(data.numEntries()))
        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('sigma',bin,value,error)
        pmap.setData('sigma2',bin,value*value,2*value*error)
    pmap.save('0')


def fitZ(pmap,outputFile,dataIn,genIn,grid=False):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('b',0.00)



    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)

    genSplitter = DataSetSplitter(pmap)
    genSplitter.load(genIn)

    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        w=ROOT.RooWorkspace('w','w')
        prepareWorkspace(w)
        w.var('mass').setMin(80)
        w.var('mass').setMax(110.)
        w.var('mass').setBins(50)
        fitter = LineshapeFitter(w)
        fitter.addObservable('mass')
        dataU = data.reduce('mass>80.0&&mass<110.0')
        dataU.get().find("mass").setMin(80.0)
        dataU.get().find("mass").setMax(110.)
        dataU.get().find("mass").setBins(50)
        dataBinned=fitter.convertToBinned(dataU,'mass',50)
        lineshape = genSplitter.data[bin].reduce('mass>78&&mass<112&&massErr/mass>0.002').reduce(ROOT.RooFit.EventRange(0,10000))
        fitter.buildZModel('model',w.var('mass'),lineshape)
        fitter.importData(dataBinned)
        fitter.save()
        w.writeToFile("fit_"+str(bin)+".root")
        sched.declareJob("fit_"+str(bin)+".root",bin)
#    sched.submitLOCAL()
    if grid:    
        sched.submit('1nh')
        sched.wait()
    else:    
        sched.submitLOCAL()
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('b',bin,value,error)

    pmap.save('0')




def harvest(pmap,outputFile,dataIn):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('scale',0.00)
    pmap.declareData('b',0.00)


    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        sched.declareFile("fit_"+str(bin)+".root",bin)

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'scale')
        pmap.setData('scale',bin,value,error)
        value,error = sched.harvest(bin,'error1')
        pmap.setData('b',bin,value,error)

    pmap.save('0')


def harvestParam(pmap,outputFile,dataIn):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('mean',0.00)
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)

    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn)


    sched = Scheduler([dataIn.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue
        sched.declareFile("fit_"+str(bin)+".root",bin)

    for bin,data in dataSplitter.data.iteritems():
        if data.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'mean')
        pmap.setData('mean',bin,value,error)
        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')


def harvestSim(pmap,outputFile,dataIn1,dataIn2):
    pmap.addLog(outputFile)
    pmap.log.cd()
    pmap.declareData('a',0.00)
    pmap.declareData('b',0.00)
    pmap.declareData('c',0.00)
    pmap.declareData('d',0.00)

    dataSplitter = DataSetSplitter(pmap)
    dataSplitter.load(dataIn1)

    dataSplitter2 = DataSetSplitter(pmap)
    dataSplitter2.load(dataIn2)


    sched = Scheduler([outputFile.split('.')[0]],1)
    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]
        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue

        sched.declareFile("fit_"+str(bin)+".root",bin)

    for bin,data in dataSplitter.data.iteritems():
        data2=dataSplitter2.data[bin]

        if data.numEntries()<50:
            continue
        if data2.numEntries()<50:
            continue

        value,error = sched.harvest(bin,'a')
        pmap.setData('a',bin,value,error)
        value,error = sched.harvest(bin,'b')
        pmap.setData('b',bin,value,error)
        value,error = sched.harvest(bin,'c')
        pmap.setData('c',bin,value,error)
        value,error = sched.harvest(bin,'d')
        pmap.setData('d',bin,value,error)

    pmap.save('0')
    
