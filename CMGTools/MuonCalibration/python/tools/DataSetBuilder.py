import ROOT
import os
import array


class DataSetBuilder (object):
    def __init__(self,pmap,workspace,file,dataset,entries,unfold = False,smear = False):
        self.map=pmap
        self.w=workspace
        f = ROOT.TFile(file)
        datas=f.Get(dataset)
        self.cache = ROOT.TFile('__cacheE__.root','RECREATE')
        limits = self.map.limits()
        cut1="curvRaw1>{curvDown}&&curvRaw1<{curvUp}&&etaRaw1>{etaDown}&&etaRaw1<{etaUp}&&phiRaw1>{phiDown}&&phiRaw1<{phiUp}".format(curvDown=limits['curv'][0],curvUp=limits['curv'][1],etaDown=limits['eta'][0],etaUp=limits['eta'][1],phiDown=limits['phi'][0],phiUp=limits['phi'][1])
        cut2="curvRaw2>{curvDown}&&curvRaw2<{curvUp}&&etaRaw2>{etaDown}&&etaRaw2<{etaUp}&&phiRaw2>{phiDown}&&phiRaw2<{phiUp}".format(curvDown=limits['curv'][0],curvUp=limits['curv'][1],etaDown=limits['eta'][0],etaUp=limits['eta'][1],phiDown=limits['phi'][0],phiUp=limits['phi'][1])
        self.tree = datas.reduce(ROOT.RooFit.EventRange(0,entries),ROOT.RooFit.Cut(cut1+'&&'+cut2))
#        self.tree = datas.reduce(ROOT.RooFit.Cut(cut1+'&&'+cut2))
        self.inputEntries = entries
        self.positiveSamples  = {}
        self.negativeSamples  = {}

        self.pairSamples  = {}


        self.w.var('curvRaw1').setMin(limits['curv'][0])
        self.w.var('curvRaw1').setMax(limits['curv'][1])
        self.w.var('curvRaw2').setMin(limits['curv'][0])
        self.w.var('curvRaw2').setMax(limits['curv'][1])
        self.w.var('etaRaw1').setMin(limits['eta'][0])
        self.w.var('etaRaw1').setMax(limits['eta'][1])
        self.w.var('etaRaw2').setMin(limits['eta'][0])
        self.w.var('etaRaw2').setMax(limits['eta'][1])
        self.w.var('phiRaw1').setMin(limits['phi'][0])
        self.w.var('phiRaw1').setMax(limits['phi'][1])
        self.w.var('phiRaw2').setMin(limits['phi'][0])
        self.w.var('phiRaw2').setMax(limits['phi'][1])

        self.tree.get().find('curvRaw1').setMin(limits['curv'][0])
        self.tree.get().find('curvRaw1').setMax(limits['curv'][1])
        self.tree.get().find('curvRaw2').setMin(limits['curv'][0])
        self.tree.get().find('curvRaw2').setMax(limits['curv'][1])
        self.tree.get().find('etaRaw1').setMin(limits['eta'][0])
        self.tree.get().find('etaRaw1').setMax(limits['eta'][1])
        self.tree.get().find('etaRaw2').setMin(limits['eta'][0])
        self.tree.get().find('etaRaw2').setMax(limits['eta'][1])
        self.tree.get().find('phiRaw1').setMin(limits['phi'][0])
        self.tree.get().find('phiRaw1').setMax(limits['phi'][1])
        self.tree.get().find('phiRaw2').setMin(limits['phi'][0])
        self.tree.get().find('phiRaw2').setMax(limits['phi'][1])

        if smear:
            self.tree=pmap.smearEbE2D(self.tree,self.w)
#            self.unsmearedTree = self.tree
#            self.tree=pmap.smear(self.tree,self.w,0.01)

        if unfold:
            self.foldedTree =self.tree
            self.tree = self.unfoldDataSet2D(self.tree)


    def getDerivativeOverContent(self,h,x):
        b = h.GetXaxis().FindBin(x)
        content = h.GetBinContent(b)
        if content==0:
            return 0.0

        if b==1:
            return (h.GetBinContent(b+1)-h.GetBinContent(b))/(content*(h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b)))
        elif b==h.GetNbinsX():
            return (h.GetBinContent(b)-h.GetBinContent(b-1))/(content*(h.GetXaxis().GetBinCenter(b)-h.GetXaxis().GetBinCenter(b-1)))
        else:
            return (h.GetBinContent(b+1)-h.GetBinContent(b-1))/(content*(h.GetXaxis().GetBinCenter(b+1)-h.GetXaxis().GetBinCenter(b-1)))


    def getDerivative2D(self,h,x,y,direction = 'X'):
        b={}
        bnext={}
        bprev={}
            
        otherdirection='X'
        if direction =='X':
            otherdirection='Y'


        b['X'] = h.GetXaxis().FindBin(x)
        b['Y'] = h.GetYaxis().FindBin(y)
        content = h.GetBinContent(h.GetBin(b['X'],b['Y']))
        
        if b[direction]==getattr(h,'GetNbins'+direction)():
            bnext[direction]=b[direction]
            bprev[direction]=b[direction]-1
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]
        elif b[direction]==1:
            bnext[direction]=b[direction]+1
            bprev[direction]=b[direction]
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]
        else:
            bnext[direction]=b[direction]+1
            bprev[direction]=b[direction]-1
            bnext[otherdirection]=b[otherdirection]
            bprev[otherdirection]=b[otherdirection]

        bnext2D=h.GetBin(bnext['X'],bnext['Y'])
        bprev2D=h.GetBin(bprev['X'],bprev['Y'])
       
        deriv=(h.GetBinContent(bnext2D)-h.GetBinContent(bprev2D))/(getattr(h,'Get'+direction+'axis')().GetBinCenter(bnext[direction])-getattr(h,'Get'+direction+'axis')().GetBinCenter(bprev[direction]))
        if content>0:
            return deriv/content
        else:
            return 0.0

    def unfoldDataSet(self,data):
        spectrum = self.convertToBinned(data,'massRaw',100).createHistogram('massRaw',100)
        newData = ROOT.RooDataSet(data.GetName(),data.GetTitle(),data.get())
        for i in range(0,data.numEntries()):
            line = data.get(i)
            x = line.find('massRaw').getVal()
            resolution = line.find('massErrRaw').getVal()
            deriv = self.getDerivativeOverContent(spectrum,x)
            line.find('massRaw').setVal(x+resolution*resolution*deriv)
            newData.add(line)
        return newData

    def unfoldDataSet2D(self,data):
        spectrum = self.convertToBinned2D(data,['curvRaw1','curvRaw2'],30,30).createHistogram('curvRaw1,curvRaw2',30,30)
        newData = ROOT.RooDataSet(data.GetName(),data.GetTitle(),data.get())
        for i in range(0,data.numEntries()):
            line = data.get(i)
            x = line.find('curvRaw1').getVal()
            y = line.find('curvRaw2').getVal()
            resolutionx = line.find('massErrRaw1').getVal()*line.find('curvRaw1').getVal()/line.find('massRaw').getVal()
            resolutiony = line.find('massErrRaw2').getVal()*line.find('curvRaw2').getVal()/line.find('massRaw').getVal()
            derivx = self.getDerivative2D(spectrum,x,y,'X')
            derivy = self.getDerivative2D(spectrum,x,y,'Y')

#            print 'before',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal(),line.find('massRaw').getVal()

            
            line.find('curvRaw1').setVal(x+resolutionx*resolutionx*derivx)
            line.find('curvRaw2').setVal(y+resolutiony*resolutiony*derivy)

            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                                   line.find('etaRaw1').getVal(),
                                   line.find('phiRaw1').getVal(),
                                   self.w.var('muMass').getVal())
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                                   line.find('etaRaw2').getVal(),
                                   line.find('phiRaw2').getVal(),
                                   self.w.var('muMass').getVal())

            line.find('massRaw').setVal((v1+v2).M())
#            print 'after',line.find('curvRaw1').getVal(),line.find('curvRaw2').getVal(),line.find('massRaw').getVal()

            newData.add(line)
        return newData
     

            

        
    def build(self,max=-1):
        for i in range(1,self.map.bins_curv()+1):
            for j in range(1,self.map.bins_eta()+1):
                for k in range(1,self.map.bins_phi()+1):
                    bin = self.map.bin(i,j,k)
                    curvDown,curvUp = self.map.boundaries_curv(bin)
                    etaDown,etaUp = self.map.boundaries_eta(bin)
                    phiDown,phiUp = self.map.boundaries_phi(bin)
                    setPos = self.tree.reduce("curvRaw1>{curvDown}&&curvRaw1<{curvUp}&&etaRaw1>{etaDown}&&etaRaw1<{etaUp}&&phiRaw1>{phiDown}&&phiRaw1<{phiUp}".format(curvDown=curvDown,curvUp=curvUp,etaDown=etaDown,etaUp=etaUp,phiDown=phiDown,phiUp=phiUp))
                    setPos.SetName("pos_"+str(bin))

                    if max>-1:
                        setPos=setPos.reduce(ROOT.RooFit.EventRange(0,max))
                    self.positiveSamples[bin]=setPos
                        
                    setNeg = self.tree.reduce("curvRaw2>{curvDown}&&curvRaw2<{curvUp}&&etaRaw2>{etaDown}&&etaRaw2<{etaUp}&&phiRaw2>{phiDown}&&phiRaw2<{phiUp}".format(curvDown=curvDown,curvUp=curvUp,etaDown=etaDown,etaUp=etaUp,phiDown=phiDown,phiUp=phiUp))
                    setNeg.SetName("neg_"+str(bin))


                    if max>-1:
                        setNeg=setNeg.reduce(ROOT.RooFit.EventRange(0,max))

                    self.negativeSamples[bin]=setNeg
        self.cache.Close()
        self.statistics()



    def buildPairs(self,minEntries = 50 , isGen = False ):
        #first create the datasets
        pairSamples={}
        for i in range(1,self.map.bins_curv()+1):
            for j in range(1,self.map.bins_eta()+1):
                for k in range(1,self.map.bins_phi()+1):
                    for l in range(1,self.map.bins_curv()+1):
                        for m in range(1,self.map.bins_eta()+1):
                            for n in range(1,self.map.bins_phi()+1):

                                bin1 = self.map.bin(i,j,k)
                                bin2 = self.map.bin(l,m,n)
                                pairSamples[str(bin1)+'_'+str(bin2)] = ROOT.RooDataSet(str(bin1)+'_'+str(bin2),
                                                                                    str(bin1)+'_'+str(bin2),
                                                                                    self.tree.get())
        print 'Datasets created now filling'


        for N in range(0,self.tree.numEntries()):
            line=self.tree.get(N)
            curvRaw1 =line.find("curvRaw1").getVal() 
            etaRaw1  =line.find("etaRaw1").getVal() 
            phiRaw1  =line.find("phiRaw1").getVal() 

            curvRaw2 =line.find("curvRaw2").getVal() 
            etaRaw2  =line.find("etaRaw2").getVal() 
            phiRaw2  =line.find("phiRaw2").getVal() 

            bin1= self.map.binFromVals(curvRaw1,etaRaw1,phiRaw1)
            bin2= self.map.binFromVals(curvRaw2,etaRaw2,phiRaw2)
            if str(bin1)+'_'+str(bin2) in pairSamples:
                pairSamples[str(bin1)+'_'+str(bin2)].add(line)

            #if it is MC add also the symmetric sample
            if isGen and  str(bin2)+'_'+str(bin1) in pairSamples:
                pairSamples[str(bin1)+'_'+str(bin2)].add(line)

            
        print 'OK now pruning the samples'

        lowStatSamples=[]
        for title,data in pairSamples.iteritems():
            if data.numEntries()<minEntries:
                lowStatSamples.append(title)
            else:
                print 'Found high stat',data.numEntries()
        for entry in lowStatSamples:
            del pairSamples[entry]

        #convert to standard format    
        print 'CONVERTING'
        for i,(title,data) in enumerate(pairSamples.iteritems()):
            bins = title.split('_')
            self.pairSamples[i] = {'bin1':int(bins[0]),'bin2':int(bins[1]),'data':data}
            
        
        self.statisticsPairs()                            


    def statistics(self):
        print 'Number of positive bins',len(self.positiveSamples)
        print 'Number of negative bins',len(self.negativeSamples)
        entries=[]
        print '----POSITIVE-----'
        for key,data in self.positiveSamples.iteritems():
            entries.append(data.numEntries())
            print self.map.boundaries_eta(key),self.map.boundaries_phi(key),self.map.boundaries_curv(key),data.numEntries()
        print '----NEGATIVE-----'
        for key,data in self.negativeSamples.iteritems():
            entries.append(data.numEntries())
            print self.map.boundaries_eta(key),self.map.boundaries_phi(key),self.map.boundaries_curv(key),data.numEntries()
        print 'Minimum number of entries',str(min(entries)), ' out of ',str(self.inputEntries),'input entries'

    def statisticsPairs(self):
        print 'Number of  bins',len(self.pairSamples)

        entries=[]

        for entry,sample in self.pairSamples.iteritems():
            entries.append(sample['data'].numEntries())
        print 'Minimum number of entries',str(min(entries)), ' out of ',str(self.tree.numEntries()),'input entries'


    def posData(self):
        return self.positiveSamples

    def negData(self):
        return self.negativeSamples


    def data(self,sign):
        if sign == 'pos':
            return self.posData()
        elif sign == 'neg':
            return self.negData()
        else:
            return self.pairSamples

    def convertToBinned(self,data,variable,bins=100):
        self.w.var(variable).setBins(bins)
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variable)),data)

    def convertToBinned2D(self,data,variables,binsx,binsy):
        self.w.var(variables[0]).setBins(binsx);
        self.w.var(variables[1]).setBins(binsy);
        
        return ROOT.RooDataHist(data.GetName(),'',ROOT.RooArgSet(self.w.var(variables[0]),self.w.var(variables[1])),data)
        

    
