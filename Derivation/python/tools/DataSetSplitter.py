import ROOT



class DataSetSplitter(object):
    def __init__(self,partitioner,fx = lambda x: 0.0,fy = lambda x: 0.0,fz = lambda x: 0.0):
        self.pmap = partitioner
        self.fx=fx
        self.fy=fy
        self.fz=fz
        self.data={}
    
    def split(self,fin,fout,modifiers,selectors=[]):
        FIN=ROOT.TFile(fin)
        data=FIN.Get('data')
        self.FOUT=ROOT.TFile(fout,'RECREATE')
        self.FOUT.cd()
        for i in range(1,self.pmap.bins_x()+1):
            for j in range(1,self.pmap.bins_y()+1):
                for k in range(1,self.pmap.bins_z()+1):
                    bin = self.pmap.bin(i,j,k)
                    self.data[bin] = ROOT.RooDataSet(data.GetName()+'_'+str(bin),data.GetName(),data.get())

        for i in range(0,data.numEntries()):
            line=data.get(i)
            for m in modifiers:
                m(line)

            passed=True    
            for s in selectors:
                if not s(line):
                    passed=False
                    break
            
            bin=self.pmap.binFromVals(self.fx(line),self.fy(line),self.fz(line))
            if bin<0 or (not passed):
                continue
            self.data[bin].add(line)
        self.FOUT.cd()
        self.statistics()    
        for bin,data  in self.data.iteritems():
            data.Write()


    def splitPairs(self,fin,fout,modifiers,selectors=[]):
        import itertools
        FIN=ROOT.TFile(fin)
        data=FIN.Get('data')
        self.FOUT=ROOT.TFile(fout,'RECREATE')
        self.FOUT.cd()

        allBins=[]
        for i in range(1,self.pmap.bins_x()+1):
            for j in range(1,self.pmap.bins_y()+1):
                for k in range(1,self.pmap.bins_z()+1):
                    allBins.append(self.pmap.bin(i,j,k))
        for bin1,bin2 in itertools.combinations(allBins,2):    
            self.data[str(bin1)+"_"+str(bin2)] = ROOT.RooDataSet(data.GetName()+'_'+str(bin1)+"_"+str(bin2),data.GetName(),data.get())

        for i in range(0,data.numEntries()):
            line=data.get(i)
            for m in modifiers:
                m(line)

            passed=True    
            for s in selectors:
                if not s(line):
                    passed=False
                    break
            
            (x1,x2) =self.fx(line)     
            (y1,y2) =self.fy(line)     
            (z1,z2) =self.fz(line)     
            bin1=self.pmap.binFromVals(x1,y1,z1)
            bin2=self.pmap.binFromVals(x2,y2,z2)
            if bin1<0 or bin2<0 or (not passed):
                continue
            if str(bin1)+'_'+str(bin2) in self.data.keys():
                self.data[str(bin1)+"_"+str(bin2)].add(line)
            elif str(bin2)+'_'+str(bin1) in self.data.keys():             
                self.data[str(bin2)+"_"+str(bin1)].add(line)
            self.FOUT.cd()
        for label,data  in self.data.iteritems():
            data.Write()



    def load(self,filename):
        self.FOUT = ROOT.TFile(filename)
        for i in range(1,self.pmap.bins_x()+1):
            for j in range(1,self.pmap.bins_y()+1):
                for k in range(1,self.pmap.bins_z()+1):
                    bin = self.pmap.bin(i,j,k)
                    self.data[bin]=self.FOUT.Get("data_"+str(bin))

    def loadPairs(self,filename):
        import itertools
        self.FOUT = ROOT.TFile(filename)
        allBins=[]
        for i in range(1,self.pmap.bins_x()+1):
            for j in range(1,self.pmap.bins_y()+1):
                for k in range(1,self.pmap.bins_z()+1):
                    allBins.append(self.pmap.bin(i,j,k))

        for bin1,bin2 in itertools.combinations(allBins,2):    
            self.data[str(bin1)+"_"+str(bin2)]=self.FOUT.Get('data_'+str(bin1)+"_"+str(bin2))

    
    def statistics(self):
        print 'Number of  bins',len(self.data)
        entries=[]
        for key,data in self.data.iteritems():
            entries.append(data.numEntries())
            print self.pmap.boundaries_x(key),self.pmap.boundaries_y(key),self.pmap.boundaries_z(key),data.numEntries()
        print 'Minimum number of entries=', (min(entries))
        
