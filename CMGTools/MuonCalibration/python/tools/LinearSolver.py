

import ROOT
class LinearSolver(object):
    def __init__(self,dataPairs):

        binsPosArr=[]
        binsNegArr=[]

        for N,inp in dataPairs.iteritems():
            binsPosArr.append(inp['bin1'])
            binsNegArr.append(inp['bin2'])

        self.binsPos = set(binsPosArr)
        self.binsNeg = set(binsNegArr)

        self.iToBin = {}
        self.binToI = {}


        for i,bin in enumerate(self.binsPos):
            self.iToBin[i] = 'pos_'+str(bin)
            self.binToI['pos_'+str(bin)] = i

        for i,bin in enumerate(self.binsNeg):
            self.iToBin[i+len(self.binsPos)] = 'neg_'+str(bin)
            self.binToI['neg_'+str(bin)] =i+len(self.binsPos)



        self.H = ROOT.TMatrixD(len(dataPairs),len(self.iToBin))
        self.H.Zero()


        for y,(N,inp) in enumerate(dataPairs.iteritems()):
            xx =self.binToI['pos_'+str(inp['bin1'])]
            self.H[y][xx] = 1./inp['curv1']
            xx =self.binToI['neg_'+str(inp['bin2'])]
            self.H[y][xx] = 1./inp['curv2']

        print 'MATRICE H'
#        self.H.Print()
        
        

    def solve(self,fitInputs,pmap):
        z=ROOT.TVectorD(self.H.GetNrows())
        z.Zero()
        e=ROOT.TVectorD(self.H.GetNrows())
        e.Zero()

        for i,(N,data) in enumerate(fitInputs.iteritems()):
            error = 2*data['error']/(data['value']*data['value']*data['value'])
            z[i] =(1./(data['value']*data['value'])-1)
#            for j in range(0,self.H.GetNcols()):
#                self.H[i][j] = self.H[i][j]/(error*error)


        svd = ROOT.TDecompSVD(self.H)
        ok = True
        solution = svd.Solve(z)

        print 'RESULT = ',solution
        z.Print()
        for i,text in self.iToBin.iteritems():
            line = text.split('_')
            pmap.setData('C'+line[0],int(line[1]),z[i],0.0)
            
        
            
            
              
              
                          
            
            
            
        

            
            
