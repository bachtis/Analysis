from CMGTools.MuonCalibration.tools.SimpleKalmanFilter import SimpleKalmanFilter
import ROOT
import math 
class KalmanCalibrator(SimpleKalmanFilter):
    def __init__(self,dataPairs,initialError = 0.4e-3,filename = "",iteration=-1 ):

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


        x = ROOT.TVectorD(len(self.binsPos)+len(self.binsNeg))
        x.Zero()
        

        H = ROOT.TMatrixD(len(dataPairs),len(self.iToBin))
        H.Zero()


        P = ROOT.TMatrixDSym(len(self.binsPos)+len(self.binsNeg))
        P.Zero()

        for y,inp in dataPairs.iteritems():
            xx =self.binToI['pos_'+str(inp['bin1'])]
            H[y][xx] = 1./inp['data'].mean(inp['data'].get().find('curvRaw1'))
            xx =self.binToI['neg_'+str(inp['bin2'])]
            
            H[y][xx] = 1./inp['data'].mean(inp['data'].get().find('curvRaw2'))


        for i in range(0,len(self.binToI)):
            P[i][i]=initialError*initialError

        if filename =="":   
            super(KalmanCalibrator, self).create(x,H,P)
        else:
            super(KalmanCalibrator, self).createFromFile(filename,iteration)
            
        

    def update(self,fitInputs,pmap):
        z=ROOT.TVectorD(self.H.GetNrows())
        z.Zero()
        R=ROOT.TMatrixDSym(self.H.GetNrows())
        R.Zero()
        for i,data in fitInputs.iteritems():
            z[i] =1./(data['value']*data['value'])-1
            error = 2*data['error']/(data['value']*data['value']*data['value'])
            R[i][i] = (error*error)
            ###ADD the correlation####
            
        super(KalmanCalibrator, self).iterate(z,R)

        for i,text in self.iToBin.iteritems():
            line = text.split('_')
            pmap.setData('C'+line[0],int(line[1]),self.x[i],math.sqrt(abs(self.P[i][i])))
            
        
            
            
              
              
                          
            
            
            
        

            
            
