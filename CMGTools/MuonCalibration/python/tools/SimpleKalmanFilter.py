import ROOT


#Kalman filter is modeled as follows
# x(k) =  x(k-1) 
# z(k) = Hx(k) +v(k) [R]  <We do the same measurement all time so H is constant>
#


class SimpleKalmanFilter(object):
    def create(self,x,P):
        self.x = x
        self.P = P
        self.i=0
        self.history = {0:x}

    def save(self,filename):
        f = ROOT.TFile(filename,"UPDATE")
        f.cd()
        self.x.Write("x"+str(self.i))
        self.P.Write("P"+str(self.i))
        self.H.Write("H"+str(self.i))
        self.HT.Write("HT"+str(self.i))
        f.Close()


    def createFromFile(self,filename,i):
        f = ROOT.TFile(filename)
        self.x = f.Get("x"+str(i))
        self.P = f.Get("P"+str(i))
        self.H = f.Get("H"+str(i))
        self.HT = f.Get("HT"+str(i))
        self.i=i+1

        
    def iterate(self,z,R,H,h,history=False):    
        self.H=H
        self.HT = ROOT.TMatrixD(H)
        self.HT.T()

        residual = z-h

        PH = ROOT.TMatrixD(self.P,ROOT.TMatrixD.kMult,self.HT)
        HPH = ROOT.TMatrixD(self.H,ROOT.TMatrixD.kMult,PH)
        S = HPH+R

#        import pdb;pdb.set_trace()
        #Be a bit picky on inversion and try different algorithms
#        SVD = ROOT.TDecompSVD(S)
#        SVDinv = SVD.Invert()
#        S = SVDinv
        S.Invert()
#        print 'SVD inversion succeeded '
        HSinv = ROOT.TMatrixD(self.HT)
        HSinv*=S
        K=ROOT.TMatrixD(self.P,ROOT.TMatrixD.kMult,HSinv)
        self.x = self.x+K*residual

        HP=ROOT.TMatrixD(self.H,ROOT.TMatrixD.kMult,self.P)

        KHP=ROOT.TMatrixD(K,ROOT.TMatrixD.kMult,HP)

        self.P =self.P-KHP
        self.i=self.i+1
        if history:
            self.history[self.i]=self.x

            
        
    
        
