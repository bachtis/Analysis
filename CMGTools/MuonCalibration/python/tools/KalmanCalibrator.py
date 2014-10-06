from CMGTools.MuonCalibration.tools.SimpleKalmanFilter import SimpleKalmanFilter
import ROOT
import math 

class KalmanMatrix(object):
    def __init__(self,infos):
        self.NToI={}
        self.infos=infos
        N=0
        self.L =0
        for t,data in infos.iteritems():
                for i in range(1,data['map'].bins_curv()+1):
                    for j in range(1,data['map'].bins_eta()+1):
                        for k in range(1,data['map'].bins_phi()+1):
                            bin = data['map'].bin(i,j,k)
                            self.NToI[N] = {'name':t,'bin':bin}
                            N=N+1
                            self.L = self.L+1
    def size(self):
        return self.L

                       
    def x(self):
        V = ROOT.TVectorD(self.L)
        V.Zero()

        for N,data in self.NToI.iteritems():
            V[N] = self.infos[data['name']]['map'].getData(data['name'],data['bin'])
        return V


    def P(self):
        P = ROOT.TMatrixDSym(self.L)
        P.Zero()
        i=0
        for N,data in self.NToI.iteritems():
            P[i][i] = self.infos[data['name']]['error']*self.infos[data['name']]['error']
            i=i+1
        return P        

    
    def update(self,x,P):
        for N,data in self.NToI.iteritems():
            self.infos[data['name']]['map'].setData(data['name'],data['bin'],x[N],math.sqrt(P[N][N]))
           


class KalmanCalibrator(SimpleKalmanFilter):
    def __init__(self,infos,filename='results.root'):
        self.state = KalmanMatrix(infos)
        P = self.state.P()
        self.h=None
        self.jacobian=None
        self.fileout = ROOT.TFile(filename,'RECREATE')
        super(KalmanCalibrator, self).create(self.state.x(),P)
            
    def setModel(self,h,jacobian):
        self.h=h
        self.jacobian=jacobian

    def getStateCovariance(self):
        h=ROOT.TH2F('correlationMatrix','correlationMatrix',self.state.size(),0,self.state.size(),self.state.size(),0,self.state.size())
        h2=ROOT.TH2F('covarianceMatrix','covarianceMatrix',self.state.size(),0,self.state.size(),self.state.size(),0,self.state.size())
        matrixD = ROOT.TMatrixD(self.state.size(),self.state.size())
        matrixD.Zero()

        for N,data in self.state.NToI.iteritems():
            h.GetXaxis().SetBinLabel(N+1,data['name']+'_'+str(data['bin']))
            h.GetYaxis().SetBinLabel(N+1,data['name']+'_'+str(data['bin']))
            h2.GetXaxis().SetBinLabel(N+1,data['name']+'_'+str(data['bin']))
            h2.GetYaxis().SetBinLabel(N+1,data['name']+'_'+str(data['bin']))
            matrixD[N][N] = 1/math.sqrt(self.P[N][N])

        PD = ROOT.TMatrixD(self.P,ROOT.TMatrixD.kMult,matrixD)
        correlation = ROOT.TMatrixD(matrixD,ROOT.TMatrixD.kMult,PD)
        for N1 in range(1,self.state.size()+1):
            for N2 in range(1,self.state.size()+1):
                bin=h.GetBin(N1,N2)
                h.SetBinContent(bin,correlation[N1-1][N2-1])
                h2.SetBinContent(bin,self.P[N1-1][N2-1])
        return h,h2
    
    def loadJPsiMatrix(self,filename):
        self.f=ROOT.TFile(filename)
        self.JMATRIX=self.f.Get('mass_fit')
        self.jpsiMatrix = self.JMATRIX.ProjectionY("jpsimatrix",1,1,1,1)

    def loadZMatrix(self,filename):
        self.f2=ROOT.TFile(filename)
        self.ZMATRIX=self.f2.Get('mass_fit')
        self.zMatrix = self.ZMATRIX.ProjectionY("zmatrix",1,1,1,1)
        self.ZMATRIXRMS=self.f2.Get('width_fit')
        self.zMatrixRMS = self.ZMATRIXRMS.ProjectionY("zmatrixRMS",1,1,1,1)

    def save(self,tag):
        self.fileout.cd()
        corr,cov = self.getStateCovariance()
        corr.Write('CORRELATION_'+tag)
        cov.Write('COVARIANCE_'+tag)
        for t,data in self.state.infos.iteritems():
            data['map'].data[t].Write(t+'_'+tag)
            
        


    def updateJPSI(self,data,EV):
        if self.jpsiMatrix is None:
            print 'MATRIX FOR J/psi not found'
            return

           
        #loop in all events
        for i in range(0,data.numEntries()):
            line = data.get(i)

            #calculate J/psi mass from matrix
            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                            line.find('etaRaw1').getVal(),
                            line.find('phiRaw1').getVal(),
                            0.1056583715)
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                            line.find('etaRaw2').getVal(),
                            line.find('phiRaw2').getVal(),
                            0.1056583715)
            
            v = v1+v2

            matrix_bin = self.jpsiMatrix.GetXaxis().FindBin(v.Eta())
            mean = self.jpsiMatrix.GetBinContent(matrix_bin)



            z=ROOT.TVectorD(1)
            z[0] = line.find('massRaw').getVal()/mean
            R=ROOT.TMatrixD(1,1)
            error = line.find('massErrRaw').getVal()/mean
            R[0][0] =error*error 

            H= self.jacobian(line,self.state)
            h =self.h(line,self.state) 
            
            super(KalmanCalibrator, self).iterate(z,R,H,h)

            if i % 20000 ==0:
                print 'EVENT ---------',i
                self.save(str(EV)+'_'+str(i))    


            self.state.update(self.x,self.P)

        self.save(str(EV))    
    


    def updateZ(self,data,EV):

        if self.zMatrix is None:
            print 'MATRIX FOR Z not found'
            return
        if self.zMatrixRMS is None:
            print 'RMS MATRIX FOR Z not found'
            return

            
        #loop in all events
        for i in range(0,data.numEntries()):
            line = data.get(i)

            #calculate J/psi mass from matrix
            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                            line.find('etaRaw1').getVal(),
                            line.find('phiRaw1').getVal(),
                            0.1056583715)
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                            line.find('etaRaw2').getVal(),
                            line.find('phiRaw2').getVal(),
                            0.1056583715)
            
            v = v1+v2

            matrix_bin = self.zMatrix.GetXaxis().FindBin(v.Eta())
            mean = self.zMatrix.GetBinContent(matrix_bin)
            width = self.zMatrixRMS.GetBinContent(matrix_bin)


            z=ROOT.TVectorD(1)
            z[0] = line.find('massRaw').getVal()/mean
            R=ROOT.TMatrixD(1,1)
            R[0][0] = (line.find('massErrRaw').getVal()*line.find('massErrRaw').getVal())/(mean*mean) + line.find('massRaw').getVal()*line.find('massRaw').getVal()*width*width/(mean*mean*mean*mean)


            H= self.jacobian(line,self.state)
            h =self.h(line,self.state) 
            
            super(KalmanCalibrator, self).iterate(z,R,H,h)


            if i % 20000 ==0:
                print 'EVENT ---------',i
                self.save(str(EV)+'_'+str(i))    


            self.state.update(self.x,self.P)

        self.save(str(EV))    

            

            
              
              
                          
            
            
            
        

            
            
