import ROOT
import math 

class KalmanMatrix(object):
    def __init__(self,infos):
        self.NToI={}
        self.IToN={}
        self.infos=infos
        N=0
        self.L =0
        for t,data in infos.iteritems():
                for i in range(1,data['map'].bins_x()+1):
                    for j in range(1,data['map'].bins_y()+1):
                        for k in range(1,data['map'].bins_z()+1):
                            bin = data['map'].bin(i,j,k)
                            self.NToI[N] = {'name':t,'bin':bin}
                            self.IToN[t+"_"+str(bin)] = N
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
            P[i,i] = self.infos[data['name']]['error']*self.infos[data['name']]['error']
            i=i+1
        return P        

    
    def update(self,x,P):
        for N,data in self.NToI.iteritems():
            self.infos[data['name']]['map'].setData(data['name'],data['bin'],x[N],math.sqrt(P[N][N]))
           


class KalmanCalibrator(object):
    def __init__(self,infos,filename='results.root'):
        self.state = KalmanMatrix(infos)
        self.jacobian=None
        self.fileout = ROOT.TFile(filename,'RECREATE')
        self.calib = ROOT.SimpleKalmanCalculator(self.state.x(),self.state.P())
            
    def setModel(self,jacobian):
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
            matrixD[N,N] = 1/math.sqrt(self.calib.P()[N][N])

        PD = ROOT.TMatrixD(self.calib.P(),ROOT.TMatrixD.kMult,matrixD)
        correlation = ROOT.TMatrixD(matrixD,ROOT.TMatrixD.kMult,PD)
        for N1 in range(1,self.state.size()+1):
            for N2 in range(1,self.state.size()+1):
                bin=h.GetBin(N1,N2)
                h.SetBinContent(bin,correlation[N1-1][N2-1])
                h2.SetBinContent(bin,self.calib.P()[N1-1][N2-1])
        return h,h2
    

    def loadTarget(self,filename):
        self.f=ROOT.TFile(filename)
        self.target=self.f.Get('target_0')



    def save(self,tag):
        self.fileout.cd()
        corr,cov = self.getStateCovariance()
        corr.Write('CORRELATION_'+tag)
        cov.Write('COVARIANCE_'+tag)
        for t,data in self.state.infos.iteritems():
            data['map'].data[t].Write(t+'_'+tag)
            
        




    def update(self,data,EV,meanFunc,width=0.0,weight=1.0):
        #loop in all events
        for i in range(0,data.numEntries()):
            line = data.get(i)

            #calculate J/psi mass from matrix
            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('c1').getVal(),
                            line.find('eta1').getVal(),
                            line.find('phi1').getVal(),
                            0.1056583715)
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('c2').getVal(),
                            line.find('eta2').getVal(),
                            line.find('phi2').getVal(),
                            0.1056583715)
            
            v = v1+v2


            mass=line.find("mass").getVal()
            massErr=line.find("massErr").getVal()

            rapidity = v.Rapidity()
            mean = meanFunc(rapidity,massErr*massErr/(mass*mass))

            z = mass/mean
            resolution2 = (massErr*massErr+width*width)/(mean *mean*weight*weight)

            H,h= self.jacobian(line,self.state)

            self.calib.iterate(z-h,resolution2,H)
            self.state.update(self.calib.x(),self.calib.P())

            if i % 200000 ==0:
                print 'EVENT ---------',i
                self.save(str(EV)+'_'+str(i))    


        self.save(str(EV))    



            
            
