from CMGTools.MuonCalibration.tools.SimpleKalmanFilter import SimpleKalmanFilter
import ROOT
import math 

class KalmanMatrix(object):
    def __init__(self,pmap,infos):
        self.NToI={}
        self.map = pmap
        N=0
        self.L =0
        for t in infos:
            for i in range(1,self.map.bins_curv()+1):
                for j in range(1,self.map.bins_eta()+1):
                    for k in range(1,self.map.bins_phi()+1):
                        bin = self.map.bin(i,j,k)
                        self.NToI[N] = {'name':t,'bin':bin}
                        N=N+1
                        self.L = self.L+1

        
    def size(self):
        return self.L
                        
    def x(self):
        V = ROOT.TVectorD(self.L)
        V.Zero()

        for N,data in self.NToI.iteritems():
            V[N] = self.map.getData(data['name'],data['bin'])
        return V


    def P(self):
        P = ROOT.TMatrixDSym(self.L)
        P.Zero()
        i=0
        for N,data in self.NToI.iteritems():
            if data['name']=='A':
                P[i][i] = 0.005*0.005
            elif data['name']=='B':
                P[i][i] = 100e-6*100e-6
            elif data['name']=='M':
                P[i][i] = 0.02*0.02
            i=i+1
        return P        

    def H(self,line):
        V = ROOT.TMatrixD(1,self.L)
        V.Zero()

        bin1 = self.map.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
        bin2 = self.map.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())
        
        c1 = line.find('curvRaw1').getVal()
        c2 = line.find('curvRaw2').getVal()
        a1 = self.map.getData('A',bin1)
        a2 = self.map.getData('A',bin2)
        b1 = self.map.getData('B',bin1)
        b2 = self.map.getData('B',bin2)
        m1 = self.map.getData('M',bin1)
        m2 = self.map.getData('M',bin2)

        for N,data in self.NToI.iteritems():
            if ((data['bin'] ==bin1 or data['bin']==bin2) and (bin1==bin2)):
                if data['name'] =='A':
                    V[0][N] = (a2+1.0/(1+m2*c2)-b2/c2)+(a1+1./(1+m1*c1)+b1/c1)
                if data['name'] =='B':
                    V[0][N] = (a2+1.0/(1+m2*c2)-b2/c2)/c1-(a1+1./(1+m1*c1)+b1/c1)/c2
                if data['name'] =='M':
                    V[0][N] = -c1*(a2+1.0/(1+m2*c2)-b2/c2)/((1+m1*c1)*(1+m1*c1))-c2*(a1+1./(1+m1*c1)+b1/c1)/((1+m2*c2)*(1+m2*c2))
            elif data['bin'] ==bin1:
                if data['name'] =='A':
                    V[0][N] = a2-1.0+1.0/(1.0+m2*c2)-b2/c2
                if data['name'] =='B':
                    V[0][N] = (a2-1.0+1.0/(1.0+m2*c2)-b2/c2)/c1
                if data['name'] =='M':
                    V[0][N] = -c1*(a2-1.0 +1.0/(1.0+m2*c2)-b2/c2)/((1+m1*c1)*(1+m1*c1))
            elif data['bin'] ==bin2:
                if data['name'] =='A':
                    V[0][N] = a1-1.0+1.0/(1.0+m1*c1)+b1/c1
                if data['name'] =='B':
                    V[0][N] = -(a1-1.0+1.0/(1.0+m1*c1)+b1/c1)/c2
                if data['name'] =='M':
                    V[0][N] = -c2*(a1-1.0+1.0/(1.0+m1*c1)+b1/c1)/((1+m2*c2)*(1+m2*c2))
            else:
                V[0][N]=0.0
                
        return V



    def HNOMAT(self,line):
        V = ROOT.TMatrixD(1,self.L)
        V.Zero()

        bin1 = self.map.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
        bin2 = self.map.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())
        
        c1 = line.find('curvRaw1').getVal()
        c2 = line.find('curvRaw2').getVal()
        a1 = self.map.getData('A',bin1)
        a2 = self.map.getData('A',bin2)
        b1 = self.map.getData('B',bin1)
        b2 = self.map.getData('B',bin2)


        for N,data in self.NToI.iteritems():
            if ((data['bin'] ==bin1 or data['bin']==bin2) and (bin1==bin2)):
                if data['name'] =='A':
                    V[0][N] = 2*a1+(1./c1-1./c2)*b1
                if data['name'] =='B':
                    V[0][N] = (1./c1-1./c2)*a1-2*b1/(c1*c2)
            elif data['bin'] ==bin1:
                if data['name'] =='A':
                    V[0][N] = a2-b2/c2
                if data['name'] =='B':
                    V[0][N] = (a2-b2/c2)/c1
            elif data['bin'] ==bin2:
                if data['name'] =='A':
                    V[0][N] = a1+b1/c1
                if data['name'] =='B':
                    V[0][N] = -(a1+b1/c1)/c2
            else:
                V[0][N]=0.0
                
        return V
    



    def hNOMAT(self,line):
        V = ROOT.TVectorD(1)

        c1 = line.find('curvRaw1').getVal()
        c2 = line.find('curvRaw2').getVal()
        bin1 = self.map.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
        bin2 = self.map.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())


        a1 = self.map.getData('A',bin1)
        b1 = self.map.getData('B',bin1)

        a2= self.map.getData('A',bin2)
        b2 = self.map.getData('B',bin2)

        V[0]= (a1+b1/c1)*(a2-b2/c2)

        return V


    def h(self,line):
        V = ROOT.TVectorD(1)

        c1 = line.find('curvRaw1').getVal()
        c2 = line.find('curvRaw2').getVal()
        bin1 = self.map.binFromVals(line.find('curvRaw1').getVal(),
                                    line.find('etaRaw1').getVal(),
                                    line.find('phiRaw1').getVal())
        bin2 = self.map.binFromVals(line.find('curvRaw2').getVal(),
                                    line.find('etaRaw2').getVal(),
                                    line.find('phiRaw2').getVal())


        a1 = self.map.getData('A',bin1)
        b1 = self.map.getData('B',bin1)
        m1 = self.map.getData('M',bin1)

        a2= self.map.getData('A',bin2)
        b2 = self.map.getData('B',bin2)
        m2 = self.map.getData('M',bin2)


        V[0]= (a1-1+1/(1+m1*c1)+b1/c1)*(a2-1.+1.0/(1+m2*c2)-b2/c2)

        return V


    def update(self,x,P):
        for N,data in self.NToI.iteritems():
            self.map.setData(data['name'],data['bin'],x[N],math.sqrt(P[N][N]))
            



class KalmanCalibrator(SimpleKalmanFilter):
    def __init__(self,data,pmap,infos,filename = "",iteration=-1 ):
        self.data = data
        self.state = KalmanMatrix(pmap,infos)
        self.map = pmap
        P = self.state.P()





        if filename =="":   
            super(KalmanCalibrator, self).create(self.state.x(),P)
        else:
            super(KalmanCalibrator, self).createFromFile(filename,iteration)
            
        

    def update(self,EV,mean= 3.0950064614389796,width = 0.0):
        #loop in all events
        for i in range(0,self.data.numEntries()):
#        for i in range(0,10000):
            line = self.data.get(i)
            z=ROOT.TVectorD(1)
            z[0] = (line.find('massRaw').getVal()*line.find('massRaw').getVal())/(mean*mean)
            R=ROOT.TMatrixD(1,1)
            error = 2*line.find('massRaw').getVal()*line.find('massErrRaw').getVal()/(mean*mean)
            R[0][0] =error*error +  2*line.find('massRaw').getVal()*line.find('massRaw').getVal()*width*width/(mean*mean*mean)


            bin1 = self.state.map.binFromVals(line.find('curvRaw1').getVal(),
                                        line.find('etaRaw1').getVal(),
                                        line.find('phiRaw1').getVal())
            bin2 = self.state.map.binFromVals(line.find('curvRaw2').getVal(),
                                        line.find('etaRaw2').getVal(),
                                        line.find('phiRaw2').getVal())

            H= self.state.H(line)
            h =self.state.h(line) 
            
            super(KalmanCalibrator, self).iterate(z,R,H,h)

            if i % 20000 ==0:
                print 'EVENT ---------',i
#                for i in range(0, self.state.size()):
#                               print self.x[i],'+-',math.sqrt(self.P[i][i]),'H=',H[0][i]               

                self.map.save(str(EV)+'_'+str(i))               
            self.state.update(self.x,self.P)

            
        

    def updateZ(self,EV,mean = 90.73183839410241,width = 1.9065751749239424):
        #loop in all events
        for i in range(0,self.data.numEntries()):
#        for i in range(0,10000):
            line = self.data.get(i)
            z=ROOT.TVectorD(1)
            z[0] = (line.find('massRaw').getVal()*line.find('massRaw').getVal())/(mean*mean)
            R=ROOT.TMatrixD(1,1)
            error = 2*line.find('massRaw').getVal()*line.find('massErrRaw').getVal()/(mean*mean)
            R[0][0] =error*error +  2*line.find('massRaw').getVal()*line.find('massRaw').getVal()*width*width/(mean*mean*mean)


            bin1 = self.state.map.binFromVals(line.find('curvRaw1').getVal(),
                                        line.find('etaRaw1').getVal(),
                                        line.find('phiRaw1').getVal())
            bin2 = self.state.map.binFromVals(line.find('curvRaw2').getVal(),
                                        line.find('etaRaw2').getVal(),
                                        line.find('phiRaw2').getVal())

            H= self.state.HNOMAT(line)
            h =self.state.hNOMAT(line) 
            
            super(KalmanCalibrator, self).iterate(z,R,H,h)

            if i % 20000 ==0:
                print 'EVENT ---------',i
#                for i in range(0, self.state.size()):
#                               print self.x[i],'+-',math.sqrt(self.P[i][i]),'H=',H[0][i]               

                self.map.save(str(EV)+'_'+str(i))               
            self.state.update(self.x,self.P)

            
        

            
              
              
                          
            
            
            
        

            
            
