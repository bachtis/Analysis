import ROOT
ROOT.gSystem.Load("libKaMuCaDerivation.so")
from math import sin,cos,pi,pow
from KaMuCa.Derivation.tools.KalmanCalibrator  import *
from KaMuCa.Derivation.tools.Partitioner  import *

############################DATA SETS#################################

fJDATA=ROOT.TFile("JDATA.root")
jdata = fJDATA.Get("data_13")
print 'Jpsi data',jdata.numEntries()

fYDATA=ROOT.TFile("YDATA.root")
ydata = fYDATA.Get("data_13")
print 'U data',ydata.numEntries()


fJMC=ROOT.TFile("JMC.root")
jmc = fJMC.Get("data_13")
print 'Jpsi MC',jmc.numEntries()

fYMC=ROOT.TFile("YMC.root")
ymc = fYMC.Get("data_13")
print 'U MC',ymc.numEntries()


jmc=jmc.reduce("abs(eta1)<2.5&&abs(eta2)<2.5")
ymc=ymc.reduce("abs(eta1)<2.5&&abs(eta2)<2.5")





###########################INPUTS####################################


phiArr=[]
for i in range(0,11):
    phiArr.append(-math.pi+i*2*math.pi/10.0)
#1x10 partition map for  B
pmap1_8 = Partitioner([1./150000.,1.0/1.0],[-2.5,-2.2,-2.0,-1.7,-1.4,-1.1,-0.9,-0.7,-0.5,-0.25,0.0,0.25,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.2,2.5],phiArr,"")
pmap1_8.declareData('B0',0)



#1x10 partition map for  M
pmap1_10 = Partitioner([1.0/150000.,1.0/1.0],[-2.5,-2.2,-2.0,-1.7,-1.4,-1.1,-0.9,-0.7,-0.5,-0.25,0.0,0.25,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.2,2.5]
,[-pi,pi],"")
pmap1_10.declareData('M',0.00)



#1x1 partition map for A,K
pmap1_1= Partitioner([1.0/150000.,1.0/1.0],[-2.6,2.6],[-pi,pi],"")
pmap1_1.declareData('A',1.0)
pmap1_1.declareData('K',0)



infos = {}
infos['A']={'map':pmap1_1,'error':0.005}
infos['K']={'map':pmap1_1,'error':0.005}
infos['B0']={'map':pmap1_8,'error':0.001}
infos['M']={'map':pmap1_10,'error':0.01}





################################MODEL################################







def jacobian(line,state):
    V = ROOT.TVectorD(state.L)
    V.Zero()


    c1 = line.find('c1').getVal()
    c2 = line.find('c2').getVal()
    e1 = line.find('eta1').getVal()
    e2 = line.find('eta2').getVal()
    p1 = line.find('phi1').getVal()
    p2 = line.find('phi2').getVal()
    st1 =math.sin(2*math.atan(math.exp(-e1))) 
    st2 =math.sin(2*math.atan(math.exp(-e2))) 
    tt1 =math.tan(2*math.atan(math.exp(-e1))) 
    tt2 =math.tan(2*math.atan(math.exp(-e2))) 


    a = state.infos['A']['map'].getData('A',13)
    k = state.infos['K']['map'].getData('K',13)


    b01_bin = state.infos['B0']['map'].binFromVals(c1,e1,p1)
    b01 = state.infos['B0']['map'].getData('B0',b01_bin)
    b02_bin = state.infos['B0']['map'].binFromVals(c2,e2,p2)
    b02 = state.infos['B0']['map'].getData('B0',b02_bin)

    m1_bin = state.infos['M']['map'].binFromVals(c1,e1,p1)
    m1 = state.infos['M']['map'].getData('M',m1_bin)
    m2_bin = state.infos['M']['map'].binFromVals(c2,e2,p2)
    m2 = state.infos['M']['map'].getData('M',m2_bin)

#    print m1_bin,m1,m2_bin,m2,c1,c2,e1,e2,p1,p2
    term1 = a-1+k/(tt1*tt1)+1.0/(1+st1*m1*c1)+b01/c1
    term2 = a-1+k/(tt2*tt2)+1.0/(1+st2*m2*c2)-b02/c2




    h= math.sqrt(term1*term2)   





    N = state.IToN['A_13']
    V[N] =(0.5/h) * (term1+term2)

    N = state.IToN['K_13']
    V[N] = (0.5/h) * (term1/(tt2*tt2)+term2/(tt1*tt1))
    
    if m1_bin != m2_bin:
        N = state.IToN['M_'+str(m1_bin)]
        V[N] = (0.5/h) * (-st1*c1)*term2/((1+m1*st1*c1)*(1+m1*st1*c1))
        N = state.IToN['M_'+str(m2_bin)]
        V[N] = (0.5/h) * (-st2*c2)*term1/((1+m2*st2*c2)*(1+m2*st2*c2))
    else:
        N = state.IToN['M_'+str(m1_bin)]
        V[N] = (0.5/h) * (-st1*c1*term2/((1+m1*st1*c1)*(1+m1*st1*c1))-st2*c2*term1/((1+m2*st2*c2)*(1+m2*st2*c2)))

    if b01_bin != b02_bin:
        N = state.IToN['B0_'+str(b01_bin)]
        V[N] = (0.5/h)*term2/c1
        N = state.IToN['B0_'+str(b02_bin)]
        V[N] = -(0.5/h)*term1/c2
    else:
        N = state.IToN['B0_'+str(b01_bin)]
        V[N] = (0.5/h)*term2/c1-(0.5/h)*term1/c2
    return V,h


calibrator = KalmanCalibrator(infos,'calibrationBinnedDATA.root')
calibrator.setModel(jacobian)
calibrator.update(jdata,1,lambda x,y: 3.09468)
calibrator.update(ydata,0,lambda x,y: 9.451862)










    
    
