import ROOT
ROOT.gSystem.Load("libKaMuCaDerivation.so")
from math import sin,cos,pi,pow
from KaMuCa.Derivation.tools.KalmanCalibrator  import *
from KaMuCa.Derivation.tools.Partitioner  import *

############################DATA SETS#################################

fJDATA=ROOT.TFile("JDATA8TeV.root")
jdata = fJDATA.Get("data_13")
print 'Jpsi data',jdata.numEntries()

#fYDATA=ROOT.TFile("YDATA13TeV.root")
#ydata = fYDATA.Get("data_13")
#print 'U data',ydata.numEntries()


fJMC=ROOT.TFile("JMC8TeV.root")
jmc = fJMC.Get("data_13")
print 'Jpsi MC',jmc.numEntries()

#fYMC=ROOT.TFile("YMC13TeV.root")
#ymc = fYMC.Get("data_13")
#print 'U MC',ymc.numEntries()


jmc=jmc.reduce("abs(eta1)<2.5&&abs(eta2)<2.5")
#ymc=ymc.reduce("abs(eta1)<2.5&&abs(eta2)<2.5")





###########################INPUTS####################################

#1x10 partition map for  B
#pmap1_8 = Partitioner([1./150000.,1.0/1.0],[-2.6,-2.2,-2.0,-1.7,-1.4,-1.1,-0.9,-0.7,-0.5,-0.25,0.0,0.25,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.2,2.6],[-pi,pi],"")
pmap1_8 = Partitioner([1./150000.,1.0/1.0],[-2.6,-2.0,-1.7,-1.4,-1.1,-0.9,-0.4,0.0,0.4,0.9,1.1,1.4,1.7,2.0,2.6],[-pi,pi],"")
pmap1_8.declareData('B0',0)
pmap1_8.declareData('B1',0)
pmap1_8.declareData('B2',0)
pmap1_8.declareData('C1',0)
pmap1_8.declareData('C2',0)



#1x10 partition map for  M
pmap1_10 = Partitioner([1.0/150000.,1.0/1.0],[-2.6,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6],[-pi,pi],"")
pmap1_10.declareData('M',0.00)



#1x1 partition map for A,K
pmap1_1= Partitioner([1.0/150000.,1.0/1.0],[-2.6,2.6],[-pi,pi],"")
pmap1_1.declareData('A',1.0)
pmap1_1.declareData('K',0)



infos = {}
infos['A']={'map':pmap1_1,'error':0.005}
infos['K']={'map':pmap1_1,'error':0.005}
infos['B0']={'map':pmap1_8,'error':500e-6}
infos['B1']={'map':pmap1_8,'error':500e-6}
infos['B2']={'map':pmap1_8,'error':500e-6}
infos['C1']={'map':pmap1_8,'error':500e-6}
infos['C2']={'map':pmap1_8,'error':500e-6}
infos['M']={'map':pmap1_10,'error':0.05}





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

    b11_bin = state.infos['B1']['map'].binFromVals(c1,e1,p1)
    b11 = state.infos['B1']['map'].getData('B1',b11_bin)
    b12_bin = state.infos['B1']['map'].binFromVals(c2,e2,p2)
    b12 = state.infos['B1']['map'].getData('B1',b12_bin)

    b21_bin = state.infos['B2']['map'].binFromVals(c1,e1,p1)
    b21 = state.infos['B2']['map'].getData('B2',b21_bin)
    b22_bin = state.infos['B2']['map'].binFromVals(c2,e2,p2)
    b22 = state.infos['B2']['map'].getData('B2',b22_bin)


    c11_bin = state.infos['C1']['map'].binFromVals(c1,e1,p1)
    c11 = state.infos['C1']['map'].getData('C1',c11_bin)
    c12_bin = state.infos['C1']['map'].binFromVals(c2,e2,p2)
    c12 = state.infos['C1']['map'].getData('C1',c12_bin)

    c21_bin = state.infos['C2']['map'].binFromVals(c1,e1,p1)
    c21 = state.infos['C2']['map'].getData('C2',c21_bin)
    c22_bin = state.infos['C2']['map'].binFromVals(c2,e2,p2)
    c22 = state.infos['C2']['map'].getData('C2',c22_bin)


    m1_bin = state.infos['M']['map'].binFromVals(c1,e1,p1)
    m1 = state.infos['M']['map'].getData('M',m1_bin)
    m2_bin = state.infos['M']['map'].binFromVals(c2,e2,p2)
    m2 = state.infos['M']['map'].getData('M',m2_bin)

    term1 = a+k*e1*e1-st1*m1*c1+(b01+b11*sin(p1)+b21*sin(2*p1) +c11*cos(p1)+c21*cos(2*p1))/c1
    term2 = a+k*e2*e2-st2*m2*c2-(b02+b12*sin(p2)+b22*sin(2*p2) +c12*cos(p2)+c22*cos(2*p2))/c2



    h= math.sqrt(term1*term2)   


    N = state.IToN['A_13']
    V[N] =(0.5/h) * (term1+term2)

    N = state.IToN['K_13']
    V[N] = (0.5/h) * (term1*e2*e2+term2*e1*e1)
    
    if m1_bin != m2_bin:
        N = state.IToN['M_'+str(m1_bin)]
        V[N] = (0.5/h) * (-st1*c1)*term2
        N = state.IToN['M_'+str(m2_bin)]
        V[N] = (0.5/h) * (-st2*c2)*term1
    else:
        N = state.IToN['M_'+str(m1_bin)]
        V[N] = (0.5/h) * (-st1*c1*term2-st2*c2*term1)

    if b01_bin != b02_bin:
        N = state.IToN['B0_'+str(b01_bin)]
        V[N] = (0.5/h)*term2/c1
        N = state.IToN['B0_'+str(b02_bin)]
        V[N] = -(0.5/h)*term1/c2
    else:
        N = state.IToN['B0_'+str(b01_bin)]
        V[N] = (0.5/h)*term2/c1-(0.5/h)*term1/c2


    if b11_bin != b12_bin:
        N = state.IToN['B1_'+str(b11_bin)]
        V[N] = (0.5/h)*sin(p1)*term2/c1
        
        N = state.IToN['B1_'+str(b12_bin)]
        V[N] = -(0.5/h)*sin(p2)*term1/c2
    else:
        N = state.IToN['B1_'+str(b11_bin)]
        V[N] = (0.5/h)*term2*sin(p1)/c1-(0.5/h)*sin(p2)*term1/c2


    if b21_bin != b22_bin:
        N = state.IToN['B2_'+str(b21_bin)]
        V[N] = (0.5/h)*sin(2*p1)*term2/c1
        N = state.IToN['B2_'+str(b22_bin)]
        V[N] = -(0.5/h)*sin(2*p2)*term1/c2
    else:
        N = state.IToN['B2_'+str(b21_bin)]
        V[N] = (0.5/h)*term2*sin(2*p1)/c1-(0.5/h)*sin(2*p2)*term1/c2


    if c11_bin != c12_bin:
        N = state.IToN['C1_'+str(c11_bin)]
        V[N] = (0.5/h)*cos(p1)*term2/c1
        N = state.IToN['C1_'+str(c12_bin)]
        V[N] = -(0.5/h)*cos(p2)*term1/c2

    else:
        N = state.IToN['C1_'+str(c11_bin)]
        V[N] = (0.5/h)*term2*cos(p1)/c1-(0.5/h)*cos(p2)*term1/c2

    if c21_bin != c22_bin:
        N = state.IToN['C2_'+str(c21_bin)]
        V[N] = (0.5/h)*cos(2*p1)*term2/c1
        N = state.IToN['C2_'+str(c22_bin)]
        V[N] = -(0.5/h)*cos(2*p2)*term1/c2
    else:
        N = state.IToN['C2_'+str(c21_bin)]
        V[N] = (0.5/h)*term2*cos(2*p1)/c1-(0.5/h)*cos(2*p2)*term1/c2

    return V,h


calibrator = KalmanCalibrator(infos,'calibrationDATA8TeVETA.root')
calibrator.setModel(jacobian)

calibrator.update(jdata,0,lambda x,y: 3.09522+8.4272e-5*x*x*x*x)
#calibrator.update(jmc,0,lambda x,y: 3.09519+5.9435e-5*x*x*x*x)











    
    
