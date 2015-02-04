from defs import *
from CMGTools.MuonCalibration.tools.KalmanCalibrator  import *

############################DATA SETS#################################

#builder = DataSetBuilder(pmap,w,'JData_Input.root','data',1000000000)
#jdata = builder.tree.reduce('massRaw>3.0&&massRaw<3.2&&abs(rapidity)<1.2')
#print 'Jpsi data',jdata.numEntries()


#builder = DataSetBuilder(pmap,w,'JMC_Input.root','data',1000000000)
#jmc = builder.tree.reduce('massRaw>3.0&&massRaw<3.2&&abs(rapidity)<1.2')
#print 'Jpsi MC',jmc.numEntries()


#builder = DataSetBuilder(pmap,w,'ZMC_Input.root','data',1000000000)
#zmc = builder.tree.reduce('massRaw>85&&massRaw<95&&abs(rapidity)<3.2')
#print 'ZMC ',zmc.numEntries()


builder = DataSetBuilder(pmap,w,'ZData_Input.root','data',1000000000)
zdata = builder.tree.reduce('massRaw>85&&massRaw<95')
print 'Z data ',zdata.numEntries()


############################DATA SETS#################################



###########################INPUTS####################################

#10x10 partition map for  B
pmap10_10 = PartitionMap(curvArr,etaArr,phiArr,"")
pmap10_10.declareData('B',0.0)


#1x10 partition map for  M
pmap1_10 = PartitionMap(curvArr,etaArr,[-math.pi,math.pi],"")
pmap1_10.declareData('M',0.0)


#1x1 partition map for A,K
pmap1_1= PartitionMap(curvArr,[-0.9,0.9],[-math.pi,math.pi],"")
pmap1_1.declareData('A',1.0)
pmap1_1.declareData('K',0.0)

infos = {}
infos['A']={'map':pmap1_1,'error':0.005}
infos['K']={'map':pmap1_1,'error':0.005}
infos['B']={'map':pmap10_10,'error':200e-6}
infos['M']={'map':pmap1_10,'error':0.015}





################################MODEL################################

def h(line,state):
    V = ROOT.TVectorD(1)

    c1 = line.find('curvRaw1').getVal()
    c2 = line.find('curvRaw2').getVal()
    e1 = line.find('etaRaw1').getVal()
    e2 = line.find('etaRaw2').getVal()
    p1 = line.find('phiRaw1').getVal()
    p2 = line.find('phiRaw2').getVal()
    st1 =math.sin(2*math.atan(math.exp(-e1))) 
    st2 =math.sin(2*math.atan(math.exp(-e2))) 


    a1_bin = state.infos['A']['map'].binFromVals(c1,e1,p1)
    a1 = state.infos['A']['map'].getData('A',a1_bin)
    a2_bin = state.infos['A']['map'].binFromVals(c2,e2,p2)
    a2 = state.infos['A']['map'].getData('A',a2_bin)

    k1_bin = state.infos['K']['map'].binFromVals(c1,e1,p1)
    k1 = state.infos['K']['map'].getData('K',k1_bin)
    k2_bin = state.infos['K']['map'].binFromVals(c2,e2,p2)
    k2 = state.infos['K']['map'].getData('K',k2_bin)


    b1_bin = state.infos['B']['map'].binFromVals(c1,e1,p1)
    b1 = state.infos['B']['map'].getData('B',b1_bin)
    b2_bin = state.infos['B']['map'].binFromVals(c2,e2,p2)
    b2 = state.infos['B']['map'].getData('B',b2_bin)

    m1_bin = state.infos['M']['map'].binFromVals(c1,e1,p1)
    m1 = state.infos['M']['map'].getData('M',m1_bin)
    m2_bin = state.infos['M']['map'].binFromVals(c2,e2,p2)
    m2 = state.infos['M']['map'].getData('M',m2_bin)


    term1 = a1-1.0+k1*e1*e1+1.0/(1+st1*m1*c1)+b1/c1
    term2 = a2-1.0+k2*e2*e2+1.0/(1+st2*m2*c2)-b2/c2


    V[0]= math.sqrt(term1*term2)   

    return V



def jacobian(line,state):
    V = ROOT.TVectorD(1)
    V = ROOT.TMatrixD(1,state.L)
    V.Zero()




    c1 = line.find('curvRaw1').getVal()
    c2 = line.find('curvRaw2').getVal()
    e1 = line.find('etaRaw1').getVal()
    e2 = line.find('etaRaw2').getVal()
    p1 = line.find('phiRaw1').getVal()
    p2 = line.find('phiRaw2').getVal()
    st1 =math.sin(2*math.atan(math.exp(-e1))) 
    st2 =math.sin(2*math.atan(math.exp(-e2))) 


    c1 = line.find('curvRaw1').getVal()
    c2 = line.find('curvRaw2').getVal()
    e1 = line.find('etaRaw1').getVal()
    e2 = line.find('etaRaw2').getVal()
    p1 = line.find('phiRaw1').getVal()
    p2 = line.find('phiRaw2').getVal()
    st1 =math.sin(2*math.atan(math.exp(-e1))) 
    st2 =math.sin(2*math.atan(math.exp(-e2))) 


    a1_bin = state.infos['A']['map'].binFromVals(c1,e1,p1)
    a1 = state.infos['A']['map'].getData('A',a1_bin)
    a2_bin = state.infos['A']['map'].binFromVals(c2,e2,p2)
    a2 = state.infos['A']['map'].getData('A',a2_bin)

    k1_bin = state.infos['K']['map'].binFromVals(c1,e1,p1)
    k1 = state.infos['K']['map'].getData('K',k1_bin)
    k2_bin = state.infos['K']['map'].binFromVals(c2,e2,p2)
    k2 = state.infos['K']['map'].getData('K',k2_bin)


    b1_bin = state.infos['B']['map'].binFromVals(c1,e1,p1)
    b1 = state.infos['B']['map'].getData('B',b1_bin)
    b2_bin = state.infos['B']['map'].binFromVals(c2,e2,p2)
    b2 = state.infos['B']['map'].getData('B',b2_bin)

    m1_bin = state.infos['M']['map'].binFromVals(c1,e1,p1)
    m1 = state.infos['M']['map'].getData('M',m1_bin)
    m2_bin = state.infos['M']['map'].binFromVals(c2,e2,p2)
    m2 = state.infos['M']['map'].getData('M',m2_bin)


    term1 = a1-1.0+k1*e1*e1+1.0/(1+st1*m1*c1)+b1/c1
    term2 = a2-1.0+k2*e2*e2+1.0/(1+st2*m2*c2)-b2/c2

    h= math.sqrt(term1*term2)   

    for N,data in state.NToI.iteritems():
        if  data['name'] =='A':
            if data['bin'] == a1_bin and (a1_bin !=a2_bin):
                V[0][N] = (0.5/h) * (term2)
            elif  data['bin'] == a2_bin and (a1_bin !=a2_bin):  
                V[0][N] = (0.5/h) * (term1)
            elif  data['bin'] == a1_bin and (a1_bin ==a2_bin):  
                V[0][N] = (0.5/h) * (term1+term2)
            else:    
                V[0][N]=0.0
        elif  data['name'] =='K':
            if data['bin'] == k1_bin and (k1_bin !=k2_bin):
                V[0][N] = (0.5/h) * (term2)*e1*e1
            elif  data['bin'] == k2_bin and (k1_bin !=k2_bin):  
                V[0][N] = (0.5/h) * (term1)*e2*e2
            elif  data['bin'] == k1_bin and (k1_bin ==k2_bin):  
                V[0][N] = (0.5/h) * (e2*e2*term1+e1*e1*term2)
            else:    
                V[0][N]=0.0
                
        elif  data['name'] =='B':
            if data['bin']==b1_bin and (b1_bin !=b2_bin):
                V[0][N] = (0.5/h)*term2/c1
            elif data['bin']==b2_bin and (b1_bin !=b2_bin):   
                V[0][N] = -(0.5/h)*term1/c2
            elif data['bin']==b1_bin and (b1_bin ==b2_bin):   
                V[0][N] = (0.5/h)*term2/c1-(0.5/h)*term1/c2
            else:    
                V[0][N]=0.0
                
        elif  data['name'] =='M':
            if data['bin']==m1_bin and (m1_bin !=m2_bin):
                V[0][N] = (0.5/h) * (-st1*c1/((1+m1*st1*c1)*(1+m1*st1*c1)))*term2
            elif data['bin']==m2_bin and (m1_bin !=m2_bin):
                V[0][N] = (0.5/h) * (-st2*c2/((1+m2*st2*c2)*(1+m2*st2*c2)))*term1
            elif data['bin']==m1_bin and (m1_bin ==m2_bin):
                V[0][N] = (0.5/h) * (-st1*c1/((1+m1*st1*c1)*(1+m1*st1*c1)))*term2 + (0.5/h) *(-st2*c2/((1+m2*st2*c2)*(1+m2*st2*c2)))*term1
            else:    
                V[0][N]=0.0


    return V


calibrator = KalmanCalibrator(infos,'kalmanScale_data_nounf.root')
calibrator.loadJPsiMatrix('../step1_KalmanInput_Jpsi/kalmanTargetJpsi_bkg.root')
calibrator.loadZMatrix('kalmanTargetZ.root')
calibrator.setModel(h,jacobian)

calibrator.updateZ(zdata,1)
calibrator.updateJPSI(jdata,0)





    
    
