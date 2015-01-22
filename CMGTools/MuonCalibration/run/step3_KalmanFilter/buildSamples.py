from defs import *

print 'This step is needed only for data to apply the bfield map corrections'


print 'J /psi data'
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = builder.tree.reduce('massErrRaw1<0.03&&massErrRaw2<0.03&&curvRaw1>0.04&&curvRaw2>0.04') 
builder.tree = correctDataSet(builder.tree,True,False,True,False)
unfolding=AnalyticalUnfolding(builder.tree,25,0.04,1./5.5)
builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('JData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'J /psi MC'

builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.tree = builder.tree.reduce('massErrRaw1<0.03&&massErrRaw2<0.03&&curvRaw1>0.04&&curvRaw2>0.04') 
builder.tree = correctDataSet(builder.tree,False,False,True,False)
unfolding=AnalyticalUnfolding(builder.tree,25,0.04,1./5.5)
builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('JMC_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()



print 'Z data'

builder = DataSetBuilder(pmap,w,'../../data/ZDATA.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5&&massErrRaw1<1&&massErrRaw2<1&curvRaw1<0.033&&curvRaw2<0.033') 
builder.tree = correctDataSet(builder.tree,True,False,True,False)
unfolding=AnalyticalUnfolding(builder.tree,50,1./100.,0.033)
builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('ZData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'Z MC1'

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5&&massErrRaw1<1&&massErrRaw2<1&&curvRaw1<0.033&&curvRaw2<0.033') 
builder.tree = correctDataSet(builder.tree,False,False,True,False)
unfolding=AnalyticalUnfolding(builder.tree,50,1./100.,0.033)
builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('ZMC_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()

print 'Z MC2'

builder = DataSetBuilder(pmap,w,'../../data/ZMC2.root','data',10000000)
builder.tree = builder.tree.reduce('abs(eta)<5&&massErrRaw1<1&&massErrRaw2<1&&curvRaw1<0.033&&curvRaw2<0.033') 
builder.tree = correctDataSet(builder.tree,False,False,True,False)
unfolding=AnalyticalUnfolding(builder.tree,50,1./100.,0.033)
builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('ZMC2_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()



