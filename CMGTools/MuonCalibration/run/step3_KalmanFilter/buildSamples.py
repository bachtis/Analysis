from defs import *

print 'This step is needed only for data to apply the bfield map corrections'


print 'J /psi data'
builder = DataSetBuilder(pmap,w,'../../data/JDATAPruned.root','data',10000000)
builder.tree = builder.tree.reduce('massRaw>3.0&&massRaw<3.2') 
builder.tree = correctDataSet(builder.tree,True,False,True,True)
#unfolding=AnalyticalUnfolding(builder.tree,50,1./25,1./5.5)
#builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('JData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()

print 'J /psi MC'

builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.tree = builder.tree.reduce('massRaw>3&&massRaw<3.2') 
builder.tree = correctDataSet(builder.tree,False,False,True,True)
#unfolding=AnalyticalUnfolding(builder.tree,25,0.04,1./5.5)
#builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('JMC_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()



print 'Z data'

builder = DataSetBuilder(pmap,w,'../../data/ZDATA.root','data',10000000)
builder.tree = correctDataSet(builder.tree,True,False,True,True)
#unfolding=AnalyticalUnfolding(builder.tree,50,1./80.,1./25)
#builder.tree = unfolding.unfoldDataSet(builder.tree)
builder.tree = builder.tree.reduce('massRaw>85.0&&massRaw<95.') 
f2=ROOT.TFile('ZData_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()


print 'Z MC1'

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
builder.tree = builder.tree.reduce('massRaw>85.&&massRaw<95.') 
builder.tree = correctDataSet(builder.tree,False,False,True,True)
#unfolding=AnalyticalUnfolding(builder.tree,50,1./100.,0.033)
#builder.tree = unfolding.unfoldDataSet(builder.tree)
f2=ROOT.TFile('ZMC_Input.root','RECREATE')
f2.cd()
builder.tree.Write('data')
f2.Close()




#print 'Z GEN'

#builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',10000000)
#builder.tree = smearEbE2D(builder.tree,1.0,1.0,True)
#builder.tree = builder.tree.reduce('massRaw>85.0&&massRaw<95.&&1./curvRaw1<65.&&1./curvRaw2<65.') 
#f2=ROOT.TFile('ZGEN_Input.root','RECREATE')
#f2.cd()
#builder.tree.Write('data')
#f2.Close()


#print 'J /psi GEN'

#builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
#builder.tree = smearEbE2D(builder.tree,1.0,1.0,True)
#builder.tree = smearFlat(builder.tree,1.0,0.00133,True)
#builder.tree = builder.tree.reduce('massErrRaw1<0.03&&massErrRaw2<0.03&&curvRaw1>0.04&&curvRaw2>0.04') 
#builder.tree = correctDataSet(builder.tree,False,False,True,False)
#unfolding=AnalyticalUnfolding(builder.tree,25,0.04,1./5.5)
#builder.tree = unfolding.unfoldDataSet(builder.tree)
#f2=ROOT.TFile('JGEN_Input.root','RECREATE')
#f2.cd()
#builder.tree.Write('data')
#f2.Close()






