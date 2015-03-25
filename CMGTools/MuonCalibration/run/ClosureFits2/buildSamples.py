
from defs import *
pmapZ = PartitionMap(curvArrZ,etaArr,phiArr,"")
pmapJ = PartitionMap(curvArrJ,etaArr,phiArr,"")
pmapY = PartitionMap(curvArrY,etaArr,phiArr,"")


print 'J /psi data'
builder = DataSetBuilder(pmapJ,w,'../../data/JDATAVertex.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>3.0&&(1.0/curvRaw1)<15.0&&(1.0/curvRaw2)>3.0&&(1.0/curvRaw2)<15.0")
builder.tree = correctDataSet(builder.tree,True,True,True)
builder.buildPtAbsEta()
builder.save("JData.root")


print 'J /psi MC'
builder = DataSetBuilder(pmapJ,w,'../../data/JMCVertex.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>3.0&&(1.0/curvRaw1)<15.0&&(1.0/curvRaw2)>3.0&&(1.0/curvRaw2)<15.0")
builder.tree=smearRelative(builder.tree,False)
builder.tree = correctDataSet(builder.tree,False,True,True)
builder.buildPtAbsEta()
builder.save("JMC.root")

print 'Y data'
builder = DataSetBuilder(pmapY,w,'../../data/YDATA.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>3.0&&(1.0/curvRaw1)<20.0&&(1.0/curvRaw2)>3.0&&(1.0/curvRaw2)<20.0")
builder.tree = correctDataSet(builder.tree,True,True,True)
builder.buildPtAbsEta()
builder.save("YData.root")


print 'Y MC'
builder = DataSetBuilder(pmapY,w,'../../data/YMC.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>3.0&&(1.0/curvRaw1)<20.0&&(1.0/curvRaw2)>3.0&&(1.0/curvRaw2)<20.0")
builder.tree=smearRelative(builder.tree,False)
builder.tree = correctDataSet(builder.tree,False,True,True)
builder.buildPtAbsEta()
builder.save("YMC.root")


print 'Z data'
builder = DataSetBuilder(pmapZ,w,'../../data/ZDATA.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>20.0&&(1.0/curvRaw1)<100.0&&(1.0/curvRaw2)>20.0&&(1.0/curvRaw2)<100.0")
builder.tree = correctDataSet(builder.tree,True,True,True)
#builder.tree=smearRelative(builder.tree,False)
builder.buildPtAbsEta()
builder.save("ZData.root")

print 'Z MC'
builder = DataSetBuilder(pmapZ,w,'../../data/ZMC1.root','data',10000000,True)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&(1.0/curvRaw1)>20.0&&(1.0/curvRaw1)<100.0&&(1.0/curvRaw2)>20.0&&(1.0/curvRaw2)<100.0")
builder.tree=smearRelative(builder.tree,False)
builder.tree = correctDataSet(builder.tree,False,True,True)
builder.buildPtAbsEta()
builder.save("ZMC.root")


