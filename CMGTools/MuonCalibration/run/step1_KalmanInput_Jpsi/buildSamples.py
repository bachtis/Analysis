

from defs import *
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.buildRapidity()
builder.save("JDATA.root")

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.buildRapidity()
builder.save("JGEN.root")


builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.tree = smearAbsolute(builder.tree)
builder.buildRapidity()
builder.save("JGENDataSmear.root")

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.tree = correctDataSet(builder.tree,False,False,True)
builder.tree = smearEbE2D(builder.tree,1,1,True)
builder.buildRapidity()
builder.save("JGENMCSmear.root")

