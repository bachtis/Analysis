from defs import *
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.buildRapidity()
builder.save("JDATA.root")

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.buildRapidity()
builder.save("JGEN.root")

