from defs import *
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = recalibrateMap(builder.tree)
builder.tree = builder.tree.reduce('1.0/curvRaw1>4.&& 1.0/curvRaw2>4.') 

builder.buildVsDiLepton()
builder.save("JDATA_Input.root")

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000,False,False)
builder.tree = builder.tree.reduce('1.0/curvRaw1>4.&& 1.0/curvRaw2>4.') 

builder.buildVsDiLepton()
builder.save("JGEN_Input.root")

