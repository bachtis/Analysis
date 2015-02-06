from defs import *



#builder = DataSetBuilder(pmap,w,'../../data/ZDATA.root','data',10000000)
#builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
#builder.build()
#builder.save("ZData_Input.root")

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
builder.build()
builder.save("ZMC_Input.root")


#builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',10000000)
#builder.tree = smearEbE2D(builder.tree,1.0,1.0,False)
#builder.build()
#builder.save("ZGEN_Input.root")

