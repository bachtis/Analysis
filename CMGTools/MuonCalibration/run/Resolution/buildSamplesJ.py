from defs import *
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
#builder.tree = correctDataSet(builder.tree,True,False,True,False,True)
#builder.build(-1,2)
builder.build()
builder.save("JData_Input.root")


builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
#builder.tree = correctDataSet(builder.tree,False,False,True,False,True)
#builder.tree = correctDataSet(builder.tree,False,False,True,False,True)
builder.build()
builder.save("JMC_Input.root")



builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
#builder.tree = correctDataSet(builder.tree,False,False,True,False,True)
#builder.tree = correctDataSet(builder.tree,False,False,True,False,True)
builder.build()
builder.save("JGEN_Input.root")




