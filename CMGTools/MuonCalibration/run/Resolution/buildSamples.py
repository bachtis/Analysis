from defs import *
builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = correctDataSet(builder.tree,False,False,True)
builder.build(-1,2)
builder.save("Data_Input.root")


builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.build(-1,2)
builder.save("MC_Input.root")


