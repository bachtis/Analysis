from defs import *



builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',10000000,False,True)
builder.tree=builder.tree.reduce('etaRaw1>-1.1&&etaRaw1<1.1&&etaRaw2>-1.1&&etaRaw2<1.1&&eta>-5&&eta<5')

builder.buildVsDiLepton()
builder.save("ZGEN_InputSmeared.root")


builder = DataSetBuilder(pmap,w,'../../data/ZGEN.root','data',10000000,False,False)
builder.tree=builder.tree.reduce('etaRaw1>-1.1&&etaRaw1<1.1&&etaRaw2>-1.1&&etaRaw2<1.1&&eta>-5&&eta<5')

builder.buildVsDiLepton()
builder.save("ZGEN_Input.root")

