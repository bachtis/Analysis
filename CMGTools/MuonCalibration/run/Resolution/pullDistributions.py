from defs import *


def pull(data):
    pull = ROOT.TH1F("pull","pull",100,-7,7)
    for i in range(0,data.numEntries()):
        line=data.get(i)
        p1=1.0/line.find("curvRaw1").getVal()
        p2=1.0/line.find("curvRaw2").getVal()
        gp1=1.0/line.find("curvGenRaw1").getVal()
        gp2=1.0/line.find("curvGenRaw2").getVal()
        mass=line.find("massRaw").getVal()
        e1=2*line.find("massErrRaw1").getVal()*p1/mass
        e2=2*line.find("massErrRaw2").getVal()*p2/mass
        
        pull.Fill((p1-gp1)/e1)
        pull.Fill((p2-gp2)/e2)
    return  pull



def pullMass(data):
    pullMass = ROOT.TH1F("pullMass","pull",1000,-20,20)
    for i in range(0,data.numEntries()):
        line=data.get(i)
        mass=line.find("massRaw").getVal()
        e=line.find("massErrRaw").getVal()
        pullMass.Fill((mass-3.09563)/e)
    return  pullMass


builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0&&curvGenRaw1>0&&curvGenRaw2>0")
builder.tree = correctDataSet(builder.tree,False,False,True,False,True)
p=pull(builder.tree)
pm=pullMass(builder.tree)


#builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
#builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
#builder.tree = correctDataSet(builder.tree,False,True,True,False,True)
#pZ=pull(builder.tree)
#p=pull(builder.tree)





    


