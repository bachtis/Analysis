import math
from KaMuCa.Derivation.tools.Partitioner import Partitioner
ptArr=[0,10000000000]
etaArr=[-2.5,-2.0,-1.7,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.7,2.0,2.5]
#etaArr=[-2.5,-2.0,-1.1,-0.9,0.0,0.9,1.1,2.0,2.5]


phiArr=[-math.pi,math.pi]
pmap = Partitioner(ptArr,etaArr,phiArr)







