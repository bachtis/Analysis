import ROOT
from array import array
from math import pi
ROOT.gSystem.Load("libKaMuCaDerivation")
pts=[3,4,5,6,8,9,10,12,14,16,18,20,25,30,35,40,45,50,55,60,70,80,90,100,120]
etas=[-2.5,-2.1,-1.7,-0.9,-0.5,-0.3,0.3,0.5,0.9,1.7,2.1,2.5]

res=[]

for i in range(0,101):
    res.append(0.8+0.4*i/100.0)


hmap = ROOT.TH3D("hmap","hmap",len(pts)-1,array('f',pts),len(etas)-1,array('f',etas),len(res)-1,array('f',res))

