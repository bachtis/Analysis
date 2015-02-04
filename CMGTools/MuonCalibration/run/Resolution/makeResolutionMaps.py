import ROOT

from defs import *


def fillResolution(histo,data):
    for i in range(0,data.numEntries()):
        line=data.get(i)
        pt1=1./line.find('curvRaw1').getVal()
        pt2=1./line.find('curvRaw2').getVal()
        eta1=line.find('etaRaw1').getVal()
        eta2=line.find('etaRaw2').getVal()
        err1=2*line.find('massErrRaw1').getVal()/line.find('massRaw').getVal()
        err2=2*line.find('massErrRaw2').getVal()/line.find('massRaw').getVal()

        histo.Fill(pt1,eta1,err1)
        histo.Fill(pt2,eta2,err2)


def subtract(data,mc,diff):
    for i in range(0,diff.GetNbinsX()+1):
        for j in range(0,diff.GetNbinsY()+1):
            bin=diff.GetBin(i,j)
            c1=data.GetBinContent(bin)
            c2=mc.GetBinContent(bin)
            e1=data.GetBinError(bin)
            e2=mc.GetBinError(bin)
            
            d=c1*c1-c2*c2
            e=math.sqrt(4*c1*c1*e1*e1+4*c2*c2*e2*e2)
            if d==0:
                diff.SetBinContent(bin,0.0)
            else:    
                diff.SetBinContent(bin,(d/abs(d))*math.sqrt(abs(d)))
            diff.SetBinError(bin,e)

ptArr=[5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.,55.,60.,80.,100.]
etaArr=[-1.1,-0.9,-0.7,-0.5,-0.3,-0.15,0,0.15,0.3,0.5,0.7,0.9,1.1]


from array import array


mcMap = ROOT.TProfile2D("mcMap","mcMap",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr),"")
dataMap = ROOT.TProfile2D("dataMap","dataMap",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr),"")
diff = ROOT.TH2D("diff","diff",len(ptArr)-1,array('d',ptArr),len(etaArr)-1,array('d',etaArr))

builder = DataSetBuilder(pmap,w,'../../data/ZMC1.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
builder.tree = correctDataSet(builder.tree,False,False,True)

fillResolution(mcMap,builder.tree)


builder = DataSetBuilder(pmap,w,'../../data/JMC.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
builder.tree = correctDataSet(builder.tree,False,False,True)


fillResolution(mcMap,builder.tree)


builder = DataSetBuilder(pmap,w,'../../data/ZDATA.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
builder.tree = correctDataSet(builder.tree,True,False,True)

fillResolution(dataMap,builder.tree)


builder = DataSetBuilder(pmap,w,'../../data/JDATA.root','data',10000000)
builder.tree = builder.tree.reduce("massErrRaw1>0&&massErrRaw2>0")
builder.tree = correctDataSet(builder.tree,True,False,True)
fillResolution(dataMap,builder.tree)
subtract(dataMap,mcMap,diff)




f=ROOT.TFile('resolutionMaps.root','RECREATE')
f.cd()
mcMap.Write()
dataMap.Write()
diff.Write()
f.Close()


