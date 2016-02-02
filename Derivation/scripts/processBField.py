import ROOT
from array import array

runArr=[]
for run in range(246908,260629):
    runArr.append(run)


h=ROOT.TH1F("current","current",len(runArr)-1,array('f',runArr))

f=open('bfield.txt')
for line in f:
    comps=line.split(', ')
    print comps
    bin=h.GetXaxis().FindBin(float(comps[0]))
    print 3.8/float(comps[1])
    h.SetBinContent(bin,3.8/float(comps[1]))

f.close()


F=ROOT.TFile("current.root","RECREATE")
F.cd()
h.Write()
F.Close()
