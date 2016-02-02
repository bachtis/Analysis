import ROOT
f1=ROOT.TFile('map.root')
map1=f1.Get('map')
f2=ROOT.TFile('map3D.root')
map2=f2.Get('map')

map1.Sumw2()
map2.Sumw2()


f=ROOT.TFile('mapCalibration.root','RECREATE')
f.cd()
map1.Write('approximateMap')
map2.Write('preciseMap')

calib = map1.Clone()
calib.SetName('mapCalibration')
calib.Divide(map2)
calib.Write('mapCorrection')
f.Close()
f1.Close()
f2.Close()



