import ROOT
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
import os


def convert(muon1,muon2,file,tree,preselection):
        w=ROOT.RooWorkspace("w","w")
        f = ROOT.TFile(file)
        t = f.Get(tree)
        w = ROOT.RooWorkspace(w)
        prepareWorkspace(w)
        cache = ROOT.TFile('__cache__.root','RECREATE')
        tree = t.CopyTree(preselection)
        data = ROOT.RooDataSet("data","DATA",ROOT.RooArgSet(w.var('curvRaw1'),
                                                     w.var('etaRaw1'),
                                                     w.var('phiRaw1'),
                                                     w.var('curvRaw2'),
                                                     w.var('etaRaw2'),
                                                     w.var('phiRaw2'),
                                                     w.var('massRaw')
                                                     ))
        for event in tree:
            pt1 = getattr(event,muon1+'_pt')
            eta1 = getattr(event,muon1+'_eta')
            phi1 = getattr(event,muon1+'_phi')
            pt2 = getattr(event,muon2+'_pt')
            eta2 = getattr(event,muon2+'_eta')
            phi2 = getattr(event,muon2+'_phi')
            
            w.var('curvRaw1').setVal(1./pt1)
            w.var('etaRaw1').setVal(eta1)
            w.var('phiRaw1').setVal(phi1)
            w.var('curvRaw2').setVal(1./pt2)
            w.var('etaRaw2').setVal(eta2)
            w.var('phiRaw2').setVal(phi2)

            v1=ROOT.TLorentzVector()
            v2=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(pt1,eta1,phi1,w.var('muMass').getVal())
            v2.SetPtEtaPhiM(pt2,eta2,phi2,w.var('muMass').getVal())
            w.var('massRaw').setVal((v1+v2).M())
            data.add(ROOT.RooArgSet(w.var('curvRaw1'),
                                    w.var('etaRaw1'),
                                    w.var('phiRaw1'),
                                    w.var('curvRaw2'),
                                    w.var('etaRaw2'),
                                    w.var('phiRaw2'),
                                    w.var('massRaw')
                                    ))
        filename=file.split('.')
        newFileName=filename[0]+'_data.root'

        f2= ROOT.TFile(newFileName,"RECREATE")
        f2.cd()
        data.Write()
        f2.Close()
        cache.Close()
        f.Close()
        os.remove('__cache__.root')
