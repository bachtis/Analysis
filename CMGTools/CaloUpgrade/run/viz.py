from CMGTools.CaloUpgrade.tools.Visualizer import Visualizer
import ROOT
ROOT.gROOT.ProcessLine(".x tdrstyle.C")
visualizer = Visualizer('visInput.root')
visualizer.process(0,0)

