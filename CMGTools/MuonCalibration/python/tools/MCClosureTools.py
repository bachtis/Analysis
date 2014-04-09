import ROOT
import numpy
import math

class MCClosureTools (object):


    def loop(self,sign,data):
        h = ROOT.TH1F("res","res",500,-0.003,0.003)
        
        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            if sign == 'pos':
                curv = 'curvRaw1'
                gencurv = 'curvGenRaw1'
            else:    
                curv = 'curvRaw2'
                gencurv = 'curvGenRaw2'
                
            h.Fill(line.find(gencurv).getVal()-line.find(curv).getVal())
        return h
    

    def getMeanAndSpreadArithmetic(self,sign,data):
        h = self.loop(sign,data)
        hist=h
        return h.GetMean(),h.GetMeanError(),h.GetRMS(),h.GetRMSError()

    def getMeanAndSpreadFit(self,sign,data):
        h= self.loop(sign,data)

        w=ROOT.RooWorkspace("w","w")
        w.factory("res[-0.003,0.003]")
        datahist = ROOT.RooDataHist("histD","histD",ROOT.RooArgList(w.var("res")),h)
        w.factory("RooGaussian::gaus(res,mean[-0.003,0.003],sigma[0,0.01])")
        w.pdf("gaus").fitTo(datahist)
        
        return w.var('mean').getVal(),w.var('mean').getError(),w.var('sigma').getVal(),w.var('sigma').getError()
    
    
    


