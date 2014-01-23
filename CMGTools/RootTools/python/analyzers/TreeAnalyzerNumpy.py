from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from ROOT import TFile,TTree
import numpy

class TreeAnalyzerNumpy( Analyzer ):
    """Base TreeAnalyzerNumpy, to create flat TTrees.

    Check out TestTreeAnalyzer for a concrete example.
    IMPORTANT: FOR NOW, CANNOT RUN SEVERAL TreeAnalyzers AT THE SAME TIME!
    Anyway, you want only one TTree, don't you?"""

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(TreeAnalyzerNumpy,self).__init__(cfg_ana, cfg_comp, looperName)
        fileName = '/'.join([self.dirName,
                             self.name+'_tree.root'])
        self.file = TFile( fileName, 'recreate' )
        self.tree = TTree('tree',self.name)
        self.vars={}

        self.declareVariables()


    def branch(self, varName,type=float ):
        
        self.vars[varName]=numpy.zeros(1,type)
        if type is float  : 
            self.tree.Branch(varName,self.vars[varName],varName+'/D')
        elif type is int    : 
            self.tree.Branch(varName,self.vars[varName],varName+'/I')


    def set(self, varName, value ):
        self.vars[varName][0]=value

    def fill(self):
        self.tree.Fill()
        
    def declareVariables(self):
        print 'TreeAnalyzerNumpy.declareVariables : overload this function.'
        pass

    def reset(self):
        for name,value in self.vars.iteritems():
            value[0]=-99


    def write(self):
        super(TreeAnalyzerNumpy, self).write()
        self.file.Write() 

