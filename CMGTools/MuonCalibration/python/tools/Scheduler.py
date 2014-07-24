import os
import shutil
import subprocess
import time
import ROOT
import socket

class Scheduler:
    def __init__(self,levels,filesPerJob = 1):
        self.levels=levels
        self.inputFiles={}
        self.filesPerJob = filesPerJob

        for i,level in enumerate(self.levels):
            subset='/'.join(self.levels[0:i+1])
            if not os.path.exists(subset):
                os.mkdir(subset)
        self.dir =os.getcwd()+'/'+ '/'.join(self.levels)
        
    def declareJob(self,file,bin):
        self.inputFiles[bin]=file
        shutil.move(file,self.dir+'/')

    def declareFile(self,file,bin):
        self.inputFiles[bin]=file

    def makePyScript(self,filename,index):
        pyScript = "import ROOT\n"
        pyScript +="ROOT.gSystem.Load('libCMGToolsMuonCalibration')\n"
        pyScript +="f=ROOT.TFile('"+filename+"','UPDATE')\n"
        pyScript +="w= f.Get('w')\n"
        pyScript+="""
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
fitter = LineshapeFitter(w)
fitter.load()
(hintResult,fitResult) = fitter.fit('model',w.data('data'),False,0)
f.cd()
#frame=w.var('massRaw').frame()
#w.data('data').plotOn(frame)
#w.pdf('model').plotOn(frame)
#c = ROOT.TCanvas('fit','fit')
#c.cd()
#frame.Draw()
#c.Write()
w.Write("",ROOT.TObject.kOverwrite)
f.Close()
"""
        return pyScript

    def makeExecScript(self,chunk):
        execScript = 'cd {cwd} \n'.format(cwd=self.dir)
        execScript += 'eval `scramv1 runtime -sh` \n'
        for item in chunk:
            execScript += "python job{p}.py >> log{p}.txt\n".format(p = item)
        return execScript


    def chunks(self,l,n):
        return [l[i:i+n] for i in range(0, len(l), n)]

    def submit(self,queue):
        cwd = os.getcwd()
        #first make all python files!
        for i,file in self.inputFiles.iteritems():
            pyScript = self.makePyScript(file,i)
            f2 = open(self.dir+'/job{p}.py'.format(p=i),'w')
            f2.write(pyScript)
            f2.close()
        #next make all exec files    
        for chunk in self.chunks(list(self.inputFiles),self.filesPerJob):
            f1 = open(self.dir+'/job{p}_{n}.sh'.format(p=chunk[0],n=chunk[-1]),'w')
            execScript = self.makeExecScript(chunk)
            f1.write(execScript)
            f1.close()
            os.chdir(self.dir)
            os.system('chmod +x job{p}_{n}.sh'.format(p=chunk[0],n=chunk[-1]))
            os.system('bsub -q {queue} -oo joblog{p}_{n}.txt job{p}_{n}.sh '.format(queue=queue,p=chunk[0],n=chunk[-1]))

            os.chdir(cwd)


    def submitLOCAL(self):
        cwd = os.getcwd()
        #first make all python files!
        for i,file in self.inputFiles.iteritems():
            pyScript = self.makePyScript(file,i)
            f2 = open(self.dir+'/job{p}.py'.format(p=i),'w')
            f2.write(pyScript)
            f2.close()
        #next make all exec files    
        for chunk in self.chunks(list(self.inputFiles),self.filesPerJob):
            f1 = open(self.dir+'/job{p}_{n}.sh'.format(p=chunk[0],n=chunk[-1]),'w')
            execScript = self.makeExecScript(chunk)
            f1.write(execScript)
            f1.close()
            os.chdir(self.dir)
            os.system('chmod +x job{p}_{n}.sh'.format(p=chunk[0],n=chunk[-1]))
            print 'RUNNING job' 
            os.system('./job{p}_{n}.sh '.format(p=chunk[0],n=chunk[-1]))

            os.chdir(cwd)

    def fakeSubmit(self):
        cwd = os.getcwd()
        #first make all python files!
        for i,file in self.inputFiles.iteritems():
            pyScript = self.makePyScript(file,i)
            f2 = open(self.dir+'/job{p}.py'.format(p=i),'w')
            f2.write(pyScript)
            f2.close()
        #next make all exec files    
        for chunk in self.chunks(list(self.inputFiles),self.filesPerJob):
            f1 = open(self.dir+'/job{p}_{n}.sh'.format(p=chunk[0],n=chunk[-1]),'w')
            execScript = self.makeExecScript(chunk)
            f1.write(execScript)
            f1.close()



    
    
    def wait(self):
        done = False

        while (not done):
            comOut=subprocess.check_output('bjobs',stderr=subprocess.STDOUT)
#            comOut=subprocess.check_output('ls',stderr=subprocess.STDOUT)

            if comOut.find(socket.gethostname().split('.')[0]) != -1:
                print 'Waiting 1 more min to see if jobs finish'
                time.sleep(30)
            else:
                done=True
        return    


    def harvest(self,bin,var = 'scale'):
        f = ROOT.TFile(self.dir+'/'+self.inputFiles[bin])
        w = f.Get('w')
        w.loadSnapshot("result")
        return ( w.var(var).getVal() , w.var(var).getError()) 
