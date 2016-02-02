import math

def prepareWorkspace(w):
    w.factory('c1[0,1]')
    w.factory('c2[0,1]')
    w.factory('N1[0,200]')
    w.factory('N2[0,200]')
    w.factory('eta1[-2.4,2.4]')
    w.factory('eta2[-2.4,2.4]')
    w.factory('massErr[0.01,0,3.]')
    w.factory('cErr1[0.01,0.000001,3.]')
    w.factory('cErr2[0.01,0.000001,3.]')
    w.factory('weight[1,0.0,1000000.]')
    w.factory('phi1[-'+str(math.pi)+','+str(math.pi)+ ']')
    w.factory('phi2[-'+str(math.pi)+','+str(math.pi)+ ']')
    w.factory('muMass[0.1056583715]')
    w.factory('mass[-1000,1000]')



