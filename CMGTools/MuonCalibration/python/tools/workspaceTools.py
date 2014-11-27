import math

def prepareWorkspace(w):
    w.factory('curvRaw1[0,1000]')
    w.factory('curvRaw2[0,1000]')
    w.factory('etaRaw1[-2.4,2.4]')
    w.factory('etaRaw2[-2.4,2.4]')
    w.factory('massErrRaw[1.0,0,3.]')
    w.factory('massErrRaw1[1.0,0.000001,10.]')
    w.factory('massErrRaw2[1.0,0.000001,10.]')
    w.factory('weight[1,0.0,1000000.]')
    w.factory('phiRaw1[-'+str(math.pi)+','+str(math.pi)+ ']')
    w.factory('phiRaw2[-'+str(math.pi)+','+str(math.pi)+ ']')
    w.factory('muMass[0.1056583715]')
    w.factory('massRaw[0,1000]')



