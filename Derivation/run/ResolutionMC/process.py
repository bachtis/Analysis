from defs import *
Era='13TEV_Rereco'

estimator=ROOT.MCEBEEstimator("results.root",hmap,"MC_Moriond17_13TeV")
estimator.processTree("/scratch3/Kalman/data/{era}/JPsiMC.root".format(era=Era),"pt1>3&&pt2>3&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0")
estimator.processTree("/scratch3/Kalman/data/{era}/JPsiMC20.root".format(era=Era),"pt1>3.0&&pt2>3.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0")
estimator.processTree("/scratch3/Kalman/data/{era}/ZMC.root".format(era=Era),"pt1>22.0&&pt2>22.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0")
estimator.processTree("/scratch3/Kalman/data/{era}/ZJETSMC.root".format(era=Era),"pt1>22.0&&pt2>22.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0")
estimator.close()



