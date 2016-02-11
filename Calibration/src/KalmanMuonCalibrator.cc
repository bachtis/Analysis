#include "../interface/KalmanMuonCalibrator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <math.h>
KalmanMuonCalibrator::KalmanMuonCalibrator() {
}

KalmanMuonCalibrator::KalmanMuonCalibrator(const std::string& filename) {

  random_ = new TRandom3(10101982);


  edm::FileInPath path("KaMuCa/Calibration/data/"+filename+".root");
  file_ = new TFile(path.fullPath().c_str());
  isData_ = (filename.find("DATA")!=std::string::npos);


  //inputs


  //magnetic map correction from Solenoid mapping 
  magnetic = (TH2F*)file_->Get("magnetic");


  //Magnetic correction :Up vs downstairs
  scale_A1 =(TH3F*)file_->Get("A1"); 
  scale_A2 =(TH3F*)file_->Get("A2"); 


  //De/Dx
  scale_e = (TH3F*)file_->Get("e") ; 


  //Misalignment fourier terms
  scale_B0 = (TH3F*)file_->Get("B0") ; 
  scale_B1 = (TH3F*)file_->Get("B1") ; 
  scale_B2 = (TH3F*)file_->Get("B2") ; 
  scale_C1 = (TH3F*)file_->Get("C1") ; 
  scale_C2 = (TH3F*)file_->Get("C2") ; 





  //Shifted versions for systematic errors
  shifted_A1 =(TH3F*)scale_A1->Clone();
  shifted_A1->SetName("shifted_A1");
  shifted_A2 =(TH3F*)scale_A2->Clone();
  shifted_A2->SetName("shifted_A2");
  shifted_e = (TH3F*)scale_e->Clone();
  shifted_e->SetName("shifted_e");
  shifted_B0 = (TH3F*)scale_B0->Clone();
  shifted_B0->SetName("shifted_B0"); 
  shifted_B1 = (TH3F*)scale_B1->Clone();
  shifted_B1->SetName("shifted_B1"); 
  shifted_B2 = (TH3F*)scale_B2->Clone();
  shifted_B2->SetName("shifted_B2"); 
  shifted_C1 = (TH3F*)scale_C1->Clone();
  shifted_C1->SetName("shifted_C1"); 
  shifted_C2 = (TH3F*)scale_C2->Clone();
  shifted_C2->SetName("shifted_C2"); 



  //Resolution Histograms
  aSRC_ = (TH1D*)file_->Get("aSRC");
  bSRC_ = (TH1D*)file_->Get("bSRC");
  cSRC_ = (TH1D*)file_->Get("cSRC");
  dSRC_ = (TH1D*)file_->Get("dSRC");

  aTARGET_ = (TH1D*)file_->Get("aTARGET");
  bTARGET_ = (TH1D*)file_->Get("bTARGET");
  cTARGET_ = (TH1D*)file_->Get("cTARGET");
  dTARGET_ = (TH1D*)file_->Get("dTARGET");

  aEBE_ = (TH1D*)file_->Get("aEBE");
  bEBE_ = (TH1D*)file_->Get("bEBE");
  cEBE_ = (TH1D*)file_->Get("cEBE");
  dEBE_ = (TH1D*)file_->Get("dEBE");




  //Closure plot
  closure_ = (TH3F*)file_->Get("closure"); 

  //Mapping between covariance Matrix and correction factors
  covHistoMap_ = (TH1I*)file_->Get("covHistoMap");
  covBinMap_ = (TH1I*)file_->Get("covBinMap");

  //Eigen vectors and eigen values
  eigenvalues_ = (TVectorD*)file_->Get("eigenvalues");
  eigenvectors_ = (TMatrixD*)file_->Get("eigenvectors");

  varyClosure_=0;
 }


void KalmanMuonCalibrator::reset() {
  varyClosure_=0;
  resetHisto(shifted_A1,scale_A1);
  resetHisto(shifted_A2,scale_A2);
  resetHisto(shifted_e,scale_e);
  resetHisto(shifted_B0,scale_B0);
  resetHisto(shifted_B1,scale_B2);
  resetHisto(shifted_B2,scale_B2);
  resetHisto(shifted_C1,scale_C1);
  resetHisto(shifted_C2,scale_C2);
}


void KalmanMuonCalibrator::resetHisto(TH1* histo,const TH1* ref) {
  for (int i=0;i<histo->GetNbinsX()+1;++i)
    for (int j=0;j<histo->GetNbinsY()+1;++j)
      for (int k=0;k<histo->GetNbinsZ()+1;++k) {
	int bin = histo->GetBin(i,j,k);
	histo->SetBinContent(bin,ref->GetBinContent(bin));
      }
}


KalmanMuonCalibrator::~KalmanMuonCalibrator() {
}


double KalmanMuonCalibrator::closure(double pt,double eta) {
  Int_t bin = closure_->GetBin(
			       closure_->GetXaxis()->FindBin(pt),	       
			       closure_->GetYaxis()->FindBin(fabs(eta)),
			       1);
  return closure_->GetBinContent(bin)-1.0;
}


double KalmanMuonCalibrator::getCorrectedPt(double pt,double eta,double phi,int charge) {
    double magneticMapFactor=1.0;
    if (isData_)
      magneticMapFactor = magnetic->GetBinContent(magnetic->GetBin(
								   magnetic->GetXaxis()->FindBin(phi),
								   magnetic->GetYaxis()->FindBin(eta)));
    double curvature = (magneticMapFactor)/pt;
    double e = shifted_e->GetBinContent(scale_e->GetBin(1,scale_e->GetYaxis()->FindBin(eta),1));
    double sinTheta  = sin(2*atan(exp(-eta))); 

    double A1 = shifted_A1->GetBinContent(13);
    //    double A1 = shifted_A1->GetBinContent(scale_A1->GetBin(1,scale_A1->GetYaxis()->FindBin(eta),scale_A1->GetZaxis()->FindBin(phi)));
    double A2 = shifted_A2->GetBinContent(13);


    double B0 = shifted_B0->GetBinContent(scale_B0->GetBin(1,scale_B0->GetYaxis()->FindBin(eta),scale_B0->GetZaxis()->FindBin(phi)));
    double B1 = shifted_B1->GetBinContent(scale_B1->GetBin(1,scale_B1->GetYaxis()->FindBin(eta),scale_B1->GetZaxis()->FindBin(phi)));
    double B2 = shifted_B2->GetBinContent(scale_B2->GetBin(1,scale_B2->GetYaxis()->FindBin(eta),scale_B2->GetZaxis()->FindBin(phi)));
    double C1 = shifted_C1->GetBinContent(scale_C1->GetBin(1,scale_C1->GetYaxis()->FindBin(eta),scale_C1->GetZaxis()->FindBin(phi)));
    double C2 = shifted_C2->GetBinContent(scale_C2->GetBin(1,scale_C2->GetYaxis()->FindBin(eta),scale_C2->GetZaxis()->FindBin(phi)));
    
    double B = B0+B1*sin(phi)+B2*sin(2*phi)+C1*cos(phi)+C2*cos(2*phi);

    double tanTheta = 2.0/(exp(eta)-exp(-eta));

    double mag=A1+A2/(tanTheta*tanTheta);
    //    double mag=A1+A2*eta*eta;
   curvature = mag*curvature -e*sinTheta*curvature*curvature+charge*B;
    return (1.0/curvature)*(1.0+varyClosure_*closure(pt,eta));
}
double KalmanMuonCalibrator::getCorrectedPtMag(double pt,double eta,double phi) {
    double magneticMapFactor=1.0;
       if (isData_)
         magneticMapFactor = magnetic->GetBinContent(magnetic->GetBin(
								      magnetic->GetXaxis()->FindBin(phi),
								      magnetic->GetYaxis()->FindBin(eta)));
    double curvature = (magneticMapFactor)/pt;
  return 1.0/curvature;
}



double KalmanMuonCalibrator::getCorrectedError(double pt,double eta,double error) {
  Int_t bin = aEBE_->GetXaxis()->FindBin(eta);
  
  double a2 = aEBE_->GetBinContent(bin);
  double b2 = bEBE_->GetBinContent(bin);
  double c2 = cEBE_->GetBinContent(bin);
  double d2 = dEBE_->GetBinContent(bin);


  double aSRC = aSRC_->GetBinContent(bin);
  double bSRC = bSRC_->GetBinContent(bin);
  double cSRC = cSRC_->GetBinContent(bin);
  double dSRC = dSRC_->GetBinContent(bin);

  double pt2=pt*pt;
  
  //new ebe^2 = ebe^2 + sigma^2-ebe_avg^2

  double error2=error*error + aSRC+bSRC/(1+dSRC/pt2)+cSRC*pt2 -(a2+b2/(1+d2/pt2)+c2*pt2);

  if (error2<0) {
    //    printf("Got Negative EbE !! Will ignore and not correct the ebe of this muon pt=%f , eta=%f ,error=%f ,residual2=%f\n",pt,eta,error,error2-error*error); 
    return error;
  }

  return sqrt(error2);
}


double KalmanMuonCalibrator::getCorrectedErrorAfterSmearing(double pt,double eta,double error) {
  Int_t bin = aEBE_->GetXaxis()->FindBin(eta);
  
  double a2 = aEBE_->GetBinContent(bin);
  double b2 = bEBE_->GetBinContent(bin);
  double c2 = cEBE_->GetBinContent(bin);
  double d2 = dEBE_->GetBinContent(bin);
  double aTARGET = aTARGET_->GetBinContent(bin);
  double bTARGET = bTARGET_->GetBinContent(bin);
  double cTARGET = cTARGET_->GetBinContent(bin);
  double dTARGET = dTARGET_->GetBinContent(bin);

  double pt2=pt*pt;
  
  //new ebe^2 = ebe^2 + sigma^2-ebe_avg^2

  double error2=error*error + aTARGET+bTARGET/(1+dTARGET/pt2)+cTARGET*pt2 -(a2+b2/(1+d2/pt2)+c2*pt2);

  if (error2<0) {
    //    printf("Got Negative EbE !! Will ignore and not correct the ebe of this muon pt=%f , eta=%f ,error=%f ,residual2=%f\n",pt,eta,error,error2-error*error); 
    return error;
  }

  return sqrt(error2);
}




int KalmanMuonCalibrator::getN() {
  return eigenvalues_->GetNoElements();
}


void KalmanMuonCalibrator::varyClosure(int sigmas) {
    varyClosure_ = sigmas;

}

void KalmanMuonCalibrator::vary(int ii,int sigmas) {
  reset();
  int N = getN();

  TVectorD * v = new TVectorD(N);
  v->Zero();

  if (ii>N) {
    printf("Hey you are trying to vary outside the number of elements,which is pretty stupid \n");
    printf("You ask for element %d out of %d\n",ii,N);   
  }

  (*v)[ii] = sqrt((*eigenvalues_)[ii])*sigmas;


  TVectorD correlated = (*eigenvectors_)*(*v);

  for (int i=0;i<N;++i) {
    int histo = covHistoMap_->GetBinContent(i+1);
    int bin = covBinMap_->GetBinContent(i+1);
    float value = correlated[i];
    switch (histo) {
    case 1:
      shifted_A1->SetBinContent(bin,scale_A1->GetBinContent(bin)+value);
      break;
    case 2:
      shifted_A2->SetBinContent(bin,scale_A2->GetBinContent(bin)+value);
      break;
    case 3:
      shifted_e->SetBinContent(bin,scale_e->GetBinContent(bin)+value);
      break;

    case 4:
      shifted_B0->SetBinContent(bin,scale_B0->GetBinContent(bin)+value);
      break;
    case 5:
      shifted_B1->SetBinContent(bin,scale_B1->GetBinContent(bin)+value);
      break;
    case 6:
      shifted_B2->SetBinContent(bin,scale_B2->GetBinContent(bin)+value);
      break;
    case 7:
      shifted_C1->SetBinContent(bin,scale_C1->GetBinContent(bin)+value);
      break;
    case 8:
      shifted_C2->SetBinContent(bin,scale_C2->GetBinContent(bin)+value);
      break;

    default:
      printf("UNKNOWN HISTO-> problem while varying (That is important\n");
    }

  }

  delete v;

}


double KalmanMuonCalibrator::smear(double pt,double eta) {
  Int_t bin = aSRC_->GetXaxis()->FindBin(eta);

  double aSRC = aSRC_->GetBinContent(bin);
  double bSRC = bSRC_->GetBinContent(bin);
  double cSRC = cSRC_->GetBinContent(bin);
  double dSRC = dSRC_->GetBinContent(bin);

  double aTARGET = aTARGET_->GetBinContent(bin);
  double bTARGET = bTARGET_->GetBinContent(bin);
  double cTARGET = cTARGET_->GetBinContent(bin);
  double dTARGET = dTARGET_->GetBinContent(bin);


  double resolutionSRC2 = aSRC+bSRC/(1+dSRC/(pt*pt))+cSRC*pt*pt;
  double resolutionTARGET2 = aTARGET+bTARGET/(1+dTARGET/(pt*pt))+cTARGET*pt*pt;

  //  printf ("Resolution SRC=%f,Resolution TARGET=%f a=%f b=%f c=%f d=%f ratio=%f\n",sqrt(resolutionSRC2),sqrt(resolutionTARGET2),sqrt(fabs(aTARGET-aSRC)),sqrt(fabs(bTARGET/(1+dTARGET/(pt*pt))-bSRC/(1+dSRC/(pt*pt)))),sqrt(fabs(cTARGET*pt*pt-cSRC*pt*pt)),sqrt(fabs(dTARGET-dSRC)),sqrt(resolutionTARGET2/resolutionSRC2));
  Double_t factor = resolutionTARGET2-resolutionSRC2;
  if (factor<0) {
    //    printf("target has better resolution than source-not smearing pt=%f,eta=%f RSRC=%f RTGR=%f\n",pt,eta,resolutionTARGET2,resolutionSRC2);
    return pt;
  }
  factor = sqrt(fabs(factor));
  float smeared = random_->Gaus(1.0,factor)/pt;
  return 1.0/smeared;
}


double KalmanMuonCalibrator::smearForSync(double pt,double eta) {
  Int_t bin = aSRC_->GetXaxis()->FindBin(eta);

  double aSRC = aSRC_->GetBinContent(bin);
  double bSRC = bSRC_->GetBinContent(bin);
  double cSRC = cSRC_->GetBinContent(bin);
  double dSRC = dSRC_->GetBinContent(bin);

  double aTARGET = aTARGET_->GetBinContent(bin);
  double bTARGET = bTARGET_->GetBinContent(bin);
  double cTARGET = cTARGET_->GetBinContent(bin);
  double dTARGET = dTARGET_->GetBinContent(bin);


  double resolutionSRC2 = aSRC+bSRC/(1+dSRC/(pt*pt))+cSRC*pt*pt;
  double resolutionTARGET2 = aTARGET+bTARGET/(1+dTARGET/(pt*pt))+cTARGET*pt*pt;
  Double_t factor = resolutionTARGET2-resolutionSRC2;
  if (factor<0) {
    //    printf("target has better resolution than source-not smearing pt=%f,eta=%f RSRC=%f RTGR=%f\n",pt,eta,resolutionTARGET2,resolutionSRC2);
    return pt;
  }
  factor = sqrt(fabs(factor))/pt;
  float smeared = 1.0/pt+factor;
  return 1.0/smeared;
}


double KalmanMuonCalibrator::smearUsingEbE(double pt,double eta,double relErr) {
  Int_t bin = aSRC_->GetXaxis()->FindBin(eta);

  double aSRC = aSRC_->GetBinContent(bin);
  double bSRC = bSRC_->GetBinContent(bin);
  double cSRC = cSRC_->GetBinContent(bin);
  double dSRC = dSRC_->GetBinContent(bin);

  double aEBE = aEBE_->GetBinContent(bin);
  double bEBE = bEBE_->GetBinContent(bin);
  double cEBE = cEBE_->GetBinContent(bin);
  double dEBE = dEBE_->GetBinContent(bin);


  double aTARGET = aTARGET_->GetBinContent(bin);
  double bTARGET = bTARGET_->GetBinContent(bin);
  double cTARGET = cTARGET_->GetBinContent(bin);
  double dTARGET = dTARGET_->GetBinContent(bin);


  double resolutionSRC2 = aSRC+bSRC/(1+dSRC/(pt*pt))+cSRC*pt*pt;
  double resolutionTARGET2 = aTARGET+bTARGET/(1+dTARGET/(pt*pt))+cTARGET*pt*pt;
  double ebeAvg2 = aEBE+bEBE/(1+dEBE/(pt*pt))+cEBE*pt*pt;

  Double_t factor = resolutionTARGET2-resolutionSRC2+relErr*relErr-ebeAvg2;
  if (factor<=0) {
    //    printf("target has better resolution than source-not smearing pt=%f,eta=%f RSRC=%f RTGR=%f\n",pt,eta,resolutionTARGET2,resolutionSRC2);
    return pt;
  }
  factor = sqrt(fabs(factor))/pt;
  float smeared = random_->Gaus(1.0/pt,factor);
  if (1.0/smeared<0.0)
    return pt;
  return 1.0/smeared;
}
