#include "CMGTools/KalmanCalibration/interface/KalmanCalibrator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <math.h>


KalmanCalibrator::KalmanCalibrator(bool isData) {

  random_ = new TRandom3(10101982);

  if (isData) {
    edm::FileInPath path("CMGTools/KalmanCalibration/data/dataInputs.root");
    file_ = new TFile(path.fullPath().c_str());
    magnetic = (TH2F*)file_->Get("magnetic");
  }
  else {
    edm::FileInPath path("CMGTools/KalmanCalibration/data/mcInputs.root");
    file_ = new TFile(path.fullPath().c_str());
  }
  isData_ = isData;


  scale_A1 =(TH3F*)file_->Get("A1"); 
  scale_A2 =(TH3F*)file_->Get("A2"); 
  scale_e = (TH3F*)file_->Get("e") ; 
  scale_B = (TH3F*)file_->Get("B") ; 


  shifted_A1 =(TH3F*)file_->Get("A1"); 
  shifted_A2 =(TH3F*)file_->Get("A2"); 
  shifted_e = (TH3F*)file_->Get("e") ; 
  shifted_B = (TH3F*)file_->Get("B") ; 

  sigma_A_target = (TH3F*)file_->Get("sigma_A_target"); 
  sigma_B_target = (TH3F*)file_->Get("sigma_B_target"); 
  sigma_C_target = (TH3F*)file_->Get("sigma_C_target"); 

  sigma_A_ref = (TH3F*)file_->Get("sigma_A_ref"); 
  sigma_B_ref = (TH3F*)file_->Get("sigma_B_ref"); 
  sigma_C_ref = (TH3F*)file_->Get("sigma_C_ref"); 

  ebe_A = (TH3F*)file_->Get("ebe_A"); 
  ebe_B = (TH3F*)file_->Get("ebe_B"); 
  ebe_C = (TH3F*)file_->Get("ebe_C"); 

 }


void KalmanCalibrator::reset() {
  resetHisto(shifted_A1,scale_A1);
  resetHisto(shifted_A2,scale_A2);
  resetHisto(shifted_e,scale_e);
  resetHisto(shifted_B,scale_B);

}

void KalmanCalibrator::randomize() {
  reset();
  randomizeHisto(shifted_A1);
  randomizeHisto(shifted_A2);
  randomizeHisto(shifted_e);
  randomizeHisto(shifted_B);
}


void KalmanCalibrator::randomizeHisto(TH1* histo) {
  for (int i=0;i<histo->GetNbinsX()+1;++i)
    for (int j=0;j<histo->GetNbinsY()+1;++j)
      for (int k=0;k<histo->GetNbinsZ()+1;++k) {
	int bin = histo->GetBin(i,j,k);
	histo->SetBinContent(bin,random_->Gaus(histo->GetBinContent(bin),histo->GetBinError(bin)));
      }
}

void KalmanCalibrator::resetHisto(TH1* histo,const TH1* ref) {
  for (int i=0;i<histo->GetNbinsX()+1;++i)
    for (int j=0;j<histo->GetNbinsY()+1;++j)
      for (int k=0;k<histo->GetNbinsZ()+1;++k) {
	int bin = histo->GetBin(i,j,k);
	histo->SetBinContent(bin,ref->GetBinContent(bin));
      }
}

KalmanCalibrator::~KalmanCalibrator() {

  if(shifted_A1) delete shifted_A1;
  if(shifted_A2) delete shifted_A2;
  if(shifted_e) delete shifted_e;
  if(shifted_B) delete shifted_B;
  file_->Close();
}


double KalmanCalibrator::getCorrectedPt(double pt,double eta,double phi,int charge) {
    double magneticMapFactor=1.0;
    if (isData_)
      magneticMapFactor = magnetic->GetBinContent(magnetic->GetBin(
								   magnetic->GetXaxis()->FindBin(phi),
								   magnetic->GetYaxis()->FindBin(eta)));
  double curvature = magneticMapFactor/pt;
  double sinTheta  = sin(2*atan(exp(-eta))); 
  double e = shifted_e->GetBinContent(scale_e->GetBin(1,
						    scale_e->GetYaxis()->FindBin(eta),
						    scale_e->GetZaxis()->FindBin(phi)));
  double A1 = shifted_A1->GetBinContent(scale_A1->GetBin(1,1,1));
  double A2 = shifted_A2->GetBinContent(scale_A2->GetBin(1,1,1));
  double B = shifted_B->GetBinContent(scale_B->GetBin(1,
						    scale_B->GetYaxis()->FindBin(eta),
						    scale_B->GetZaxis()->FindBin(phi)));
  curvature = (A2*eta*eta+A1-1)*curvature +curvature/(1+sinTheta*e*curvature) +charge*B;
  return 1.0/curvature;
}
double KalmanCalibrator::getCorrectedPtMag(double pt,double eta,double phi) {
    double magneticMapFactor=1.0;
    if (isData_)
      magneticMapFactor = magnetic->GetBinContent(magnetic->GetBin(
								   magnetic->GetXaxis()->FindBin(phi),
								   magnetic->GetYaxis()->FindBin(eta)));
  double curvature = magneticMapFactor/pt;
  return 1.0/curvature;
}



double KalmanCalibrator::getCorrectedError(double pt,double eta,double error) {
  double eA  = ebe_A->GetBinContent(ebe_A->GetBin(1,ebe_A->GetYaxis()->FindBin(eta),1));
  double eB  = ebe_B->GetBinContent(ebe_B->GetBin(1,ebe_B->GetYaxis()->FindBin(eta),1));
  double eC  = ebe_C->GetBinContent(ebe_C->GetBin(1,ebe_C->GetYaxis()->FindBin(eta),1));

  double diff = error*error+eA+eB*pt+eC*pt*pt;
 
  if (diff<0.0) {
    printf("Negative square root !!Ignoring\n");
    return error;
  }
  return sqrt(diff);
}



double KalmanCalibrator::smearGEN(double pt,double eta) {
  double A  = sigma_A_target->GetBinContent(sigma_A_target->GetBin(1,sigma_A_target->GetYaxis()->FindBin(eta),1));
  double B  = sigma_B_target->GetBinContent(sigma_B_target->GetBin(1,sigma_B_target->GetYaxis()->FindBin(eta),1));
  double C  = sigma_C_target->GetBinContent(sigma_C_target->GetBin(1,sigma_C_target->GetYaxis()->FindBin(eta),1));
  
  double resolution = sqrt( 4*A +4*B*pt+4*C*pt*pt )*pt;

  return pt+random_->Gaus(0.0,resolution);
}

double KalmanCalibrator::smear(double pt,double eta,bool reverse) {
  double Aref  = sigma_A_ref->GetBinContent(sigma_A_ref->GetBin(1,sigma_A_ref->GetYaxis()->FindBin(eta),1));
  double Bref  = sigma_B_ref->GetBinContent(sigma_B_ref->GetBin(1,sigma_B_ref->GetYaxis()->FindBin(eta),1));
  double Cref  = sigma_C_ref->GetBinContent(sigma_C_ref->GetBin(1,sigma_C_ref->GetYaxis()->FindBin(eta),1));
  
  double resolutionRef2 =  4*Aref +4*Bref*pt+4*Cref*pt*pt ;

  double Atarget  = sigma_A_target->GetBinContent(sigma_A_target->GetBin(1,sigma_A_target->GetYaxis()->FindBin(eta),1));
  double Btarget  = sigma_B_target->GetBinContent(sigma_B_target->GetBin(1,sigma_B_target->GetYaxis()->FindBin(eta),1));
  double Ctarget  = sigma_C_target->GetBinContent(sigma_C_target->GetBin(1,sigma_C_target->GetYaxis()->FindBin(eta),1));

  double resolutionTarget2 = 4*Atarget +4*Btarget*pt+4*Ctarget*pt*pt;


  if (resolutionTarget2<=resolutionRef2) {
    //assume we smear MC
    double factor = sqrt(resolutionRef2-resolutionTarget2)*pt;
    if(!reverse)
      return pt+random_->Gaus(0,factor);
    else
      return pt;
  }
  else {
    //assume we smear Data
    printf("You are smearing something strange!.. \n");
    double factor = sqrt(resolutionTarget2-resolutionRef2)*pt;
    if(reverse)
      return pt+random_->Gaus(0,factor);
    else
      return pt;

  }
}
