#include "KalmanMuonCalibrator.h"
#include "SimpleKalmanCalculator.h"
KalmanMuonCalibrator::KalmanMuonCalibrator() {
 
  double etaBins[] ={-2.5,-2.1,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.1,2.5};
  nEtaBins_ = 22;
  etaAxis_ = new TAxis(nEtaBins_,etaBins);
  double etaMaterialBins[] ={-2.5,-2.1,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.1,2.5};
  nMaterialEtaBins_ = 22;
  etaMaterialAxis_ = new TAxis(nMaterialEtaBins_,etaMaterialBins);
  totalBins_ = 1+1+1+nEtaBins_*13+nMaterialEtaBins_;
  data_ = new double[totalBins_];
  derivative_ = new double[totalBins_];

  //Reset all zeros 
  for (unsigned int i=0;i<totalBins_;++i) {
    data_[i]=0.0;
  }
  printf("Initialized Kalman Filter calibration with %d parameters\n",totalBins_); 
  calculator_ = new SimpleKalmanCalculator(totalBins_,data_);
}

KalmanMuonCalibrator::~KalmanMuonCalibrator() {
  if (etaAxis_)
    delete etaAxis_;
  if (etaMaterialAxis_)
    delete etaMaterialAxis_;
  if (data_)
    delete data_;
  if (derivative_)
    delete derivative_;
  if (calculator_)
    delete calculator_;
}


unsigned int KalmanMuonCalibrator::getBin(Measurement histo,float eta) {
  switch(histo) {
  case A:
    return 0;
  case K:
    return 1;
  case L:
    return 2;
  case A11:
    return 2+etaAxis_->FindBin(eta);
  case A12:
    return 2+nEtaBins_+etaAxis_->FindBin(eta);
  case A21:
    return 2+2*nEtaBins_+etaAxis_->FindBin(eta);
  case A22:
    return 2+3*nEtaBins_+etaAxis_->FindBin(eta);
  case A31:
    return 2+4*nEtaBins_+etaAxis_->FindBin(eta);
  case A32:
    return 2+5*nEtaBins_+etaAxis_->FindBin(eta);
  case B0:
    return 2+6*nEtaBins_+etaAxis_->FindBin(eta);
  case B11:
    return 2+7*nEtaBins_+etaAxis_->FindBin(eta);
  case B12:
    return 2+8*nEtaBins_+etaAxis_->FindBin(eta);
  case B21:
    return 2+9*nEtaBins_+etaAxis_->FindBin(eta);
  case B22:
    return 2+10*nEtaBins_+etaAxis_->FindBin(eta);
  case B31:
    return 2+11*nEtaBins_+etaAxis_->FindBin(eta);
  case B32:
    return 2+12*nEtaBins_+etaAxis_->FindBin(eta);
  case e:
    return 2+13*nEtaBins_+etaMaterialAxis_->FindBin(eta);
  }


  return -1;
}

void KalmanMuonCalibrator::resetDerivative() {
  for (unsigned int i=0;i<totalBins_;++i) {
    derivative_[i]=0.0;
  }
}

  
void KalmanMuonCalibrator::processLine(const RooArgSet* line) {
  double c1 = line->getRealValue("c1");
  double eta1 = line->getRealValue("eta1");
  double phi1 = line->getRealValue("phi1");
  double c2 = line->getRealValue("c2");
  double eta2 = line->getRealValue("eta2");
  double phi2 = line->getRealValue("phi2");
  double scale = line->getRealValue("scale");
  double resolution = line->getRealValue("resolution");
  

  double a_1 = getData(A,eta1);
  double k_1 = getData(K,eta1);
  double l_1 = getData(L,eta1);
  double a11_1 = getData(A11,eta1);
  double a12_1 = getData(A12,eta1);
  double a21_1 = getData(A21,eta1);
  double a22_1 = getData(A22,eta1);
  double a31_1 = getData(A31,eta1);
  double a32_1 = getData(A32,eta1);
  double b0_1 = getData(B0,eta1);
  double b11_1 = getData(B11,eta1);
  double b12_1 = getData(B12,eta1);
  double b21_1 = getData(B21,eta1);
  double b22_1 = getData(B22,eta1);
  double b31_1 = getData(B31,eta1);
  double b32_1 = getData(B32,eta1);
  double e_1 = getData(e,eta1);
  double st1 =sin(2*atan(exp(-eta1))); 
  double a_2 = getData(A,eta2);
  double k_2 = getData(K,eta2);
  double l_2 = getData(L,eta2);
  double a11_2 = getData(A11,eta2);
  double a12_2 = getData(A12,eta2);
  double a21_2 = getData(A21,eta2);
  double a22_2 = getData(A22,eta2);
  double a31_2 = getData(A31,eta2);
  double a32_2 = getData(A32,eta2);
  double b0_2 = getData(B0,eta2);
  double b11_2 = getData(B11,eta2);
  double b12_2 = getData(B12,eta2);
  double b21_2 = getData(B21,eta2);
  double b22_2 = getData(B22,eta2);
  double b31_2 = getData(B31,eta2);
  double b32_2 = getData(B32,eta2);
  double e_2 = getData(e,eta2);
  double st2 =sin(2*atan(exp(-eta2))); 

  
  double magnetic1 = a_1+k_1*eta1*eta1+l_1*eta1*eta1*eta1;
  double scaleTracker1 = a11_1*sin(phi1)+ a12_1*cos(phi1)+a21_1*sin(2*phi1)+ a22_1*cos(2*phi1)+a31_1*sin(3*phi1)+ a32_1*cos(3*phi1); 
  double material1 = 1.0/(1+e_1*c1*st1);
  double alignment1 = b0_1+b11_1*sin(phi1)+ b12_1*cos(phi1)+b21_1*sin(2*phi1)+ b22_1*cos(2*phi1)+b31_1*sin(3*phi1)+ b32_1*cos(3*phi1);
  double term1 = magnetic1+scaleTracker1+material1+alignment1/c1;

  
  double magnetic2 = a_2+k_2*eta2*eta2+l_2*eta2*eta2*eta2;
  double scaleTracker2 = a11_2*sin(phi2)+ a12_2*cos(phi2)+a21_2*sin(2*phi2)+ a22_2*cos(2*phi2)+a31_2*sin(3*phi2)+ a32_2*cos(3*phi2); 
  double material2 = 1.0/(1+e_2*c2*st2);
  double alignment2 = b0_2+b11_2*sin(phi2)+ b12_2*cos(phi2)+b21_2*sin(2*phi2)+ b22_2*cos(2*phi2)+b31_2*sin(3*phi2)+ b32_2*cos(3*phi2);
  double term2 = magnetic2+scaleTracker2+material2-alignment2/c2;
  double h= 1.0/sqrt(term1*term2);

  //now the derivative
  resetDerivative();
  addDerivative(A,eta1,0.5*h*term2);
  addDerivative(A,eta2,0.5*h*term1);
  addDerivative(K,eta1,0.5*h*term2*eta1*eta1);
  addDerivative(K,eta2,0.5*h*term1*eta2*eta2);
  addDerivative(L,eta1,0.5*h*term2*eta1*eta1*eta1);
  addDerivative(L,eta2,0.5*h*term1*eta2*eta2*eta2);
  addDerivative(A11,eta1,0.5*h*term2*sin(phi1));
  addDerivative(A11,eta2,0.5*h*term1*sin(phi2));
  addDerivative(A12,eta1,0.5*h*term2*cos(phi1));
  addDerivative(A12,eta2,0.5*h*term1*cos(phi2));
  addDerivative(A21,eta1,0.5*h*term2*sin(2*phi1));
  addDerivative(A21,eta2,0.5*h*term1*sin(2*phi2));
  addDerivative(A22,eta1,0.5*h*term2*cos(2*phi1));
  addDerivative(A22,eta2,0.5*h*term1*cos(2*phi2));
  addDerivative(A31,eta1,0.5*h*term2*sin(3*phi1));
  addDerivative(A31,eta2,0.5*h*term1*sin(3*phi2));
  addDerivative(A32,eta1,0.5*h*term2*cos(3*phi1));
  addDerivative(A32,eta2,0.5*h*term1*cos(3*phi2));
  addDerivative(B0,eta1,0.5*h*term2);
  addDerivative(B0,eta2,-0.5*h*term1);
  addDerivative(B11,eta1,0.5*h*term2*sin(phi1));
  addDerivative(B11,eta2,-0.5*h*term1*sin(phi2));
  addDerivative(B12,eta1,0.5*h*term2*cos(phi1));
  addDerivative(B12,eta2,-0.5*h*term1*cos(phi2));
  addDerivative(B21,eta1,0.5*h*term2*sin(2*phi1));
  addDerivative(B21,eta2,-0.5*h*term1*sin(2*phi2));
  addDerivative(B22,eta1,0.5*h*term2*cos(2*phi1));
  addDerivative(B22,eta2,-0.5*h*term1*cos(2*phi2));
  addDerivative(B31,eta1,0.5*h*term2*sin(3*phi1));
  addDerivative(B31,eta2,-0.5*h*term1*sin(3*phi2));
  addDerivative(B32,eta1,0.5*h*term2*cos(3*phi1));
  addDerivative(B32,eta2,-0.5*h*term1*cos(3*phi2));
  addDerivative(e,eta1,0.5*h*term2*(-st1*c1)/((1+e_1*c1*st1)*(1+e_1*c1*st1)));
  addDerivative(e,eta2,0.5*h*term1*(-st2*c2)/((1+e_2*c2*st2)*(1+e_2*c2*st2)));



  calculator_->iterate(scale-h,resolution,derivative_);
    
}

 
double KalmanMuonCalibrator::getData(Measurement histo,float eta) {
  unsigned int bin = getBin(histo,eta);
  return data_[bin];
}

 void KalmanMuonCalibrator::setData(Measurement histo,float eta,double data) {
  unsigned int bin = getBin(histo,eta);
  data_[bin]=data;
}


 void KalmanMuonCalibrator::addDerivative(Measurement histo,float eta,double data) {
  unsigned int bin = getBin(histo,eta);
  derivative_[bin]=derivative_[bin]+data;
 }
 


 void KalmanMuonCalibrator::processFile(const char* file,const char* datasetName) {
   TFile *f = new TFile(file);
   RooDataSet *dataset = (RooDataSet*)f->Get(datasetName);
   unsigned int entries = dataset->numEntries();
   for (unsigned int i=0;i<entries;++i) {
     const RooArgSet* line = dataset->get(i);
     processLine(line);
     if (i % 100000 == 0) {
       printf("Processed %d/%d entries\n",i,entries);
     }

   }
   f->Close();
 }

 
