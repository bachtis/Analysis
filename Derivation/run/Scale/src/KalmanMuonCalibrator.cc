#include "KalmanMuonCalibrator.h"
#include "SimpleKalmanCalculator.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include <math.h>
#include <unistd.h>
KalmanMuonCalibrator::KalmanMuonCalibrator() {
  //  double etaBins[] ={-2.5,-2.1,-1.5,-1.1,-0.9,0.0,0.9,1.1,1.5,2.1,2.5};
  //  nEtaBins_ = 10;

  //  double etaMaterialBins[] ={-2.5,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.15,0.0,0.15,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.5};
  double etaMaterialBins[] ={-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-0.9,-0.5,0.0,0.5,0.9,1.2,1.4,1.6,1.8,2.0,2.2,2.4};

  //  nEtaMaterialBins_ = 46;
  nEtaMaterialBins_ = 18;
  etaMaterialAxis_ = new TAxis(nEtaMaterialBins_,etaMaterialBins);
  double etaBins[] ={-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-0.9,-0.5,0.0,0.5,0.9,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
  nEtaBins_ = 18;
  etaAxis_ = new TAxis(nEtaBins_,etaBins);

  double phiBins[] ={-M_PI,-5.0*M_PI/6.0,-4.0*M_PI/6.0,-3.0*M_PI/6.0,-2*M_PI/6.0,-1.0*M_PI/6.0,0.0,M_PI/6.0,2*M_PI/6.0,3*M_PI/6.0,4*M_PI/6.0,5*M_PI/6.0,M_PI};
    nPhiBins_ = 12;
    //double phiBins[] ={-M_PI,M_PI};
    //  nPhiBins_ = 1;

  phiAxis_ = new TAxis(nPhiBins_,phiBins);
  totalBins_ = 2+nEtaBins_+nEtaMaterialBins_+nEtaBins_*nPhiBins_;
  data_ = new double[totalBins_];
  derivative_ = new double[totalBins_];
  //Initiate errors
  double *errors = new double[totalBins_];
  
  //Reset all zeros 

  for (unsigned int i=0;i<totalBins_;++i) {
    data_[i]=0.0;
    errors[i]=0.0;
  }

  data_[0] = 1.0;
  errors[0] = 0.5;
  errors[1] = 0.5;


  for (unsigned int i=0;i<nEtaBins_;++i)
    errors[2+i]=1.0/5.0;
  for (unsigned int i=0;i<nEtaMaterialBins_;++i)
    errors[2+nEtaBins_+i]=5;
  for (unsigned int i=0;i<nEtaBins_*nPhiBins_;++i)
    errors[2+i+nEtaBins_+nEtaMaterialBins_]=0.001;
  printf("Initializing Kalman Filter calibration with %d parameters\n",totalBins_); 
  calculator_ = new SimpleKalmanCalculator(totalBins_,data_,errors);
  printf("Initialized Kalman Filter calibration with %d parameters\n",totalBins_); 



}

KalmanMuonCalibrator::~KalmanMuonCalibrator() {
  if (etaAxis_)
    delete etaAxis_;
  if (etaMaterialAxis_)
    delete etaMaterialAxis_;
  if (phiAxis_)
    delete phiAxis_;
  if (data_)
    delete data_;
  if (derivative_)
    delete derivative_;
  if (calculator_)
    delete calculator_;
}


unsigned int KalmanMuonCalibrator::getBin(Measurement histo,float eta,float phi) {

  switch(histo) {
  case A:
    return 0;
  case K:
    return 1;
  case L:
    return 2+etaAxis_->FindBin(eta)-1;
  case e:
    return 2+nEtaBins_+etaMaterialAxis_->FindBin(eta)-1;
  case B:
    return 2+nEtaBins_+nEtaMaterialBins_+etaAxis_->FindBin(eta)+nEtaBins_*(phiAxis_->FindBin(phi)-1)-1;
  }
  return -1;
}

void KalmanMuonCalibrator::resetDerivative() {
  for (unsigned int i=0;i<totalBins_;++i) {
    derivative_[i]=0.0;
  }
}

  
void KalmanMuonCalibrator::processLine(double c1,double eta1,double phi1,double c2,double eta2,double phi2,double scale,double resolution,int Z) {
  double a_1 = getData(A,eta1,phi1);
  double k_1 = getData(K,eta1,phi1);
  double l_1 = getData(L,eta1,phi1);
  double b_1 = getData(B,eta1,phi1);
  double e_1 = getData(e,eta1,phi1);
  double st1 = sin(2*atan(exp(-eta1))); 
  //  double tt1 = tan(2*atan(exp(-eta1))); 

  double a_2 = getData(A,eta2,phi2);
  double k_2 = getData(K,eta2,phi2);
  double l_2 = getData(L,eta2,phi2);
  double b_2 = getData(B,eta2,phi2);
  double e_2 = getData(e,eta2,phi2);
  double st2 = sin(2*atan(exp(-eta2))); 
  //  double tt2 = tan(2*atan(exp(-eta2))); 


  double magnetic1 = a_1+k_1*eta1*eta1;
  double material1= -e_1*st1*c1;
  double term1 = magnetic1+material1+b_1/c1+l_1/c1;
  double magnetic2 = a_2+k_2*eta2*eta2;
  double material2= -e_2*st2*c2;
  double term2 = magnetic2+material2-b_2/c2+l_2/c2;

  if (term1*term2<0) {
    printf("ERROR NEGATIVE ROOT!!!!! \n");
    printf("%f+%f (%f)^2 -%f %f +%f/%f \n",a_1,k_1,eta1,e_1,c1,b_1,c1);
    printf("%f+%f (%f)^2 -%f %f +%f/%f \n",a_2,k_2,eta2,e_2,c2,b_2,c2);
    printf(" %f %f\n",term1,term2);
    pause();

    return;
  }

  double h= sqrt(term1*term2);
  double q=1.0/(2*h);
  resetDerivative();

  addDerivative(A,eta1,phi1,q*term2);
  addDerivative(A,eta2,phi2,q*term1);

  addDerivative(K,eta1,phi1,q*term2*eta1*eta1);
  addDerivative(K,eta2,phi2,q*term1*eta2*eta2);

  addDerivative(e,eta1,phi1,q*term2*(-st1*c1));
  addDerivative(e,eta2,phi2,q*term1*(-st2*c2));

  addDerivative(B,eta1,phi1,q*term2/c1);
  addDerivative(B,eta2,phi2,-q*term1/c2);

  //addDerivative(L,eta1,phi1,q*term2/c1);
  //addDerivative(L,eta2,phi2,q*term1/c2);

  calculator_->iterate(scale-h,resolution,derivative_);  
}









 
double KalmanMuonCalibrator::getData(Measurement histo,float eta,float phi) {
  unsigned int bin = getBin(histo,eta,phi);
  return data_[bin];
}

void KalmanMuonCalibrator::setData(Measurement histo,float eta,float phi,double data) {
  unsigned int bin = getBin(histo,eta,phi);
  data_[bin]=data;
}


void KalmanMuonCalibrator::addDerivative(Measurement histo,float eta,float phi,double data) {
  unsigned int bin = getBin(histo,eta,phi);
  derivative_[bin]=derivative_[bin]+data;
 }
 


 void KalmanMuonCalibrator::processFile(const char* file,const char* treeName) {
   TFile *f = new TFile(file);
   TTree *dataset = (TTree*)f->Get(treeName);
   unsigned int entries = dataset->GetEntries();
   float c1,eta1,phi1,c2,eta2,phi2,scale,resolution;

   int Z=1;
   char *str = (char*)file;
   char *pch = strstr(str, "LowPt");

   if (pch) {
     Z=0;
     printf ("Fitting Low PT for material"); 
   }
   else {
     printf ("Fitting High PT\n"); 
   }

   dataset->SetBranchAddress("c1",&c1);
   dataset->SetBranchAddress("c2",&c2);
   dataset->SetBranchAddress("eta1",&eta1);
   dataset->SetBranchAddress("eta2",&eta2);
   dataset->SetBranchAddress("phi1",&phi1);
   dataset->SetBranchAddress("phi2",&phi2);
   dataset->SetBranchAddress("scale",&scale);
   dataset->SetBranchAddress("resolution",&resolution);

   for (unsigned int i=0;i<entries;++i) {
     dataset->GetEntry(i);
     processLine(c1,eta1,phi1,c2,eta2,phi2,scale,resolution,Z);
     if (i % 5000000 == 0) {
       printf("Processed %d/%d entries\n",i,entries);
       char* newFile;
       if(asprintf(&newFile,"previewMM_%d.root",i)<0)
      	 continue;
       save(newFile);
      }

   }
   f->Close();
 }






TMatrixDSym KalmanMuonCalibrator::correlationMatrix() {
  TMatrixDSym P = *(calculator_->P()); 
  TMatrixDSym diagonal(P.GetNrows());
  diagonal.Zero();
  for (unsigned int i=0;i<P.GetNrows();++i) {
    diagonal[i][i]=sqrt(P[i][i]);
  }
  diagonal.Invert();
  return P.Similarity(diagonal);
}


void KalmanMuonCalibrator::load(const char* file) {
  TFile *f = new TFile(file);
  TH1D* data = (TH1D*)f->Get("data");
  for (unsigned int i=1;i<data->GetNbinsX()+1;++i)
    data_[i-1] = data->GetBinContent(i);
}

 
 void KalmanMuonCalibrator::save(const char* file) {
   TFile *f = new TFile(file,"RECREATE");
   f->cd();
   const TMatrixDSym* covariance = calculator_->P();
   TMatrixDSym correlation = correlationMatrix();
   TH1D * data = new TH1D("data","data",totalBins_,0,totalBins_);

   for (unsigned int i=0;i<totalBins_;++i) {
     data->SetBinContent(i+1,data_[i]);
     double sigma2 = (*covariance)[i][i];
     if (sigma2>0.0)
       data->SetBinError(i+1,sqrt(sigma2));
   }
   data->Write();  
   covariance->Write("covariance");
   correlation.Write("correlation");
   //Diagonalize the matrix and save eigenvalues and eigenvectors
   TMatrixDSymEigen eigen(*covariance);
   TVectorD eigenvalues = eigen.GetEigenValues();
   TMatrixD eigenvectors = eigen.GetEigenVectors();
   eigenvalues.Write("eigenvalues");
   eigenvectors.Write("eigenvectors");
  
   TH1D * magnetic = new TH1D("magnetic","A",2,0,2);
   magnetic->SetBinContent(1,data_[0]);
   magnetic->SetBinError(1,data->GetBinError(1));
   magnetic->SetBinContent(2,data_[1]);
   magnetic->SetBinError(2,data->GetBinError(2));
   magnetic->Write();

    TH1D * K = new TH1D("L","K",nEtaBins_,etaAxis_->GetXbins()->GetArray());
    for (unsigned int i=1;i<=nEtaBins_;++i) {
      K->SetBinContent(i,data_[2+i-1]);
      K->SetBinError(i,data->GetBinError(2+i));
    }
    K->Write();


    TH1D * e = new TH1D("e","e",nEtaMaterialBins_,etaMaterialAxis_->GetXbins()->GetArray());
    for (unsigned int i=1;i<=nEtaMaterialBins_;++i) {
      e->SetBinContent(i,data_[2+nEtaBins_+i-1]);
      e->SetBinError(i,data->GetBinError(2+nEtaBins_+i));
    }
    e->Write();



    TH2D * b = new TH2D("B","B",nEtaBins_,etaAxis_->GetXbins()->GetArray(),nPhiBins_,phiAxis_->GetXbins()->GetArray());
    for (unsigned int i=1;i<=nEtaBins_;++i) {
      for (unsigned int j=1;j<=nPhiBins_;++j) {
	b->SetBinContent(i,j,data_[2+nEtaBins_+nEtaMaterialBins_+i+(j-1)*nEtaBins_-1]);
	b->SetBinError(i,j,data->GetBinError(2+nEtaBins_+nEtaMaterialBins_+i+(j-1)*nEtaBins_));
      }
    }
    b->Write();


   f->Close();

   
 }




 void KalmanMuonCalibrator::processFileMC(const char* file,const char* treeName) {
   TFile *f = new TFile(file);
   TTree *dataset = (TTree*)f->Get(treeName);
   unsigned int entries = dataset->GetEntries();
   float c1,eta1,phi1,c2,eta2,phi2,mass,gc1,gc2,cErr1,cErr2,genMass,resolution;

   int Z=0;
   char *str = (char*)file;
   char *pch = strstr(str, "Z");
   if (pch) {
     Z=1;
     printf ("Fitting Z,ignore the material\n"); 
   }
   dataset->SetBranchAddress("c1",&c1);
   dataset->SetBranchAddress("c2",&c2);
   dataset->SetBranchAddress("cErr1",&cErr1);
   dataset->SetBranchAddress("cErr2",&cErr2);
   dataset->SetBranchAddress("gc1",&gc1);
   dataset->SetBranchAddress("gc2",&gc2);
   dataset->SetBranchAddress("eta1",&eta1);
   dataset->SetBranchAddress("eta2",&eta2);
   dataset->SetBranchAddress("phi1",&phi1);
   dataset->SetBranchAddress("phi2",&phi2);
   dataset->SetBranchAddress("mass",&mass);
   dataset->SetBranchAddress("genMass",&genMass);
   dataset->SetBranchAddress("resolution",&resolution);
   


   for (unsigned int i=0;i<entries;++i) {
     dataset->GetEntry(i);
     processLine(c1,eta1,phi1,c2,eta2,phi2,genMass/mass,resolution,Z);
     if (i % 5000000 == 0) {
       printf("Processed %d/%d entries\n",i,entries);
       char* newFile;
       if(asprintf(&newFile,"previewMC_%d.root",i)<0)
	 continue;
       save(newFile);

     }

   }
   f->Close();
 }



