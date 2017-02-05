#include "KalmanMuonCalibrator.h"
#include "SimpleKalmanCalculator.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"

KalmanMuonCalibrator::KalmanMuonCalibrator() {
  //  double etaBins[] ={-2.5,-2.1,-1.5,-1.1,-0.9,0.0,0.9,1.1,1.5,2.1,2.5};
  //  nEtaBins_ = 10;
  double etaBins[] ={-2.5,-2.25,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.5,-0.3,0.0,0.3,0.5,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.25,2.5};
  nEtaBins_ = 24;


  etaAxis_ = new TAxis(nEtaBins_,etaBins);
  double etaMaterialBins[] ={-2.5,-2.25,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.5,-0.3,0.0,0.3,0.5,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.25,2.5};
  nMaterialEtaBins_ = 24;
  etaMaterialAxis_ = new TAxis(nMaterialEtaBins_,etaMaterialBins);
  totalBins_ = nEtaBins_*14+2*nMaterialEtaBins_;
  data_ = new double[totalBins_];
  derivative_ = new double[totalBins_];

  //Reset all zeros 
  data_[0]=1.0;
  for (unsigned int i=0;i<nEtaBins_;++i) {
    data_[i]=1.0;
  }

  for (unsigned int i=nEtaBins_;i<totalBins_-nEtaBins_;++i) {
    data_[i]=0.0;
  }
  printf("Initialized Kalman Filter calibration with %d parameters\n",totalBins_); 


  //Initiate errors
  double *errors = new double[totalBins_];
  for (unsigned int i=0;i<7*nEtaBins_;++i)
    errors[i]=0.01;
  for (unsigned int i=7*nEtaBins_;i<14*nEtaBins_;++i)
    errors[i]=0.001;
  for (unsigned int i=14*nEtaBins_;i<totalBins_-nMaterialEtaBins_;++i)
    errors[i]=0.5;
  for (unsigned int i=totalBins_-nMaterialEtaBins_;i<totalBins_;++i)
    errors[i]=2.0;
  calculator_ = new SimpleKalmanCalculator(totalBins_,data_,errors);





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
    return etaAxis_->FindBin(eta)-1;
  case A11:
    return nEtaBins_+etaAxis_->FindBin(eta)-1;
  case A12:
    return 2*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case A21:
    return 3*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case A22:
    return 4*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case A31:
    return 5*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case A32:
    return 6*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B0:
    return 7*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B11:
    return 8*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B12:
    return 9*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B21:
    return 10*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B22:
    return 11*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B31:
    return 12*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case B32:
    return 13*nEtaBins_+etaAxis_->FindBin(eta)-1;
  case e:
    return 14*nEtaBins_+etaMaterialAxis_->FindBin(eta)-1;
  case e2:
    return 14*nEtaBins_+nMaterialEtaBins_+etaMaterialAxis_->FindBin(eta)-1;
  }


  return -1;
}

void KalmanMuonCalibrator::resetDerivative() {
  for (unsigned int i=0;i<totalBins_;++i) {
    derivative_[i]=0.0;
  }
}

  
void KalmanMuonCalibrator::processLine(double c1,double eta1,double phi1,double c2,double eta2,double phi2,double scale,double resolution,int Z) {
  double a_1 = getData(A,eta1);
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
  double e2_1 = getData(e2,eta1);
  double st1 =sin(2*atan(exp(-eta1))); 
  double tt1 =tan(2*atan(exp(-eta1))); 
  double a_2 = getData(A,eta2);
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
  double e2_2 = getData(e2,eta2);
  double st2 =sin(2*atan(exp(-eta2))); 
  double tt2 =tan(2*atan(exp(-eta2))); 

  
  double magnetic1 = a_1;
  double scaleTracker1 = a11_1*sin(phi1)+ a12_1*cos(phi1)+a21_1*sin(2*phi1)+ a22_1*cos(2*phi1)+a31_1*sin(3*phi1)+ a32_1*cos(3*phi1); 
  double material1= -e_1*st1*c1-e2_1*st1*c1*c1*c1;
  double alignment1 = b0_1+b11_1*sin(phi1)+ b12_1*cos(phi1)+b21_1*sin(2*phi1)+ b22_1*cos(2*phi1)+b31_1*sin(3*phi1)+ b32_1*cos(3*phi1);
  double term1 = magnetic1+scaleTracker1+material1+alignment1/c1;
  
  double magnetic2 = a_2;
  double scaleTracker2 = a11_2*sin(phi2)+ a12_2*cos(phi2)+a21_2*sin(2*phi2)+ a22_2*cos(2*phi2)+a31_2*sin(3*phi2)+ a32_2*cos(3*phi2); 
  double material2= -e_2*st2*c2-e2_2*st2*c2*c2*c2;
  double alignment2 = b0_2+b11_2*sin(phi2)+ b12_2*cos(phi2)+b21_2*sin(2*phi2)+ b22_2*cos(2*phi2)+b31_2*sin(3*phi2)+ b32_2*cos(3*phi2);
  double term2 = magnetic2+scaleTracker2+material2-alignment2/c2;
  double h= sqrt(term1*term2);
  if (term1*term2<0)
    printf("ERROR NEGATIVE ROOT!!!!! \n");

  //now the derivative
  resetDerivative();
  addDerivative(A,eta1,(0.5/h)*term2);
  addDerivative(A,eta2,(0.5/h)*term1);
  //addDerivative(K,eta1,(0.5/h)*term2/(tt1*tt1));
  //  addDerivative(K,eta2,(0.5/h)*term1/(tt2*tt2));
  //  addDerivative(L,eta1,(0.5/h)*term2*eta1*eta1*eta1);
  //  addDerivative(L,eta2,(0.5/h)*term1*eta2*eta2*eta2);

  addDerivative(A11,eta1,(0.5/h)*term2*sin(phi1));
  addDerivative(A11,eta2,(0.5/h)*term1*sin(phi2));
  addDerivative(A12,eta1,(0.5/h)*term2*cos(phi1));
  addDerivative(A12,eta2,(0.5/h)*term1*cos(phi2));
  addDerivative(A21,eta1,(0.5/h)*term2*sin(2*phi1));
  addDerivative(A21,eta2,(0.5/h)*term1*sin(2*phi2));
  addDerivative(A22,eta1,(0.5/h)*term2*cos(2*phi1));
  addDerivative(A22,eta2,(0.5/h)*term1*cos(2*phi2));
  addDerivative(A31,eta1,(0.5/h)*term2*sin(3*phi1));
  addDerivative(A31,eta2,(0.5/h)*term1*sin(3*phi2));
  addDerivative(A32,eta1,(0.5/h)*term2*cos(3*phi1));
  addDerivative(A32,eta2,(0.5/h)*term1*cos(3*phi2));
  addDerivative(B0,eta1,(0.5/h)*term2/c1);
  addDerivative(B0,eta2,-(0.5/h)*term1/c2);
  addDerivative(B11,eta1,(0.5/h)*term2*sin(phi1)/c1);
  addDerivative(B11,eta2,-(0.5/h)*term1*sin(phi2)/c2);
  addDerivative(B12,eta1,(0.5/h)*term2*cos(phi1)/c1);
  addDerivative(B12,eta2,-(0.5/h)*term1*cos(phi2)/c2);
  addDerivative(B21,eta1,(0.5/h)*term2*sin(2*phi1)/c1);
  addDerivative(B21,eta2,-(0.5/h)*term1*sin(2*phi2)/c2);
  addDerivative(B22,eta1,(0.5/h)*term2*cos(2*phi1)/c1);
  addDerivative(B22,eta2,-(0.5/h)*term1*cos(2*phi2)/c2);
  addDerivative(B31,eta1,(0.5/h)*term2*sin(3*phi1)/c1);
  addDerivative(B31,eta2,-(0.5/h)*term1*sin(3*phi2)/c2);
  addDerivative(B32,eta1,(0.5/h)*term2*cos(3*phi1)/c1);
  addDerivative(B32,eta2,-(0.5/h)*term1*cos(3*phi2)/c2);
  if (Z==0) {
    addDerivative(e,eta1,(0.5/h)*term2*(-st1*c1));
    addDerivative(e,eta2,(0.5/h)*term1*(-st2*c2));
    //    addDerivative(e2,eta1,(0.5/h)*term2*(-st1*c1*c1*c1));
    //    addDerivative(e2,eta2,(0.5/h)*term1*(-st2*c2*c2*c2));

  }
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
 


 void KalmanMuonCalibrator::processFile(const char* file,const char* treeName) {
   TFile *f = new TFile(file);
   TTree *dataset = (TTree*)f->Get(treeName);
   unsigned int entries = dataset->GetEntries();
   float c1,eta1,phi1,c2,eta2,phi2,scale,resolution;

   int Z=0;
   char *str = (char*)file;
   char *pch = strstr(str, "Z");

   if (pch) {
     Z=1;
     printf ("Fitting Z,ignore the material\n"); 
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
  
   //save the magnetic function
   // TH1D * magnetic = new TH1D("magnetic","magnetic",3,0,3);
   // magnetic->SetBinContent(1,data_[0]);
   // magnetic->SetBinError(1,data->GetBinError(1));
   // magnetic->SetBinContent(2,data_[1]);
   // magnetic->SetBinError(2,data->GetBinError(2));
   // magnetic->SetBinContent(3,data_[2]);
   // magnetic->SetBinError(3,data->GetBinError(3));
   // magnetic->Write();
   //Aii

   TH1D * a = new TH1D("A","A",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a->SetBinContent(i,data_[i-1]);
     a->SetBinError(i,data->GetBinError(i));
   }
   a->Write();


   TH1D * a11 = new TH1D("A11","A11",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a11->SetBinContent(i,data_[nEtaBins_+i-1]);
     a11->SetBinError(i,data->GetBinError(nEtaBins_+i));
   }
   a11->Write();

   TH1D * a12 = new TH1D("A12","A12",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a12->SetBinContent(i,data_[2*nEtaBins_+i-1]);
     a12->SetBinError(i,data->GetBinError(2*nEtaBins_+i));

   }
   a12->Write();

   TH1D * a21 = new TH1D("A21","A21",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a21->SetBinContent(i,data_[3*nEtaBins_+i-1]);
     a21->SetBinError(i,data->GetBinError(3*nEtaBins_+i));
   }
   a21->Write();

   TH1D * a22 = new TH1D("A22","A22",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a22->SetBinContent(i,data_[4*nEtaBins_+i-1]);
     a22->SetBinError(i,data->GetBinError(4*nEtaBins_+i));

   }
   a22->Write();

   TH1D * a31 = new TH1D("A31","A31",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a31->SetBinContent(i,data_[5*nEtaBins_+i-1]);
     a31->SetBinError(i,data->GetBinError(5*nEtaBins_+i));

   }
   a31->Write();

   TH1D * a32 = new TH1D("A32","A32",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     a32->SetBinContent(i,data_[6*nEtaBins_+i-1]);
     a32->SetBinError(i,data->GetBinError(6*nEtaBins_+i));

   }
   a32->Write();

   TH1D * b0 = new TH1D("B0","B0",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b0->SetBinContent(i,data_[7*nEtaBins_+i-1]);
     b0->SetBinError(i,data->GetBinError(7*nEtaBins_+i));

   }
   b0->Write();

   TH1D * b11 = new TH1D("B11","B11",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b11->SetBinContent(i,data_[8*nEtaBins_+i-1]);
     b11->SetBinError(i,data->GetBinError(8*nEtaBins_+i));

   }
   b11->Write();

   TH1D * b12 = new TH1D("B12","B12",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b12->SetBinContent(i,data_[9*nEtaBins_+i-1]);
     b12->SetBinError(i,data->GetBinError(9*nEtaBins_+i));

   }
   b12->Write();

   TH1D * b21 = new TH1D("B21","B21",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b21->SetBinContent(i,data_[10*nEtaBins_+i-1]);
     b21->SetBinError(i,data->GetBinError(10*nEtaBins_+i));

   }
   b21->Write();

   TH1D * b22 = new TH1D("B22","B22",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b22->SetBinContent(i,data_[11*nEtaBins_+i-1]);
     b22->SetBinError(i,data->GetBinError(11*nEtaBins_+i));
   }
   b22->Write();

   TH1D * b31 = new TH1D("B31","B31",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b31->SetBinContent(i,data_[12*nEtaBins_+i-1]);
     b31->SetBinError(i,data->GetBinError(12*nEtaBins_+i));

   }
   b31->Write();

   TH1D * b32 = new TH1D("B32","B32",nEtaBins_,etaAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nEtaBins_;++i) {
     b32->SetBinContent(i,data_[13*nEtaBins_+i-1]);
     b32->SetBinError(i,data->GetBinError(13*nEtaBins_+i));

   }
   b32->Write();


   TH1D * e1 = new TH1D("e1","e1",nMaterialEtaBins_,etaMaterialAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nMaterialEtaBins_;++i) {
     e1->SetBinContent(i,data_[14*nEtaBins_+i-1]);
     e1->SetBinError(i,data->GetBinError(14*nEtaBins_+i));

   }
   e1->Write();

   TH1D * e2 = new TH1D("e2","e2",nMaterialEtaBins_,etaMaterialAxis_->GetXbins()->GetArray());
   for (unsigned int i=1;i<=nMaterialEtaBins_;++i) {
     e2->SetBinContent(i,data_[14*nEtaBins_+nMaterialEtaBins_+i-1]);
     e2->SetBinError(i,data->GetBinError(14*nEtaBins_+nMaterialEtaBins_+i));

   }
   e2->Write();
   f->Close();
 }




 void KalmanMuonCalibrator::processFileMC(const char* file,const char* treeName) {
   TFile *f = new TFile(file);
   TTree *dataset = (TTree*)f->Get(treeName);
   unsigned int entries = dataset->GetEntries();
   float c1,eta1,phi1,c2,eta2,phi2,mass,genMass,resolution;

   int Z=0;
   char *str = (char*)file;
   char *pch = strstr(str, "Z");
   if (pch) {
     Z=1;
     printf ("Fitting Z,ignore the material\n"); 
   }
   dataset->SetBranchAddress("c1",&c1);
   dataset->SetBranchAddress("c2",&c2);
   dataset->SetBranchAddress("eta1",&eta1);
   dataset->SetBranchAddress("eta2",&eta2);
   dataset->SetBranchAddress("phi1",&phi1);
   dataset->SetBranchAddress("phi2",&phi2);
   dataset->SetBranchAddress("mass",&mass);
   dataset->SetBranchAddress("genMass",&genMass);
   dataset->SetBranchAddress("resolution",&resolution);


   for (unsigned int i=0;i<entries;++i) {
     dataset->GetEntry(i);
     processLine(c1,eta1,phi1,c2,eta2,phi2,mass/genMass,resolution,Z);
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



