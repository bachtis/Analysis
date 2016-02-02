#ifndef KaMuCa_Derivation_SimpleKalmanCalculator_h
#define KaMuCa_Derivation_SimpleKalmanCalculator_h


#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TFile.h"

class SimpleKalmanCalculator  {
 public:
  SimpleKalmanCalculator(const TVectorD& , const TMatrixDSym& );

  void iterate(double residual,double resolution,const TVectorD& H);
  const TVectorD* x() {
    return x_;
  }

  const TMatrixDSym* P() {
    return P_;
  }

 private:
  TVectorD *x_;
  TMatrixDSym *P_; 
  

};



#endif
