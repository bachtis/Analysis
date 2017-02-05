#ifndef SimpleKalmanCalculator_h
#define SimpleKalmanCalculator_h


#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TFile.h"
class SimpleKalmanCalculator  {
 public:
  SimpleKalmanCalculator(unsigned int,double*,double*);
  ~SimpleKalmanCalculator();
  void iterate(double residual,double resolution,const double* H);
  const TMatrixDSym* P() {
    return P_;
  }



 private:
  double* x_;
  unsigned int N_;
  TVectorD *H_;
  TMatrixDSym *P_;
  
};



#endif
