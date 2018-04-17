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
  //  void iterate2D(const TMatrixD& residual,const TMatrixDSym& resolution,const TMatrixD& H);
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
