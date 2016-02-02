#include "../interface/SimpleKalmanCalculator.h"
SimpleKalmanCalculator::SimpleKalmanCalculator(const TVectorD& x, const TMatrixDSym& P):
    x_(new TVectorD(x)),
    P_( new TMatrixDSym(P))
{
  
}


void 
SimpleKalmanCalculator::iterate(double residual,double resolution, const TVectorD& H) {
  //In OUR 1D measurement cases create the scalar hph+R and invert it!
  //  printf("Initial P\n");
  //  P_->Print();
    
  double sinv = 1.0/(P_->Similarity(H)+resolution);
  


  //  printf("Similarity=%f, Resolution=%f ,Inverse = %f residual=%f\n",P_->Similarity(H),resolution,sinv,residual);
  //Then calculate the PH transpose
  TVectorD PH = (*P_)*H;


  x_->Add(residual*sinv*PH);

  P_->Rank1Update(PH,-sinv);


}


