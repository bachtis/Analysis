#include "SimpleKalmanCalculator.h"
SimpleKalmanCalculator::SimpleKalmanCalculator(unsigned int nElements,double* elements):
  x_(elements), 
  N_(nElements),
  H_(new TVectorD(nElements)),
  P_( new TMatrixDSym(nElements))
{
  //X  is the state vector 1xN
  //H_ is a derivative vector 1x N (not needed in GPU implementation )
  //P_ is a NXN symmetric matrix (covariance matrix)

  
}


SimpleKalmanCalculator::~SimpleKalmanCalculator(){
  if (H_)
    delete H_;
  if (P_)
    delete P_;
}


void 
SimpleKalmanCalculator::iterate(double residual,double resolution, const double* derivative) {
  //Here I put the elements in a TVector to use the ROOT Matrix package. This should not be done in GPU
  //implementation but rather use directly the vector 
  H_->SetElements(derivative);
     
  //First calculation  S  = HPH^T  is a number
  //Sinv= 1/(S+resolution)
  double sinv = 1.0/(P_->Similarity((*H_))+resolution); 


  //second calculation PH = P x H  (NxN) matrix times (Nx1) matrix 
  TVectorD PH = (*P_)*(*H_);

  
  //Then set the values X = X+Sinv*PH*residual
  for (unsigned int i=0;i<N_;++i)
    (x_)[i] = (x_)[i]+residual*sinv*PH[i];


  //Finally update the covariance matrix
  //P=P-sinv*PH*(PH)^T
  P_->Rank1Update(PH,-sinv);
}


