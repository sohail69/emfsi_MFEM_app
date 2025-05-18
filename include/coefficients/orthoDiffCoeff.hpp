#pragma once
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! produces a orthotropic matrix coefficient
! based on a given set of orthotropic fibre 
! vector fields and a given vector field of
! scalar coefficients
!
! The diffusion tensor is given as
! D_ij = a_p (f_i f_j)_p
!
\*****************************************/
class orthoDiffGFuncCoeff : public MatrixCoefficient
{
private:
   Array<GridFunction*> & fibreBasis; //Fibre GFuncs
   GridFunction         & diff;       //Diffusion GFun
   unsigned dim;

public:
   //Define the Orthotropic diffusion coefficient
   //matrix, using a vector of scalar diff coeffs
   //and an array of vector fibre fields
   orthoDiffGFuncCoeff(Array<GridFunction*> & fibreBasis_, GridFunction & diff_, unsigned dim_)
      : dim(dim_), MatrixCoefficient(dim), fibreBasis(fibreBasis_), diff(diff_){};

   /// Evaluate the matrix coefficient at @a ip.
   void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~orthoDiffGFuncCoeff(){};
};


/*****************************************\
!
! Evaluate the diffusion tensor at the 
! integration points
!
! The diffusion tensor is given as
! D_ij = a_p (f_i f_j)_p
!
\*****************************************/
void orthoDiffGFuncCoeff::Eval(DenseMatrix &DiffMat, ElementTransformation &T,const IntegrationPoint &ip){
  const int nDIM = fibreBasis.Size();
  real_t x[nDIM];
  Vector D_i, f_i, transip(x, nDIM);
  T.Transform(ip, transip);

  diff.GetVectorValue(T, ip, D_i);
  for(int I=0; I<nDIM; I++){
    fibreBasis[I]->GetVectorValue(T, ip, f_i);
    for(int J=0; J<nDIM; J++){
      for(int K=0; K<nDIM; K++){
        DiffMat(J,K) += D_i(I) * f_i(J) * f_i(K);
      }
    }
  }
};
