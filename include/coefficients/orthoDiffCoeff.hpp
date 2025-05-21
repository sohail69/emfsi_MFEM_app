#pragma once
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!  Author: Sohail Rathore
!
! Constant Diff-Coeff Version
!
! produces a orthotropic matrix coefficient
! based on a given set of orthotropic fibre 
! vector fields and a given vector of
! scalar coefficients. 
!
! The diffusion tensor is given as
! D_ij = a_p (f_i f_j)_p
!
\*****************************************/
class orthoDiffCoeff : public MatrixCoefficient
{
private:
    Array<ParGridFunction*> *fibreBasis; //Fibre GFuncs
   const Vector  & diff_v;                  //Diffusion Vector
   const int & dim;
public:
   //Define the Orthotropic diffusion coefficient
   //matrix, using a vector of scalar diff coeffs
   //and an array of vector fibre fields
   orthoDiffCoeff(Array<ParGridFunction*> *fibreBasis_, const Vector & diff_, const int & dim_)
      :  MatrixCoefficient(dim_), dim(dim_), fibreBasis(fibreBasis_), diff_v(diff_){};

   using MatrixCoefficient::Eval;

   /// Evaluate the matrix coefficient at @a ip.l
   void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~orthoDiffCoeff(){};
};


// Evaluate the diffusion tensor at the 
// integration points
// The diffusion tensor is given as
// D_ij = a_p (f_i f_j)_p
void orthoDiffCoeff::Eval(DenseMatrix &DiffMat, ElementTransformation &T,const IntegrationPoint &ip){
  if(DiffMat.Size() != dim)  DiffMat.SetSize(dim);

  //Zero out the matrix
  for(int J=0; J<dim; J++)
    for(int K=0; K<dim; K++)
      DiffMat(J,K) = 0.00;


  //Add the Orthotropic Components
  for(int I=0; I<dim; I++){
    Vector f_i;
    (*fibreBasis)[I]->GetVectorValue(T, ip, f_i);
    for(int J=0; J<dim; J++){
      for(int K=0; K<dim; K++){
        if(J==K) DiffMat(J,K) += diff_v(I) * f_i(J) * f_i(K);
      }
    }
  }
};
