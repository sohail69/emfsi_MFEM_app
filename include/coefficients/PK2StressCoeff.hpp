/*****************************************\
! Author: Sohail Rathore
! Date  : 2025-05-18
!
! Produces the 2nd Piola Kirchoff
! stress-tensor multiplied  by a series 
! of strain matrices for Passive strain
! case, used for the evaluation of the
! following form:
!
! rho_ij = (S_ik + S_ki)*(du_k/dx_j)+S_ik*d_kj 
! (1)
!
! This is from the solid mechanics energy
! form, which is given as the following
! tensor contraction:
!
! PI = S_ij*E_ij
! (2)
!
! The E_ij term is the Green-Lagrange
! strain measure, it is given by the 
! following expression:
!
! (du_k/dx_i)*(du_k/dx_j)+2*d_ik*(du_k/dx_j)
! (3)
!
! The reason for defining the coefficient
! as rho_ij is so that the following
! equivalence can be made:
!
! rho_ij*(du_i/dx_j) = S_ij*E_ij
! (4)
!
! This significantly simplifies the
! solid-residual-jacobian integration
! procedure
\*****************************************/
#include <functional>
#include "mfem.hpp"




class PK2StressCoeff : public MatrixCoefficient
{
private:
  Array<ParGridFunction*> *fibreBasis; //Fibre GFuncs
  ParGridFunction *u;                  //Displacement field
  ParGridFunction *p;                  //Pressure field

  const int & dim;
public:
   //Define the Orthotropic diffusion coefficient
   //matrix, using a vector of scalar diff coeffs
   //and an array of vector fibre fields
   PK2StressCoeff(Array<ParGridFunction*> *fibreBasis_, const Vector & diff_, const int & dim_)
      :  MatrixCoefficient(dim_), dim(dim_), fibreBasis(fibreBasis_), diff_v(diff_){};

   using MatrixCoefficient::Eval;

   /// Evaluate the matrix coefficient at @a ip.l
   void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~PK2StressCoeff(){};
};


// Evaluate the diffusion tensor at the 
// integration points
// The diffusion tensor is given as
// D_ij = a_p (f_i f_j)_p
void PK2StressCoeff::Eval(DenseMatrix &DiffMat, ElementTransformation &T,const IntegrationPoint &ip){
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
};

