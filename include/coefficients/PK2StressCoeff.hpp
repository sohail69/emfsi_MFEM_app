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

real_t kDelta(I,J){return ( (I==J)? 1.00:0.00);};

class NeoHookeanPK2StressCoeff : public MatrixCoefficient
{
private:
  //Sources of Data
  Array<ParGridFunction*> *fibreBasis; //Fibre GFuncs
  ParGridFunction *gama;               //Active strain field
  ParGridFunction *u;                  //Displacement field
  ParGridFunction *p;                  //Pressure field
  const int & dim;
  const real_t Y=1.00, nu=0.499999;
  real_t mu, lmbda;
public:
   //Define the Orthotropic diffusion coefficient
   //matrix, using a vector of scalar diff coeffs
   //and an array of vector fibre fields
   NeoHookeanPK2StressCoeff(Array<ParGridFunction*> *fibreBasis_, const Vector & diff_, const int & dim_)
      :  MatrixCoefficient(dim_), dim(dim_), fibreBasis(fibreBasis_), diff_v(diff_)
  {
    mu    = (4.60/2.20)/(Y/(2.00+2.00*nu));
    lmbda = Y*nu/((1.00+nu)*(1.00-2.00*nu));
  };

   using MatrixCoefficient::Eval;

   /// Evaluate the matrix coefficient at @a ip.l
   void Eval(DenseMatrix &rho_ij, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~NeoHookeanPK2StressCoeff(){};
};


// Evaluate the Stress tensor at the 
// integration points
// The diffusion tensor is given as
void NeoHookeanPK2StressCoeff::Eval(DenseMatrix &rho_ij, ElementTransformation &T,const IntegrationPoint &ip){
  if(rho_ij.Size() != dim)  rho_ij.SetSize(dim);
  DenseMatrix S_ij(dim), gradU(dim)

  //Calculate deformation measures
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
    }
  }

  //Calculate S_ij
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
    }
  }  

  //Calculate rho_ij
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      for(int K=0; K<dim; K++){
        rho_ij(I,J) = (S_ij(I,K) + S_ij(K,I))*gradU(K,J) + S_ij(I,K)*kDelta(K,J);
      }
    }
  }
};


