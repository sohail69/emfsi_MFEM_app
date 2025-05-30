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
! rho_ij = 0.5*((S_ik + S_ki)*(du_k/dx_j)+S_ik*d_kj + S_jk d_ki)
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
! 0.5*(du_k/dx_i)*(du_k/dx_j)+d_ik*(du_k/dx_j)+d_jk*(du_k/dx_i)
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
!
! This hard codes the deviatoric PK2 stress
! for the Neo-Hookean Model
!
\*****************************************/
#pragma once
#include "mfem.hpp"
#include <functional>
#include "../utilityFuncs.hpp"

using namespace std;
using namespace mfem;

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
   NeoHookeanPK2StressCoeff(Array<ParGridFunction*> *fibreBasis_, ParGridFunction *u_, ParGridFunction *p_
                          , ParGridFunction *gama_, const int & dim_)
      :  MatrixCoefficient(dim_), dim(dim_), fibreBasis(fibreBasis_), u(u_), p(p_), gama(gama_)
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
  real_t gs0, gf0, gn0, press, S_kk, logJ;
  DenseMatrix S_ij(dim), gradU(dim), CInv(dim), FsInv(dim), F0(dim), Fe(dim), FeT(dim);
  Array<Vector*> f0(dim);

  //Get values at local Ip and calc
  //the active strain
  #pragma unroll
  for(int I=0; I<dim; I++) (*fibreBasis)[I]->GetVectorValue(T, ip, *f0[I]);
  gf0 = -(gama->GetValue(T, ip));
  gs0 =  4.00*gf0;
  gn0 =  (1.00/((1.00 + gf0)*(1.00 + gn0))) - 1.00;
  u->GetVectorGradient(T, gradU);
  press = p->GetValue(T, ip);

  //Calculate deformation measures
  #pragma unroll
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      FsInv(I,J) = kDelta(I,J) + gs0*(*f0[0])(I)*(*f0[0])(J);
      if(dim > 1) FsInv(I,J) = FsInv(I,J) + gf0*(*f0[1])(I)*(*f0[1])(J);
      if(dim > 2) FsInv(I,J) = FsInv(I,J) + gn0*(*f0[2])(I)*(*f0[2])(J);
      F0(I,J)= + gradU(I,J) + kDelta(I,J);
    }
  }
  FsInv.Invert();
  FeT.Transpose(Fe);
  Mult(F0, FsInv, Fe);
  Mult(FeT, Fe, CInv);
  CInv.Invert();
  logJ = log(Fe.Det());

  //Calculate S_ij for standard formulation
  S_kk = 0.00;
  #pragma unroll
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      S_ij(I,J) = mu*(kDelta(I,J) - CInv(I,J)) + lmbda*logJ*CInv(I,J);
      S_kk += ((I==J) ? S_ij(I,J) : 0.00);
    }
  }

  //Augment S_ij for the mixed formulation
  #pragma unroll
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      S_ij(I,J) = S_ij(I,J) + press*CInv(I,J) - S_kk*kDelta(I,J)/3.00;
    }
  }

  //Calculate rho_ij
  #pragma unroll
  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      for(int K=0; K<dim; K++){
        rho_ij(I,J) = 0.50*( (S_ij(I,K) + S_ij(K,I))*gradU(K,J) + S_ij(I,K)*kDelta(K,J) + S_ij(J,K)*kDelta(K,I));
      }
    }
  }
};
