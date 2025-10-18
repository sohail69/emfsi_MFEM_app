#pragma once

#include "mfem.hpp"
#include <functional>
#include "../utilityFuncs.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Generates the totalFluidCoeff coefficient
!
!
!  totalFluidCoeff : uu^T - mu Grad(u) + pI
!
\*****************************************/
class totalFluidCoeff : public MatrixCoefficient
{
private:
   int dim;
   real_t mu=1.0;
   GridFunction *veloc; //Velocity field
   GridFunction *press; //Pressure field
public:
   /// Define a time-independent vector coefficient from a std function
   totalFluidCoeff(int dim_, GridFunction *vel_, GridFunction *pres_) : 
                 MatrixCoefficient(dim_), dim(dim_), veloc(vel_), press(pres_) {}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~totalFluidCoeff(){}
};

//
// Evaluate the convection at
// a integration point
//
void totalFluidCoeff::Eval(DenseMatrix &Fmat, ElementTransformation &T, const IntegrationPoint &ip){
  if(Fmat.Size() != dim) Fmat.SetSize(dim);
  real_t p;
  Vector U_i(dim);
  DenseMatrix gradU_ij(dim);

  Mesh *gf_mesh = veloc->FESpace()->GetMesh();
  if (T.mesh->GetNE() == gf_mesh->GetNE())
  {
    p = press->GetValue(T, ip);
    veloc->GetVectorValue(T, ip, U_i);
    veloc->GetVectorGradient(T, gradU_ij);
  }
  else
  {
    IntegrationPoint coarse_ip;
    ElementTransformation *coarse_T = RefinedToCoarse(*gf_mesh, T, ip, coarse_ip);
    p = press->GetValue(*coarse_T, ip);
    veloc->GetVectorValue(*coarse_T, coarse_ip, U_i);
    veloc->GetVectorGradient(*coarse_T, gradU_ij);
  }

  for(int I=0; I<dim; I++){
    for(int J=0; J<dim; J++){
      Fmat(I,J) = U_i(I)*U_i(J) - mu*gradU_ij(I,J) + p*kDelta(I,J);
    }
  }
};

/*****************************************\
!
! Generates the vector gradient coefficient
!
!
!  vector Gradient : Grad(u)
!
\*****************************************/
class vectorGradientCoeff : public MatrixCoefficient
{
private:
   int dim;
   GridFunction *vel; //Velocity field

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   vectorGradientCoeff(int dim_, GridFunction *vel_) : MatrixCoefficient(dim_), dim(dim_), vel(vel_){}

   /// Evaluate the matrix coefficient at @a ip.l
   void Eval(DenseMatrix &gradU_ij, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~vectorGradientCoeff() {}
};

//
// Evaluate the convection at
// a integration point
//
void vectorGradientCoeff::Eval(DenseMatrix &gradU_ij, ElementTransformation &T, const IntegrationPoint &ip){
  if(gradU_ij.Size() != dim) gradU_ij.SetSize(dim);
  Mesh *gf_mesh = vel->FESpace()->GetMesh();
  if (T.mesh->GetNE() == gf_mesh->GetNE())
  {
    vel->GetVectorGradient(T, gradU_ij);
  }
  else
  {
    IntegrationPoint coarse_ip;
    ElementTransformation *coarse_T = RefinedToCoarse(*gf_mesh, T, ip, coarse_ip);
    vel->GetVectorGradient(*coarse_T, gradU_ij);
  }
};

/*****************************************\
!
! Generates the convection coefficient
!
!
!  convectiveCoeff : Convection u.Grad(u)
!
\*****************************************/
class convectiveCoeff : public VectorCoefficient
{
private:
   int dim;
   GridFunction *vel; //Velocity field

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   convectiveCoeff(int dim_, GridFunction *vel_) : VectorCoefficient(dim_), dim(dim_), vel(vel_){}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~convectiveCoeff() {}
};

//
// Evaluate the convection at
// a integration point
//
void convectiveCoeff::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){

  if(V.Size() != dim) V.SetSize(dim);
  Vector U(dim);
  DenseMatrix gradU(dim);
  Mesh *gf_mesh = vel->FESpace()->GetMesh();
  if (T.mesh->GetNE() == gf_mesh->GetNE())
  {
    vel->GetVectorValue(T, ip, U);
    vel->GetVectorGradient(T, gradU);
  }
  else
  {
    IntegrationPoint coarse_ip;
    ElementTransformation *coarse_T = RefinedToCoarse(*gf_mesh, T, ip, coarse_ip);
    vel->GetVectorValue(*coarse_T, coarse_ip, U);
    vel->GetVectorGradient(*coarse_T, gradU);
  }
  gradU.Mult(U,V);
};

