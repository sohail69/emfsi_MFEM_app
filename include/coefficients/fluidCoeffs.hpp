#pragma once

#include "mfem.hpp"
#include "../utilityFuncs.hpp"

using namespace std;
using namespace mfem;

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


/*****************************************\
!
! Generates the vector Laplacian coefficient
!
!  vector Laplacian : Div(Grad(u))
!
\*****************************************/
class vectorLaplacianCoeff : public VectorCoefficient
{
private:
   int dim;
   GridFunction *vel; //Velocity field

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   vectorLaplacianCoeff(int dim_, GridFunction *vel_) : VectorCoefficient(dim_), dim(dim_), vel(vel_){}

   /// Evaluate the matrix coefficient at @a ip.l
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~vectorLaplacianCoeff() {}
};

//
// Evaluate the convection at
// a integration point
//
void vectorLaplacianCoeff::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  DenseMatrix gradU_ij(dim);  
  if(V.Size() != dim) V.SetSize(dim);
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
