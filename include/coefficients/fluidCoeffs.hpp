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
   convectiveCoeff(int dim_, GridFunction *vel_) : dim(dim_), VectorCoefficient(dim), vel(vel_){}

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
  Vector      U(dim);
  DenseMatrix gradU(dim,dim);
  u->GetVectorValue(T, ip, U);
  u->GetVectorGradient(T, gradU);
  gradU.Mult(U,V);
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
   vectorGradientCoeff(int dim_, GridFunction *vel_) : dim(dim_), MatrixCoefficient(dim_), vel(vel_){}

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
  u->GetVectorGradient(T, gradU_ij);
};
