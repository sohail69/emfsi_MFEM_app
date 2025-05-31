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
