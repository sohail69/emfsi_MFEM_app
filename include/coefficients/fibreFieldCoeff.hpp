#pragma once
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Generates the fibre field map from
! the potential field gradient and a
! central direction. specifically used
! for axisymmetric geometries
!
!
!
\*****************************************/
class fibreFieldCoeff1 : public VectorCoefficient
{
private:
   GridFunction & phi; //Potential field GFun

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff1(int dim, GridFunction & phi_) : VectorCoefficient(dim), phi(phi_){}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~fibreFieldCoeff1() {}
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
void fibreFieldCoeff::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  phi.GetGradient(T,V);
};

/*
//   Array<GridFunction*> & fibreBasis; //Fibre GFuncs
  // GridFunction         & diff;       //Diffusion GFun

   real_t x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   V.SetSize(vdim);
   if (Function)
   {
      Function(transip, V);
   }
   else
   {
      TDFunction(transip, GetTime(), V);
   }
   if (Q)
   {
      V *= Q->Eval(T, ip, GetTime());
   }


*/
