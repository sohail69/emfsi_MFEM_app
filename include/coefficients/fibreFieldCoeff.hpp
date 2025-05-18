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
!  fibreFieldCoeff1 : Sheetlet-direction
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

//
// Evaluate the cross-sheet-direction
// using the gradient of the potenial
// field
//
void fibreFieldCoeff1::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  real_t v_norm = 0.00;
  phi.GetGradient(T,V);
  for(int I=0; I<V.Size();I++) v_norm += V(I)*V(I);
  v_norm = sqrt(v_norm);
  for(int I=0; I<V.Size();I++) V(I) = V(I)/v_norm;
};

/*****************************************\
!
! Generates the fibre field map from
! the potential field gradient and a
! central direction. specifically used
! for axisymmetric geometries
!
!  fibreFieldCoeff2 : fibre-direction
!
\*****************************************/
class fibreFieldCoeff2 : public VectorCoefficient
{
private:
   GridFunction & phi; //Potential field GFun
   GridFunction & S0;  //The sheetlet direction

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff2(int dim, GridFunction & phi_, GridFunction & S0_)
                 : VectorCoefficient(dim), phi(phi_), S0(S0_){}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~fibreFieldCoeff2() {}
};

//
// Evaluate the cross-sheet-direction
// using the gradient of the potenial
// field
//
void fibreFieldCoeff2::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  real_t v_norm = 0.00;
  phi.GetGradient(T,V);
  for(int I=0; I<V.Size();I++) v_norm += V(I)*V(I);

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
