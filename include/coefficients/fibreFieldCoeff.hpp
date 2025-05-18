#pragma once
#include "mfem.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Cross product of two 3D vectors
!
\*****************************************/
void Crossproduct3D(const Vector & A, const Vector & B, Vector & AxB){
  AxB(0) = A(1)*B(2) - A(2)*B(1);
  AxB(1) = A(2)*B(0) - A(0)*B(2);
  AxB(2) = A(0)*B(1) - A(1)*B(0);
};

/*****************************************\
!
! Calculate a rotation matrix using the
! Rodrigues' rotation formula 
!
\*****************************************/
void RodriguesRotMat(const Vector & S0, const real_t & theta, DenseMatrix & R0){

/*
  R0 = 0._iwp
  CALL CROSS_PD_MATRIX(s0x, s0, ndim)
  DO I = 1,ndim
    DO J = 1,ndim
      R0(I,J) = kdelta(I,J) + sin(theta)*s0x(I,J) &
              + 2._iwp*sin(theta/2._iwp)          &
              *sin(theta/2._iwp)*(s0(I)*s0(J) - kdelta(I,J))
    ENDDO
  ENDDO
*/

};

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
   int dim;
   GridFunction & phi; //Potential field GFun
   GridFunction & S0;  //The sheetlet direction
   Vector & C0;        //The azimuthal vector
   real_t & theta_epi; //Epicardial fibre angle
   real_t & theta_end; //Endocardial fibre angle
public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff2(int dim_, GridFunction & phi_, GridFunction & S0_
                  , Vector C0_, real_t & theta_epi_, real_t & theta_end_)
                 : dim(dim_), VectorCoefficient(dim), phi(phi_), S0(S0_), C0(C0_)
                 , theta_epi(theta_epi_), theta_end(theta_end_) {}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~fibreFieldCoeff2() {}
};

//
// Evaluate the fibre-direction
// using the Rodrigues' rotation
// formula 
//
void fibreFieldCoeff2::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  real_t theta, phi_i, kp_norm = 0.00;
  DenseMatrix R0(dim,dim);
  Vector Kp(dim), f0(dim), s0;
  phi_i = phi.GetValue(T, ip);
  S0.GetVectorValue(T, ip, s0);
  theta = (theta_epi - theta_end)*phi_i + theta_end;
  Crossproduct3D(s0, Kp, f0);
  RodriguesRotMat(s0, theta, R0);
  R0.Mult(f0, V);
};
/*
   real_t & theta_epi; //Epicardial fibre angle
   real_t & theta_end; //Endocardial fibre angle

  // Array<GridFunction*> & fibreBasis; //Fibre GFuncs
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
