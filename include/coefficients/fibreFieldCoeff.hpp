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

void OrthoVec2D(const Vector & A, Vector & B){
  B(0) =  A(1);
  B(1) = -A(0);
};

/*****************************************\
!
! Calculate a rotation matrix using the
! Rodrigues' rotation formula 
!
\*****************************************/
// Kronecker delta function
real_t kDelta(int I, int J){
  return ((I==J) ? 1.00 : 0.00);
}

void RodriguesRotMat(const Vector & S0, const real_t & theta, DenseMatrix & R0){
  int DIM = S0.Size();
  DenseMatrix S0x(DIM,DIM);
  real_t sinT  = sin(theta);
  real_t sinT2 = sin(theta/2.00);

  S0x(0,0)= 0.00;   S0x(0,0)=-S0(2);  S0x(0,0)= S0(1);
  S0x(1,0)= S0(2);  S0x(1,1)= 0.00;   S0x(0,0)=-S0(0);
  S0x(2,0)=-S0(1);  S0x(0,0)= S0(0);  S0x(2,2)= 0.00;

  for(int I=0; I<DIM; I++){
    for(int J=0; J<DIM; J++){
      R0(I,J) = kDelta(I,J) 
              + sinT*S0x(I,J)
              + 2.00*sinT2*sinT2*(S0(I)*S0(J) - kDelta(I,J) );
    }
  }
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
   int dim;
   GridFunction *phi; //Potential field GFun

public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff1(int dim_, GridFunction *phi_) : dim(dim_), VectorCoefficient(dim), phi(phi_){}

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
  if(V.Size() != dim) V.SetSize(dim);
  real_t v_normInv = 0.00, tol = 1.0E-15;
  phi->GetGradient(T,V);
  for(int I=0; I<dim; I++) v_normInv += V(I)*V(I);
  v_normInv = 1.00/sqrt(v_normInv);
  V *= v_normInv;
};

/*****************************************\
!
! Generates the fibre field map from
! the potential field, sheetlet vector and
! a central direction vector. specifically
! used for axisymmetric geometries
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
  real_t theta, phi_i, kp_norm = 0.00, alpha = 0.00;
  DenseMatrix R0(dim,dim);
  Vector Kp(dim), f0(dim), s0;
  phi_i = phi.GetValue(T, ip);
  S0.GetVectorValue(T, ip, s0);
  alpha = InnerProduct(s0,C0);
  Kp = s0;
  Kp *= -alpha;
  Kp += C0;
  kp_norm = sqrt(InnerProduct(Kp,Kp));
  Kp *= kp_norm;
  theta = (theta_epi - theta_end)*phi_i + theta_end;
  Crossproduct3D(s0, Kp, f0);
  RodriguesRotMat(s0, theta, R0);
  R0.Mult(f0, V);
};

/*****************************************\
!
! Generates the fibre field map from
! the fibre vector and the sheetlet vector 
! specifically used for axisymmetric
! geometries
!
!  fibreFieldCoeff3 : Cross-sheet-direction
!
\*****************************************/
class fibreFieldCoeff3 : public VectorCoefficient
{
private:
   int dim;
   GridFunction & F0; //Potential field GFun
   GridFunction & S0;  //The sheetlet direction
public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff3(int dim_, GridFunction & F0_, GridFunction & S0_)
                 : dim(dim_), VectorCoefficient(dim), F0(F0_), S0(S0_){}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~fibreFieldCoeff3() {}
};

//
// Evaluate the Cross-Sheet-Direction
// using the cross product of the
// fibre and the sheetlet vectors
//
void fibreFieldCoeff3::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  Vector f0, s0;
  F0.GetVectorValue(T, ip, f0);
  S0.GetVectorValue(T, ip, s0);
  Crossproduct3D(s0, f0, V);
};

/*****************************************\
!
! Generates the fibre field map from
! the fibre vector specifically used for
! 2D geometries
!
!  fibreFieldCoeff4 : 2D Cross-sheet-direction
!
\*****************************************/
class fibreFieldCoeff4 : public VectorCoefficient
{
private:
   int dim;
   GridFunction *F0; //Potential field GFun
public:
   /// Define a time-independent vector coefficient from a std function
   /** vector coefficient **/
   fibreFieldCoeff4(int dim_, GridFunction *F0_)
                 : dim(dim_), VectorCoefficient(dim), F0(F0_){}

   /// Evaluate the vector coefficient at @a ip.
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   virtual ~fibreFieldCoeff4() {}
};

//
// Evaluate the Cross-Sheet-Direction
// using the cross product of the
// fibre and the sheetlet vectors
//
void fibreFieldCoeff4::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip){
  if(V.Size() != dim) V.SetSize(dim);
  Vector f0(dim);
  F0->GetVectorValue(T, ip, f0);
  OrthoVec2D(f0, V);
};










