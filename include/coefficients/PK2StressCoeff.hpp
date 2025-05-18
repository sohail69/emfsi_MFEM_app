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


class PK2StressCoeff : MatrixCoefficient
{
protected:
   int height, width;
   std::function<void(const DenseMatrix &K)> S_IJ;

public:
   /// Construct a h x w matrix coefficient.
   PK2StressCoeff(int h, int w, bool symm=false) :
      height(h), width(w), time(0.), symmetric(symm) { }

   /// Construct the matrix coefficients
   PK2StressCoeff(GridFunction & GFuncs)

   /** @brief Evaluate the matrix coefficient in the element described by @a T
       at the point @a ip, storing the result in @a K. */
   /** @note When this method is called, the caller must make sure that the
       IntegrationPoint associated with @a T is the same as @a ip. This can be
       achieved by calling T.SetIntPoint(&ip). */
   virtual void Eval(DenseMatrix &K, ElementTransformation &T,
                     const IntegrationPoint &ip) = 0;

   /// @brief Fill the QuadratureFunction @a qf by evaluating the coefficient at
   /// the quadrature points. The matrix will be transposed or not according to
   /// the boolean argument @a transpose.
   ///
   /// The @a vdim of the QuadratureFunction should be equal to the height times
   /// the width of the matrix.
   virtual void Project(QuadratureFunction &qf, bool transpose=false);
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
void Eval(DenseMatrix &K, ElementTransformation &T,const IntegrationPoint &ip){

 

};

