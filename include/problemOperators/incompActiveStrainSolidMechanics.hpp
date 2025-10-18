#pragma once
#include "mfem.hpp"
#include "../coefficients/PK2StressCoeff.hpp"
#include "../integrators/gradContractionIntegrator.hpp"

using namespace std;
using namespace mfem;


/*****************************************\
!
! Solves the incompressible solid-mechanics
! problem using a mixed U-P formulation for
! a material subject to active strain with
! applied traction stress BC's and external 
! forcing vector b. The equation being
! solved minimised is:
!
! find (u,p) such that:
! MIN(int [(S_ij + p*Cinv_ij) Eij - L*p^2] dV - int [ft_i u_i] dS - int [fb_i u_i] dV )
!
! Traction force BC's are applied using the
! following form:
! ft_i = sigma_ij n_j
!
! where sigma_ij is the surface traction
! stress given by another model
! (Traction is handled this way for FSI
!  applications)
!
/*****************************************/
class solidMechanicsIASOper : public Operator
{
protected:
   ParFiniteElementSpace *fesU, *fesP;
   ParMixedBilinearForm  *K_up, *K_pu;
   ParBilinearForm       *K_uu, *K_pp;
   ParLinearForm         *u_Res, *p_Res;
   HypreParMatrix _Kmat;

   //Coefficients
   NeoHookeanPK2StressCoeff *PK2;

   //Internal gridfunctions used for coefficient
   //update and boundary conditions
   mutable ParGridFunction *u_gf, *p_gf, *phi_gf;

   //Dirchelet boundary condition
   Array<int>    _ess_bcs_markers;
   Array<double> _BC_Vals;
public:
   //Constructs fibre map operator
   solidMechanicsIASOper(ParFiniteElementSpace *u_space, ParFiniteElementSpace *p_space);

   //Updates the values of the active strain
   void UpdatePhi(const Vector &x) const;

   //Calculates and returns the Residual
   //and recalculates the Jacobian
   virtual void Mult(const Vector &x, Vector &Residual) const;

   //Returns a handle to the Jacobian
   mfem::Operator & GetGradient(const mfem::Vector &x) const override;

   //Destroys the fibre map operator
   ~solidMechanicsIASOper();

};


/*****************************************\
!
! Construct the problem Operator
!
\*****************************************/
solidMechanicsIASOper::solidMechanicsIASOper(ParFiniteElementSpace *u_space, ParFiniteElementSpace *p_space):
                                             Operator(fesU->GetTrueVSize()+fesP->GetTrueVSize())
                                           , fesU(u_space), fesP(p_space)
{
  // Build the reference grid
  // functions
  u_gf   = new ParGridFunction(fesU);
  p_gf   = new ParGridFunction(fesP);
  phi_gf = new ParGridFunction(fesP);

  // Build the coffeicients
  // necessary for the Linear
  // and Bilinear Forms
//  PK2 = new NeoHookeanPK2StressCoeff(Array<ParGridFunction*> *fibreBasis_, u_gf, p_gf, phi_gf, const int & dim_)

  // Build the linear/residual 
  // forms and add integrators
  u_Res = new ParLinearForm(fesU);
//  u_Res->AddDomainIntegrator(new gradContractionIntegrator(*NSM_coeff,dim) ); //-div(rho_ij . grad(u) )

  p_Res = new ParLinearForm(fesP);
//  p_Res->AddDomainIntegrator(new DomainLFIntegrator() ); //div(p . u)
}


/*****************************************\
!
! Updates the values of the active strain
!
\*****************************************/
void solidMechanicsIASOper::UpdatePhi(const Vector &x) const
{

};


/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
void solidMechanicsIASOper::Mult(const Vector &x, Vector &Residual) const
{

};


/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
 mfem::Operator & solidMechanicsIASOper::GetGradient(const mfem::Vector &x) const
{
  return;
};


/*****************************************\
!
   //Destroys the fibre map operator
!
\*****************************************/
solidMechanicsIASOper::~solidMechanicsIASOper()
{
};

















