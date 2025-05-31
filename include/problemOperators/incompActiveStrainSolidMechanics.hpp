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
   Array<int>    & _ess_bcs_markers;
   Array<double> & _BC_Vals;
   ParFiniteElementSpace *fesU, *fesP;
   ParMixedBilinearForm *K_uu, *K_up, *K_pu, *K_pp;
   ParBilinearForm *K_uu_uncon, *K_pp_uncon;
   HypreParMatrix _Kmat;

   HypreSmoother _K_prec; // Preconditioner for the diffusion matrix K

   //The Preconditioning objects
   HypreBoomerAMG *invM=NULL;


public:
   //Constructs fibre map operator
   solidMechanicsIASOper(ParFiniteElementSpace &f, ParFiniteElementSpace &f
                       , Array<int> & ess_bcs_markers, Array<double> & BC_Vals);

   //Calculates and returns the Residual
   //and recalculates the Jacobian
   virtual void Mult(const Vector &x, Vector &Residual) const;

   //Returns a handle to the Jacobian
   mfem::Operator & GetGradient(const mfem::Vector &x) const override;

   //Destroys the fibre map operator
   virtual ~solidMechanicsIASOper();

};
