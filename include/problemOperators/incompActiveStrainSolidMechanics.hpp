#pragma once
#include "mfem.hpp"
#include "../coefficients/fibreFieldCoeff.hpp"

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
   ParFiniteElementSpace & _fespace;
   ParMixedBilinearForm *K_uu, *K_up, *K_pu, *K_pp;
   ParBilinearForm *_K_uncon;
   HypreParMatrix _Kmat, _Kmat_uncon;

   CGSolver _K_solver;    // Krylov solver for inverting the mass matrix K
   HypreSmoother _K_prec; // Preconditioner for the diffusion matrix K

   //The Preconditioning objects
   HypreBoomerAMG *invM=NULL;
   mutable ParGridFunction *z; // auxiliary GFunc used for BC's

public:
   //Constructs fibre map operator
   solidMechanicsIASOper(ParFiniteElementSpace &f, Array<int> & ess_bcs_markers, Array<double> & BC_Vals);

   //Solves the Poisson equation
   virtual void Mult(const Vector &x, Vector &y) const;

   //Destroys the fibre map operator
   virtual ~solidMechanicsIASOper();

};
