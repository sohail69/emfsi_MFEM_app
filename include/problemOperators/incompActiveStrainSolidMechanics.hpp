#pragma once
#include "mfem.hpp"
#include "../coefficients/PK2StressCoeff.hpp"
#include "../integrators/gradContractionIntegrator.hpp"
#include "schurrComplementPrecon.hpp"

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
   unsigned nDIM;
   Array<int> btoffs;
   ParFiniteElementSpace *fesU, *fesP;
   ParLinearForm         *u_Residual, *p_Residual;
   ParMixedBilinearForm  *K_up, *K_pu;
   ParBilinearForm       *K_uu, *K_pp;

   // Jacobian matrix storage
   mutable OperatorPtr opUU, opUP, opPU, opPP;
   mutable BlockOperator *Jacobian=NULL;

   // Preconditioning matrix storage
   mutable schurrComplementPrecon *scPrecon=NULL;

   //Coefficients
   NeoHookeanPK2StressCoeff *PK2;

   //Internal gridfunctions used for coefficient
   //update and boundary conditions
   mutable Array<ParGridFunction*> *fibreBasis;
   mutable ParGridFunction *u_gf, *p_gf, *phi_gf;
   mutable BlockVector xBlock, rhsBlock;       //auxillary vector(s)

   //Dirchelet boundary condition
   Array<int>    _ess_bcs_markers;
   Array<double> _BC_Vals;


   // Boundary conditions
   bool UfromCoeffs=false, PfromCoeffs=false;
   mutable ParGridFunction *uBDR_gf, *pBDR_gf; //Dirch boundary condition values
   Array<int> empty_tdofs;
   Array<int> U_ess_BCTags,  P_ess_BCTags;     //Dirch boundary condition tags
   Array<int> U_ess_BCDofs,  P_ess_BCDofs;     //Dirch boundary condition dofs
   Array<VectorCoefficient*> U_ess_BCs;        //Velocity boundary conditions
   Array<Coefficient*>       P_ess_BCs;        //Pressure boundary conditions

   //Build the operators for the Jacobian
   void reassembleJacobian() const;
public:
   //Constructs fibre map operator
   solidMechanicsIASOper(ParFiniteElementSpace *u_space, ParFiniteElementSpace *p_space
                       , Array<ParGridFunction*> *fibreBasis_);

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
solidMechanicsIASOper::solidMechanicsIASOper(ParFiniteElementSpace *u_space, ParFiniteElementSpace *p_space
                                           , Array<ParGridFunction*> *fibreBasis_):
                                             Operator(fesU->GetTrueVSize()+fesP->GetTrueVSize())
                                           , fesU(u_space), fesP(p_space), btoffs(3)
{
   const real_t rel_tol = 1e-8;
   nDIM = fibreBasis_->Size();;

   //Set the Offsets and block Vectors
   btoffs[0] = 0;
   btoffs[1] = fesU->GetTrueVSize();
   btoffs[2] = fesP->GetTrueVSize();
   btoffs.PartialSum();
   xBlock.Update(btoffs);
   rhsBlock.Update(btoffs);

  // Build the reference grid
  // functions
  fibreBasis = fibreBasis_;
  u_gf   = new ParGridFunction(fesU);
  p_gf   = new ParGridFunction(fesP);
  phi_gf = new ParGridFunction(fesP);

  // Build the coffeicients
  // necessary for the Linear
  // and Bilinear Forms
  PK2 = new NeoHookeanPK2StressCoeff(fibreBasis, u_gf, p_gf, phi_gf, nDIM);

  // Build the linear/residual 
  // forms and add integrators
  u_Residual = new ParLinearForm(fesU);
  u_Residual->AddDomainIntegrator(new gradContractionIntegrator(*PK2,nDIM) ); //-div(rho_ij . grad(u) )

  p_Residual = new ParLinearForm(fesP);
//  p_Residual->AddDomainIntegrator(new DomainLFIntegrator() ); // ln(J3(u))
//  p_Residual->AddDomainIntegrator(new DomainLFIntegrator() ); // -p.(1.0/lmbda)
}


/*****************************************\
!
! Updates the values of the active strain
!
\*****************************************/
void solidMechanicsIASOper::UpdatePhi(const Vector &x) const
{
  phi_gf->Distribute(x);
};


/*****************************************\
!
  //Build the operators for the Jacobian
!
\*****************************************/
void solidMechanicsIASOper::reassembleJacobian() const
{
  //Update and re-assemble the matrices
  K_uu->Update();
  K_uu->Assemble();

  K_up->Update();
  K_up->Assemble();

  K_pu->Update();
  K_pu->Assemble();

  K_pp->Update();
  K_pp->Assemble();

  //Form the matrix operators
  K_uu->FormSystemMatrix(U_ess_BCDofs, opUU);
  K_up->FormRectangularSystemMatrix(U_ess_BCDofs, P_ess_BCDofs, opUP);
  K_pu->FormRectangularSystemMatrix(P_ess_BCDofs, U_ess_BCDofs, opPU);
  K_pp->FormSystemMatrix(P_ess_BCDofs, opPP);

  //Set the block matrix operator
  Jacobian->SetBlock(0,0, opUU.Ptr());
  Jacobian->SetBlock(1,0, opUP.Ptr(), -1.0);
  Jacobian->SetBlock(0,1, opPU.Ptr(), -1.0);
  Jacobian->SetBlock(1,1, opPP.Ptr());
};


/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
void solidMechanicsIASOper::Mult(const Vector &x, Vector &Residual) const
{
  //Copy in the current solution
  copyVec(x, xBlock);

  //Apply Dirchelet BC values
  if(UfromCoeffs) ApplyDircheletBCs<VectorCoefficient>(U_ess_BCTags, U_ess_BCs, uBDR_gf);
  if(PfromCoeffs) ApplyDircheletBCs<Coefficient>(P_ess_BCTags, P_ess_BCs, uBDR_gf);
  applyDirchValues(*uBDR_gf, xBlock.GetBlock(0), U_ess_BCDofs);
  applyDirchValues(*pBDR_gf, xBlock.GetBlock(1), P_ess_BCDofs);

  //Update the GridFunctions
  u_gf->Distribute(xBlock.GetBlock(0));
  p_gf->Distribute(xBlock.GetBlock(1));

  //Reassemble the residuals
  u_Residual->Assemble();
  p_Residual->Assemble();
  u_Residual->ParallelAssemble(rhsBlock.GetBlock(0));
  p_Residual->ParallelAssemble(rhsBlock.GetBlock(1));

  //eliminate Dirchelet BC residuals
  if(U_ess_BCDofs.Size() != 0) rhsBlock.GetBlock(0).SetSubVector(U_ess_BCDofs,0.00);
  if(P_ess_BCDofs.Size() != 0) rhsBlock.GetBlock(1).SetSubVector(P_ess_BCDofs,0.00);

  //Copy out the residual
  rhsBlock *= -1.00;
  copyVec(rhsBlock, Residual);
};

/*****************************************\
!
   //Calculates and returns the Residual
   //and recalculates the Jacobian
!
\*****************************************/
 mfem::Operator & solidMechanicsIASOper::GetGradient(const mfem::Vector &x) const
{
  //Copy in the current solution
  copyVec(x, xBlock);

  //Apply Dirchelet BC values
  if(UfromCoeffs) ApplyDircheletBCs<VectorCoefficient>(U_ess_BCTags, U_ess_BCs, uBDR_gf);
  if(PfromCoeffs) ApplyDircheletBCs<Coefficient>(P_ess_BCTags, P_ess_BCs, uBDR_gf);
  applyDirchValues(*uBDR_gf, xBlock.GetBlock(0), U_ess_BCDofs);
  applyDirchValues(*pBDR_gf, xBlock.GetBlock(1), P_ess_BCDofs);

  //Update the GridFunctions
  u_gf->Distribute(xBlock.GetBlock(0));
  p_gf->Distribute(xBlock.GetBlock(1));

  //Update the Jacobian
  reassembleJacobian();
  return *Jacobian;
};


/*****************************************\
!
   //Destroys the fibre map operator
!
\*****************************************/
solidMechanicsIASOper::~solidMechanicsIASOper()
{
  delete u_gf, p_gf, phi_gf;
  delete u_Residual, p_Residual;
};

















