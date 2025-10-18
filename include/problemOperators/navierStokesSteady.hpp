#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

#include "../utilityFuncs.hpp"
#include "../coefficients/fluidCoeffs.hpp"
#include "../coefficients/totalFluidCoeff.hpp"
#include "../integrators/gradContractionIntegrator.hpp"
#include "schurrComplementPrecon.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Solves the Incompressible steady-state 
! Navier-Stokes equation
!
/*****************************************/
class navierStokesSteadyOper : public Operator
{
protected:
   //
   // Forms, FE-spaces and operators
   //
   int dim;
   Array<int> btoffs;
   ParFiniteElementSpace *fesU, *fesP;
   ParBilinearForm *K_uu=NULL;
   ParMixedBilinearForm *K_up=NULL, *K_pu=NULL;
   ParLinearForm *uResidual=NULL, *pResidual=NULL;

   //
   // Integrator Coefficients
   //
   ConstantCoefficient *mu, *rho, *one, *Mone;
   convectiveCoeff *uGradu_c;
   vectorGradientCoeff *Gradu_c;
   GradientGridFunctionCoefficient *Gradp_c;
   VectorGridFunctionCoefficient *u_c;
   totalFluidCoeff *NSM_coeff;

   //
   // Jacobian matrix storage
   //
   mutable OperatorPtr opUU, opUP, opPU;
   mutable BlockOperator *Jacobian=NULL;
   real_t dt, time;

   //
   // Preconditioning matrix storage
   //
   mutable schurrComplementPrecon *scPrecon=NULL;

   //
   // Boundary conditions
   //
   bool UfromCoeffs=false, PfromCoeffs=false;
   Array<int> empty_tdofs;
   Array<int> U_ess_BCTags, P_ess_BCTags;      //Dirch boundary condition tags
   Array<int> U_ess_BCDofs, P_ess_BCDofs;      //Dirch boundary condition dofs
   mutable BlockVector xBlock, rhsBlock;       //auxillary vector(s)
   mutable ParGridFunction *u_gf, *p_gf;       //Coefficient GridFunctions
   mutable ParGridFunction *uBDR_gf, *pBDR_gf; //Dirch boundary condition values
   Array<VectorCoefficient*> U_ess_BCs;        //Velocity boundary conditions
   Array<Coefficient*>       P_ess_BCs;        //Pressure boundary conditions

   void reassembleJacobian() const;
public:
   //Constructor
   navierStokesSteadyOper(ParFiniteElementSpace *f_u, ParFiniteElementSpace *f_p, int &dim_, real_t &dt);

   //Destructor
   ~navierStokesSteadyOper();

   //Nonlinear Mult
   virtual void Mult(const Vector &u, Vector &du_dt) const;

   //Returns a handle to the Jacobian
   mfem::Operator & GetGradient(const mfem::Vector &x) const override;

   //Set the boundary condition tagging and TDof arrays
   void SetPressureBCTags(const Array<int> & P_ess_BCTags_);
   void SetVelocityBCTags(const Array<int> & U_ess_BCTags_);

   //Set the Boundary condition coefficient arrays
   //or the BC gridfunctions
   void SetVelocityBCCoeffs(const Array<VectorCoefficient*> & U_ess_BCs_);
   void SetPressureBCCoeffs(const Array<Coefficient*>       & P_ess_BCs_);

   void SetVelocityBCGfs(const ParGridFunction & uBDR_gf_);
   void SetPressureBCGfs(const ParGridFunction & pBDR_gf_);

   //Get the preconditioner
   mfem::Solver & GetPrecon() const;
};

/*****************************************\
!
! Constructor builds the integrators for
! the problem class
!
/*****************************************/
navierStokesSteadyOper::navierStokesSteadyOper(ParFiniteElementSpace *f_u, ParFiniteElementSpace *f_p, int &dim_, real_t &dt)
   : Operator(f_u->GetTrueVSize() + f_p->GetTrueVSize()), fesU(f_u), fesP(f_p),
     btoffs(3), dt(0.0)
{
   const real_t rel_tol = 1e-8;
   dim = dim_;

   //Set the Offsets and block Vectors
   btoffs[0] = 0;
   btoffs[1] = f_u->GetTrueVSize();
   btoffs[2] = f_p->GetTrueVSize();
   btoffs.PartialSum();
   xBlock.Update(btoffs);
   rhsBlock.Update(btoffs);

   //
   // Setup the gridfunctions
   //
   u_gf    = new ParGridFunction(fesU);
   uBDR_gf = new ParGridFunction(fesU);
   p_gf    = new ParGridFunction(fesP);
   pBDR_gf = new ParGridFunction(fesP);
   *u_gf    = 0.00;
   *uBDR_gf = 0.00;
   *p_gf    = 0.00;
   *pBDR_gf = 0.00;

   //
   // Setup the coefficients
   //
   one       = new ConstantCoefficient(1.0);
   Mone      = new ConstantCoefficient(-1.0);
   mu        = new ConstantCoefficient(1.0);
   rho       = new ConstantCoefficient(1.0);
   uGradu_c  = new convectiveCoeff(dim,u_gf);
   Gradp_c   = new GradientGridFunctionCoefficient(p_gf);
   NSM_coeff = new totalFluidCoeff(dim, u_gf, p_gf);

   u_c       = new VectorGridFunctionCoefficient(u_gf);
   Gradu_c   = new vectorGradientCoeff(dim,u_gf);


   // Within the forms u,p represent the pressure
   // and velocity fields while v,h represent the
   // test functions that correspond to those fields

   // Residual Forms
   //

   //Momentum equation (du/dt)
   uResidual = new ParLinearForm(fesU);
   uResidual->AddDomainIntegrator(new gradContractionIntegrator(*NSM_coeff,dim) ); //-div(uu^T - mu Grad(u) + pI)

   //Continuity equation (dp/dt = 0.0)
   pResidual = new ParLinearForm(fesP);
   pResidual->AddDomainIntegrator(new DomainLFGradIntegrator(*u_c)); // div(u)

   // Jacobian Forms
   //

   //Velocity Momentum
   K_uu = new ParBilinearForm(fesU);
   K_uu->AddDomainIntegrator(new VectorDiffusionIntegrator(*mu));             // -div(mu grad(u))  Diffusion
   K_uu->AddDomainIntegrator(new VectorMassIntegrator(*Gradu_c));             // [v, grad(u) v] Convection0
   K_uu->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(*u_c)); // [v, u grad(v)] Convection1

   //Pressure Momentum
   K_up = new ParMixedBilinearForm(fesU,fesP);
   K_up->AddDomainIntegrator(new GradientIntegrator(*one));         // [v,grad(h)]  Pressure gradient

   //Continuity
   K_pu = new ParMixedBilinearForm(fesP,fesU);
   K_pu->AddDomainIntegrator(new VectorDivergenceIntegrator(*Mone)); // [div(v),h]  Continuity Error

   //Build the Jacobian
   if(Jacobian != NULL ) delete Jacobian;
   Jacobian = NULL;
   Jacobian = new BlockOperator(btoffs);
   reassembleJacobian();

   //build the preconditioner
   if(scPrecon != NULL ) delete scPrecon;
   scPrecon = NULL;
   scPrecon= new schurrComplementPrecon(opUU.Ptr(), opUP.Ptr(), opPU.Ptr());
   scPrecon->preconRebuild();
};


/*****************************************\
!
! Solve the Backward-Euler equation: 

! k = f(u + dt*k, t), for the unknown k.
! This is the only requirement for 
! high-order SDIRK implicit integration.
!
/*****************************************/
navierStokesSteadyOper::~navierStokesSteadyOper(){
   delete K_uu, K_up, K_pu;
   delete uResidual, pResidual;
};


/*****************************************\
!
! Update the Jacobian block Operator
!
/*****************************************/
void navierStokesSteadyOper::reassembleJacobian() const{
  //Update and re-assemble the matrices
  K_uu->Update();
  K_uu->Assemble();

  K_up->Update();
  K_up->Assemble();

  K_pu->Update();
  K_pu->Assemble();

  //Form the matrix operators
  K_uu->FormSystemMatrix(U_ess_BCDofs, opUU);
  K_up->FormRectangularSystemMatrix(U_ess_BCDofs, P_ess_BCDofs, opUP);
  K_pu->FormRectangularSystemMatrix(P_ess_BCDofs, U_ess_BCDofs, opPU);

  //Set the block matrix operator
  Jacobian->SetBlock(0,0, opUU.Ptr());
  Jacobian->SetBlock(1,0, opUP.Ptr(), -1.0);
  Jacobian->SetBlock(0,1, opPU.Ptr(), -1.0);
};


/*****************************************\
!
!       Get the preconditioner
!
/*****************************************/
mfem::Solver & navierStokesSteadyOper::GetPrecon() const{
  scPrecon->preconRebuild();
  return *scPrecon;
};


/*****************************************\
!
! Calculates the du/dt entry by computing
! M du/dt = K u
! then inverts M to give
! 0 = Minv K u
!
/*****************************************/
void navierStokesSteadyOper::Mult(const Vector &u, Vector &du_dt) const{
  //Copy in the current solution
  copyVec(u, xBlock);

  //Apply Dirchelet BC values
  if(UfromCoeffs) ApplyDircheletBCs<VectorCoefficient>(U_ess_BCTags, U_ess_BCs, uBDR_gf);
  if(PfromCoeffs) ApplyDircheletBCs<Coefficient>(P_ess_BCTags, P_ess_BCs, uBDR_gf);
  applyDirchValues(*uBDR_gf, xBlock.GetBlock(0), U_ess_BCDofs);
  applyDirchValues(*pBDR_gf, xBlock.GetBlock(1), P_ess_BCDofs);

  //Update the GridFunctions
  u_gf->Distribute(xBlock.GetBlock(0));
  p_gf->Distribute(xBlock.GetBlock(1));

  //Reassemble the residuals
  uResidual->Assemble();
  pResidual->Assemble();
  uResidual->ParallelAssemble(rhsBlock.GetBlock(0));
  pResidual->ParallelAssemble(rhsBlock.GetBlock(1));

  //eliminate Dirchelet BC residuals
  if(U_ess_BCDofs.Size() != 0) rhsBlock.GetBlock(0).SetSubVector(U_ess_BCDofs,0.00);
  if(P_ess_BCDofs.Size() != 0) rhsBlock.GetBlock(1).SetSubVector(P_ess_BCDofs,0.00);

  //Copy out the residual
  rhsBlock *= -1.00;
  copyVec(rhsBlock, du_dt);
};


/*****************************************\
!
!  Return a reference to the Jacobian
!              matrix
!
/*****************************************/
mfem::Operator & navierStokesSteadyOper::GetGradient(const mfem::Vector &x) const
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
! Setup the BC tag arrays
!
/*****************************************/
void navierStokesSteadyOper::SetPressureBCTags(const Array<int> & P_ess_BCTags_){
  P_ess_BCTags.SetSize(P_ess_BCTags_.Size());
  for(int I=0; I<P_ess_BCTags_.Size(); I++) P_ess_BCTags[I] = P_ess_BCTags_[I];
  fesU->GetEssentialTrueDofs(P_ess_BCTags, U_ess_BCDofs);
};


void navierStokesSteadyOper::SetVelocityBCTags(const Array<int> & U_ess_BCTags_){
  U_ess_BCTags.SetSize(U_ess_BCTags_.Size());
  for(int I=0; I<U_ess_BCTags_.Size(); I++) U_ess_BCTags[I] = U_ess_BCTags_[I];
  fesU->GetEssentialTrueDofs(U_ess_BCTags, U_ess_BCDofs);
};


/*****************************************\
!
! Set the Boundary condition coefficient
!     arrays or the BC gridfunctions
!
/*****************************************/
void navierStokesSteadyOper::SetVelocityBCCoeffs(const Array<VectorCoefficient*> & U_ess_BCs_){
  if(U_ess_BCs.Size() != 0) for(int I=0; I<U_ess_BCs.Size(); I++) U_ess_BCs[I]=NULL;
  U_ess_BCs.SetSize(U_ess_BCs_.Size());
  for(int I=0; I<U_ess_BCs_.Size(); I++) U_ess_BCs[I] = U_ess_BCs_[I];
  UfromCoeffs=true;
};


void navierStokesSteadyOper::SetPressureBCCoeffs(const Array<Coefficient*> & P_ess_BCs_){
  if(P_ess_BCs.Size() != 0) for(int I=0; I<P_ess_BCs.Size(); I++) P_ess_BCs[I]=NULL;
  P_ess_BCs.SetSize(P_ess_BCs_.Size());
  for(int I=0; I<P_ess_BCs_.Size(); I++) P_ess_BCs[I] = P_ess_BCs_[I];
  PfromCoeffs=true;
};


void navierStokesSteadyOper::SetVelocityBCGfs(const ParGridFunction & uBDR_gf_){
  *uBDR_gf = uBDR_gf_;
};


void navierStokesSteadyOper::SetPressureBCGfs(const ParGridFunction & pBDR_gf_){
  *pBDR_gf = pBDR_gf_;
};


