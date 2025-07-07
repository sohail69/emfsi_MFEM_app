#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

#include "../utilityFuncs.hpp"
#include "../coefficients/fluidCoeffs.hpp"
#include "../coefficients/totalFluidCoeff.hpp"
#include "../integrators/gradContractionIntegrator.hpp"

using namespace std;
using namespace mfem;

/*****************************************\
!
! Solves the Incompressible Navier-Stokes
! equation
!
/*****************************************/
class navierStokesOper : public TimeDependentOperator
{
protected:
   //
   // Forms, FE-spaces and operators
   //
   int dim;
   Array<int> btoffs;
   ParFiniteElementSpace *fesU, *fesP;
   ParMixedBilinearForm *M_uu=NULL, *K_uu=NULL, *K_up=NULL, *K_pu=NULL;
   ParBilinearForm *K_uu_uncon=NULL, *M_uu_uncon=NULL;
   ParLinearForm *uResidual=NULL, *pResidual=NULL;

   //
   // Integrator Coefficients
   //
   ConstantCoefficient *mu, *rho, *one;
   convectiveCoeff *uGradu_c;
   vectorGradientCoeff *Gradu_c;
   GradientGridFunctionCoefficient *Gradp_c;
   VectorGridFunctionCoefficient *u_c;
   totalFluidCoeff *NSM_coeff;

   //
   // Jacobian and Mass matrix storage
   //
   mutable OperatorPtr opMUU, opUU, opUP, opPU;
   mutable BlockOperator *Jacobian=NULL, *Mass=NULL;

   HypreParMatrix Mmat, Kmat;
   HypreParMatrix *T=NULL; // T = M + dt K
   real_t current_dt;

   MINRESSolver  M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver
   
   //
   // Boundary conditions
   //
   Array<int> empty_tdofs;
   Array<int> U_ess_BCTags, P_ess_BCTags;      //Dirch boundary condition tags
   Array<int> U_ess_BCDofs, P_ess_BCDofs;      //Dirch boundary condition dofs
   mutable BlockVector xBlock, rhsBlock;       //auxillary vector(s)
   mutable ParGridFunction *u_gf, *p_gf;       //Coefficient GridFunctions
   mutable ParGridFunction *uBDR_gf, *pBDR_gf; //Dirch boundary condition values
   Array<VectorCoefficient*> U_ess_BCs;        //Velocity boundary conditions
   Array<Coefficient*>       P_ess_BCs;        //Pressure boundary conditions

   void updateBCVals() const;
   void reassembleJacobian() const;
public:
   navierStokesOper(ParFiniteElementSpace *f_u, ParFiniteElementSpace *f_p, int &dim_, real_t &dt);

   virtual void Mult(const Vector &u, Vector &du_dt) const;

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k);

   virtual ~navierStokesOper();

   //Set the boundary condition arrays
   void SetPressureBCTags(Array<int> & P_ess_BCTags_);
   void SetVelocityBCTags(Array<int> & U_ess_BCTags_);
};


/*****************************************\
!
! Setup the BC tag arrays
!
/*****************************************/
void navierStokesOper::SetVelocityBCTags(Array<int> & U_ess_BCTags_){
  U_ess_BCTags.SetSize(U_ess_BCTags_.Size());
  for(int I=0; I<U_ess_BCTags_.Size(); I++) U_ess_BCTags[I] = U_ess_BCTags_[I];
  fesU->GetEssentialTrueDofs(U_ess_BCTags, U_ess_BCDofs);
};


/*****************************************\
!
! Constructor builds the integrators for
! the problem class
!
/*****************************************/
navierStokesOper::navierStokesOper(ParFiniteElementSpace *f_u, ParFiniteElementSpace *f_p, int &dim_, real_t &dt)
   : TimeDependentOperator(f_u->GetTrueVSize() + f_p->GetTrueVSize(), (real_t) 0.0), fesU(f_u), fesP(f_p),
     btoffs(3), T(NULL), current_dt(0.0), M_solver(f_u->GetComm()), T_solver(f_u->GetComm())
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
   mu        = new ConstantCoefficient(1.0);
   rho       = new ConstantCoefficient(1.0);
   uGradu_c  = new convectiveCoeff(dim,u_gf);
   Gradu_c   = new vectorGradientCoeff(dim,u_gf);
   Gradp_c   = new GradientGridFunctionCoefficient(p_gf);
   u_c       = new VectorGridFunctionCoefficient(u_gf);
   NSM_coeff = new totalFluidCoeff(dim, u_gf, p_gf);


   // Within the forms u,p represent the pressure
   // and velocity fields while v,h represent the
   // test functions that correspond to those fields

   //
   // Residual Forms
   //

   //Momentum equation (du/dt)
   uResidual = new ParLinearForm(fesU);
   uResidual->AddDomainIntegrator(new gradContractionIntegrator(*NSM_coeff,dim) ); //-div(uu^T - mu Grad(u) + pI)

   //Continuity equation (dp/dt = 0.0)
   pResidual = new ParLinearForm(fesP);
   pResidual->AddDomainIntegrator(new DomainLFGradIntegrator(*u_c));        // div(u)          Continuity Error

   //
   // Jacobian Forms
   //

   //Velocity Momentum
   K_uu = new ParMixedBilinearForm(fesU,fesU);
   K_uu->AddDomainIntegrator(new VectorDiffusionIntegrator(*mu));  // -div(mu grad(u))  Diffusion
   K_uu->AddDomainIntegrator(new VectorMassIntegrator(*Gradu_c));  // [v, grad(u) v] Convection0
//   K_uu->AddDomainIntegrator(new VectorMassIntegrator());  // [v, u grad(v)] Convection1
//*K_uu_uncon=NULL, *M_uu_uncon=NULL;


   //Pressure Momentum
   K_up = new ParMixedBilinearForm(fesU,fesP);
   K_up->AddDomainIntegrator(new GradientIntegrator(*one));         // [v,grad(h)]  Pressure gradient

   //Continuity
   K_pu = new ParMixedBilinearForm(fesP,fesU);
   K_pu->AddDomainIntegrator(new VectorDivergenceIntegrator(*one)); // [div(v),h]  Continuity Error

   //
   // Build the Mass Matrix and forms
   //
   M_uu = new ParMixedBilinearForm(fesU,fesU);
   M_uu->AddDomainIntegrator(new VectorMassIntegrator(*rho));
   M_uu->Assemble();
   M_uu->FormRectangularSystemMatrix(empty_tdofs, empty_tdofs, opMUU);
   Mass  = new BlockOperator(btoffs);
   Mass->SetBlock(0,0, opMUU.Ptr(), 1.0);
   reassembleJacobian();

   //
   // Set the mass solver
   //
   unsigned maxIter=200;
   double atol=1.00e-10, rtol=1.00e-7;
   M_solver.SetOperator(*Mass);
   M_solver.SetAbsTol(atol);
   M_solver.SetRelTol(rtol);
   M_solver.SetMaxIter(maxIter);
};

/*****************************************\
!
! Update the Jacobian block Operator
!
!
/*****************************************/
void navierStokesOper::reassembleJacobian() const{
  //Update and re-assemble the matrices
  K_uu->Update();
  K_up->Update();
  K_pu->Update();
/*
  K_uu->Assemble();
  K_up->Assemble();
  K_pu->Assemble();*/
/*
  //Form the matrix operators
  K_uu->FormRectangularSystemMatrix(empty_tdofs, U_ess_BCDofs, opUU );
  K_up->FormRectangularSystemMatrix(empty_tdofs, U_ess_BCDofs, opUP );
  K_pu->FormRectangularSystemMatrix(empty_tdofs, P_ess_BCDofs, opPU );
*/
  //Set the block matrix operator
/*  if(Jacobian != NULL ) delete Jacobian;
  Jacobian = NULL;
  Jacobian = new BlockOperator(btoffs);
  Jacobian->SetBlock(0,0, opUU);
  Jacobian->SetBlock(0,1, opUP, -1.0);
  Jacobian->SetBlock(1,0, opPU, -1.0);*/
};

/*****************************************\
!
! Update the Dirchelet boundary condition
! GridFunctions
!
\*****************************************/
void navierStokesOper::updateBCVals() const{
  if( U_ess_BCTags.Size() != 0){
    for(int I=0; I<U_ess_BCTags.Size(); I++){
      if(U_ess_BCTags[I] == 1){
        Array<int> tmp_BDR_tags( U_ess_BCTags.Size() );
        tmp_BDR_tags = 0;
        tmp_BDR_tags[I] = 1;
        uBDR_gf->ProjectBdrCoefficient(*(U_ess_BCs[I]), tmp_BDR_tags);
      }
    }
  }
  if( P_ess_BCTags.Size() != 0){
    for(int I=0; I<P_ess_BCTags.Size(); I++){
      if(P_ess_BCTags[I] == 1){
        Array<int> tmp_BDR_tags( P_ess_BCTags.Size() );
        tmp_BDR_tags = 0;
        tmp_BDR_tags[I] = 1;
        pBDR_gf->ProjectBdrCoefficient(*(P_ess_BCs[I]), tmp_BDR_tags);
      }
    }
  }
};


/*****************************************\
!
! Calculates the du/dt entry by computing
! M du/dt = K u
! then inverts M to give
! du/dt = Minv K u
!
/*****************************************/
void navierStokesOper::Mult(const Vector &u, Vector &du_dt) const{
  //Copy in the current solution
  copyVec(u, xBlock);

  //Apply Dirchelet BC values
//  updateBCVals();
//  applyDirchValues(*uBDR_gf, xBlock.GetBlock(0), U_ess_BCDofs);
//  applyDirchValues(*pBDR_gf, xBlock.GetBlock(1), P_ess_BCDofs);


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

  //Update the Jacobian
//  reassembleJacobian();

  //Copy out the residual
  copyVec(rhsBlock, du_dt);
};

/*****************************************\
!
! Solve the Backward-Euler equation: 
! k = f(u + dt*k, t), for the unknown k.
! This is the only requirement for 
! high-order SDIRK implicit integration.
!    Solve the equation:
!     du_dt = M^{-1}*[-K(u + dt*du_dt)]
!
/*****************************************/
void navierStokesOper::ImplicitSolve(const real_t dt, const Vector &u, Vector &k){
 // copyVec(u, xBlock);
   if (!T){
//      T = Add(1.0, dynamic_cast<HypreParMatrix*>(Mass), dt, *dynamic_cast<HypreParMatrix*>(Jacobian));
    //  current_dt = dt;
   //   T_solver.SetOperator(*T);
   }
  // MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
  Mult(u, xBlock);
  xBlock *= dt;
  xBlock += u;
//  M_solver.Mult(xBlock,k);
//   xBlock.Neg();
 //  T_solver.Mult(xBlock, k);
  copyVec(xBlock, k);
};


/*****************************************\
!
! Solve the Backward-Euler equation: 
! k = f(u + dt*k, t), for the unknown k.
! This is the only requirement for 
! high-order SDIRK implicit integration.
!
/*****************************************/
navierStokesOper::~navierStokesOper(){
   delete K_uu, K_up, K_pu;
   delete uResidual, pResidual;
};









