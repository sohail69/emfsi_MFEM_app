#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

#include "../utilityFuncs.hpp"

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
   ParFiniteElementSpace *fesU, *fesP;
   ParMixedBilinearForm *K_uu, *K_up, *K_pu, *K_pp;
   ParBilinearForm *K_uu_uncon, *K_pp_uncon;

   ParLinearForm *uResidual, *pResidual;
   BlockOperator *Jacobian;

   Array<int> empty_tdofs;
   ParFiniteElementSpace &fespace;
   ParBilinearForm *M=NULL, *K=NULL;

   HypreParMatrix Mmat, Kmat;
   HypreParMatrix *T=NULL; // T = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver
   
   //
   // Boundary conditions
   //
   Array<int> U_ess_BCTags, P_ess_BCTags;      //Dirch boundary condition tags
   Array<int> U_ess_BCDofs, P_ess_BCDofs;      //Dirch boundary condition dofs
   mutable BlockVector xBlock, rhsBlock;       //auxillary vector
   mutable ParGridFunction *u_gf, *p_gf;       //Coefficient GridFunctions
   mutable ParGridFunction *uBDR_gf, *pBDR_gf; //Dirch boundary condition values
   Array<VectorCoefficient*> U_ess_BCs;        //Velocity boundary conditions
   Array<Coefficient*>       P_ess_BCs;        //Pressure boundary conditions

   void updateBCVals();
public:
   navierStokesOper(ParFiniteElementSpace &f, MatrixCoefficient &D_ij_, real_t &dt);

   virtual void Mult(const Vector &u, Vector &du_dt) const;

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k);

   virtual ~navierStokesOper();
};


/*****************************************\
!
! Constructor builds the matrices for the
! problem class
!
/*****************************************/
navierStokesOper::navierStokesOper(ParFiniteElementSpace &f)
   : TimeDependentOperator(f.GetTrueVSize(), (real_t) 0.0), fespace(f),
     M(NULL), K(NULL), T(NULL), current_dt(0.0),
     M_solver(f.GetComm()), T_solver(f.GetComm()), z(height)
{
   const real_t rel_tol = 1e-8;

   M = new ParBilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator());
   M->Assemble(0); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(empty_tdofs, Mmat);

   K = new ParBilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(D_ij));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(empty_tdofs, Kmat);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);
};


/*****************************************\
!
! Update the Dirchelet boundary condition
! GridFunctions
!
/*****************************************/
void updateBCVals(){
  if( U_ess_BCTags.Size() != 0){
    for(int I=0; I<U_ess_BCTags.Size(); I++){
      if(U_ess_BCTags[I] = 1){
        Array<int> tmp_BDR_tags( U_ess_BCTags.Size() );
        tmp_BDR_tags = 0;
        tmp_BDR_tags[I] = 1;
        uBDR_gf->ProjectBdrCoefficient(*(U_ess_BCs[I]), tmp_BDR_tags);
      }
    }
  }
  if( P_ess_BCTags.Size() != 0){
    for(int I=0; I<P_ess_BCTags.Size(); I++){
      if(P_ess_BCTags[I] = 1){
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
  //Copy the vector into blockVector
  copyVec(u, xBlock);

  //Apply Dirchelet BC values
  updateBCVals();
  applyDirchValues(*uBDR_gf, xBlock.GetBlock(0), U_ess_BCDofs);
  applyDirchValues(*pBDR_gf, xBlock.GetBlock(1), P_ess_BCDofs);

  //Update the GridFunctions
  u_gf->Distribute(xBlock.GetBlock(0));
  p_gf->Distribute(xBlock.GetBlock(1));

  //Reassemble the residuals
  uResidual->ParallelAssemble(&rhsBlock.GetBlock(0));
  pResidual->ParallelAssemble(&rhsBlock.GetBlock(1));

  //eliminate Dirchelet BC residuals
  rhsBlock.GetBlock(0).SetSubVector(U_ess_BCDofs,0.00);
  rhsBlock.GetBlock(1).SetSubVector(P_ess_BCDofs,0.00);

  //Update the Jacobian
  if(false){
  }

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
     du_dt = M^{-1}*[-K(u + dt*du_dt)]
!
/*****************************************/
void navierStokesOper::ImplicitSolve(const real_t dt, const Vector &u, Vector &k){
   if (!T){
      T = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      T_solver.SetOperator(*T);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   T_solver.Mult(z, k);
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
   delete T;
   delete M;
   delete K;
};
