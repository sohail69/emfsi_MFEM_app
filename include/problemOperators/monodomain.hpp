#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

/*****************************************\
!
! Solves the orthotropic-diffusion
! coefficient heat equation for a single
! time-step (This is for usage in the 
! mono-domain equation which is equivalent)
!
/*****************************************/
class monodomainOper : public TimeDependentOperator
{
protected:
   Array<int> empty_tdofs;
   ParFiniteElementSpace &fespace;

   MatrixCoefficient & D_ij;
   ParBilinearForm *M=NULL, *K=NULL;

   HypreParMatrix Mmat, Kmat;
   HypreParMatrix *T=NULL; // T = M + dt K
   real_t current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver
   mutable Vector z; // auxiliary vector

public:
   monodomainOper(ParFiniteElementSpace &f, MatrixCoefficient &D_ij_);

   virtual void Mult(const Vector &u, Vector &du_dt) const;

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k);

   virtual ~monodomainOper();
};


/*****************************************\
!
! Constructor builds the matrices for the
! problem class
!
/*****************************************/
monodomainOper::monodomainOper(ParFiniteElementSpace &f, MatrixCoefficient &D_ij_)
   : TimeDependentOperator(f.GetTrueVSize(), (real_t) 0.0), fespace(f),
     M(NULL), K(NULL), T(NULL), current_dt(0.00), D_ij(D_ij_),
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
! Calculates the du/dt entry by computing
! M du/dt = K u
! then inverts M to give
! du/dt = Minv K u
!
/*****************************************/
void monodomainOper::Mult(const Vector &u, Vector &du_dt) const{
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   M_solver.Mult(z, du_dt);
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
void monodomainOper::ImplicitSolve(const real_t dt, const Vector &u, Vector &k){
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
monodomainOper::~monodomainOper(){
   delete T;
   delete M;
   delete K;
};
