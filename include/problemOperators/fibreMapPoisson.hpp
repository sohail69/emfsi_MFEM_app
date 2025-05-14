#pragma once

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

/*****************************************\
!
! Solves the poisson equation and generates
! the necesssary fibre maps needed for other
! calculations (Active strain, orthotropic
! diffusion, etc...)
!
/*****************************************/
class fibreMapPoissonOper : public Operator
{
protected:
   Array<int>    & _ess_bcs;
   Array<double> & _BC_Vals
   ParFiniteElementSpace & _fespace;
   ParBilinearForm *_K;
   HypreParMatrix _Kmat;

   CGSolver _K_solver;    // Krylov solver for inverting the mass matrix K
   HypreSmoother _K_prec; // Preconditioner for the diffusion matrix K

   real_t alpha, kappa;
   mutable Vector z; // auxiliary vector

public:
   //Constructs fibre map operator
   fibreMapPoissonOper(ParFiniteElementSpace &f, Array<int> & ess_bcs, Array<double> & BC_Vals);

   //Solves the Poisson equation
   virtual void Mult(const Vector &x, Vector &y) const;

   //Destroys the fibre map operator
   virtual ~fibreMapPoissonOper();

   //Gets a grid-function array of the vector basis
   //needed to define the fibre fields
   void getFibreMapGFuncs(const Vector &x, Array<ParGridFunction> FibreFields);
};

/*****************************************\
!
! Constructs the fibre map operator
! Assumes a H1 FE-space
!
/*****************************************/
fibreMapPoissonOper::fibreMapPoissonOper(ParFiniteElementSpace &f
                                       , Array<int> & ess_bcs
                                       , Array<double> & BC_Vals): 
                                       _ess_bcs(ess_bcs), _BC_Vals(BC_Vals), _fespace(f)
{
  _K = new ParBilinearForm(&_fespace);
  _K->(new DiffusionIntegrator);
  _K->Assemble();
  _K->FormSystemMatrix(_ess_bcs, _Kmat);

   const real_t rel_tol = 1e-8;
   _K_solver.iterative_mode = false;
   _K_solver.SetRelTol(rel_tol);
   _K_solver.SetAbsTol(0.0);
   _K_solver.SetMaxIter(100);
   _K_solver.SetPrintLevel(0);
   _K_prec.SetType(HypreSmoother::Jacobi);
   _K_solver.SetPreconditioner(_K_prec);
   _K_solver.SetOperator(_Kmat);
};

/*****************************************\
!
! Solves the Poisson problem
!
/*****************************************/
void fibreMapPoissonOper::Mult(const Vector &x, Vector &y) const
{
  

  _K_solver->Mult(x,y);
};

/*****************************************\
!
! Destroys the fibre map operator
!
/*****************************************/
fibreMapPoissonOper::~fibreMapPoissonOper(){



};
