#pragma once
#include "mfem.hpp"


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
   Array<int>    & _ess_bcs_markers;
   Array<double> & _BC_Vals;
   ParFiniteElementSpace & _fespace;
   ParBilinearForm *_K;
   HypreParMatrix _Kmat;

   CGSolver _K_solver;    // Krylov solver for inverting the mass matrix K
   HypreSmoother _K_prec; // Preconditioner for the diffusion matrix K

   mutable ParGridFunction *z; // auxiliary GFunc used for BC's

public:
   //Constructs fibre map operator
   fibreMapPoissonOper(ParFiniteElementSpace &f, Array<int> & ess_bcs_markers, Array<double> & BC_Vals);

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
                                       , Array<int> & ess_bcs_markers
                                       , Array<double> & BC_Vals): 
                                       Operator(f.GetTrueVSize()), _ess_bcs_markers(ess_bcs_markers)
                                       , _BC_Vals(BC_Vals), _fespace(f)
{
  //Construct the matrix Operators
  Array<int> ess_bcs_tdofs;
  _K = new ParBilinearForm(&_fespace);
  _K->AddDomainIntegrator(new DiffusionIntegrator);
  _K->Assemble();
  _fespace.GetEssentialTrueDofs(_ess_bcs_markers, ess_bcs_tdofs);
  _K->FormSystemMatrix(ess_bcs_tdofs, _Kmat);

  //Construct the solver
  const real_t rel_tol = 1e-8;
  _K_solver.iterative_mode = false;
  _K_solver.SetRelTol(rel_tol);
  _K_solver.SetAbsTol(0.0);
  _K_solver.SetMaxIter(100);
  _K_solver.SetPrintLevel(0);
  _K_prec.SetType(HypreSmoother::Jacobi);
  _K_solver.SetPreconditioner(_K_prec);
  _K_solver.SetOperator(_Kmat);

  //Construct the temporary gridFunction
  z = new ParGridFunction(&_fespace);
};

/*****************************************\
!
! Solves the Poisson problem
!
/*****************************************/
void fibreMapPoissonOper::Mult(const Vector &x, Vector &y) const
{
  //Set the tmp Gridfunction before BCs
  z->Distribute(&x);

  //Apply the BC's
  Array<int> ess_bcs_tmp(_BC_Vals.Size());
  for(int I=0; I<_BC_Vals.Size(); I++){
    if(_ess_bcs_markers[I] != 0){
      ess_bcs_tmp = 0;
      ess_bcs_tmp[I] = 1;
      ConstantCoefficient BDR_coeff(_BC_Vals[I]);
      z->ProjectBdrCoefficient(BDR_coeff, ess_bcs_tmp);
    }
  }

  //Solve the problem
  _K_solver.Mult(*z,y);
};

/*****************************************\
!
! Destroys the fibre map operator
!
/*****************************************/
fibreMapPoissonOper::~fibreMapPoissonOper(){
  delete _K;
  delete z;
};


/*****************************************\
!
! Constructs the actual fibre maps and puts
! them into a gridfunction
!
/*****************************************/
void fibreMapPoissonOper::getFibreMapGFuncs(const Vector &x, Array<ParGridFunction> FibreFields){};
